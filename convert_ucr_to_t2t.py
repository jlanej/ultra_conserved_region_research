import os
import sys
import csv
import hashlib
import urllib.request
import subprocess
import datetime
import pandas as pd
import stat

# --- Configuration ---
# Excel file URL (Supplementary Table S1 from PMC6857462)
# Table S1 lists genes overlapping UCEs from the Stephen et al. 2008 catalogue
# (UCR IDs run to 13,736 for the ≥100 bp set; 2,189 for the ≥200 bp set).
# Rows without a chromosome assignment are excluded during coordinate extraction.
EXCEL_URL = "https://pmc.ncbi.nlm.nih.gov/articles/instance/6857462/bin/Supp_TableS1.xlsx"

# UCSC LiftOver binary and chain file URLs
LIFTOVER_URL = "https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver"
CHAIN_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHs1.over.chain.gz"

# --- Path resolution ---
# The liftOver binary is baked into the container image alongside the script.
# When running outside a container it is downloaded into OUTPUT_DIR on first run.
# Large data files (chain file, Excel) always download at runtime into OUTPUT_DIR.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.environ.get("OUTPUT_DIR", os.getcwd())

# Prefer the liftOver binary next to the script (i.e. inside the container),
# fall back to OUTPUT_DIR for bare-metal / local runs.
_container_liftover = os.path.join(SCRIPT_DIR, "liftOver")
LIFTOVER_BIN = _container_liftover if os.path.exists(_container_liftover) \
    else os.path.join(OUTPUT_DIR, "liftOver")

# Data files always land in OUTPUT_DIR (downloaded at runtime)
CHAIN_FILE = os.path.join(OUTPUT_DIR, "hg38ToHs1.over.chain.gz")

# Prefer the bundled resources file (checked into the repo) if available;
# fall back to downloading into OUTPUT_DIR when running outside the repo.
_bundled_excel = os.path.normpath(os.path.join(SCRIPT_DIR, "resources", "Supp_TableS1.xlsx"))
EXCEL_FILE = _bundled_excel if os.path.exists(_bundled_excel) \
    else os.path.join(OUTPUT_DIR, "Supp_TableS1.xlsx")

HG38_BED = os.path.join(OUTPUT_DIR, "ucr_hg38.bed")
T2T_BED = os.path.join(OUTPUT_DIR, "ucr_t2t_chm13.bed")
UNMAPPED_BED = os.path.join(OUTPUT_DIR, "ucr_unmapped.bed")
AUDIT_REPORT = os.path.join(OUTPUT_DIR, "ucr_liftover_audit.tsv")


def sha256sum(filepath):
    """Return the SHA-256 hex digest of a file."""
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def download_file(url, filename, make_executable=False):
    """Download a file if it does not already exist."""
    if not os.path.exists(filename):
        print(f"Downloading {filename}...")
        urllib.request.urlretrieve(url, filename)
        if make_executable:
            st = os.stat(filename)
            os.chmod(filename, st.st_mode | stat.S_IEXEC)
    else:
        print(f"{filename} already exists, skipping download.")


def normalize_chrom(x):
    """Normalise a chromosome name to the UCSC 'chrN' convention.

    Handles:
    * Plain numbers/letters: 1 → chr1, X → chrX
    * Already-prefixed (any case): chr1 → chr1, CHR1 → chr1, chrx → chrX
    * Mitochondrial: MT → chrM, chrMT → chrM
    * Strips surrounding whitespace
    * NaN / None → returned unchanged (dropped by drop_duplicates / validation)
    """
    if pd.isna(x):
        return x
    s = str(x).strip()
    # Strip any existing 'chr' prefix, case-insensitively, to re-add it cleanly
    if s.lower().startswith('chr'):
        s = s[3:]
    # Normalise mitochondrial chromosome to UCSC convention (M, not MT)
    if s.upper() == 'MT':
        s = 'M'
    # Uppercase non-numeric identifiers so X, Y, M are consistently capitalised
    if not s.isdigit():
        s = s.upper()
    return f'chr{s}'


def extract_coordinates():
    """Read Table S1, extract unique UCR coordinates, and write hg38 BED."""
    print(f"Reading {EXCEL_FILE}...")
    df = pd.read_excel(EXCEL_FILE, skiprows=1, engine='openpyxl')

    # Table S1 may list the same UCR multiple times if it overlaps multiple genes
    ucr_df = df[['Chr.', 'UCR start (bp)', 'UCR end (bp)', 'UCR ID']].drop_duplicates()

    # Drop rows where the chromosome is missing (gene-level rows without a
    # chromosome assignment produce NaN in the Chr. column)
    ucr_df = ucr_df.dropna(subset=['Chr.'])

    # Normalise chromosome names to UCSC 'chrN' convention
    ucr_df['Chr.'] = ucr_df['Chr.'].apply(normalize_chrom)

    # BED format: 0-based start, 1-based end
    ucr_df['BED_start'] = ucr_df['UCR start (bp)'].astype(int) - 1
    ucr_df['BED_end'] = ucr_df['UCR end (bp)'].astype(int)

    bed_df = ucr_df[['Chr.', 'BED_start', 'BED_end', 'UCR ID']]

    print(f"Extracted {len(bed_df)} unique UCR regions.")
    bed_df.to_csv(HG38_BED, sep='\t', header=False, index=False)
    print(f"Saved hg38 coordinates to {HG38_BED}")

    return bed_df


def run_liftover():
    """Execute the UCSC liftOver binary."""
    print("Running UCSC liftOver...")
    cmd = [LIFTOVER_BIN, HG38_BED, CHAIN_FILE, T2T_BED, UNMAPPED_BED]
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"liftOver stderr:\n{result.stderr}", file=sys.stderr)
        result.check_returncode()

    if result.stderr:
        print(f"liftOver warnings:\n{result.stderr}")

    print("LiftOver complete!")
    print(f"  Mapped   → {T2T_BED}")
    print(f"  Unmapped → {UNMAPPED_BED}")


def parse_bed(filepath):
    """Parse a 4-column BED file into a list of dicts."""
    regions = []
    if not os.path.exists(filepath):
        return regions
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) >= 4:
                regions.append({
                    'chrom': parts[0],
                    'start': int(parts[1]),
                    'end': int(parts[2]),
                    'name': parts[3],
                })
    return regions


def parse_unmapped(filepath):
    """Parse the liftOver unmapped file and return a dict keyed by region name.

    The unmapped file has comment lines starting with '#' that explain why the
    region was not mapped, followed by the original BED line.
    """
    unmapped = {}
    if not os.path.exists(filepath):
        return unmapped
    reason = ""
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('#'):
                reason = line.lstrip('#').strip()
            elif line:
                parts = line.split('\t')
                if len(parts) >= 4:
                    unmapped[parts[3]] = reason
                reason = ""
    return unmapped


def generate_audit_report(hg38_df):
    """Produce a comprehensive audit TSV that pairs every input UCR with its
    liftover result, including mapping status, coordinate shifts, size changes,
    and chromosome concordance so the pipeline can be fully QC'd.
    """
    print(f"Generating audit report → {AUDIT_REPORT}")

    # Build lookup tables from output BEDs
    t2t_regions = {r['name']: r for r in parse_bed(T2T_BED)}
    unmapped_reasons = parse_unmapped(UNMAPPED_BED)

    # Prepare input rows from the original hg38 BED
    hg38_regions = parse_bed(HG38_BED)

    mapped_count = 0
    unmapped_count = 0
    chrom_change_count = 0
    size_change_count = 0

    with open(AUDIT_REPORT, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t')

        # --- Header block with provenance metadata ---
        fh.write(f"## UCR Liftover Audit Report\n")
        fh.write(f"## Generated: {datetime.datetime.now(datetime.timezone.utc).isoformat()}\n")
        fh.write(f"## Source: Giacopuzzi et al. 2020 (PMC6857462) Supplementary Table S1\n")
        fh.write(f"## Source assembly: GRCh38 / hg38\n")
        fh.write(f"## Target assembly: T2T-CHM13v2.0 / Hs1\n")
        fh.write(f"## Chain file: hg38ToHs1.over.chain.gz\n")
        fh.write(f"## Excel SHA-256: {sha256sum(EXCEL_FILE)}\n")
        fh.write(f"## Chain SHA-256: {sha256sum(CHAIN_FILE)}\n")
        fh.write(f"##\n")

        # Column header
        writer.writerow([
            'ucr_id',
            'hg38_chrom', 'hg38_start', 'hg38_end', 'hg38_length',
            'status',
            't2t_chrom', 't2t_start', 't2t_end', 't2t_length',
            'delta_length', 'delta_start',
            'chrom_changed', 'reason',
        ])

        for region in hg38_regions:
            name = region['name']
            hg38_len = region['end'] - region['start']

            if name in t2t_regions:
                t = t2t_regions[name]
                t2t_len = t['end'] - t['start']
                delta_len = t2t_len - hg38_len
                delta_start = t['start'] - region['start']
                chrom_changed = (region['chrom'] != t['chrom'])

                if chrom_changed:
                    chrom_change_count += 1
                if delta_len != 0:
                    size_change_count += 1

                writer.writerow([
                    name,
                    region['chrom'], region['start'], region['end'], hg38_len,
                    'MAPPED',
                    t['chrom'], t['start'], t['end'], t2t_len,
                    delta_len, delta_start,
                    'YES' if chrom_changed else 'NO',
                    '',
                ])
                mapped_count += 1
            else:
                reason = unmapped_reasons.get(name, 'Unknown')
                writer.writerow([
                    name,
                    region['chrom'], region['start'], region['end'], hg38_len,
                    'UNMAPPED',
                    '', '', '', '',
                    '', '',
                    '',
                    reason,
                ])
                unmapped_count += 1

        # --- Summary block ---
        total = mapped_count + unmapped_count
        fh.write(f"##\n")
        fh.write(f"## === SUMMARY ===\n")
        fh.write(f"## Total UCRs:           {total}\n")
        fh.write(f"## Mapped:               {mapped_count}\n")
        fh.write(f"## Unmapped:             {unmapped_count}\n")
        fh.write(f"## Chromosome changed:   {chrom_change_count}\n")
        fh.write(f"## Size changed:         {size_change_count}\n")
        fh.write(f"## Mapping rate:         {mapped_count/total*100:.1f}%\n" if total > 0 else "")

    print(f"Audit report written ({mapped_count} mapped, {unmapped_count} unmapped)")

    # Print summary to stdout as well
    print(f"\n{'='*50}")
    print(f"  LIFTOVER SUMMARY")
    print(f"{'='*50}")
    print(f"  Total UCRs:           {total}")
    print(f"  Mapped:               {mapped_count}")
    print(f"  Unmapped:             {unmapped_count}")
    print(f"  Chromosome changed:   {chrom_change_count}")
    print(f"  Size changed:         {size_change_count}")
    if total > 0:
        print(f"  Mapping rate:         {mapped_count/total*100:.1f}%")
    print(f"{'='*50}")

    if unmapped_count > 0:
        print("\n  ⚠  Unmapped regions — review ucr_unmapped.bed and the audit report")
    if chrom_change_count > 0:
        print("  ⚠  Some regions mapped to a different chromosome — review the audit report")
    if size_change_count > 0:
        print(f"  ⚠  {size_change_count} regions changed size during liftover")


if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # 1. Download data files.
    # Excel: use the bundled resources file when available; download otherwise.
    if EXCEL_FILE == os.path.join(OUTPUT_DIR, "Supp_TableS1.xlsx"):
        download_file(EXCEL_URL, EXCEL_FILE)
    else:
        print(f"Using bundled Excel file: {EXCEL_FILE}")
    download_file(CHAIN_URL, CHAIN_FILE)

    # liftOver binary is baked into the container; download only for local runs
    download_file(LIFTOVER_URL, LIFTOVER_BIN, make_executable=True)

    # 2. Extract and format hg38 coordinates
    hg38_df = extract_coordinates()

    # 3. Convert hg38 to T2T-CHM13
    run_liftover()

    # 4. Generate comprehensive audit report
    generate_audit_report(hg38_df)

    print(f"\nAll output files are in: {OUTPUT_DIR}")
