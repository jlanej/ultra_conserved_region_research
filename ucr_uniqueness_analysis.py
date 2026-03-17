#!/usr/bin/env python3
"""UCR Uniqueness Analysis: end-to-end liftover → validation → k-mer overlap.

Orchestrates the full UCR analysis pipeline:
  1. Lifts UCR regions from hg38 to T2T-CHM13v2.0 (via convert_ucr_to_t2t.py)
  2. Validates liftover with per-base alignment (via validate_liftover.py)
  3. Intersects each lifted UCR with T2T unique-mappability tracks
  4. Reports per-region uniqueness with full and non-unique sequences

Produces:
    ucr_uniqueness_report.tsv       – per-region, per-kmer uniqueness metrics
    ucr_uniqueness_summary.json     – overall pipeline summary

Usage (Apptainer one-liner):
    apptainer exec --bind "$(pwd):/output" \\
        docker://ghcr.io/jlanej/ultra_conserved_region_research:latest \\
        python /app/ucr_uniqueness_analysis.py

Usage (local):
    python ucr_uniqueness_analysis.py
    python ucr_uniqueness_analysis.py --kmers 24,50,100
"""

import argparse
import bisect
import csv
import datetime
import json
import os
import stat
import subprocess
import sys
import urllib.request

from Bio import SeqIO

# ──────────────────────────── Constants ────────────────────────────────

ASSEMBLY = "hs1"
DEFAULT_KMERS = (24, 36, 50, 100, 150, 250)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

HS1_2BIT_URL = (
    "https://hgdownload.cse.ucsc.edu/goldenPath/hs1/bigZips/hs1.2bit"
)
TWOBITTOFA_URL = (
    "https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa"
)
BIGBEDTOBED_URL = (
    "https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed"
)

CONVERT_SCRIPT = os.path.join(SCRIPT_DIR, "convert_ucr_to_t2t.py")
VALIDATE_SCRIPT = os.path.join(SCRIPT_DIR, "validate_liftover.py")


# ──────────────────────────── Helpers ─────────────────────────────────

def download_file(url, filepath, make_executable=False):
    """Download *url* to *filepath* unless it already exists."""
    if os.path.exists(filepath):
        print(f"  {os.path.basename(filepath)} already present – skipping.")
        return
    print(f"  Downloading {os.path.basename(filepath)} …")
    urllib.request.urlretrieve(url, filepath)
    if make_executable:
        st = os.stat(filepath)
        os.chmod(filepath, st.st_mode | stat.S_IEXEC)
    size_mb = os.path.getsize(filepath) / 1_048_576
    print(f"  ✓ {os.path.basename(filepath)} ({size_mb:.1f} MB)")


def kmer_bb_url(k):
    """UCSC download URL for a k-mer mappability bigBed on hs1."""
    return (
        f"https://hgdownload.soe.ucsc.edu/gbdb/{ASSEMBLY}/"
        f"hoffmanMappability/k{k}.Unique.Mappability.bb"
    )


def bb_path(k, outdir):
    return os.path.join(outdir, f"k{k}.Unique.Mappability.bb")


def bed_path(k, outdir):
    return os.path.join(outdir, f"k{k}.unique.bed")


def bigbed_to_bed(bigbedtobed_bin, bb_file, bed_file):
    """Convert a UCSC bigBed to BED."""
    result = subprocess.run(
        [bigbedtobed_bin, bb_file, bed_file],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print(f"bigBedToBed stderr:\n{result.stderr}", file=sys.stderr)
        result.check_returncode()


# ──────────────────────────── BED I/O ─────────────────────────────────

def parse_bed4(filepath):
    """Parse a BED4 file → list of dicts with chrom, start, end, name."""
    regions = []
    if not os.path.exists(filepath):
        return regions
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 4:
                regions.append({
                    "chrom": parts[0],
                    "start": int(parts[1]),
                    "end": int(parts[2]),
                    "name": parts[3],
                })
    return regions


def load_fasta(path):
    """Return ``{record_id: uppercase_sequence}`` from a FASTA file."""
    return {
        rec.id: str(rec.seq).upper()
        for rec in SeqIO.parse(path, "fasta")
    }


def extract_sequences(twobittofa_bin, twobit_file, bed_file, fasta_out):
    """Run ``twoBitToFa`` to extract BED regions from a 2bit genome."""
    print(f"  Extracting sequences from {os.path.basename(twobit_file)} …")
    cmd = [twobittofa_bin, f"-bed={bed_file}", twobit_file, fasta_out]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"twoBitToFa error:\n{result.stderr}", file=sys.stderr)
        sys.exit(1)
    n = sum(1 for _ in SeqIO.parse(fasta_out, "fasta"))
    print(f"  ✓ {n} sequences → {os.path.basename(fasta_out)}")


# ──────────────── Unique interval loading ─────────────────────────────

def load_unique_intervals(bed_file):
    """Load unique intervals from a k-mer BED file.

    Applies the same strict-score filtering as ``compute_unique_fraction.py``:
    when a track has mixed scores, only intervals at the **maximum** score
    are kept as truly unique.

    Returns
    -------
    by_chrom : dict
        ``{chrom: [(start, end), ...]}`` with merged, sorted intervals.
    strict_filter : bool
        Whether non-maximal-score intervals were excluded.
    """
    records = []
    numeric_scores = []
    with open(bed_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            if end <= start:
                continue
            score = None
            if len(fields) >= 5:
                try:
                    score = float(fields[4])
                    numeric_scores.append(score)
                except ValueError:
                    pass
            elif len(fields) >= 4:
                try:
                    score = float(fields[3])
                    numeric_scores.append(score)
                except ValueError:
                    pass
            records.append((chrom, start, end, score))

    strict_filter = False
    threshold = None
    if numeric_scores:
        min_s, max_s = min(numeric_scores), max(numeric_scores)
        if max_s > min_s:
            strict_filter = True
            threshold = max_s

    by_chrom = {}
    for chrom, start, end, score in records:
        if strict_filter and (score is None or score < threshold):
            continue
        by_chrom.setdefault(chrom, []).append((start, end))

    for chrom in by_chrom:
        by_chrom[chrom] = merge_intervals(sorted(by_chrom[chrom]))

    return by_chrom, strict_filter


# ──────────────── Interval utilities ──────────────────────────────────

def merge_intervals(sorted_intervals):
    """Merge overlapping/adjacent intervals.  Input must be sorted by start."""
    if not sorted_intervals:
        return []
    merged = [list(sorted_intervals[0])]
    for s, e in sorted_intervals[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [tuple(iv) for iv in merged]


def compute_overlap_bp(ucr_start, ucr_end, merged_intervals):
    """Count bp in ``[ucr_start, ucr_end)`` covered by *merged_intervals*.

    Parameters
    ----------
    ucr_start, ucr_end : int
        Half-open UCR interval.
    merged_intervals : list of (int, int)
        Sorted, non-overlapping intervals (e.g. unique regions).

    Returns
    -------
    int
        Number of base pairs covered.
    """
    if not merged_intervals:
        return 0
    starts = [iv[0] for iv in merged_intervals]
    right = bisect.bisect_left(starts, ucr_end)
    left = bisect.bisect_right(starts, ucr_start) - 1
    if left < 0:
        left = 0
    total = 0
    for i in range(left, right):
        s, e = merged_intervals[i]
        os_ = max(s, ucr_start)
        oe = min(e, ucr_end)
        if os_ < oe:
            total += oe - os_
    return total


def get_non_unique_intervals(ucr_start, ucr_end, merged_intervals):
    """Return sub-intervals of ``[ucr_start, ucr_end)`` NOT covered by unique.

    Returns a list of ``(start, end)`` tuples in absolute coordinates.
    """
    if not merged_intervals:
        return [(ucr_start, ucr_end)]

    starts = [iv[0] for iv in merged_intervals]
    right = bisect.bisect_left(starts, ucr_end)
    left = bisect.bisect_right(starts, ucr_start) - 1
    if left < 0:
        left = 0

    clipped = []
    for i in range(left, right):
        s, e = merged_intervals[i]
        cs = max(s, ucr_start)
        ce = min(e, ucr_end)
        if cs < ce:
            clipped.append((cs, ce))

    non_unique = []
    prev = ucr_start
    for s, e in clipped:
        if prev < s:
            non_unique.append((prev, s))
        prev = max(prev, e)
    if prev < ucr_end:
        non_unique.append((prev, ucr_end))

    return non_unique


def extract_non_unique_seq(ucr_start, ucr_sequence, non_unique_intervals):
    """Extract bases from *ucr_sequence* at non-unique positions.

    Parameters
    ----------
    ucr_start : int
        Absolute genomic start coordinate of the UCR.
    ucr_sequence : str
        Full UCR sequence.
    non_unique_intervals : list of (int, int)
        Absolute genomic coordinates of non-unique sub-intervals.

    Returns
    -------
    str
        Concatenated non-unique bases, or ``""`` if all unique.
    """
    parts = []
    for s, e in non_unique_intervals:
        rel_s = s - ucr_start
        rel_e = e - ucr_start
        parts.append(ucr_sequence[rel_s:rel_e])
    return "".join(parts)


# ────────────── Validation results loader ─────────────────────────────

def load_validation_results(report_tsv):
    """Parse ``ucr_alignment_report.tsv`` → ``{ucr_id: {field: value}}``."""
    results = {}
    if not os.path.exists(report_tsv):
        return results
    with open(report_tsv) as fh:
        reader = csv.DictReader(
            (line for line in fh if not line.startswith("##")),
            delimiter="\t",
        )
        for row in reader:
            results[row["ucr_id"]] = row
    return results


# ──────────────────── Report writers ──────────────────────────────────

def write_per_region_report(path, rows):
    """Write per-region, per-kmer uniqueness report as a TSV."""
    fields = [
        "ucr_id",
        "t2t_chrom", "t2t_start", "t2t_end", "ucr_length",
        "identity_pct",
        "kmer",
        "unique_bp", "unique_fraction",
        "non_unique_bp", "non_unique_fraction",
        "ucr_sequence", "non_unique_sequence",
    ]
    with open(path, "w", newline="") as fh:
        fh.write("## UCR Uniqueness Report\n")
        fh.write(
            f"## Generated: "
            f"{datetime.datetime.now(datetime.timezone.utc).isoformat()}\n"
        )
        fh.write("## Assembly: T2T-CHM13v2.0 (hs1)\n")
        fh.write("## Unique mappability source: UCSC hoffmanMappability\n")
        fh.write("##\n")
        writer = csv.DictWriter(
            fh, fieldnames=fields, delimiter="\t", extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            out = dict(row)
            out["unique_fraction"] = f"{row['unique_fraction']:.6f}"
            out["non_unique_fraction"] = f"{row['non_unique_fraction']:.6f}"
            writer.writerow(out)

    print(f"  ✓ {os.path.basename(path)}")


def write_summary(path, per_region_rows, validation, kmers,
                   total_ucrs_lifted):
    """Write overall uniqueness summary as JSON."""
    per_kmer = {}
    for k in kmers:
        k_rows = [r for r in per_region_rows if r["kmer"] == k]
        total_ucr_bp = sum(r["ucr_length"] for r in k_rows)
        total_unique_bp = sum(r["unique_bp"] for r in k_rows)
        fully_unique = sum(
            1 for r in k_rows if r["unique_bp"] == r["ucr_length"]
        )
        partially_unique = sum(
            1 for r in k_rows
            if 0 < r["unique_bp"] < r["ucr_length"]
        )
        not_unique = sum(1 for r in k_rows if r["unique_bp"] == 0)
        per_kmer[str(k)] = {
            "total_ucr_bp": total_ucr_bp,
            "total_unique_bp": total_unique_bp,
            "overall_unique_fraction": (
                round(total_unique_bp / total_ucr_bp, 6)
                if total_ucr_bp else 0.0
            ),
            "fully_unique_ucrs": fully_unique,
            "partially_unique_ucrs": partially_unique,
            "not_unique_ucrs": not_unique,
        }

    val_summary = {}
    if validation:
        total = len(validation)
        identical = sum(
            1 for v in validation.values()
            if v.get("identical", "").upper() == "YES"
        )
        pcts = []
        for v in validation.values():
            try:
                pcts.append(float(v.get("identity_pct", 0)))
            except (ValueError, TypeError):
                pass
        val_summary = {
            "paired_ucrs": total,
            "identical_sequences": identical,
            "average_identity_pct": round(
                sum(pcts) / len(pcts), 2
            ) if pcts else 0.0,
        }

    payload = {
        "generated_at_utc": datetime.datetime.now(
            datetime.timezone.utc
        ).isoformat(),
        "assembly": f"{ASSEMBLY} (T2T-CHM13v2.0)",
        "kmers_analyzed": list(kmers),
        "total_ucrs_lifted_over": total_ucrs_lifted,
        "per_kmer_summary": per_kmer,
        "validation": val_summary,
    }
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2)
    print(f"  ✓ {os.path.basename(path)}")


# ─────────────────────────── Main ─────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=(
            "End-to-end UCR uniqueness analysis: liftover → validation → "
            "k-mer overlap on T2T-CHM13v2.0."
        ),
    )
    parser.add_argument(
        "--output-dir",
        default=os.environ.get("OUTPUT_DIR", os.getcwd()),
        help="Output directory (default: OUTPUT_DIR env var, or cwd).",
    )
    parser.add_argument(
        "--kmers",
        default=",".join(str(k) for k in DEFAULT_KMERS),
        help=(
            "Comma-separated k-mer sizes to evaluate "
            f"(default: {','.join(str(k) for k in DEFAULT_KMERS)})."
        ),
    )
    parser.add_argument(
        "--skip-validation", action="store_true",
        help="Skip the validate_liftover.py step.",
    )
    args = parser.parse_args()
    outdir = args.output_dir
    kmers = tuple(
        int(k.strip()) for k in args.kmers.split(",") if k.strip()
    )
    os.makedirs(outdir, exist_ok=True)

    # Derived paths
    t2t_bed = os.path.join(outdir, "ucr_t2t_chm13.bed")
    t2t_fasta = os.path.join(outdir, "ucr_sequences_t2t.fa")
    hs1_2bit = os.path.join(outdir, "hs1.2bit")
    report_tsv = os.path.join(outdir, "ucr_uniqueness_report.tsv")
    summary_json = os.path.join(outdir, "ucr_uniqueness_summary.json")
    alignment_report = os.path.join(outdir, "ucr_alignment_report.tsv")

    # Tool paths – prefer container-baked binaries next to the script
    _container_bb2bed = os.path.join(SCRIPT_DIR, "bigBedToBed")
    bigbedtobed_bin = (
        _container_bb2bed
        if os.path.exists(_container_bb2bed)
        else os.path.join(outdir, "bigBedToBed")
    )
    _container_2bitfa = os.path.join(SCRIPT_DIR, "twoBitToFa")
    twobittofa_bin = (
        _container_2bitfa
        if os.path.exists(_container_2bitfa)
        else os.path.join(outdir, "twoBitToFa")
    )

    # ── Step 1: Liftover ──
    print("\n═══ Step 1: Liftover (hg38 → T2T-CHM13) ═══")
    if os.path.exists(t2t_bed):
        print(
            f"  {os.path.basename(t2t_bed)} already exists – "
            "skipping liftover."
        )
    else:
        print("  Running convert_ucr_to_t2t.py …")
        env = {**os.environ, "OUTPUT_DIR": outdir}
        subprocess.run(
            [sys.executable, CONVERT_SCRIPT], check=True, env=env,
        )

    # ── Step 2: Validation ──
    print("\n═══ Step 2: Sequence validation ═══")
    if args.skip_validation:
        print("  --skip-validation: skipping.")
    else:
        print("  Running validate_liftover.py …")
        subprocess.run(
            [sys.executable, VALIDATE_SCRIPT, "--output-dir", outdir],
            check=True,
        )

    # ── Step 3: Download tools & data ──
    print("\n═══ Step 3: Download tools and data ═══")
    download_file(BIGBEDTOBED_URL, bigbedtobed_bin, make_executable=True)
    download_file(TWOBITTOFA_URL, twobittofa_bin, make_executable=True)
    download_file(HS1_2BIT_URL, hs1_2bit)

    # ── Step 4: Extract T2T UCR sequences ──
    print("\n═══ Step 4: Extract T2T UCR sequences ═══")
    if os.path.exists(t2t_fasta):
        print(
            f"  {os.path.basename(t2t_fasta)} already exists – reusing."
        )
    else:
        extract_sequences(twobittofa_bin, hs1_2bit, t2t_bed, t2t_fasta)

    ucr_seqs = load_fasta(t2t_fasta)
    ucr_coords = parse_bed4(t2t_bed)
    print(f"  {len(ucr_coords)} UCR regions loaded.")

    # ── Step 5: Load validation results ──
    validation = load_validation_results(alignment_report)
    if validation:
        print(f"  Loaded {len(validation)} validation results.")

    # ── Step 6: K-mer uniqueness analysis ──
    print(
        f"\n═══ Step 5: K-mer uniqueness analysis "
        f"(k={','.join(str(k) for k in kmers)}) ═══"
    )
    per_region_rows = []
    for k in kmers:
        print(f"\n  ── k={k} ──")
        bb = bb_path(k, outdir)
        bed = bed_path(k, outdir)
        download_file(kmer_bb_url(k), bb)
        if not os.path.exists(bed):
            print(f"  Converting bigBed → BED for k={k} …")
            bigbed_to_bed(bigbedtobed_bin, bb, bed)
        else:
            print(f"  {os.path.basename(bed)} already exists – reusing.")

        print(f"  Loading unique intervals for k={k} …")
        unique_by_chrom, strict_filter = load_unique_intervals(bed)
        print(f"    strict score filter: {strict_filter}")

        for ucr in ucr_coords:
            chrom = ucr["chrom"]
            ucr_start = ucr["start"]
            ucr_end = ucr["end"]
            ucr_len = ucr_end - ucr_start
            ucr_id = ucr["name"]

            merged = unique_by_chrom.get(chrom, [])
            unique_bp = compute_overlap_bp(ucr_start, ucr_end, merged)
            non_unique_bp = ucr_len - unique_bp

            ucr_seq = ucr_seqs.get(ucr_id, "")
            non_unique_seq = ""
            if non_unique_bp > 0 and ucr_seq:
                nu_intervals = get_non_unique_intervals(
                    ucr_start, ucr_end, merged,
                )
                non_unique_seq = extract_non_unique_seq(
                    ucr_start, ucr_seq, nu_intervals,
                )

            val = validation.get(ucr_id, {})
            identity_pct = val.get("identity_pct", "")

            per_region_rows.append({
                "ucr_id": ucr_id,
                "t2t_chrom": chrom,
                "t2t_start": ucr_start,
                "t2t_end": ucr_end,
                "ucr_length": ucr_len,
                "identity_pct": identity_pct,
                "kmer": k,
                "unique_bp": unique_bp,
                "unique_fraction": (
                    unique_bp / ucr_len if ucr_len else 0.0
                ),
                "non_unique_bp": non_unique_bp,
                "non_unique_fraction": (
                    non_unique_bp / ucr_len if ucr_len else 0.0
                ),
                "ucr_sequence": ucr_seq,
                "non_unique_sequence": non_unique_seq,
            })

    # ── Step 7: Write reports ──
    print("\n═══ Step 6: Write reports ═══")
    write_per_region_report(report_tsv, per_region_rows)
    write_summary(
        summary_json, per_region_rows, validation, kmers,
        total_ucrs_lifted=len(ucr_coords),
    )

    # ── Print summary to stdout ──
    print(f"\n{'=' * 60}")
    print("  UCR UNIQUENESS ANALYSIS COMPLETE")
    print(f"{'=' * 60}")
    print(f"  UCRs analyzed:          {len(ucr_coords)}")
    print(f"  K-mer sizes:            {', '.join(str(k) for k in kmers)}")
    for k in kmers:
        k_rows = [r for r in per_region_rows if r["kmer"] == k]
        fully_unique = sum(
            1 for r in k_rows if r["unique_bp"] == r["ucr_length"]
        )
        print(
            f"    k={k}: {fully_unique}/{len(k_rows)} UCRs fully unique"
        )
    if validation:
        ident = sum(
            1 for v in validation.values()
            if v.get("identical", "").upper() == "YES"
        )
        print(f"  Validation:             {ident}/{len(validation)} identical")
    print(f"{'=' * 60}")
    print(f"\n  Reports:")
    print(f"    {report_tsv}")
    print(f"    {summary_json}")


if __name__ == "__main__":
    main()
