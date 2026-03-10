#!/usr/bin/env python3
"""Validate UCR liftover by extracting and aligning genomic sequences.

Downloads the hg38 and T2T-CHM13v2.0 reference genomes in UCSC 2bit format,
extracts UCR sequences from both assemblies using the BED coordinates produced
by the liftover pipeline, and performs per-base Needleman–Wunsch global
pairwise alignment for every paired region.

Produces four output files:
    ucr_sequences_hg38.fa       – hg38 UCR sequences (FASTA)
    ucr_sequences_t2t.fa        – T2T-CHM13 UCR sequences (FASTA)
    ucr_alignment_report.tsv    – Per-region alignment statistics
    ucr_alignment_details.txt   – Visual alignments for non-identical pairs

Prerequisites:
    Run convert_ucr_to_t2t.py first so that ucr_hg38.bed and
    ucr_t2t_chm13.bed exist in the output directory.

Usage:
    python validate_liftover.py
    python validate_liftover.py --output-dir /path/to/results
"""

import argparse
import csv
import datetime
import os
import stat
import subprocess
import sys
import urllib.request

from Bio import SeqIO
from Bio.Align import PairwiseAligner

# ──────────────────────────────── URLs ────────────────────────────────
HG38_2BIT_URL = "https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit"
HS1_2BIT_URL = "https://hgdownload.cse.ucsc.edu/goldenPath/hs1/bigZips/hs1.2bit"
TWOBITTOFA_URL = (
    "https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa"
)

# Container-baked twoBitToFa lives next to the script; for bare-metal runs
# it is downloaded into OUTPUT_DIR on first use.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_CONTAINER_TWOBITTOFA = os.path.join(SCRIPT_DIR, "twoBitToFa")


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


def extract_sequences(twobittofa_bin, twobit_file, bed_file, fasta_out):
    """Run ``twoBitToFa`` to pull BED regions out of a 2bit genome file.

    Parameters
    ----------
    twobittofa_bin : str
        Path to the ``twoBitToFa`` executable.
    twobit_file : str
        Path to the ``.2bit`` genome file.
    bed_file : str
        BED4 file with regions to extract (chrom, start, end, name).
    fasta_out : str
        Destination FASTA file for the extracted sequences.
    """
    print(f"  Extracting sequences from {os.path.basename(twobit_file)} …")
    cmd = [twobittofa_bin, f"-bed={bed_file}", twobit_file, fasta_out]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"twoBitToFa error:\n{result.stderr}", file=sys.stderr)
        sys.exit(1)
    n = sum(1 for _ in SeqIO.parse(fasta_out, "fasta"))
    print(f"  ✓ {n} sequences → {os.path.basename(fasta_out)}")


def load_fasta(path):
    """Return ``{record_id: uppercase_sequence}`` from a FASTA file."""
    return {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(path, "fasta")}


def load_bed(path):
    """Return ``{name: (chrom, start, end)}`` from a BED4 file."""
    regions = {}
    with open(path) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                regions[parts[3]] = (parts[0], int(parts[1]), int(parts[2]))
    return regions


# ──────────────────── Alignment / comparison ─────────────────────────

def _get_aligned_sequences(alignment):
    """Reconstruct gap-padded aligned strings from an ``Alignment`` object.

    Walks the coordinate path produced by BioPython's ``PairwiseAligner``
    and inserts gap characters (``-``) where one sequence has an insertion
    relative to the other.
    """
    coords = alignment.coordinates
    target_seq = str(alignment.target)
    query_seq = str(alignment.query)
    aligned_target, aligned_query = [], []

    for i in range(len(coords[0]) - 1):
        ts, te = int(coords[0][i]), int(coords[0][i + 1])
        qs, qe = int(coords[1][i]), int(coords[1][i + 1])
        tlen = abs(te - ts)
        qlen = abs(qe - qs)

        if tlen and qlen:
            aligned_target.append(target_seq[ts:te])
            aligned_query.append(query_seq[qs:qe])
        elif tlen:
            aligned_target.append(target_seq[ts:te])
            aligned_query.append("-" * tlen)
        elif qlen:
            aligned_target.append("-" * qlen)
            aligned_query.append(query_seq[qs:qe])

    return "".join(aligned_target), "".join(aligned_query)


def _format_visual(top, mid, bot, width=80):
    """Wrap a three-line alignment into fixed-width blocks."""
    lines = []
    for i in range(0, len(top), width):
        lines.append(f"  hg38  {top[i:i + width]}")
        lines.append(f"        {mid[i:i + width]}")
        lines.append(f"  t2t   {bot[i:i + width]}")
        lines.append("")
    return "\n".join(lines)


def compare_sequences(seq1, seq2):
    """Per-base comparison of two sequences.

    * If the sequences are identical, returns immediately.
    * If they are the same length, performs a direct character comparison
      (no gaps possible).
    * If they differ in length, runs a Needleman–Wunsch global alignment
      (BioPython ``PairwiseAligner``) and counts matches / mismatches / gaps.

    Returns a dict with keys: ``identical``, ``matches``, ``mismatches``,
    ``gaps``, ``aligned_length``, ``identity_pct``, ``visual``.
    """
    if seq1 == seq2:
        n = len(seq1)
        return {
            "identical": True,
            "matches": n,
            "mismatches": 0,
            "gaps": 0,
            "aligned_length": n,
            "identity_pct": 100.0,
            "visual": None,
        }

    if len(seq1) == len(seq2):
        # Same length → direct per-base comparison (no gaps)
        matches = sum(a == b for a, b in zip(seq1, seq2))
        mismatches = len(seq1) - matches
        mid = "".join("|" if a == b else "*" for a, b in zip(seq1, seq2))
        return {
            "identical": False,
            "matches": matches,
            "mismatches": mismatches,
            "gaps": 0,
            "aligned_length": len(seq1),
            "identity_pct": matches / len(seq1) * 100,
            "visual": _format_visual(seq1, mid, seq2),
        }

    # Different lengths → Needleman–Wunsch global alignment
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -3
    aligner.extend_gap_score = -0.5

    aln = aligner.align(seq1, seq2)[0]
    t_aln, q_aln = _get_aligned_sequences(aln)

    matches = mismatches = gaps = 0
    mid_chars = []
    for a, b in zip(t_aln, q_aln):
        if a == "-" or b == "-":
            gaps += 1
            mid_chars.append(" ")
        elif a == b:
            matches += 1
            mid_chars.append("|")
        else:
            mismatches += 1
            mid_chars.append("*")

    mid = "".join(mid_chars)
    aligned_length = matches + mismatches + gaps
    return {
        "identical": False,
        "matches": matches,
        "mismatches": mismatches,
        "gaps": gaps,
        "aligned_length": aligned_length,
        "identity_pct": (
            matches / aligned_length * 100 if aligned_length else 0
        ),
        "visual": _format_visual(t_aln, mid, q_aln),
    }


# ──────────────────────── Report writers ──────────────────────────────

def _write_report_tsv(path, results):
    """Write per-region alignment statistics as a TSV with provenance header."""
    with open(path, "w", newline="") as fh:
        fh.write("## UCR Sequence Alignment Report\n")
        fh.write(
            f"## Generated: "
            f"{datetime.datetime.now(datetime.timezone.utc).isoformat()}\n"
        )
        fh.write("## hg38 genome: hg38.2bit (GRCh38)\n")
        fh.write("## T2T genome: hs1.2bit (T2T-CHM13v2.0)\n")
        fh.write("##\n")

        fields = [
            "ucr_id",
            "hg38_chrom", "hg38_start", "hg38_end", "hg38_len",
            "t2t_chrom", "t2t_start", "t2t_end", "t2t_len",
            "identical", "identity_pct",
            "matches", "mismatches", "gaps", "aligned_length",
        ]
        writer = csv.DictWriter(
            fh, fieldnames=fields, delimiter="\t", extrasaction="ignore",
        )
        writer.writeheader()
        for r in results:
            row = dict(r)
            row["identical"] = "YES" if row["identical"] else "NO"
            row["identity_pct"] = f"{row['identity_pct']:.2f}"
            writer.writerow(row)

    print(f"  ✓ {os.path.basename(path)}")


def _write_details(path, results):
    """Write visual alignments for every non-identical pair."""
    divergent = [r for r in results if not r["identical"]]
    with open(path, "w") as fh:
        fh.write("UCR Sequence Alignment Details\n")
        fh.write(
            f"Generated: "
            f"{datetime.datetime.now(datetime.timezone.utc).isoformat()}\n"
        )
        fh.write(f"{'=' * 70}\n\n")

        if not divergent:
            fh.write(
                "All paired UCR sequences are identical between hg38 "
                "and T2T-CHM13.\n"
            )
        else:
            fh.write(f"{len(divergent)} region(s) with differences:\n\n")
            for r in divergent:
                fh.write(f"─── UCR {r['ucr_id']} ───\n")
                fh.write(
                    f"  hg38: {r['hg38_chrom']}:{r['hg38_start']}"
                    f"-{r['hg38_end']} ({r['hg38_len']} bp)\n"
                )
                fh.write(
                    f"  T2T:  {r['t2t_chrom']}:{r['t2t_start']}"
                    f"-{r['t2t_end']} ({r['t2t_len']} bp)\n"
                )
                fh.write(
                    f"  Identity: {r['identity_pct']:.2f}%  "
                    f"Matches: {r['matches']}  "
                    f"Mismatches: {r['mismatches']}  "
                    f"Gaps: {r['gaps']}\n\n"
                )
                if r["visual"]:
                    fh.write(r["visual"])
                    fh.write("\n")
                fh.write("\n")

    print(f"  ✓ {os.path.basename(path)}")


# ──────────────────────────── Main ────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Extract UCR sequences from hg38 and T2T-CHM13v2.0 reference "
            "genomes and perform per-base pairwise alignment validation."
        ),
    )
    parser.add_argument(
        "--output-dir",
        default=os.environ.get("OUTPUT_DIR", os.getcwd()),
        help=(
            "Directory containing liftover BED files and where results are "
            "written (default: OUTPUT_DIR env var, or current directory)."
        ),
    )
    args = parser.parse_args()
    outdir = args.output_dir
    os.makedirs(outdir, exist_ok=True)

    # ── Derived paths ──
    hg38_bed = os.path.join(outdir, "ucr_hg38.bed")
    t2t_bed = os.path.join(outdir, "ucr_t2t_chm13.bed")
    hg38_2bit = os.path.join(outdir, "hg38.2bit")
    hs1_2bit = os.path.join(outdir, "hs1.2bit")
    hg38_fasta = os.path.join(outdir, "ucr_sequences_hg38.fa")
    t2t_fasta = os.path.join(outdir, "ucr_sequences_t2t.fa")
    report_tsv = os.path.join(outdir, "ucr_alignment_report.tsv")
    details_txt = os.path.join(outdir, "ucr_alignment_details.txt")

    # twoBitToFa: prefer container-baked binary, else download
    twobittofa = (
        _CONTAINER_TWOBITTOFA
        if os.path.exists(_CONTAINER_TWOBITTOFA)
        else os.path.join(outdir, "twoBitToFa")
    )

    # ── Validate prerequisites ──
    for bed, label in [(hg38_bed, "hg38 BED"), (t2t_bed, "T2T BED")]:
        if not os.path.exists(bed):
            print(f"Error: {label} not found at {bed}", file=sys.stderr)
            print(
                "Run convert_ucr_to_t2t.py first to produce the BED files.",
                file=sys.stderr,
            )
            sys.exit(1)

    # ── 1. Download reference genomes and twoBitToFa ──
    print("\n── Step 1: Download reference genomes ──")
    print("  (hg38.2bit ≈ 800 MB, hs1.2bit ≈ 780 MB – first run only)")
    download_file(HG38_2BIT_URL, hg38_2bit)
    download_file(HS1_2BIT_URL, hs1_2bit)
    download_file(TWOBITTOFA_URL, twobittofa, make_executable=True)

    # ── 2. Extract UCR sequences ──
    print("\n── Step 2: Extract UCR sequences ──")
    extract_sequences(twobittofa, hg38_2bit, hg38_bed, hg38_fasta)
    extract_sequences(twobittofa, hs1_2bit, t2t_bed, t2t_fasta)

    # ── 3. Load and pair sequences ──
    print("\n── Step 3: Per-base pairwise alignment ──")
    hg38_seqs = load_fasta(hg38_fasta)
    t2t_seqs = load_fasta(t2t_fasta)
    hg38_coords = load_bed(hg38_bed)
    t2t_coords = load_bed(t2t_bed)

    paired_ids = sorted(
        set(hg38_seqs) & set(t2t_seqs),
        key=lambda x: int(x) if x.isdigit() else x,
    )
    skipped = set(hg38_seqs) - set(t2t_seqs)
    print(f"  {len(paired_ids)} UCRs paired for comparison")
    if skipped:
        print(f"  {len(skipped)} unmapped UCRs skipped")

    # ── 4. Align each pair ──
    results = []
    non_identical = 0
    for ucr_id in paired_ids:
        stats = compare_sequences(hg38_seqs[ucr_id], t2t_seqs[ucr_id])
        hc = hg38_coords.get(ucr_id, ("", 0, 0))
        tc = t2t_coords.get(ucr_id, ("", 0, 0))
        row = {
            "ucr_id": ucr_id,
            "hg38_chrom": hc[0], "hg38_start": hc[1], "hg38_end": hc[2],
            "t2t_chrom": tc[0], "t2t_start": tc[1], "t2t_end": tc[2],
            "hg38_len": len(hg38_seqs[ucr_id]),
            "t2t_len": len(t2t_seqs[ucr_id]),
            **stats,
        }
        results.append(row)
        if not stats["identical"]:
            non_identical += 1

    # ── 5. Write reports ──
    print("\n── Step 4: Write reports ──")
    _write_report_tsv(report_tsv, results)
    _write_details(details_txt, results)

    # ── Summary ──
    total = len(results)
    identical = total - non_identical
    print(f"\n{'=' * 55}")
    print("  SEQUENCE VALIDATION SUMMARY")
    print(f"{'=' * 55}")
    print(f"  Paired UCRs compared:   {total}")
    print(f"  Identical sequences:    {identical}")
    print(f"  Non-identical:          {non_identical}")
    if total:
        avg_id = sum(r["identity_pct"] for r in results) / total
        print(f"  Average identity:       {avg_id:.2f} %")
    print(f"{'=' * 55}")
    print(f"\n  Reports:  {report_tsv}")
    print(f"            {details_txt}")
    print(f"  FASTA:    {hg38_fasta}")
    print(f"            {t2t_fasta}")


if __name__ == "__main__":
    main()
