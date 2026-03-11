import csv
import os
import re
import stat
import subprocess
import sys
import urllib.request

ASSEMBLY = "hs1"
DEFAULT_KMERS = (24, 36, 50, 100, 150, 250)

CHROM_SIZES_URL = f"https://hgdownload.soe.ucsc.edu/goldenPath/{ASSEMBLY}/bigZips/{ASSEMBLY}.chrom.sizes"
BIGBEDTOBED_URL = "https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed"
BIGBEDINFO_URL = "https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigBedInfo"

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.environ.get("OUTPUT_DIR", os.getcwd())

_container_bigbedtobed = os.path.join(SCRIPT_DIR, "bigBedToBed")
BIGBEDTOBED_BIN = (
    _container_bigbedtobed
    if os.path.exists(_container_bigbedtobed)
    else os.path.join(OUTPUT_DIR, "bigBedToBed")
)
_container_bigbedinfo = os.path.join(SCRIPT_DIR, "bigBedInfo")
BIGBEDINFO_BIN = (
    _container_bigbedinfo
    if os.path.exists(_container_bigbedinfo)
    else os.path.join(OUTPUT_DIR, "bigBedInfo")
)

CHROM_SIZES_FILE = os.path.join(OUTPUT_DIR, f"{ASSEMBLY}.chrom.sizes")
SUMMARY_TSV = os.path.join(OUTPUT_DIR, f"{ASSEMBLY}.unique_fraction_by_kmer.tsv")
CHROM_SUMMARY_TSV = os.path.join(OUTPUT_DIR, f"{ASSEMBLY}.unique_fraction_by_kmer_by_chrom.tsv")
COMPARISON_TSV = os.path.join(OUTPUT_DIR, f"{ASSEMBLY}.unique_fraction_comparison.tsv")
SUMMARY_TXT = os.path.join(OUTPUT_DIR, f"{ASSEMBLY}.unique_fraction_summary.txt")

# "Primary" here includes chrM by default; --exclude-chrM removes it.
PRIMARY_RE = re.compile(r"^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$")


def kmer_bb_url(k):
    return f"https://hgdownload.soe.ucsc.edu/gbdb/{ASSEMBLY}/hoffmanMappability/k{k}.Unique.Mappability.bb"


def bb_file_path(k):
    return os.path.join(OUTPUT_DIR, f"k{k}.Unique.Mappability.bb")


def bed_file_path(k):
    return os.path.join(OUTPUT_DIR, f"k{k}.unique.bed")


def download_file(url, filename, make_executable=False):
    if os.path.exists(filename):
        print(f"{filename} already exists, skipping download.")
        return
    print(f"Downloading {filename}...")
    urllib.request.urlretrieve(url, filename)
    if make_executable:
        st = os.stat(filename)
        os.chmod(filename, st.st_mode | stat.S_IEXEC)


def parse_kmers(argv):
    for i, arg in enumerate(argv):
        if arg == "--kmers" and i + 1 < len(argv):
            return tuple(int(x.strip()) for x in argv[i + 1].split(",") if x.strip())
        if arg.startswith("--kmers="):
            val = arg.split("=", 1)[1]
            return tuple(int(x.strip()) for x in val.split(",") if x.strip())
    return DEFAULT_KMERS


def parse_args(argv):
    primary_only = "--primary-only" in argv
    exclude_chrm = "--exclude-chrM" in argv
    kmers = parse_kmers(argv)
    return primary_only, exclude_chrm, tuple(sorted(set(kmers)))


def chrom_allowed(chrom, primary_only, exclude_chrm):
    if primary_only and not PRIMARY_RE.match(chrom):
        return False
    if exclude_chrm and chrom == "chrM":
        return False
    return True


def genome_size_bp(chrom_sizes_file, primary_only=False, exclude_chrm=False):
    return sum(genome_size_by_chrom(chrom_sizes_file, primary_only=primary_only, exclude_chrm=exclude_chrm).values())


def genome_size_by_chrom(chrom_sizes_file, primary_only=False, exclude_chrm=False):
    by_chrom = {}
    with open(chrom_sizes_file) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            chrom, size = line.split("\t")[:2]
            if chrom_allowed(chrom, primary_only, exclude_chrm):
                by_chrom[chrom] = int(size)
    return by_chrom


def unique_covered_bp(bed_file, primary_only=False, exclude_chrm=False):
    covered_by_chrom, strict_filter = unique_covered_bp_by_chrom(
        bed_file,
        primary_only=primary_only,
        exclude_chrm=exclude_chrm,
    )
    return sum(covered_by_chrom.values()), strict_filter


def unique_covered_bp_by_chrom(bed_file, primary_only=False, exclude_chrm=False):
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
            if not chrom_allowed(chrom, primary_only, exclude_chrm):
                continue
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
                    score = None
            elif len(fields) >= 4:
                try:
                    score = float(fields[3])
                    numeric_scores.append(score)
                except ValueError:
                    score = None
            records.append((chrom, start, end, score))

    if not records:
        return {}, False

    filter_non_maximal_scores = False
    unique_score_threshold = None
    if numeric_scores:
        min_score = min(numeric_scores)
        max_score = max(numeric_scores)
        filter_non_maximal_scores = max_score > min_score
        unique_score_threshold = max_score

    by_chrom = {}
    for chrom, start, end, score in records:
        if filter_non_maximal_scores:
            if score is None or score < unique_score_threshold:
                continue
        by_chrom.setdefault(chrom, []).append((start, end))

    covered_by_chrom = {}
    for chrom in by_chrom:
        intervals = sorted(by_chrom[chrom])
        if not intervals:
            continue
        cur_start, cur_end = intervals[0]
        covered = 0
        for start, end in intervals[1:]:
            if start <= cur_end:
                if end > cur_end:
                    cur_end = end
            else:
                covered += cur_end - cur_start
                cur_start, cur_end = start, end
        covered += cur_end - cur_start
        covered_by_chrom[chrom] = covered
    return covered_by_chrom, filter_non_maximal_scores


def bigbed_to_bed(bb_file, bed_file):
    result = subprocess.run(
        [BIGBEDTOBED_BIN, bb_file, bed_file],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print(f"bigBedToBed stderr:\n{result.stderr}", file=sys.stderr)
        result.check_returncode()


def inspect_bigbed(bb_file):
    result = subprocess.run(
        [BIGBEDINFO_BIN, bb_file],
        capture_output=True,
        text=True,
    )
    if result.returncode == 0 and result.stdout:
        print("bigBedInfo output:")
        print(result.stdout.strip())
    elif result.stderr:
        print(f"bigBedInfo warning:\n{result.stderr}", file=sys.stderr)


def build_comparison_rows(rows):
    comparisons = []
    previous = None
    for row in rows:
        if previous is None:
            comparisons.append({
                "kmer": row["kmer"],
                "fraction_unique": row["fraction_unique"],
                "delta_fraction_from_previous_k": "",
                "delta_percent_from_previous_k": "",
            })
        else:
            delta_fraction = row["fraction_unique"] - previous["fraction_unique"]
            comparisons.append({
                "kmer": row["kmer"],
                "fraction_unique": row["fraction_unique"],
                "delta_fraction_from_previous_k": delta_fraction,
                "delta_percent_from_previous_k": 100.0 * delta_fraction,
            })
        previous = row
    return comparisons


def write_summary(rows, chrom_rows=None):
    with open(SUMMARY_TSV, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "assembly",
            "kmer",
            "unique_bp",
            "genome_bp",
            "fraction_unique",
            "percent_unique",
            "strict_unique_filter_applied",
        ])
        for row in rows:
            writer.writerow([
                row["assembly"],
                row["kmer"],
                row["unique_bp"],
                row["genome_bp"],
                row["fraction_unique"],
                row["percent_unique"],
                row["strict_unique_filter_applied"],
            ])

    comparisons = build_comparison_rows(rows)
    with open(COMPARISON_TSV, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "kmer",
            "fraction_unique",
            "delta_fraction_from_previous_k",
            "delta_percent_from_previous_k",
        ])
        for row in comparisons:
            writer.writerow([
                row["kmer"],
                row["fraction_unique"],
                row["delta_fraction_from_previous_k"],
                row["delta_percent_from_previous_k"],
            ])

    if chrom_rows:
        with open(CHROM_SUMMARY_TSV, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow([
                "assembly",
                "kmer",
                "chromosome",
                "unique_bp",
                "genome_bp",
                "fraction_unique",
                "percent_unique",
                "strict_unique_filter_applied",
            ])
            for row in chrom_rows:
                writer.writerow([
                    row["assembly"],
                    row["kmer"],
                    row["chromosome"],
                    row["unique_bp"],
                    row["genome_bp"],
                    row["fraction_unique"],
                    row["percent_unique"],
                    row["strict_unique_filter_applied"],
                ])

    with open(SUMMARY_TXT, "w") as fh:
        fh.write(f"Assembly: {ASSEMBLY}\n")
        fh.write("Metric: strict unique mappability fraction by k-mer\n")
        fh.write("Numerator: union bp at maximal score per track when scores vary; otherwise binary interval union.\n")
        fh.write("Denominator: total bp in selected chrom.sizes set.\n\n")
        fh.write("k-mer summary:\n")
        for row in rows:
            fh.write(
                f"  k={row['kmer']}: {row['percent_unique']:.6f}% "
                f"(unique_bp={row['unique_bp']}, genome_bp={row['genome_bp']})\n"
            )
        fh.write("\nPairwise change vs previous k:\n")
        for row in comparisons[1:]:
            fh.write(
                f"  k={row['kmer']}: delta={row['delta_percent_from_previous_k']:.6f} percentage points\n"
            )


def main(argv=None):
    argv = argv or sys.argv[1:]
    primary_only, exclude_chrm, kmers = parse_args(argv)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    download_file(BIGBEDTOBED_URL, BIGBEDTOBED_BIN, make_executable=True)
    download_file(BIGBEDINFO_URL, BIGBEDINFO_BIN, make_executable=True)
    download_file(CHROM_SIZES_URL, CHROM_SIZES_FILE)

    genome_bp = genome_size_bp(
        CHROM_SIZES_FILE,
        primary_only=primary_only,
        exclude_chrm=exclude_chrm,
    )
    if genome_bp == 0:
        raise ValueError("Genome size is zero after chromosome filters.")

    rows = []
    chrom_rows = []
    genome_bp_by_chrom = genome_size_by_chrom(
        CHROM_SIZES_FILE,
        primary_only=primary_only,
        exclude_chrm=exclude_chrm,
    )
    for k in kmers:
        bb_file = bb_file_path(k)
        bed_file = bed_file_path(k)
        download_file(kmer_bb_url(k), bb_file)
        print(f"\n=== k={k} ===")
        inspect_bigbed(bb_file)
        print("Converting bigBed to BED...")
        bigbed_to_bed(bb_file, bed_file)

        unique_bp_by_chrom, strict_filter = unique_covered_bp_by_chrom(
            bed_file,
            primary_only=primary_only,
            exclude_chrm=exclude_chrm,
        )
        unique_bp = sum(unique_bp_by_chrom.values())
        fraction = unique_bp / genome_bp
        percent = 100.0 * fraction
        row = {
            "assembly": ASSEMBLY,
            "kmer": k,
            "unique_bp": unique_bp,
            "genome_bp": genome_bp,
            "fraction_unique": fraction,
            "percent_unique": percent,
            "strict_unique_filter_applied": strict_filter,
        }
        rows.append(row)
        for chrom, chrom_genome_bp in genome_bp_by_chrom.items():
            chrom_unique_bp = unique_bp_by_chrom.get(chrom, 0)
            chrom_fraction = chrom_unique_bp / chrom_genome_bp if chrom_genome_bp else 0.0
            chrom_rows.append({
                "assembly": ASSEMBLY,
                "kmer": k,
                "chromosome": chrom,
                "unique_bp": chrom_unique_bp,
                "genome_bp": chrom_genome_bp,
                "fraction_unique": chrom_fraction,
                "percent_unique": 100.0 * chrom_fraction,
                "strict_unique_filter_applied": strict_filter,
            })
        print(
            f"k={k}\tunique_bp={unique_bp}\tgenome_bp={genome_bp}"
            f"\tfraction_unique={fraction}\tpercent_unique={percent}"
        )

    write_summary(rows, chrom_rows=chrom_rows)
    print(f"Saved k-mer summary to {SUMMARY_TSV}")
    if chrom_rows:
        print(f"Saved per-chromosome k-mer summary to {CHROM_SUMMARY_TSV}")
    print(f"Saved comparison table to {COMPARISON_TSV}")
    print(f"Saved text summary to {SUMMARY_TXT}")


if __name__ == "__main__":
    main()
