import csv
import os
import re
import stat
import subprocess
import sys
import urllib.request

K24_UNIQUE_BB_URL = (
    "https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/"
    "k24.Unique.Mappability.bb"
)
HG38_CHROM_SIZES_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"
)
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

K24_UNIQUE_BB_FILE = os.path.join(OUTPUT_DIR, "k24.Unique.Mappability.bb")
HG38_CHROM_SIZES_FILE = os.path.join(OUTPUT_DIR, "hg38.chrom.sizes")
K24_UNIQUE_BED_FILE = os.path.join(OUTPUT_DIR, "k24.unique.bed")
SUMMARY_TSV = os.path.join(OUTPUT_DIR, "k24.unique_fraction.tsv")

PRIMARY_RE = re.compile(r"^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$")


def download_file(url, filename, make_executable=False):
    if os.path.exists(filename):
        print(f"{filename} already exists, skipping download.")
        return
    print(f"Downloading {filename}...")
    urllib.request.urlretrieve(url, filename)
    if make_executable:
        st = os.stat(filename)
        os.chmod(filename, st.st_mode | stat.S_IEXEC)


def parse_args(argv):
    primary_only = "--primary-only" in argv
    exclude_chrm = "--exclude-chrM" in argv
    return primary_only, exclude_chrm


def chrom_allowed(chrom, primary_only, exclude_chrm):
    if primary_only and not PRIMARY_RE.match(chrom):
        return False
    if exclude_chrm and chrom == "chrM":
        return False
    return True


def genome_size_bp(chrom_sizes_file, primary_only=False, exclude_chrm=False):
    total = 0
    with open(chrom_sizes_file) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            chrom, size = line.split("\t")[:2]
            if chrom_allowed(chrom, primary_only, exclude_chrm):
                total += int(size)
    return total


def unique_covered_bp(bed_file, primary_only=False, exclude_chrm=False):
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
        return 0

    keep_only_max_scored = False
    unique_score_threshold = None
    if numeric_scores:
        min_score = min(numeric_scores)
        max_score = max(numeric_scores)
        keep_only_max_scored = max_score > min_score
        unique_score_threshold = max_score

    by_chrom = {}
    for chrom, start, end, score in records:
        if keep_only_max_scored:
            if score is None or score < unique_score_threshold:
                continue
        by_chrom.setdefault(chrom, []).append((start, end))

    covered = 0
    for chrom in by_chrom:
        intervals = sorted(by_chrom[chrom])
        if not intervals:
            continue
        cur_start, cur_end = intervals[0]
        for start, end in intervals[1:]:
            if start <= cur_end:
                if end > cur_end:
                    cur_end = end
            else:
                covered += cur_end - cur_start
                cur_start, cur_end = start, end
        covered += cur_end - cur_start
    return covered


def bigbed_to_bed():
    result = subprocess.run(
        [BIGBEDTOBED_BIN, K24_UNIQUE_BB_FILE, K24_UNIQUE_BED_FILE],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print(f"bigBedToBed stderr:\n{result.stderr}", file=sys.stderr)
        result.check_returncode()


def inspect_bigbed():
    """Print bigBed metadata so users can verify whether values are binary or scored."""
    result = subprocess.run(
        [BIGBEDINFO_BIN, K24_UNIQUE_BB_FILE],
        capture_output=True,
        text=True,
    )
    if result.returncode == 0 and result.stdout:
        print("bigBedInfo output:")
        print(result.stdout.strip())
    elif result.stderr:
        print(f"bigBedInfo warning:\n{result.stderr}", file=sys.stderr)


def write_summary(unique_bp, genome_bp, fraction, percent):
    with open(SUMMARY_TSV, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["unique_bp", "genome_bp", "fraction_unique_24", "percent_unique_24"])
        writer.writerow([unique_bp, genome_bp, fraction, percent])


def main(argv=None):
    argv = argv or sys.argv[1:]
    primary_only, exclude_chrm = parse_args(argv)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    download_file(BIGBEDTOBED_URL, BIGBEDTOBED_BIN, make_executable=True)
    download_file(BIGBEDINFO_URL, BIGBEDINFO_BIN, make_executable=True)
    download_file(K24_UNIQUE_BB_URL, K24_UNIQUE_BB_FILE)
    download_file(HG38_CHROM_SIZES_URL, HG38_CHROM_SIZES_FILE)

    inspect_bigbed()
    print("Converting bigBed to BED...")
    bigbed_to_bed()

    genome_bp = genome_size_bp(
        HG38_CHROM_SIZES_FILE,
        primary_only=primary_only,
        exclude_chrm=exclude_chrm,
    )
    unique_bp = unique_covered_bp(
        K24_UNIQUE_BED_FILE,
        primary_only=primary_only,
        exclude_chrm=exclude_chrm,
    )

    if genome_bp == 0:
        raise ValueError("Genome size is zero after chromosome filters.")

    fraction = unique_bp / genome_bp
    percent = 100.0 * fraction

    print(f"unique_bp\t{unique_bp}")
    print(f"genome_bp\t{genome_bp}")
    print(f"fraction_unique_24\t{fraction}")
    print(f"percent_unique_24\t{percent}")

    write_summary(unique_bp, genome_bp, fraction, percent)
    print(f"Saved summary to {SUMMARY_TSV}")


if __name__ == "__main__":
    main()
