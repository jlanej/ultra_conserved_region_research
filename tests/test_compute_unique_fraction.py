import os
import sys
import pytest

REPO_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, REPO_ROOT)

import compute_unique_fraction as cuf  # noqa: E402


def test_chrom_allowed_primary_filter():
    assert cuf.chrom_allowed("chr1", primary_only=True, exclude_chrm=False)
    assert cuf.chrom_allowed("chrX", primary_only=True, exclude_chrm=False)
    assert cuf.chrom_allowed("chrM", primary_only=True, exclude_chrm=False)
    assert not cuf.chrom_allowed("chrUn_KI270521v1", primary_only=True, exclude_chrm=False)


def test_chrom_allowed_exclude_chrm():
    assert not cuf.chrom_allowed("chrM", primary_only=False, exclude_chrm=True)
    assert cuf.chrom_allowed("chr1", primary_only=False, exclude_chrm=True)


def test_genome_size_bp_filters(tmp_path):
    chrom_sizes = tmp_path / "hg38.chrom.sizes"
    chrom_sizes.write_text(
        "chr1\t100\n"
        "chr2\t200\n"
        "chrM\t10\n"
        "chrUn_KI270521v1\t50\n"
    )

    assert cuf.genome_size_bp(str(chrom_sizes)) == 360
    assert cuf.genome_size_bp(str(chrom_sizes), primary_only=True) == 310
    assert cuf.genome_size_bp(str(chrom_sizes), primary_only=True, exclude_chrm=True) == 300


def test_unique_covered_bp_binary_intervals(tmp_path):
    bed = tmp_path / "k24.unique.bed"
    bed.write_text(
        "chr1\t0\t10\n"
        "chr1\t5\t20\n"
        "chr1\t20\t30\n"
    )

    # Union coverage for chr1 is [0,30): 30 bp
    covered, strict_filter = cuf.unique_covered_bp(str(bed))
    assert covered == 30
    assert strict_filter is False


def test_unique_covered_bp_filters_non_maximal_scores(tmp_path):
    bed = tmp_path / "k24.unique.bed"
    bed.write_text(
        # score in BED5; interpreted as 0..1000 UCSC score and normalized
        "chr1\t0\t10\tnameA\t500\n"
        "chr1\t5\t15\tnameB\t1000\n"
    )

    # Mixed score intervals should keep only maximum score as "truly unique".
    # Here score 1000 is kept and score 500 is excluded -> [5,15): 10 bp.
    covered, strict_filter = cuf.unique_covered_bp(str(bed))
    assert covered == 10
    assert strict_filter is True


def test_parse_kmers_defaults_and_custom():
    assert cuf.parse_kmers([]) == cuf.DEFAULT_KMERS
    assert cuf.parse_kmers(["--kmers", "24,50,150"]) == (24, 50, 150)
    assert cuf.parse_kmers(["--kmers=36,100"]) == (36, 100)


def test_build_comparison_rows():
    rows = [
        {"kmer": 24, "fraction_unique": 0.5},
        {"kmer": 36, "fraction_unique": 0.6},
        {"kmer": 50, "fraction_unique": 0.61},
    ]
    out = cuf.build_comparison_rows(rows)
    assert out[0]["delta_fraction_from_previous_k"] == ""
    assert out[1]["delta_percent_from_previous_k"] == pytest.approx(10.0)
    assert out[1]["delta_fraction_from_previous_k"] == pytest.approx(0.1)
    assert out[2]["delta_percent_from_previous_k"] == pytest.approx(1.0)
