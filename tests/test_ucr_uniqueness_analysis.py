"""
Unit tests for ucr_uniqueness_analysis.py.

Covers:
  - merge_intervals() with various interval configurations
  - compute_overlap_bp() for overlap computation with merged intervals
  - get_non_unique_intervals() for gap detection in unique coverage
  - extract_non_unique_seq() for sequence extraction at non-unique positions
  - load_unique_intervals() with score filtering
  - load_validation_results() with synthetic alignment report
  - parse_bed4() with synthetic BED files
  - write_per_region_report() and write_summary() output format
"""

import csv
import json
import os
import sys

import pytest

REPO_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, REPO_ROOT)

import ucr_uniqueness_analysis as ua  # noqa: E402


# ===========================================================================
# 1. merge_intervals()
# ===========================================================================


class TestMergeIntervals:
    """Unit tests for interval merging."""

    def test_empty(self):
        assert ua.merge_intervals([]) == []

    def test_single_interval(self):
        assert ua.merge_intervals([(10, 20)]) == [(10, 20)]

    def test_non_overlapping(self):
        result = ua.merge_intervals([(10, 20), (30, 40)])
        assert result == [(10, 20), (30, 40)]

    def test_overlapping(self):
        result = ua.merge_intervals([(10, 25), (20, 35)])
        assert result == [(10, 35)]

    def test_adjacent(self):
        result = ua.merge_intervals([(10, 20), (20, 30)])
        assert result == [(10, 30)]

    def test_contained(self):
        result = ua.merge_intervals([(10, 50), (20, 30)])
        assert result == [(10, 50)]

    def test_multiple_merges(self):
        result = ua.merge_intervals([(1, 5), (3, 8), (10, 15), (12, 20)])
        assert result == [(1, 8), (10, 20)]


# ===========================================================================
# 2. compute_overlap_bp()
# ===========================================================================


class TestComputeOverlapBp:
    """Unit tests for bp overlap computation."""

    def test_no_intervals(self):
        assert ua.compute_overlap_bp(100, 200, []) == 0

    def test_no_overlap(self):
        merged = [(10, 20), (300, 400)]
        assert ua.compute_overlap_bp(100, 200, merged) == 0

    def test_full_overlap(self):
        merged = [(50, 300)]
        assert ua.compute_overlap_bp(100, 200, merged) == 100

    def test_partial_overlap_left(self):
        merged = [(80, 150)]
        assert ua.compute_overlap_bp(100, 200, merged) == 50

    def test_partial_overlap_right(self):
        merged = [(150, 250)]
        assert ua.compute_overlap_bp(100, 200, merged) == 50

    def test_internal_overlap(self):
        merged = [(120, 160)]
        assert ua.compute_overlap_bp(100, 200, merged) == 40

    def test_multiple_overlaps(self):
        merged = [(110, 130), (150, 170)]
        assert ua.compute_overlap_bp(100, 200, merged) == 40

    def test_exact_boundary(self):
        merged = [(100, 200)]
        assert ua.compute_overlap_bp(100, 200, merged) == 100

    def test_adjacent_no_overlap(self):
        """Interval ending at ucr_start has no overlap."""
        merged = [(90, 100)]
        assert ua.compute_overlap_bp(100, 200, merged) == 0


# ===========================================================================
# 3. get_non_unique_intervals()
# ===========================================================================


class TestGetNonUniqueIntervals:
    """Unit tests for non-unique interval detection."""

    def test_no_unique_intervals(self):
        result = ua.get_non_unique_intervals(100, 200, [])
        assert result == [(100, 200)]

    def test_fully_unique(self):
        merged = [(50, 300)]
        result = ua.get_non_unique_intervals(100, 200, merged)
        assert result == []

    def test_gap_in_middle(self):
        merged = [(100, 130), (170, 200)]
        result = ua.get_non_unique_intervals(100, 200, merged)
        assert result == [(130, 170)]

    def test_gap_at_start(self):
        merged = [(150, 200)]
        result = ua.get_non_unique_intervals(100, 200, merged)
        assert result == [(100, 150)]

    def test_gap_at_end(self):
        merged = [(100, 150)]
        result = ua.get_non_unique_intervals(100, 200, merged)
        assert result == [(150, 200)]

    def test_multiple_gaps(self):
        merged = [(110, 120), (140, 150), (180, 200)]
        result = ua.get_non_unique_intervals(100, 200, merged)
        assert result == [(100, 110), (120, 140), (150, 180)]


# ===========================================================================
# 4. extract_non_unique_seq()
# ===========================================================================


class TestExtractNonUniqueSeq:
    """Unit tests for non-unique sequence extraction."""

    def test_empty_when_all_unique(self):
        result = ua.extract_non_unique_seq(100, "ACGTACGT", [])
        assert result == ""

    def test_full_sequence_when_none_unique(self):
        result = ua.extract_non_unique_seq(100, "ACGTACGT", [(100, 108)])
        assert result == "ACGTACGT"

    def test_partial_non_unique(self):
        # UCR at position 100-108, non-unique at 102-105
        result = ua.extract_non_unique_seq(100, "ACGTACGT", [(102, 105)])
        assert result == "GTA"

    def test_multiple_non_unique_regions(self):
        # UCR at 100-110, non-unique at 100-102 and 106-108
        seq = "ACGTACGTAC"
        result = ua.extract_non_unique_seq(100, seq, [(100, 102), (106, 108)])
        assert result == "AC" + "GT"


# ===========================================================================
# 5. load_unique_intervals() with score filtering
# ===========================================================================


class TestLoadUniqueIntervals:
    """Unit tests for unique interval loading with score filtering."""

    def test_binary_intervals_no_filter(self, tmp_path):
        bed = tmp_path / "k24.unique.bed"
        bed.write_text(
            "chr1\t0\t10\n"
            "chr1\t5\t20\n"
            "chr2\t100\t110\n"
        )
        by_chrom, strict = ua.load_unique_intervals(str(bed))
        assert strict is False
        # chr1 intervals should be merged: (0, 20)
        assert by_chrom["chr1"] == [(0, 20)]
        assert by_chrom["chr2"] == [(100, 110)]

    def test_score_filter_keeps_max_only(self, tmp_path):
        bed = tmp_path / "k24.unique.bed"
        bed.write_text(
            "chr1\t0\t10\tnameA\t500\n"
            "chr1\t10\t20\tnameB\t1000\n"
            "chr1\t30\t40\tnameC\t1000\n"
        )
        by_chrom, strict = ua.load_unique_intervals(str(bed))
        assert strict is True
        # Only score-1000 intervals kept
        assert by_chrom["chr1"] == [(10, 20), (30, 40)]

    def test_empty_file(self, tmp_path):
        bed = tmp_path / "empty.bed"
        bed.write_text("")
        by_chrom, strict = ua.load_unique_intervals(str(bed))
        assert by_chrom == {}
        assert strict is False

    def test_comment_lines_skipped(self, tmp_path):
        bed = tmp_path / "k24.unique.bed"
        bed.write_text(
            "# header\n"
            "chr1\t10\t20\n"
        )
        by_chrom, strict = ua.load_unique_intervals(str(bed))
        assert by_chrom["chr1"] == [(10, 20)]


# ===========================================================================
# 6. load_validation_results()
# ===========================================================================


class TestLoadValidationResults:
    """Unit tests for parsing the alignment report TSV."""

    def test_basic_parsing(self, tmp_path):
        tsv = tmp_path / "ucr_alignment_report.tsv"
        tsv.write_text(
            "## UCR Sequence Alignment Report\n"
            "## Generated: 2025-01-01T00:00:00\n"
            "##\n"
            "ucr_id\tidentical\tidentity_pct\n"
            "42\tYES\t100.00\n"
            "99\tNO\t99.50\n"
        )
        results = ua.load_validation_results(str(tsv))
        assert len(results) == 2
        assert results["42"]["identical"] == "YES"
        assert results["99"]["identity_pct"] == "99.50"

    def test_missing_file_returns_empty(self, tmp_path):
        results = ua.load_validation_results(
            str(tmp_path / "nonexistent.tsv")
        )
        assert results == {}


# ===========================================================================
# 7. parse_bed4()
# ===========================================================================


class TestParseBed4:
    """Unit tests for BED4 file parsing."""

    def test_basic_parsing(self, tmp_path):
        bed = tmp_path / "test.bed"
        bed.write_text("chr1\t100\t200\tucr1\nchr2\t300\t400\tucr2\n")
        regions = ua.parse_bed4(str(bed))
        assert len(regions) == 2
        assert regions[0] == {
            "chrom": "chr1", "start": 100, "end": 200, "name": "ucr1",
        }

    def test_comment_lines_skipped(self, tmp_path):
        bed = tmp_path / "test.bed"
        bed.write_text("# comment\nchr1\t100\t200\tucr1\n")
        assert len(ua.parse_bed4(str(bed))) == 1

    def test_missing_file_returns_empty(self, tmp_path):
        assert ua.parse_bed4(str(tmp_path / "missing.bed")) == []


# ===========================================================================
# 8. write_per_region_report()
# ===========================================================================


class TestWritePerRegionReport:
    """Integration test for the per-region TSV writer."""

    def test_output_format(self, tmp_path):
        path = tmp_path / "report.tsv"
        rows = [
            {
                "ucr_id": "42",
                "t2t_chrom": "chr1", "t2t_start": 100, "t2t_end": 300,
                "ucr_length": 200,
                "identity_pct": "100.00",
                "kmer": 24,
                "unique_bp": 150,
                "unique_fraction": 0.75,
                "non_unique_bp": 50,
                "non_unique_fraction": 0.25,
                "ucr_sequence": "ACGT" * 50,
                "non_unique_sequence": "ACGT" * 12 + "AC",
            },
        ]
        ua.write_per_region_report(str(path), rows)

        text = path.read_text()
        assert "## UCR Uniqueness Report" in text

        data_rows = list(csv.DictReader(
            (l for l in text.splitlines() if not l.startswith("##")),
            delimiter="\t",
        ))
        assert len(data_rows) == 1
        assert data_rows[0]["ucr_id"] == "42"
        assert data_rows[0]["unique_fraction"] == "0.750000"
        assert data_rows[0]["non_unique_fraction"] == "0.250000"


# ===========================================================================
# 9. write_summary()
# ===========================================================================


class TestWriteSummary:
    """Integration test for the JSON summary writer."""

    def test_output_structure(self, tmp_path):
        path = tmp_path / "summary.json"
        rows = [
            {
                "ucr_id": "42", "kmer": 24,
                "ucr_length": 200, "unique_bp": 150,
            },
            {
                "ucr_id": "99", "kmer": 24,
                "ucr_length": 300, "unique_bp": 300,
            },
        ]
        validation = {
            "42": {"identical": "YES", "identity_pct": "100.00"},
            "99": {"identical": "NO", "identity_pct": "99.50"},
        }
        ua.write_summary(
            str(path), rows, validation, kmers=(24,),
            total_ucrs_lifted=2,
        )

        payload = json.loads(path.read_text())
        assert payload["total_ucrs_lifted_over"] == 2
        assert "24" in payload["per_kmer_summary"]
        k24 = payload["per_kmer_summary"]["24"]
        assert k24["total_ucr_bp"] == 500
        assert k24["total_unique_bp"] == 450
        assert k24["fully_unique_ucrs"] == 1
        assert k24["partially_unique_ucrs"] == 1
        assert payload["validation"]["paired_ucrs"] == 2
        assert payload["validation"]["identical_sequences"] == 1
