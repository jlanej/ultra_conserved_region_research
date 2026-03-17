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


# ===========================================================================
# 10. Litmus control
# ===========================================================================


class TestLitmusControl:
    """Tests for the built-in litmus control region and validation logic."""

    def test_litmus_control_constants_defined(self):
        """LITMUS_CONTROL_REGION has valid BED4-compatible fields."""
        r = ua.LITMUS_CONTROL_REGION
        assert isinstance(r["chrom"], str) and r["chrom"]
        assert isinstance(r["start"], int)
        assert isinstance(r["end"], int)
        assert r["end"] > r["start"], "end must be > start"
        assert r["name"] == ua.LITMUS_CONTROL_ID
        assert 0.0 < ua.LITMUS_MAX_UNIQUE_FRACTION < 1.0

    def test_litmus_detected_as_non_unique_with_no_intervals(self):
        """If no unique intervals exist the control has unique_fraction 0."""
        region = ua.LITMUS_CONTROL_REGION
        unique_bp = ua.compute_overlap_bp(region["start"], region["end"], [])
        ucr_len = region["end"] - region["start"]
        unique_fraction = unique_bp / ucr_len
        assert unique_fraction == 0.0
        assert unique_fraction <= ua.LITMUS_MAX_UNIQUE_FRACTION

    def test_litmus_passes_when_non_unique(self):
        """Litmus logic correctly reports PASSED when unique_fraction is 0."""
        rows = [
            {
                "ucr_id": ua.LITMUS_CONTROL_ID, "kmer": 24,
                "ucr_length": 1000, "unique_bp": 0,
                "unique_fraction": 0.0,
            },
        ]
        passed = all(
            r["unique_fraction"] <= ua.LITMUS_MAX_UNIQUE_FRACTION
            for r in rows
        )
        assert passed

    def test_litmus_fails_when_all_bases_appear_unique(self):
        """Litmus logic correctly reports FAILED when unique_fraction is 1."""
        rows = [
            {
                "ucr_id": ua.LITMUS_CONTROL_ID, "kmer": 24,
                "ucr_length": 1000, "unique_bp": 1000,
                "unique_fraction": 1.0,
            },
        ]
        passed = all(
            r["unique_fraction"] <= ua.LITMUS_MAX_UNIQUE_FRACTION
            for r in rows
        )
        assert not passed

    def test_litmus_excluded_from_ucr_summary_stats(self, tmp_path):
        """LITMUS_CONTROL rows do not inflate aggregate UCR statistics."""
        path = tmp_path / "summary.json"
        rows = [
            {
                "ucr_id": "42", "kmer": 24,
                "ucr_length": 200, "unique_bp": 200,
            },
            {
                "ucr_id": ua.LITMUS_CONTROL_ID, "kmer": 24,
                "ucr_length": 1000, "unique_bp": 0,
                "unique_fraction": 0.0,
            },
        ]
        ua.write_summary(str(path), rows, {}, kmers=(24,), total_ucrs_lifted=1)
        payload = json.loads(path.read_text())

        k24 = payload["per_kmer_summary"]["24"]
        # Only the real UCR (200 bp) should be counted.
        assert k24["total_ucr_bp"] == 200
        assert k24["fully_unique_ucrs"] == 1
        assert k24["not_unique_ucrs"] == 0

    def test_litmus_section_present_in_summary_when_control_included(
        self, tmp_path
    ):
        """write_summary() adds a litmus_test section when control rows exist."""
        path = tmp_path / "summary.json"
        rows = [
            {
                "ucr_id": "42", "kmer": 24,
                "ucr_length": 200, "unique_bp": 200,
            },
            {
                "ucr_id": ua.LITMUS_CONTROL_ID, "kmer": 24,
                "ucr_length": 1000, "unique_bp": 0,
                "unique_fraction": 0.0,
            },
        ]
        ua.write_summary(str(path), rows, {}, kmers=(24,), total_ucrs_lifted=1)
        payload = json.loads(path.read_text())

        assert "litmus_test" in payload
        lt = payload["litmus_test"]
        assert lt["overall_passed"] is True
        assert "24" in lt["per_kmer"]
        assert lt["per_kmer"]["24"]["unique_fraction"] == 0.0
        assert lt["per_kmer"]["24"]["passed"] is True
        assert lt["expected_unique_fraction_max"] == ua.LITMUS_MAX_UNIQUE_FRACTION

    def test_litmus_section_absent_when_no_control_rows(self, tmp_path):
        """write_summary() omits litmus_test when no control rows are present."""
        path = tmp_path / "summary.json"
        rows = [{"ucr_id": "42", "kmer": 24, "ucr_length": 200, "unique_bp": 200}]
        ua.write_summary(str(path), rows, {}, kmers=(24,), total_ucrs_lifted=1)
        payload = json.loads(path.read_text())
        assert "litmus_test" not in payload

    def test_litmus_fails_reflected_in_summary(self, tmp_path):
        """overall_passed is False when any kmer row exceeds the threshold."""
        path = tmp_path / "summary.json"
        rows = [
            {
                "ucr_id": ua.LITMUS_CONTROL_ID, "kmer": 24,
                "ucr_length": 1000, "unique_bp": 1000,
                "unique_fraction": 1.0,
            },
        ]
        ua.write_summary(str(path), rows, {}, kmers=(24,), total_ucrs_lifted=0)
        payload = json.loads(path.read_text())
        assert payload["litmus_test"]["overall_passed"] is False
        assert payload["litmus_test"]["per_kmer"]["24"]["passed"] is False


# ===========================================================================
# 11. NON_UNIQUE_FRACTION_THRESHOLD constant
# ===========================================================================


class TestNonUniqueFractionThreshold:
    """Sanity checks for the 10 % non-unique threshold constant."""

    def test_threshold_value(self):
        assert ua.NON_UNIQUE_FRACTION_THRESHOLD == 0.10

    def test_threshold_in_range(self):
        assert 0.0 < ua.NON_UNIQUE_FRACTION_THRESHOLD < 1.0


# ===========================================================================
# 12. write_non_unique_query_fasta()
# ===========================================================================


class TestWriteNonUniqueQueryFasta:
    """Unit tests for FASTA writer used by the locus-mapping step."""

    def _make_row(self, ucr_id, kmer, nu_fraction, nu_intervals,
                  chrom="chr1", ucr_start=1000, ucr_seq="A" * 300):
        return {
            "ucr_id": ucr_id,
            "kmer": kmer,
            "t2t_chrom": chrom,
            "t2t_start": ucr_start,
            "t2t_end": ucr_start + len(ucr_seq),
            "non_unique_fraction": nu_fraction,
            "non_unique_intervals": nu_intervals,
        }

    def test_writes_above_threshold(self, tmp_path):
        fa = tmp_path / "q.fa"
        seq = "ACGT" * 75  # 300 bp
        rows = [
            self._make_row("42", 24, 0.50, [(1050, 1100)], ucr_seq=seq),
        ]
        ucr_seqs = {"42": seq}
        count = ua.write_non_unique_query_fasta(rows, ucr_seqs, str(fa))
        assert count == 1
        text = fa.read_text()
        assert "NUQUERY|42|k24|chr1|1050|1100" in text

    def test_skips_below_threshold(self, tmp_path):
        fa = tmp_path / "q.fa"
        seq = "ACGT" * 75
        rows = [
            self._make_row("42", 24, 0.05, [(1050, 1100)], ucr_seq=seq),
        ]
        ucr_seqs = {"42": seq}
        count = ua.write_non_unique_query_fasta(rows, ucr_seqs, str(fa))
        assert count == 0
        assert fa.read_text() == ""

    def test_skips_very_short_intervals(self, tmp_path):
        fa = tmp_path / "q.fa"
        seq = "ACGT" * 75
        # non_unique interval only 5 bp long – below _MIN_NU_QUERY_LEN
        rows = [
            self._make_row("42", 24, 0.50, [(1050, 1055)], ucr_seq=seq),
        ]
        ucr_seqs = {"42": seq}
        count = ua.write_non_unique_query_fasta(rows, ucr_seqs, str(fa))
        assert count == 0

    def test_skips_litmus_control(self, tmp_path):
        fa = tmp_path / "q.fa"
        seq = "ACGT" * 75
        rows = [
            self._make_row(ua.LITMUS_CONTROL_ID, 24, 0.90, [(1050, 1200)],
                           ucr_seq=seq),
        ]
        ucr_seqs = {ua.LITMUS_CONTROL_ID: seq}
        count = ua.write_non_unique_query_fasta(rows, ucr_seqs, str(fa))
        assert count == 0

    def test_multiple_intervals_and_ucrs(self, tmp_path):
        fa = tmp_path / "q.fa"
        seq = "ACGT" * 75
        rows = [
            self._make_row("42", 24, 0.50,
                           [(1050, 1100), (1150, 1200)], ucr_seq=seq),
            self._make_row("99", 36, 0.30,
                           [(1060, 1130)], ucr_seq=seq),
        ]
        ucr_seqs = {"42": seq, "99": seq}
        count = ua.write_non_unique_query_fasta(rows, ucr_seqs, str(fa))
        assert count == 3
        text = fa.read_text()
        assert "NUQUERY|42|k24|chr1|1050|1100" in text
        assert "NUQUERY|42|k24|chr1|1150|1200" in text
        assert "NUQUERY|99|k36|chr1|1060|1130" in text


# ===========================================================================
# 13. parse_paf()
# ===========================================================================


class TestParsePaf:
    """Unit tests for PAF parser."""

    def _paf_line(self, query="NUQUERY|42|k24|chr1|1000|1200",
                  qlen=200, qstart=0, qend=200,
                  strand="+", tname="chr1", tlen=248387497,
                  tstart=999, tend=1200, matches=200, alen=201, mapq=60):
        return "\t".join(str(x) for x in [
            query, qlen, qstart, qend, strand,
            tname, tlen, tstart, tend, matches, alen, mapq,
        ])

    def test_basic_parsing(self, tmp_path):
        paf = tmp_path / "test.paf"
        paf.write_text(self._paf_line() + "\n")
        records = ua.parse_paf(str(paf))
        assert len(records) == 1
        r = records[0]
        assert r["query_name"] == "NUQUERY|42|k24|chr1|1000|1200"
        assert r["query_len"] == 200
        assert r["hit_chrom"] == "chr1"
        assert r["hit_start"] == 999
        assert r["hit_end"] == 1200
        assert r["match_bp"] == 200
        assert r["align_len"] == 201
        assert r["mapq"] == 60
        assert r["strand"] == "+"

    def test_multiple_hits(self, tmp_path):
        paf = tmp_path / "test.paf"
        paf.write_text(
            self._paf_line(tname="chr1") + "\n"
            + self._paf_line(tname="chr5", tstart=5000, tend=5200, mapq=30)
            + "\n"
        )
        records = ua.parse_paf(str(paf))
        assert len(records) == 2
        assert records[1]["hit_chrom"] == "chr5"
        assert records[1]["mapq"] == 30

    def test_short_lines_skipped(self, tmp_path):
        paf = tmp_path / "test.paf"
        paf.write_text("col1\tcol2\tcol3\n")  # only 3 cols, needs 12
        records = ua.parse_paf(str(paf))
        assert records == []

    def test_empty_file(self, tmp_path):
        paf = tmp_path / "empty.paf"
        paf.write_text("")
        assert ua.parse_paf(str(paf)) == []

    def test_missing_file(self, tmp_path):
        assert ua.parse_paf(str(tmp_path / "nonexistent.paf")) == []


# ===========================================================================
# 14. write_non_unique_loci_report()
# ===========================================================================


class TestWriteNonUniqueLociReport:
    """Unit tests for the locus mapping TSV writer."""

    def _make_paf_record(self, query="NUQUERY|42|k24|chr1|1000|1200",
                          query_len=200, hit_chrom="chr1",
                          hit_start=999, hit_end=1200, strand="+",
                          match_bp=200, align_len=201, mapq=60):
        return {
            "query_name": query,
            "query_len": query_len,
            "hit_chrom": hit_chrom,
            "hit_start": hit_start,
            "hit_end": hit_end,
            "strand": strand,
            "match_bp": match_bp,
            "align_len": align_len,
            "mapq": mapq,
        }

    def test_basic_output_format(self, tmp_path):
        path = tmp_path / "loci.tsv"
        records = [self._make_paf_record()]
        nu_frac_map = {("42", 24): 0.50}
        ua.write_non_unique_loci_report(str(path), records, nu_frac_map)

        text = path.read_text()
        assert "## UCR Non-Unique Locus Mapping Report" in text

        data_rows = list(csv.DictReader(
            (l for l in text.splitlines() if not l.startswith("##")),
            delimiter="\t",
        ))
        assert len(data_rows) == 1
        row = data_rows[0]
        assert row["ucr_id"] == "42"
        assert row["kmer"] == "24"
        assert row["nu_interval"] == "chr1:1000-1200"
        assert row["hit_chrom"] == "chr1"
        assert row["is_self_locus"] == "YES"

    def test_off_target_locus_marked_no(self, tmp_path):
        path = tmp_path / "loci.tsv"
        # Hit on a completely different chromosome
        records = [self._make_paf_record(
            hit_chrom="chr5", hit_start=5000, hit_end=5200,
        )]
        nu_frac_map = {("42", 24): 0.50}
        ua.write_non_unique_loci_report(str(path), records, nu_frac_map)

        data_rows = list(csv.DictReader(
            (l for l in path.read_text().splitlines()
             if not l.startswith("##")),
            delimiter="\t",
        ))
        assert data_rows[0]["is_self_locus"] == "NO"
        assert data_rows[0]["hit_chrom"] == "chr5"

    def test_identity_computed_correctly(self, tmp_path):
        path = tmp_path / "loci.tsv"
        records = [self._make_paf_record(match_bp=180, align_len=200)]
        nu_frac_map = {("42", 24): 0.25}
        ua.write_non_unique_loci_report(str(path), records, nu_frac_map)
        data_rows = list(csv.DictReader(
            (l for l in path.read_text().splitlines()
             if not l.startswith("##")),
            delimiter="\t",
        ))
        assert abs(float(data_rows[0]["identity"]) - 0.9) < 1e-5

    def test_invalid_query_name_skipped(self, tmp_path):
        path = tmp_path / "loci.tsv"
        bad_rec = self._make_paf_record(query="BADFORMAT_nodelimiters")
        ua.write_non_unique_loci_report(str(path), [bad_rec], {})
        data_rows = list(csv.DictReader(
            (l for l in path.read_text().splitlines()
             if not l.startswith("##")),
            delimiter="\t",
        ))
        assert len(data_rows) == 0

    def test_empty_records(self, tmp_path):
        path = tmp_path / "loci.tsv"
        ua.write_non_unique_loci_report(str(path), [], {})
        text = path.read_text()
        assert "## UCR Non-Unique Locus Mapping Report" in text
        data_rows = list(csv.DictReader(
            (l for l in text.splitlines() if not l.startswith("##")),
            delimiter="\t",
        ))
        assert len(data_rows) == 0

