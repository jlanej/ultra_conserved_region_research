"""
Integration tests for convert_ucr_to_t2t.py.

Covers:
  - Bundled resources/Supp_TableS1.xlsx format and content
  - normalize_chrom() with various chromosome name conventions
  - extract_coordinates() end-to-end against the real Excel file
  - parse_bed() and parse_unmapped() with synthetic data
  - generate_audit_report() with mock liftOver output files
  - Full pipeline smoke-test using a stub liftOver script
"""

import csv
import os
import re
import stat
import sys
import textwrap

import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Make the repo root importable so we can import convert_ucr_to_t2t without
# installing it as a package.
# ---------------------------------------------------------------------------
REPO_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, REPO_ROOT)

import convert_ucr_to_t2t as uct  # noqa: E402  (import after sys.path tweak)

RESOURCES_DIR = os.path.join(REPO_ROOT, "resources")
BUNDLED_EXCEL = os.path.join(RESOURCES_DIR, "Supp_TableS1.xlsx")

# Expected column names (after skiprows=1) from Table S1
EXPECTED_COLUMNS = {
    "Chr.",
    "UCR start (bp)",
    "UCR end (bp)",
    "UCR ID",
}

# Minimum number of unique UCRs we expect after filtering rows with no
# chromosome assignment.  Derived from resources/Supp_TableS1.xlsx.
MIN_UCR_COUNT = 700


# ===========================================================================
# 1. Bundled Excel file – format & content
# ===========================================================================


class TestBundledExcel:
    """Validate the checked-in resources/Supp_TableS1.xlsx."""

    def test_file_exists(self):
        assert os.path.exists(BUNDLED_EXCEL), (
            f"Bundled Excel file not found at {BUNDLED_EXCEL}. "
            "Run: git lfs pull  or verify the resources/ directory."
        )

    def test_file_is_readable_excel(self):
        df = pd.read_excel(BUNDLED_EXCEL, skiprows=1, engine="openpyxl")
        assert len(df) > 0, "Excel file is empty after skipping header row."

    def test_required_columns_present(self):
        df = pd.read_excel(BUNDLED_EXCEL, skiprows=1, engine="openpyxl")
        missing = EXPECTED_COLUMNS - set(df.columns)
        assert not missing, f"Missing columns in Table S1: {missing}"

    def test_ucr_id_column_non_empty(self):
        df = pd.read_excel(BUNDLED_EXCEL, skiprows=1, engine="openpyxl")
        assert df["UCR ID"].notna().any(), "UCR ID column is entirely NaN."

    def test_ucr_coordinates_are_numeric(self):
        df = pd.read_excel(BUNDLED_EXCEL, skiprows=1, engine="openpyxl")
        # At least one row must have numeric start/end values
        assert pd.to_numeric(df["UCR start (bp)"], errors="coerce").notna().any()
        assert pd.to_numeric(df["UCR end (bp)"], errors="coerce").notna().any()

    def test_start_less_than_end(self):
        df = pd.read_excel(BUNDLED_EXCEL, skiprows=1, engine="openpyxl")
        valid = df.dropna(subset=["UCR start (bp)", "UCR end (bp)"])
        bad = valid[valid["UCR start (bp)"] >= valid["UCR end (bp)"]]
        assert len(bad) == 0, (
            f"{len(bad)} rows have UCR start >= UCR end:\n{bad.head()}"
        )

    def test_unique_ucr_count_above_minimum(self):
        df = pd.read_excel(BUNDLED_EXCEL, skiprows=1, engine="openpyxl")
        ucr_df = df[["Chr.", "UCR start (bp)", "UCR end (bp)", "UCR ID"]].drop_duplicates()
        ucr_df = ucr_df.dropna(subset=["Chr."])
        assert len(ucr_df) >= MIN_UCR_COUNT, (
            f"Expected at least {MIN_UCR_COUNT} unique UCRs, got {len(ucr_df)}."
        )


# ===========================================================================
# 2. normalize_chrom()
# ===========================================================================


class TestNormalizeChrom:
    """Unit tests for the chromosome name normalisation helper."""

    @pytest.mark.parametrize("raw, expected", [
        # Plain numbers
        ("1", "chr1"),
        ("22", "chr22"),
        # Already prefixed (exact case)
        ("chr1", "chr1"),
        ("chrX", "chrX"),
        # Case variants
        ("CHR1", "chr1"),
        ("Chr2", "chr2"),
        ("chrx", "chrX"),
        ("CHRX", "chrX"),
        # Mitochondrial
        ("MT", "chrM"),
        ("mt", "chrM"),
        ("chrMT", "chrM"),
        ("CHRMT", "chrM"),
        # Letter chromosomes
        ("X", "chrX"),
        ("Y", "chrY"),
        # Leading/trailing whitespace
        ("  3  ", "chr3"),
        (" chrY ", "chrY"),
    ])
    def test_normalisation(self, raw, expected):
        assert uct.normalize_chrom(raw) == expected

    def test_nan_passthrough(self):
        import math
        result = uct.normalize_chrom(float("nan"))
        assert result is None or (isinstance(result, float) and math.isnan(result)), (
            "NaN input should be returned unchanged."
        )


# ===========================================================================
# 3. extract_coordinates() – real Excel file
# ===========================================================================


class TestExtractCoordinates:
    """Integration test: run extract_coordinates() against the bundled file."""

    def test_produces_bed_file(self, tmp_path, monkeypatch):
        monkeypatch.setattr(uct, "EXCEL_FILE", BUNDLED_EXCEL)
        monkeypatch.setattr(uct, "HG38_BED", str(tmp_path / "ucr_hg38.bed"))

        bed_df = uct.extract_coordinates()

        assert os.path.exists(uct.HG38_BED), "BED file was not written."
        assert len(bed_df) >= MIN_UCR_COUNT

    def test_no_nan_chromosomes_in_bed(self, tmp_path, monkeypatch):
        monkeypatch.setattr(uct, "EXCEL_FILE", BUNDLED_EXCEL)
        monkeypatch.setattr(uct, "HG38_BED", str(tmp_path / "ucr_hg38.bed"))

        uct.extract_coordinates()

        regions = uct.parse_bed(uct.HG38_BED)
        bad = [r for r in regions if not r["chrom"] or r["chrom"] == "nan"]
        assert not bad, (
            f"{len(bad)} BED rows have a missing/NaN chromosome: {bad[:3]}"
        )

    def test_chromosomes_use_chr_prefix(self, tmp_path, monkeypatch):
        monkeypatch.setattr(uct, "EXCEL_FILE", BUNDLED_EXCEL)
        monkeypatch.setattr(uct, "HG38_BED", str(tmp_path / "ucr_hg38.bed"))

        uct.extract_coordinates()

        regions = uct.parse_bed(uct.HG38_BED)
        bad = [r for r in regions if not r["chrom"].startswith("chr")]
        assert not bad, (
            f"{len(bad)} BED rows have chromosomes without 'chr' prefix: {bad[:3]}"
        )

    def test_bed_coordinates_are_zero_based(self, tmp_path, monkeypatch):
        """BED start must be 0-based (≥ 0); end must be > start."""
        monkeypatch.setattr(uct, "EXCEL_FILE", BUNDLED_EXCEL)
        monkeypatch.setattr(uct, "HG38_BED", str(tmp_path / "ucr_hg38.bed"))

        uct.extract_coordinates()

        regions = uct.parse_bed(uct.HG38_BED)
        for r in regions:
            assert r["start"] >= 0, f"Negative start in {r}"
            assert r["end"] > r["start"], f"end <= start in {r}"

    def test_ucr_ids_are_unique(self, tmp_path, monkeypatch):
        monkeypatch.setattr(uct, "EXCEL_FILE", BUNDLED_EXCEL)
        monkeypatch.setattr(uct, "HG38_BED", str(tmp_path / "ucr_hg38.bed"))

        uct.extract_coordinates()

        regions = uct.parse_bed(uct.HG38_BED)
        names = [r["name"] for r in regions]
        assert len(names) == len(set(names)), "Duplicate UCR IDs found in BED file."


# ===========================================================================
# 4. parse_bed()
# ===========================================================================


class TestParseBed:
    """Unit tests for parse_bed() with synthetic BED content."""

    def _write_bed(self, path, content):
        path.write_text(textwrap.dedent(content))

    def test_basic_parsing(self, tmp_path):
        bed = tmp_path / "test.bed"
        self._write_bed(bed, """\
            chr1\t100\t200\tucr1
            chr2\t300\t400\tucr2
        """)
        regions = uct.parse_bed(str(bed))
        assert len(regions) == 2
        assert regions[0] == {"chrom": "chr1", "start": 100, "end": 200, "name": "ucr1"}
        assert regions[1] == {"chrom": "chr2", "start": 300, "end": 400, "name": "ucr2"}

    def test_comment_lines_skipped(self, tmp_path):
        bed = tmp_path / "test.bed"
        self._write_bed(bed, """\
            # this is a comment
            chr1\t100\t200\tucr1
        """)
        regions = uct.parse_bed(str(bed))
        assert len(regions) == 1

    def test_empty_lines_skipped(self, tmp_path):
        bed = tmp_path / "test.bed"
        self._write_bed(bed, """\
            chr1\t100\t200\tucr1

            chr3\t500\t600\tucr3
        """)
        regions = uct.parse_bed(str(bed))
        assert len(regions) == 2

    def test_missing_file_returns_empty(self, tmp_path):
        regions = uct.parse_bed(str(tmp_path / "nonexistent.bed"))
        assert regions == []

    def test_insufficient_columns_ignored(self, tmp_path):
        bed = tmp_path / "test.bed"
        self._write_bed(bed, "chr1\t100\n")
        regions = uct.parse_bed(str(bed))
        assert regions == []


# ===========================================================================
# 5. parse_unmapped()
# ===========================================================================


class TestParseUnmapped:
    """Unit tests for parse_unmapped() with synthetic liftOver output."""

    def test_basic_parsing(self, tmp_path):
        um = tmp_path / "unmapped.bed"
        um.write_text(
            "# Deleted in new\nchr1\t100\t200\tucr1\n"
            "# Partially deleted in new\nchr2\t300\t400\tucr2\n"
        )
        result = uct.parse_unmapped(str(um))
        assert result == {
            "ucr1": "Deleted in new",
            "ucr2": "Partially deleted in new",
        }

    def test_missing_file_returns_empty(self, tmp_path):
        result = uct.parse_unmapped(str(tmp_path / "nonexistent.bed"))
        assert result == {}

    def test_reason_is_empty_when_no_preceding_comment(self, tmp_path):
        """parse_unmapped() should store an empty string when a BED entry has no
        preceding '#' comment line (i.e. the reason resets between entries)."""
        um = tmp_path / "unmapped.bed"
        um.write_text(
            "# Reason A\nchr1\t1\t2\tucr1\n"
            "chr2\t3\t4\tucr2\n"   # BED entry with no preceding comment
        )
        result = uct.parse_unmapped(str(um))
        assert result["ucr1"] == "Reason A"
        assert result["ucr2"] == ""   # no comment → empty reason


# ===========================================================================
# 6. generate_audit_report() with mock BED files
# ===========================================================================


class TestGenerateAuditReport:
    """Integration test: generate_audit_report() with fully synthetic inputs."""

    def _setup(self, tmp_path, monkeypatch):
        """Write minimal synthetic BED files and patch module-level paths."""
        hg38_bed = tmp_path / "ucr_hg38.bed"
        t2t_bed = tmp_path / "ucr_t2t_chm13.bed"
        unmapped_bed = tmp_path / "ucr_unmapped.bed"
        audit = tmp_path / "ucr_liftover_audit.tsv"
        excel = tmp_path / "Supp_TableS1.xlsx"

        # Copy the real Excel so sha256sum() works
        import shutil
        shutil.copy(BUNDLED_EXCEL, str(excel))

        hg38_bed.write_text(
            "chr1\t100\t350\t5\n"
            "chr1\t500\t750\t32\n"
            "chrX\t1000\t1253\t100\n"
        )
        t2t_bed.write_text(
            "chr1\t101\t351\t5\n"    # mapped, shifted by +1
            "chr1\t500\t750\t32\n"   # mapped, identical
        )
        unmapped_bed.write_text(
            "# Deleted in new\nchrX\t1000\t1253\t100\n"
        )

        monkeypatch.setattr(uct, "HG38_BED", str(hg38_bed))
        monkeypatch.setattr(uct, "T2T_BED", str(t2t_bed))
        monkeypatch.setattr(uct, "UNMAPPED_BED", str(unmapped_bed))
        monkeypatch.setattr(uct, "AUDIT_REPORT", str(audit))
        monkeypatch.setattr(uct, "EXCEL_FILE", str(excel))
        monkeypatch.setattr(uct, "CHAIN_FILE", str(excel))  # any real file for sha256

        return audit

    def test_audit_file_created(self, tmp_path, monkeypatch):
        audit = self._setup(tmp_path, monkeypatch)
        uct.generate_audit_report(None)
        assert audit.exists()

    def test_audit_mapped_count(self, tmp_path, monkeypatch):
        audit = self._setup(tmp_path, monkeypatch)
        uct.generate_audit_report(None)

        rows = [r for r in csv.DictReader(
            (l for l in audit.read_text().splitlines() if not l.startswith("#")),
            delimiter="\t",
        )]
        mapped = [r for r in rows if r["status"] == "MAPPED"]
        unmapped = [r for r in rows if r["status"] == "UNMAPPED"]
        assert len(mapped) == 2
        assert len(unmapped) == 1

    def test_audit_delta_length(self, tmp_path, monkeypatch):
        audit = self._setup(tmp_path, monkeypatch)
        uct.generate_audit_report(None)

        rows = {r["ucr_id"]: r for r in csv.DictReader(
            (l for l in audit.read_text().splitlines() if not l.startswith("#")),
            delimiter="\t",
        )}
        # ucr 5: hg38 250 bp, t2t 250 bp → delta_length == 0; shifted by +1
        assert int(rows["5"]["delta_length"]) == 0
        assert int(rows["5"]["delta_start"]) == 1
        # ucr 32: identical coordinates → delta 0, no shift
        assert int(rows["32"]["delta_length"]) == 0
        assert int(rows["32"]["delta_start"]) == 0
        # ucr 100: unmapped
        assert rows["100"]["status"] == "UNMAPPED"
        assert rows["100"]["reason"] == "Deleted in new"

    def test_audit_contains_provenance_header(self, tmp_path, monkeypatch):
        audit = self._setup(tmp_path, monkeypatch)
        uct.generate_audit_report(None)

        text = audit.read_text()
        assert "## UCR Liftover Audit Report" in text
        assert "SHA-256" in text
        assert "GRCh38" in text
        assert "T2T-CHM13" in text

    def test_audit_summary_footer(self, tmp_path, monkeypatch):
        audit = self._setup(tmp_path, monkeypatch)
        uct.generate_audit_report(None)

        text = audit.read_text()
        assert "SUMMARY" in text
        assert "Mapped:" in text
        assert "Unmapped:" in text


# ===========================================================================
# 7. Full pipeline smoke-test with a stub liftOver
# ===========================================================================


class TestFullPipelineSmoke:
    """
    End-to-end smoke-test: extract_coordinates() → stub liftOver → audit.

    A tiny shell script acts as liftOver: it copies the input BED to the
    output BED unchanged (simulating 100 % mapping) and writes an empty
    unmapped file.
    """

    def test_pipeline_runs_without_error(self, tmp_path, monkeypatch):
        # --- Write stub liftOver binary ---
        stub = tmp_path / "liftOver"
        stub.write_text(
            "#!/bin/sh\n"
            "# argv: liftOver <in.bed> <chain> <out.bed> <unmapped>\n"
            "cp \"$1\" \"$3\"\n"   # copy input BED → mapped BED
            "touch \"$4\"\n"       # create empty unmapped file
        )
        stub.chmod(stub.stat().st_mode | stat.S_IEXEC)

        chain = tmp_path / "dummy.chain.gz"
        chain.write_bytes(b"")   # liftOver stub ignores it

        audit_path = tmp_path / "ucr_liftover_audit.tsv"

        monkeypatch.setattr(uct, "EXCEL_FILE", BUNDLED_EXCEL)
        monkeypatch.setattr(uct, "HG38_BED", str(tmp_path / "ucr_hg38.bed"))
        monkeypatch.setattr(uct, "T2T_BED", str(tmp_path / "ucr_t2t_chm13.bed"))
        monkeypatch.setattr(uct, "UNMAPPED_BED", str(tmp_path / "ucr_unmapped.bed"))
        monkeypatch.setattr(uct, "AUDIT_REPORT", str(audit_path))
        monkeypatch.setattr(uct, "LIFTOVER_BIN", str(stub))
        monkeypatch.setattr(uct, "CHAIN_FILE", str(chain))

        bed_df = uct.extract_coordinates()
        uct.run_liftover()
        uct.generate_audit_report(bed_df)

        assert os.path.exists(uct.HG38_BED)
        assert os.path.exists(uct.T2T_BED)
        assert audit_path.exists()

    def test_pipeline_all_regions_mapped(self, tmp_path, monkeypatch):
        """With the passthrough stub, every extracted UCR must appear as MAPPED."""
        stub = tmp_path / "liftOver"
        stub.write_text("#!/bin/sh\ncp \"$1\" \"$3\"\ntouch \"$4\"\n")
        stub.chmod(stub.stat().st_mode | stat.S_IEXEC)

        chain = tmp_path / "dummy.chain.gz"
        chain.write_bytes(b"")

        audit_path = tmp_path / "ucr_liftover_audit.tsv"

        monkeypatch.setattr(uct, "EXCEL_FILE", BUNDLED_EXCEL)
        monkeypatch.setattr(uct, "HG38_BED", str(tmp_path / "ucr_hg38.bed"))
        monkeypatch.setattr(uct, "T2T_BED", str(tmp_path / "ucr_t2t_chm13.bed"))
        monkeypatch.setattr(uct, "UNMAPPED_BED", str(tmp_path / "ucr_unmapped.bed"))
        monkeypatch.setattr(uct, "AUDIT_REPORT", str(audit_path))
        monkeypatch.setattr(uct, "LIFTOVER_BIN", str(stub))
        monkeypatch.setattr(uct, "CHAIN_FILE", str(chain))

        bed_df = uct.extract_coordinates()
        uct.run_liftover()
        uct.generate_audit_report(bed_df)

        rows = [r for r in csv.DictReader(
            (l for l in audit_path.read_text().splitlines() if not l.startswith("#")),
            delimiter="\t",
        )]
        unmapped = [r for r in rows if r["status"] == "UNMAPPED"]
        assert unmapped == [], (
            f"Expected 0 unmapped regions with passthrough stub, got {len(unmapped)}: "
            f"{unmapped[:3]}"
        )
        assert len(rows) >= MIN_UCR_COUNT


# ===========================================================================
# 8. Container definitions include liftOver runtime dependencies
# ===========================================================================


class TestContainerDefinitions:
    """Integration tests for container runtime package requirements."""

    def _apt_installs_package(self, text, package):
        apt_install = re.compile(r"\bapt-get\s+install\b")
        blocks = []
        block = []
        collecting = False
        for line in text.splitlines():
            if apt_install.search(line):
                collecting = True
                block = []
            if collecting:
                block.append(line)
                if not line.rstrip().endswith("\\"):
                    blocks.append(block)
                    collecting = False
        for lines in blocks:
            normalized = re.sub(r"\s+", " ", " ".join(lines).replace("\\", " "))
            if re.search(rf"\b{re.escape(package)}\b", normalized):
                return True
        return False

    def test_dockerfile_installs_libcurl4(self):
        dockerfile = os.path.join(REPO_ROOT, "Dockerfile")
        with open(dockerfile, encoding="utf-8") as fh:
            text = fh.read()
        assert self._apt_installs_package(text, "libcurl4"), (
            "Dockerfile must install libcurl4 so UCSC liftOver can load "
            "libcurl.so.4 at runtime."
        )

    def test_apptainer_installs_libcurl4(self):
        apptainer_def = os.path.join(REPO_ROOT, "Apptainer.def")
        with open(apptainer_def, encoding="utf-8") as fh:
            text = fh.read()
        assert self._apt_installs_package(text, "libcurl4"), (
            "Apptainer definition must install libcurl4 so UCSC liftOver can load "
            "libcurl.so.4 at runtime."
        )
