"""
Integration tests for convert_ucr_to_t2t.py.

Covers:
  - Bundled resources/hg38.ultraConserved.bb presence and format
  - normalize_chrom() with various chromosome name conventions
  - extract_coordinates() via mocked bigBedToBed conversion
  - parse_bed() and parse_unmapped() with synthetic data
  - generate_audit_report() with mock liftOver output files
  - Full pipeline smoke-test using stub bigBedToBed and liftOver scripts
"""

import csv
import os
import re
import stat
import sys
import textwrap

import pytest

# ---------------------------------------------------------------------------
# Make the repo root importable so we can import convert_ucr_to_t2t without
# installing it as a package.
# ---------------------------------------------------------------------------
REPO_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, REPO_ROOT)

import convert_ucr_to_t2t as uct  # noqa: E402  (import after sys.path tweak)

RESOURCES_DIR = os.path.join(REPO_ROOT, "resources")
BUNDLED_ULTRAS_BB = os.path.join(RESOURCES_DIR, "hg38.ultraConserved.bb")

MIN_UCR_COUNT = 3


# ===========================================================================
# 1. Bundled UCSC ultras bigBed file
# ===========================================================================


class TestBundledUltrasBigBed:
    """Validate the checked-in resources/hg38.ultraConserved.bb."""

    def test_file_exists(self):
        assert os.path.exists(BUNDLED_ULTRAS_BB), (
            f"Bundled bigBed file not found at {BUNDLED_ULTRAS_BB}. "
            "Run: git lfs pull  or verify the resources/ directory."
        )

    def test_file_is_non_empty_binary(self):
        assert os.path.getsize(BUNDLED_ULTRAS_BB) > 0

    def test_file_has_bigbed_magic(self):
        # bigBed magic number: 0x8789F2EB (little-endian bytes below)
        with open(BUNDLED_ULTRAS_BB, "rb") as fh:
            magic = fh.read(4)
        assert magic == b"\xeb\xf2\x89\x87"


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
# 3. extract_coordinates() – mocked bigBedToBed conversion
# ===========================================================================


class TestExtractCoordinates:
    """Integration tests for extract_coordinates() with mocked conversion."""

    def _mock_bigbed_to_bed(self, monkeypatch, tmp_path):
        raw_bed = tmp_path / "ultras_hg38_raw.bed"
        out_bed = tmp_path / "ucr_hg38.bed"

        def _fake_run(cmd, capture_output, text):
            assert cmd[0] == "fake-bigBedToBed"
            assert cmd[1].endswith(".bb")
            assert cmd[2] == str(raw_bed)
            raw_bed.write_text(
                "1\t100\t200\tuc.1\n"
                "chrX\t500\t650\tuc.2\n"
                "MT\t900\t1000\tuc.3\n"
            )
            class _Result:
                returncode = 0
                stderr = ""
            return _Result()

        monkeypatch.setattr(uct, "BIGBEDTOBED_BIN", "fake-bigBedToBed")
        monkeypatch.setattr(uct, "ULTRAS_BB_FILE", str(tmp_path / "input.bb"))
        monkeypatch.setattr(uct, "ULTRAS_RAW_BED", str(raw_bed))
        monkeypatch.setattr(uct, "HG38_BED", str(out_bed))
        monkeypatch.setattr(uct.subprocess, "run", _fake_run)
        return out_bed

    def test_produces_bed_file(self, tmp_path, monkeypatch):
        self._mock_bigbed_to_bed(monkeypatch, tmp_path)

        bed_df = uct.extract_coordinates()

        assert os.path.exists(uct.HG38_BED), "BED file was not written."
        assert len(bed_df) >= MIN_UCR_COUNT

    def test_no_nan_chromosomes_in_bed(self, tmp_path, monkeypatch):
        self._mock_bigbed_to_bed(monkeypatch, tmp_path)

        uct.extract_coordinates()

        regions = uct.parse_bed(uct.HG38_BED)
        bad = [r for r in regions if not r["chrom"] or r["chrom"] == "nan"]
        assert not bad, (
            f"{len(bad)} BED rows have a missing/NaN chromosome: {bad[:3]}"
        )

    def test_chromosomes_use_chr_prefix(self, tmp_path, monkeypatch):
        self._mock_bigbed_to_bed(monkeypatch, tmp_path)

        uct.extract_coordinates()

        regions = uct.parse_bed(uct.HG38_BED)
        bad = [r for r in regions if not r["chrom"].startswith("chr")]
        assert not bad, (
            f"{len(bad)} BED rows have chromosomes without 'chr' prefix: {bad[:3]}"
        )

    def test_bed_coordinates_are_zero_based(self, tmp_path, monkeypatch):
        """BED start must be 0-based (≥ 0); end must be > start."""
        self._mock_bigbed_to_bed(monkeypatch, tmp_path)

        uct.extract_coordinates()

        regions = uct.parse_bed(uct.HG38_BED)
        for r in regions:
            assert r["start"] >= 0, f"Negative start in {r}"
            assert r["end"] > r["start"], f"end <= start in {r}"

    def test_ucr_ids_are_unique(self, tmp_path, monkeypatch):
        self._mock_bigbed_to_bed(monkeypatch, tmp_path)

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
        ultras_bb = tmp_path / "hg38.ultraConserved.bb"

        # Copy the real bigBed so sha256sum() works
        import shutil
        shutil.copy(BUNDLED_ULTRAS_BB, str(ultras_bb))

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
        monkeypatch.setattr(uct, "ULTRAS_BB_FILE", str(ultras_bb))
        monkeypatch.setattr(uct, "CHAIN_FILE", str(ultras_bb))  # any real file for sha256

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
        assert "UCSC unusualcons/ultras table" in text
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

    Tiny shell scripts act as:
      - bigBedToBed: writes a small synthetic BED4 file
      - liftOver: copies input BED to output unchanged (100 % mapping)
    """

    def test_pipeline_runs_without_error(self, tmp_path, monkeypatch):
        # --- Write stub bigBedToBed binary ---
        bigbed = tmp_path / "bigBedToBed"
        bigbed.write_text(
            "#!/bin/sh\n"
            "# argv: bigBedToBed <in.bb> <out.bed>\n"
            "cat > \"$2\" <<'EOF'\n"
            "chr1\t100\t200\tuc.1\n"
            "chr2\t300\t400\tuc.2\n"
            "chrX\t500\t600\tuc.3\n"
            "EOF\n"
        )
        bigbed.chmod(bigbed.stat().st_mode | stat.S_IEXEC)

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

        ultras_bb = tmp_path / "dummy.bb"
        ultras_bb.write_bytes(b"bb")

        audit_path = tmp_path / "ucr_liftover_audit.tsv"

        monkeypatch.setattr(uct, "BIGBEDTOBED_BIN", str(bigbed))
        monkeypatch.setattr(uct, "ULTRAS_BB_FILE", str(ultras_bb))
        monkeypatch.setattr(uct, "ULTRAS_RAW_BED", str(tmp_path / "ultras_hg38_raw.bed"))
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
        bigbed = tmp_path / "bigBedToBed"
        bigbed.write_text(
            "#!/bin/sh\n"
            "cat > \"$2\" <<'EOF'\n"
            "chr1\t100\t200\tuc.1\n"
            "chr2\t300\t400\tuc.2\n"
            "chrX\t500\t600\tuc.3\n"
            "EOF\n"
        )
        bigbed.chmod(bigbed.stat().st_mode | stat.S_IEXEC)

        stub = tmp_path / "liftOver"
        stub.write_text("#!/bin/sh\ncp \"$1\" \"$3\"\ntouch \"$4\"\n")
        stub.chmod(stub.stat().st_mode | stat.S_IEXEC)

        chain = tmp_path / "dummy.chain.gz"
        chain.write_bytes(b"")
        ultras_bb = tmp_path / "dummy.bb"
        ultras_bb.write_bytes(b"bb")

        audit_path = tmp_path / "ucr_liftover_audit.tsv"

        monkeypatch.setattr(uct, "BIGBEDTOBED_BIN", str(bigbed))
        monkeypatch.setattr(uct, "ULTRAS_BB_FILE", str(ultras_bb))
        monkeypatch.setattr(uct, "ULTRAS_RAW_BED", str(tmp_path / "ultras_hg38_raw.bed"))
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
