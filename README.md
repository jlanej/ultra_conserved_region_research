# Ultra Conserved Region Research

Automated pipeline to lift over human **Ultraconserved Regions (UCRs)** from the GRCh38/hg38 reference assembly to the T2T-CHM13v2.0 (Hs1) telomere-to-telomere assembly, with optional per-base sequence validation.

## Background

Ultraconserved elements (UCEs) were first identified by
[Bejerano *et al.* 2004](https://doi.org/10.1126/science.1098119), who reported
genomic segments ≥ 200 bp that are 100 % identical between human, rat, and
mouse.

This project uses the UCSC Genome Browser `compGeno/unusualcons/ultras` track
as the hg38 source of ultraconserved regions and lifts those coordinates from
GRCh38/hg38 to the first complete, gapless human genome assembly
([T2T-CHM13v2.0](https://www.ncbi.nlm.nih.gov/assembly/GCA_009914755.4/)).

### Data source

Bundled `resources/hg38.ultraConserved.bb` is the primary source file. It was
exported from the UCSC table browser schema:
[`compGeno` → `unusualcons` → `ultras`](https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=compGeno&hgta_track=unusualcons&hgta_table=ultras&hgta_doSchema=describe+table+schema).
The pipeline converts this bigBed file to BED4 (`chrom`, `start`, `end`,
`name`) before liftover.

## How it works

### Step 1 – Liftover (`convert_ucr_to_t2t.py`)

1. **Source** – The bundled `resources/hg38.ultraConserved.bb` is used as the
   primary input. If the file is not found (e.g. a bare-metal run outside the
   repository tree), it is downloaded from UCSC (`gbdb/hg38/bbi/ultras.bb`).
   The `hg38ToHs1.over.chain.gz` chain file is downloaded at runtime.
2. **Extract** – UCSC [`bigBedToBed`](https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/)
   converts `hg38.ultraConserved.bb` to BED. The pipeline keeps BED4 columns
   (`chrom`, `start`, `end`, `name`) and normalises chromosome names to `chrN`.
3. **Lift over** – `liftOver` maps each region from hg38 to T2T-CHM13v2.0.
4. **Audit** – A comprehensive audit report (`ucr_liftover_audit.tsv`) is
   generated pairing every input UCR with its liftover result, including:
   mapping status, coordinate shifts, size deltas, chromosome concordance,
   and SHA-256 checksums of input files for reproducibility. A structured
   summary (`ucr_liftover_audit_summary.json`) is also produced for automation
   and CI diffing.
5. **Output** – Five files are produced in the output directory:

| File | Description |
|---|---|
| `ucr_hg38.bed` | Input UCR coordinates on hg38 |
| `ucr_t2t_chm13.bed` | Successfully mapped coordinates on T2T-CHM13 |
| `ucr_unmapped.bed` | Unmapped regions with liftOver failure reasons |
| `ucr_liftover_audit.tsv` | Full audit report for QC and traceability |
| `ucr_liftover_audit_summary.json` | Machine-readable rollup (counts, reason breakdown, delta-length stats, hashes) |

### Step 2 – Sequence validation (`validate_liftover.py`)

An optional validation step that extracts the actual genomic sequences
underlying each UCR from both reference assemblies and performs per-base
pairwise alignment to confirm sequence identity across the liftover.

1. **Download genomes** – The hg38 and T2T-CHM13v2.0 reference genomes are
   downloaded in UCSC 2bit format (~800 MB each, first run only).
2. **Extract sequences** – The UCSC
   [`twoBitToFa`](https://genome.ucsc.edu/goldenPath/help/twoBit.html) utility
   extracts the UCR intervals from each 2bit file using the BED coordinates
   produced by the liftover step, writing one FASTA file per assembly.
3. **Pairwise alignment** – Each hg38 UCR sequence is paired with its
   corresponding T2T-CHM13 sequence and aligned using a
   [Needleman–Wunsch](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm)
   global pairwise alignment (via
   [BioPython `PairwiseAligner`](https://biopython.org/docs/latest/api/Bio.Align.html#Bio.Align.PairwiseAligner)).
   Identical-length pairs use a direct per-base comparison for efficiency;
   pairs that differ in length receive a full gapped alignment.
4. **Reports** – Two reports are generated:

| File | Description |
|---|---|
| `ucr_sequences_hg38.fa` | Extracted hg38 UCR sequences (FASTA) |
| `ucr_sequences_t2t.fa` | Extracted T2T-CHM13 UCR sequences (FASTA) |
| `ucr_alignment_report.tsv` | Per-region alignment statistics (identity %, matches, mismatches, gaps) |
| `ucr_alignment_details.txt` | Visual base-by-base alignments for any non-identical pairs |

## Software and dependencies

### External tools (bundled in container images)

| Tool | Version | Source | Purpose |
|---|---|---|---|
| [UCSC `bigBedToBed`](https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/) | latest linux.x86_64 | [UCSC Downloads](https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/) | Convert UCSC bigBed (`ultras.bb`) to BED |
| [UCSC `bigBedInfo`](https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/) | latest linux.x86_64 | [UCSC Downloads](https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/) | Inspect bigBed schema/metadata to verify binary vs scored intervals |
| [UCSC `bigBedSummary`](https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/) | latest linux.x86_64 | [UCSC Downloads](https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/) | Summarise mappability values from UCSC bigBed tracks |
| [UCSC `liftOver`](https://genome.ucsc.edu/cgi-bin/hgLiftOver) | latest linux.x86_64 | [UCSC Downloads](https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/) | Coordinate conversion between genome assemblies |
| [UCSC `twoBitToFa`](https://genome.ucsc.edu/goldenPath/help/twoBit.html) | latest linux.x86_64 | [UCSC Downloads](https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/) | Sequence extraction from 2bit genome files |

All UCSC binaries above are baked into the Docker and Apptainer container images at
build time. For bare-metal runs, binaries required by a given script are
downloaded automatically on first execution.

### Python dependencies (`requirements.txt`)

| Package | Purpose |
|---|---|
| [pandas](https://pandas.pydata.org/) | BED file parsing, deduplication, and type conversion |
| [requests](https://docs.python-requests.org/) | HTTP client library |
| [biopython](https://biopython.org/) | FASTA I/O (`Bio.SeqIO`) and Needleman–Wunsch pairwise alignment (`Bio.Align.PairwiseAligner`) |

### Runtime data files

| File | Size | Source | Used by |
|---|---|---|---|
| `resources/hg38.ultraConserved.bb` | ~43 KB | Bundled in repository (fallback: [UCSC gbdb ultras.bb](https://hgdownload.soe.ucsc.edu/gbdb/hg38/bbi/ultras.bb)) | `convert_ucr_to_t2t.py` |
| `hg38ToHs1.over.chain.gz` | ~10 MB | [UCSC goldenPath](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHs1.over.chain.gz) | `convert_ucr_to_t2t.py` |
| `hg38.2bit` | ~800 MB | [UCSC goldenPath](https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit) | `validate_liftover.py` |
| `hs1.2bit` | ~780 MB | [UCSC goldenPath](https://hgdownload.cse.ucsc.edu/goldenPath/hs1/bigZips/hs1.2bit) | `validate_liftover.py` |
| `k24/k36/k50/k100/k150/k250.Unique.Mappability.bb` | large track files (downloaded at runtime) | [UCSC gbdb hs1/hoffmanMappability](https://hgdownload.soe.ucsc.edu/gbdb/hs1/hoffmanMappability/) | `compute_unique_fraction.py` |
| `hs1.chrom.sizes` | small text file | [UCSC bigZips](https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.chrom.sizes) | `compute_unique_fraction.py` |

What the `.bb` files represent in this project:

- `hg38.ultraConserved.bb` is a UCSC **bigBed** file for the hg38
  `compGeno/unusualcons/ultras` dataset: genomic intervals for ultraconserved
  regions (classically, 200+ bp regions with complete human/rat/mouse
  conservation in the original UCE definition). The pipeline converts this
  track to BED4 and then performs liftover.
- `k*.Unique.Mappability.bb` files are UCSC **bigBed** tracks from
  `hs1/hoffmanMappability` describing unique mappability at each k-mer size.
  Depending on track schema, intervals may be binary or score-bearing; when
  scores are present, this project reports strict uniqueness by retaining only
  intervals at the maximum score for that track.

### Container platforms

| Platform | Definition file | Notes |
|---|---|---|
| [Docker](https://www.docker.com/) | `Dockerfile` | Base image: `python:3.12-slim` |
| [Apptainer](https://apptainer.org/) (formerly Singularity) | `Apptainer.def` | Recommended for HPC; no daemon required |

## Quick start — HPC one-liner (Apptainer)

No Docker daemon required. [Apptainer](https://apptainer.org/) (formerly
Singularity) is available on most HPC clusters. Everything downloads at runtime:

```bash
apptainer run --bind "$(pwd):/output" docker://ghcr.io/jlanej/ultra_conserved_region_research:latest
```

That's it. The `--bind` flag maps your current directory to the container's
`/output` so results are written to the host filesystem. UCSC tools are bundled
in the image; data/tool files (`ultras.bb`, `hg38ToHs1` chain,
`bigBedToBed`) download at runtime when not already present. All output files
land in your current working directory.
No build step, no `pip install`, fully batteries-included.

To write output to a specific directory:

```bash
apptainer run --bind /scratch/my_project:/output docker://ghcr.io/jlanej/ultra_conserved_region_research:latest
```

To build a reusable SIF image (useful on air-gapped nodes or for repeated runs):

```bash
apptainer build ucr-liftover.sif Apptainer.def
apptainer run ucr-liftover.sif
```

### Sequence validation (Apptainer)

After the liftover completes, run the validation utility to extract sequences
and perform pairwise alignment (downloads ~1.6 GB of genome data on first run):

```bash
apptainer exec --bind "$(pwd):/output" docker://ghcr.io/jlanej/ultra_conserved_region_research:latest \
    python /app/validate_liftover.py
```

### hs1 unique-genome fraction by k-mer (Apptainer one-liner)

Compute strict unique mappability fractions on **T2T/hs1** across
`k24`, `k36`, `k50`, `k100`, `k150`, and `k250`:

```bash
apptainer exec --bind "$(pwd):/output" docker://ghcr.io/jlanej/ultra_conserved_region_research:latest \
    python /app/compute_unique_fraction.py
```

Optional denominator and k-mer filters:

```bash
# Primary chromosomes only (chr1-22, X, Y, M)
apptainer exec --bind "$(pwd):/output" docker://ghcr.io/jlanej/ultra_conserved_region_research:latest \
    python /app/compute_unique_fraction.py --primary-only

# Exclude chrM from denominator and numerator
apptainer exec --bind "$(pwd):/output" docker://ghcr.io/jlanej/ultra_conserved_region_research:latest \
    python /app/compute_unique_fraction.py --primary-only --exclude-chrM

# Run a custom subset of kmers
apptainer exec --bind "$(pwd):/output" docker://ghcr.io/jlanej/ultra_conserved_region_research:latest \
    python /app/compute_unique_fraction.py --kmers 24,50,100
```

The module writes:

- `hs1.unique_fraction_by_kmer.tsv` (per-k strict unique fraction table)
- `hs1.unique_fraction_by_kmer_by_chrom.tsv` (per-k, per-chromosome strict unique fraction table)
- `hs1.unique_fraction_comparison.tsv` (per-k deltas vs previous k)
- `hs1.unique_fraction_summary.txt` (concise scientific summary text)

`hs1.unique_fraction_by_kmer.tsv` columns:

| Column | Meaning |
|---|---|
| `assembly` | Assembly label for the mappability track (currently `hs1`, i.e. T2T-CHM13v2.0). |
| `kmer` | K-mer size used for the corresponding UCSC mappability track (`k24`, `k36`, `k50`, `k100`, `k150`, `k250`). |
| `unique_bp` | Number of base pairs in intervals counted as unique. For scored tracks, this is the union of intervals at the **maximum** score only (strict unique); for binary tracks, this is the union of all intervals. |
| `genome_bp` | Denominator in base pairs from `hs1.chrom.sizes` after applying `--primary-only` and/or `--exclude-chrM` filters. |
| `fraction_unique` | `unique_bp / genome_bp` as a fraction in `[0, 1]`. |
| `percent_unique` | `fraction_unique * 100`. |
| `strict_unique_filter_applied` | `TRUE` when the track had mixed scores and non-maximal intervals were excluded; `FALSE` when no score-based filtering was needed. |

`hs1.unique_fraction_by_kmer_by_chrom.tsv` columns:

| Column | Meaning |
|---|---|
| `assembly` | Assembly label for the mappability track (`hs1`). |
| `kmer` | K-mer size used for the track. |
| `chromosome` | Chromosome name from `hs1.chrom.sizes` after any filters. |
| `unique_bp` | Unique covered bp on that chromosome for the given k-mer. |
| `genome_bp` | Chromosome size bp for that chromosome (denominator component). |
| `fraction_unique` | `unique_bp / genome_bp` for that chromosome. |
| `percent_unique` | `fraction_unique * 100` for that chromosome. |
| `strict_unique_filter_applied` | Same strict-score filtering flag used for the corresponding k-mer track. |

It also prints `bigBedInfo` metadata per track so you can verify whether the
track behaves as binary intervals or contains scored values. If scored values
are present, the numerator keeps only intervals at the maximum score (i.e.
strictly unique regions) and excludes lower-scored/partial intervals.

## Running locally

### Prerequisites

- Python 3.10+
- Linux x86-64 (the UCSC binaries are platform-specific)

### Steps

```bash
pip install -r requirements.txt

# Step 1: Liftover
python convert_ucr_to_t2t.py

# Step 2: Sequence validation (optional – downloads ~1.6 GB of genome data)
python validate_liftover.py

# Step 3: hs1 unique fraction summary across multiple kmers (optional)
python compute_unique_fraction.py
```

Both scripts are self-contained: they download all required data files and
tools automatically on first run. Subsequent runs skip downloads if the files
are already present.

### With Docker

```bash
docker build -t ucr-liftover .

# Step 1: Liftover
docker run --rm -v "$(pwd)/results:/output" ucr-liftover

# Step 2: Sequence validation (optional)
docker run --rm -v "$(pwd)/results:/output" ucr-liftover /app/validate_liftover.py
```

Results will appear in the `results/` directory.

## Repository structure

```
.
├── convert_ucr_to_t2t.py          # Liftover pipeline (hg38 → T2T-CHM13)
├── validate_liftover.py           # Sequence extraction and pairwise alignment
├── compute_unique_fraction.py     # hs1 strict unique mappability summary across multiple kmers
├── requirements.txt               # Python dependencies
├── Dockerfile                     # Docker container definition
├── Apptainer.def                  # Apptainer/Singularity definition (HPC)
├── resources/
│   └── hg38.ultraConserved.bb     # Bundled UCSC ultras bigBed source (hg38)
├── tests/
│   └── test_convert_ucr_to_t2t.py # Integration tests
├── .github/workflows/
│   ├── docker-build-publish.yml   # Build & push image to GHCR
│   └── ucr-liftover.yml           # Run liftover and commit results
└── data/                          # (generated) BED + audit output files
```

## CI / CD

Two GitHub Actions workflows automate the full pipeline:

| Workflow | Trigger | What it does |
|---|---|---|
| **Docker Build & Publish** | Push to `main` (when pipeline files change), manual | Builds the Docker image and pushes it to GitHub Container Registry (`ghcr.io`) |
| **UCR Liftover** | After a successful Docker build, or manual | Pulls the latest image, runs the liftover, and commits the resulting BED + audit files to `data/` |

To run manually, go to **Actions → UCR Liftover → Run workflow**.

## Audit report

The audit report (`ucr_liftover_audit.tsv`) is a tab-separated file with:

- **Provenance header** — timestamp, UCSC source dataset, assembly versions,
  SHA-256 checksums of input files
- **Per-region columns** — `ucr_id`, `hg38_chrom`, `hg38_start`, `hg38_end`,
  `hg38_length`, `status` (MAPPED/UNMAPPED), `t2t_chrom`, `t2t_start`,
  `t2t_end`, `t2t_length`, `delta_length`, `delta_start`, `chrom_changed`,
  `reason`
- **Summary footer** — total/mapped/unmapped counts, chromosome-change and
  size-change counts, mapping rate, and unmapped-reason counts

Regions flagged with `chrom_changed=YES` or non-zero `delta_length` deserve
manual review. Unmapped regions include the failure reason from liftOver.

The JSON companion (`ucr_liftover_audit_summary.json`) captures the same
high-level audit metrics in machine-readable form:
total/mapped/unmapped counts, mapping rate, unmapped reason histogram,
delta-length min/max/mean, and input file hashes.

## Alignment report

The alignment report (`ucr_alignment_report.tsv`) contains one row per
successfully lifted-over UCR with:

- hg38 and T2T-CHM13 coordinates and sequence lengths
- Whether the extracted sequences are identical (`YES` / `NO`)
- Percent identity, match count, mismatch count, gap count, and aligned length

The companion file `ucr_alignment_details.txt` provides visual base-by-base
alignments for any non-identical pairs, making it easy to inspect exactly which
positions differ.

## Methods

### Coordinate liftover

Coordinates from UCSC `unusualcons/ultras` (`hg38.ultraConserved.bb`) are
converted to BED format (0-based half-open intervals) using `bigBedToBed`.
Chromosome names are normalised to UCSC `chrN` convention. The UCSC `liftOver`
tool maps each region from GRCh38/hg38 to T2T-CHM13v2.0 (Hs1) using the
`hg38ToHs1.over.chain.gz` chain file from the
[UCSC goldenPath](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/).

### Sequence extraction

Reference genome sequences are obtained in UCSC
[2bit format](https://genome.ucsc.edu/goldenPath/help/twoBit.html) — a
compact binary encoding (~800 MB per genome). The UCSC `twoBitToFa` utility
extracts the genomic intervals defined by each BED file, producing a
multi-record FASTA file per assembly. Record names correspond to UCR IDs.

### Pairwise alignment

Each hg38 UCR sequence is globally aligned to its T2T-CHM13 counterpart using
the Needleman–Wunsch algorithm as implemented by BioPython's
[`PairwiseAligner`](https://biopython.org/docs/latest/api/Bio.Align.html#Bio.Align.PairwiseAligner)
with the following scoring scheme:

| Parameter | Value |
|---|---|
| Match score | +2 |
| Mismatch score | −1 |
| Gap open penalty | −3 |
| Gap extend penalty | −0.5 |

For equal-length sequence pairs, a direct per-base character comparison is
used (no gaps possible), which is equivalent to a gap-free global alignment.
Identity is reported as the fraction of aligned positions that are exact
matches.

## References

- Bejerano, G. *et al.* (2004). "Ultraconserved elements in the human genome."
  *Science*, 304(5675), 1321–1325.
  [doi:10.1126/science.1098119](https://doi.org/10.1126/science.1098119)
- UCSC Table Browser schema for `unusualcons/ultras` (hg38):
  [schema description](https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=compGeno&hgta_track=unusualcons&hgta_table=ultras&hgta_doSchema=describe+table+schema)
- Nurk, S. *et al.* (2022). "The complete sequence of a human genome."
  *Science*, 376(6588), 44–53.
  [doi:10.1126/science.abj6987](https://doi.org/10.1126/science.abj6987)
- Hinrichs, A. S. *et al.* (2006). "The UCSC Genome Browser Database: update 2006."
  *Nucleic Acids Research*, 34(suppl_1), D590–D598.
  [doi:10.1093/nar/gkj144](https://doi.org/10.1093/nar/gkj144)
- Cock, P. J. A. *et al.* (2009). "Biopython: freely available Python tools for
  computational molecular biology and bioinformatics."
  *Bioinformatics*, 25(11), 1422–1423.
  [doi:10.1093/bioinformatics/btp163](https://doi.org/10.1093/bioinformatics/btp163)
