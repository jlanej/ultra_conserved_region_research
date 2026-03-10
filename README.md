# Ultra Conserved Region Research

Automated pipeline to lift over the 481 human **Ultraconserved Regions (UCRs)** from the GRCh38/hg38 reference assembly to the T2T-CHM13v2.0 (Hs1) telomere-to-telomere assembly.

## Background

Ultraconserved Regions are segments of the human genome (≥ 200 bp) that are 100 %
identical between human, mouse, and rat. This project takes the canonical UCR
catalogue from [Giacopuzzi *et al.* 2020 (PMC6857462)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6857462/)
and converts the coordinates to the first complete, gapless human genome assembly
([T2T-CHM13v2.0](https://www.ncbi.nlm.nih.gov/assembly/GCA_009914755.4/)).

### Data source

[Supplementary Table S1](https://pmc.ncbi.nlm.nih.gov/articles/instance/6857462/bin/Supp_TableS1.xlsx)
provides the non-redundant list of UCRs together with overlapping gene annotations.
The first data row looks like:

| Gene ID | Associated gene name | Chr. | Gene start (bp) | Gene end (bp) | Gene length (bp) | Gene type | Number of UCRs in this gene | UCR ID | UCR start (bp) | UCR end (bp) | 10⁶×(Number of UCRs/Gene Length) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ENSG00000142611 | PRDM16 | 1 | 3069168 | 3438621 | 369454 | protein_coding | 1 | 5 | 3075350 | 3075602 | 2.7 |

## How it works

1. **Download** – The script downloads the Excel file and the
   `hg38ToHs1.over.chain.gz` chain file at runtime.  The UCSC
   [`liftOver`](https://genome.ucsc.edu/cgi-bin/hgLiftOver) binary is bundled
   inside the container image; for bare-metal runs it is downloaded automatically.
2. **Extract** – Unique UCR coordinates (`Chr.`, `UCR start`, `UCR end`, `UCR ID`)
   are extracted from the spreadsheet. Chromosome names are normalised to `chrN`
   format and start positions are converted from 1-based to 0-based for BED.
3. **Lift over** – `liftOver` maps each region from hg38 to T2T-CHM13v2.0.
4. **Audit** – A comprehensive audit report (`ucr_liftover_audit.tsv`) is
   generated pairing every input UCR with its liftover result, including:
   mapping status, coordinate shifts, size deltas, chromosome concordance,
   and SHA-256 checksums of input files for full reproducibility.
5. **Output** – Four files are produced in the working directory:

| File | Description |
|---|---|
| `ucr_hg38.bed` | Input UCR coordinates on hg38 |
| `ucr_t2t_chm13.bed` | Successfully mapped coordinates on T2T-CHM13 |
| `ucr_unmapped.bed` | Unmapped regions with liftOver failure reasons |
| `ucr_liftover_audit.tsv` | Full audit report for QC and traceability |

## Quick start — HPC one-liner (Apptainer)

No Docker daemon required. [Apptainer](https://apptainer.org/) (formerly
Singularity) is available on most HPC clusters. Everything downloads at runtime:

```bash
apptainer run docker://ghcr.io/jlanej/ultra_conserved_region_research:latest
```

That's it. The `liftOver` binary is bundled in the image; only the data files
(chain file, Excel spreadsheet) download at runtime. All output files land in
your current working directory. No build step, no `pip install`,
fully batteries-included.

To write output to a specific directory:

```bash
OUTPUT_DIR=/scratch/my_project apptainer run docker://ghcr.io/jlanej/ultra_conserved_region_research:latest
```

To build a reusable SIF image (useful on air-gapped nodes or for repeated runs):

```bash
apptainer build ucr-liftover.sif Apptainer.def
apptainer run ucr-liftover.sif
```

## Running locally

### Prerequisites

- Python 3.10+
- Linux x86-64 (the `liftOver` binary is platform-specific)

### Steps

```bash
pip install -r requirements.txt
python convert_ucr_to_t2t.py
```

The script is self-contained: it downloads the Excel file, `liftOver` binary,
and chain file automatically on first run.  Subsequent runs skip downloads if
the files are already present.

### With Docker

```bash
docker build -t ucr-liftover .
docker run --rm -v "$(pwd)/results:/output" ucr-liftover
```

Results will appear in the `results/` directory.

## Repository structure

```
.
├── convert_ucr_to_t2t.py          # Main liftover script
├── requirements.txt               # Python dependencies
├── Dockerfile                     # Containerised pipeline (Docker / GitHub Actions)
├── Apptainer.def                  # Apptainer/Singularity definition (HPC)
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

- **Provenance header** — timestamp, source paper, assembly versions, SHA-256
  checksums of input files
- **Per-region columns** — `ucr_id`, `hg38_chrom`, `hg38_start`, `hg38_end`,
  `hg38_length`, `status` (MAPPED/UNMAPPED), `t2t_chrom`, `t2t_start`,
  `t2t_end`, `t2t_length`, `delta_length`, `delta_start`, `chrom_changed`,
  `reason`
- **Summary footer** — total/mapped/unmapped counts, chromosome-change and
  size-change counts, mapping rate

Regions flagged with `chrom_changed=YES` or non-zero `delta_length` deserve
manual review. Unmapped regions include the failure reason from liftOver.

## References

- Giacopuzzi, E. *et al.* (2020). "Population-scale distribution and copy number
  variation of ultraconserved elements in the human genome."
  *Human Mutation*, 41(6), 1101–1111.
  [PMC6857462](https://pmc.ncbi.nlm.nih.gov/articles/PMC6857462/)
- Nurk, S. *et al.* (2022). "The complete sequence of a human genome."
  *Science*, 376(6588), 44–53.
  [doi:10.1126/science.abj6987](https://doi.org/10.1126/science.abj6987)
