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

1. **Download** – The script fetches the Excel file and the UCSC
   [`liftOver`](https://genome.ucsc.edu/cgi-bin/hgLiftOver) binary together with
   the `hg38ToHs1.over.chain.gz` chain file.
2. **Extract** – Unique UCR coordinates (`Chr.`, `UCR start`, `UCR end`, `UCR ID`)
   are extracted from the spreadsheet. Chromosome names are normalised to `chrN`
   format and start positions are converted from 1-based to 0-based for the BED
   standard.
3. **Lift over** – `liftOver` maps each region from hg38 to T2T-CHM13v2.0.
4. **Output** – Three BED files are produced:
   - `ucr_hg38.bed` – input coordinates on hg38
   - `ucr_t2t_chm13.bed` – successfully mapped coordinates on T2T-CHM13
   - `ucr_unmapped.bed` – any regions that could not be mapped

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
└── data/                          # (generated) BED output files
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

The script is self-contained: it downloads the Excel file, the `liftOver` binary,
and the chain file automatically on first run.

### With Docker

```bash
docker build -t ucr-liftover .
docker run --rm -v "$(pwd)/results:/output" ucr-liftover
```

Results will appear in the `results/` directory.

### On HPC with Apptainer (one-liner)

No Docker daemon required. Apptainer (formerly Singularity) is available on most
HPC clusters. Just run:

```bash
apptainer run docker://ghcr.io/jlanej/ultra_conserved_region_research:latest
```

Output BED files are written to the current working directory. That's it — no
build step, no dependencies to install, fully batteries-included.

To build a reusable SIF image (useful on air-gapped nodes or for repeated runs):

```bash
# Build once (requires internet)
apptainer build ucr-liftover.sif docker://ghcr.io/jlanej/ultra_conserved_region_research:latest

# Run anywhere
apptainer run ucr-liftover.sif
```

Or build from the definition file in this repository:

```bash
apptainer build ucr-liftover.sif Apptainer.def
apptainer run ucr-liftover.sif
```

To write output to a specific directory:

```bash
OUTPUT_DIR=/scratch/my_project apptainer run ucr-liftover.sif
```

## CI / CD

Two GitHub Actions workflows automate the full pipeline:

| Workflow | Trigger | What it does |
|---|---|---|
| **Docker Build & Publish** | Push to `main` (when pipeline files change), manual | Builds the Docker image and pushes it to GitHub Container Registry (`ghcr.io`) |
| **UCR Liftover** | After a successful Docker build, or manual | Pulls the latest image, runs the liftover, and commits the resulting BED files to `data/` |

To run manually, go to **Actions → UCR Liftover → Run workflow**.

## References

- Giacopuzzi, E. *et al.* (2020). "Population-scale distribution and copy number
  variation of ultraconserved elements in the human genome."
  *Human Mutation*, 41(6), 1101–1111.
  [PMC6857462](https://pmc.ncbi.nlm.nih.gov/articles/PMC6857462/)
- Nurk, S. *et al.* (2022). "The complete sequence of a human genome."
  *Science*, 376(6588), 44–53.
  [doi:10.1126/science.abj6987](https://doi.org/10.1126/science.abj6987)
