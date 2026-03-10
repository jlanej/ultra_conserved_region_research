import os
import urllib.request
import subprocess
import pandas as pd
import stat

# --- Configuration ---
# Excel file URL (Supplementary Table S1 from PMC6857462)
# Table S1 provides the clean, non-redundant list of 481 UCRs with gene overlap info
EXCEL_URL = "https://pmc.ncbi.nlm.nih.gov/articles/instance/6857462/bin/Supp_TableS1.xlsx"

# UCSC LiftOver and Chain File URLs (used only when not already present)
LIFTOVER_URL = "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver"
CHAIN_URL = "https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/liftOver/hg38ToHs1.over.chain.gz"

# --- Path resolution ---
# Tools are baked into the container image at /app; when running locally they are
# downloaded to the current directory.  OUTPUT_DIR controls where results land
# (defaults to $PWD so Apptainer bind-mounts just work).
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.environ.get("OUTPUT_DIR", os.getcwd())

LIFTOVER_BIN = os.path.join(SCRIPT_DIR, "liftOver")
CHAIN_FILE = os.path.join(SCRIPT_DIR, "hg38ToHs1.over.chain.gz")
EXCEL_FILE = os.path.join(OUTPUT_DIR, "Supp_TableS1.xlsx")

HG38_BED = os.path.join(OUTPUT_DIR, "ucr_hg38.bed")
T2T_BED = os.path.join(OUTPUT_DIR, "ucr_t2t_chm13.bed")
UNMAPPED_BED = os.path.join(OUTPUT_DIR, "ucr_unmapped.bed")


def download_file(url, filename, make_executable=False):
    if not os.path.exists(filename):
        print(f"Downloading {filename}...")
        urllib.request.urlretrieve(url, filename)
        if make_executable:
            st = os.stat(filename)
            os.chmod(filename, st.st_mode | stat.S_IEXEC)
    else:
        print(f"{filename} already exists, skipping download.")


def extract_coordinates():
    print(f"Reading {EXCEL_FILE}...")
    # Read the Excel file, skipping the first row (the title "SUPPLEMENTARY TABLE S1...")
    df = pd.read_excel(EXCEL_FILE, skiprows=1)

    # Extract unique UCRs by their coordinates and ID
    # Table S1 may list the same UCR multiple times if it overlaps multiple genes
    ucr_df = df[['Chr.', 'UCR start (bp)', 'UCR end (bp)', 'UCR ID']].drop_duplicates()

    # Format chromosomes properly (e.g., '1' -> 'chr1')
    ucr_df['Chr.'] = ucr_df['Chr.'].apply(
        lambda x: str(x) if str(x).startswith('chr') else f"chr{x}"
    )

    # Create BED format dataframe (0-based start, 1-based end)
    # Assuming 1-based closed intervals in the Excel file, we subtract 1 from the start for BED format:
    ucr_df['BED_start'] = ucr_df['UCR start (bp)'].astype(int) - 1
    ucr_df['BED_end'] = ucr_df['UCR end (bp)'].astype(int)

    bed_df = ucr_df[['Chr.', 'BED_start', 'BED_end', 'UCR ID']]

    print(f"Extracted {len(bed_df)} unique UCR regions.")

    # Save to BED file
    bed_df.to_csv(HG38_BED, sep='\t', header=False, index=False)
    print(f"Saved hg38 coordinates to {HG38_BED}")


def run_liftover():
    print("Running UCSC liftOver...")
    cmd = [
        LIFTOVER_BIN,
        HG38_BED,
        CHAIN_FILE,
        T2T_BED,
        UNMAPPED_BED
    ]

    try:
        subprocess.run(cmd, check=True)
        print("LiftOver complete!")
        print(f"Successfully mapped coordinates saved to: {T2T_BED}")
        print(f"Unmapped coordinates saved to: {UNMAPPED_BED}")
    except subprocess.CalledProcessError as e:
        print(f"Error running liftOver: {e}")
        raise


if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # 1. Download necessary assets
    download_file(EXCEL_URL, EXCEL_FILE)
    download_file(LIFTOVER_URL, LIFTOVER_BIN, make_executable=True)
    download_file(CHAIN_URL, CHAIN_FILE)

    # 2. Extract and format hg38 coordinates
    extract_coordinates()

    # 3. Convert hg38 to T2T-CHM13
    run_liftover()

    # 4. Preview the results
    if os.path.exists(T2T_BED):
        print("\nPreview of T2T-CHM13 coordinates:")
        with open(T2T_BED, 'r') as f:
            for _ in range(5):
                print(f.readline().strip())
