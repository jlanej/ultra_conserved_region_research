FROM python:3.12-slim

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        wget \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY convert_ucr_to_t2t.py .

# Download UCSC liftOver binary and chain file at build time
RUN wget -q -O liftOver \
        "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver" && \
    chmod +x liftOver && \
    wget -q -O hg38ToHs1.over.chain.gz \
        "https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/liftOver/hg38ToHs1.over.chain.gz"

# Default output to /output so callers can bind-mount a host directory there
ENV OUTPUT_DIR=/output
VOLUME /output

ENTRYPOINT ["python", "convert_ucr_to_t2t.py"]
