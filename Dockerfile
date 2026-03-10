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

# Bake the liftOver binary into the image (it is a small tool, ~3 MB).
# Large data files (chain file, Excel) are downloaded at runtime.
# Use retries for resilience against transient UCSC server issues.
RUN wget --tries=3 --waitretry=5 -O liftOver \
        "https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver" && \
    chmod +x liftOver

ENV OUTPUT_DIR=/output
VOLUME /output

ENTRYPOINT ["python", "convert_ucr_to_t2t.py"]
