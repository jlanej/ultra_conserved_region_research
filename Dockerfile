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
COPY validate_liftover.py .

# Bake the UCSC binaries into the image (small tools, ~3 MB each).
# Large data files (chain file, Excel, genomes) are downloaded at runtime.
# Use retries for resilience against transient UCSC server issues.
RUN wget --tries=3 --waitretry=5 -O liftOver \
        "https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver" && \
    wget --tries=3 --waitretry=5 -O twoBitToFa \
        "https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa" && \
    chmod +x liftOver twoBitToFa

ENV OUTPUT_DIR=/output
VOLUME /output

ENTRYPOINT ["python"]
CMD ["/app/convert_ucr_to_t2t.py"]
