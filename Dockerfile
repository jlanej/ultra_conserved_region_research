FROM python:3.12-slim

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY convert_ucr_to_t2t.py .

# All resources (liftOver binary, chain file, Excel) are downloaded at runtime
# into OUTPUT_DIR.  Default to /output so callers can bind-mount a host directory.
ENV OUTPUT_DIR=/output
VOLUME /output

ENTRYPOINT ["python", "convert_ucr_to_t2t.py"]
