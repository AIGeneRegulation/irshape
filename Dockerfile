# irshape -- intron retention measurement engine.
# Reference data (GTF/FASTA/mappability/PolyASite/mask caches) is NOT baked
# into this image -- it is a separate download (irshape-ref-GRCh38-gencodeV44/,
# see build_reference_bundle.py) mounted at runtime via
# `-v <bundle>:/ref:ro` (Docker) or `--bind <bundle>:/ref` (Apptainer; see
# Apptainer.def, the packaging target for HPC use).
#
# Pinned to the exact dependency versions the engine was developed and
# validated against (see CHANGELOG.md).
FROM python:3.10-slim

LABEL org.opencontainers.image.title="irshape" \
      org.opencontainers.image.description="Intron retention measurement engine (IRFinder-principle, corrected masking, coverage-shape aware)" \
      org.opencontainers.image.source="."

# build-essential/zlib/bz2/lzma headers: fallback in case a platform lacks a
# prebuilt wheel for pysam/pyBigWig and pip needs to compile from sdist.
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/irshape
COPY pyproject.toml README.md SCHEMA.md /opt/irshape/
COPY irshape /opt/irshape/irshape

RUN pip install --no-cache-dir \
        numpy==2.2.6 \
        pandas==2.3.3 \
        pysam==0.22.1 \
        pyBigWig==0.3.25 \
    && pip install --no-cache-dir --no-deps -e /opt/irshape

ENV IRSHAPE_REF=/ref
VOLUME ["/ref"]

ENTRYPOINT ["irshape"]
CMD ["--help"]
