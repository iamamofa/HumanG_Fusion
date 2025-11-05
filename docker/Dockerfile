# docker/Dockerfile
# Full image that installs conda and the bioinformatics tools via conda (bioconda + conda-forge)
# This image will be large (several GB) because of STAR, STAR-Fusion CTAT libs (not included) and other tools.
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# Basic system deps
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates curl wget git bzip2 gzip unzip build-essential \
    ca-certificates default-jre-headless libssl-dev libcurl4-openssl-dev \
    procps vim locales \
  && rm -rf /var/lib/apt/lists/*

# Set locale to avoid warnings in some tools
RUN locale-gen en_US.UTF-8
ENV LANG=en_US.UTF-8
ENV LC_ALL=en_US.UTF-8

# Install Miniconda (installer) and mamba for faster installs
ENV CONDA_DIR=/opt/conda
ENV PATH=${CONDA_DIR}/bin:${PATH}
RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p ${CONDA_DIR} && rm /tmp/miniconda.sh && \
    ${CONDA_DIR}/bin/conda init bash

# Accept conda ToS
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Ensure channels and install packages via mamba (faster)
RUN conda install -y -n base -c conda-forge mamba && \
    mamba config set show_channel_urls true && \
    mamba update -n base -y conda

# Install core bioinformatics packages. Adjust packages or add versions if you want strict pinning.
# Note: STAR-Fusion CTAT lib is not installed here (too large). Mount it at runtime with --mount or -v.
RUN mamba install -y -c bioconda -c conda-forge \
    nextflow \
    python \
    pip \
    samtools \
    fastqc \
    cutadapt \
    kraken2 \
    star \
        star-fusion \
    mygene \
    pigz \
  && mamba clean -afy

# Add pipeline tools into the image (copy repository tools). Expect build context to include tools/ and bin/
COPY tools/ /opt/pipeline/tools/
COPY bin/ /opt/pipeline/bin/
RUN chmod -R a+rX /opt/pipeline && chmod +x /opt/pipeline/tools/*.py || true

# Working directory for container runtime
WORKDIR /data

# Helpful environment variables
ENV PATH=/opt/pipeline/bin:$PATH
ENV NF_HOME=/data

# Default entrypoint: nextflow (so `docker run nextflow ...` works). Users can override CMD or ENTRYPOINT.
ENTRYPOINT ["nextflow"]
CMD ["-h"]
