# syntax=docker/dockerfile:1
# ------------------------------------------------------------
# Bioinformatics tools container (from install.sh)
# ------------------------------------------------------------

FROM ubuntu:22.04

# Avoid interactive apt prompts
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

# Basic system setup
RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
        curl wget git bzip2 gzip unzip build-essential default-jre-headless \
        libssl-dev libcurl4-openssl-dev procps vim locales ca-certificates \
        python3-venv && \
    rm -rf /var/lib/apt/lists/*

# Set locale (UTF-8)
RUN locale-gen en_US.UTF-8
ENV LANG=en_US.UTF-8
ENV LC_ALL=en_US.UTF-8

# ------------------------------------------------------------
# Install Miniconda
# ------------------------------------------------------------
ENV CONDA_DIR=/opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH

RUN curl -fsSL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p $CONDA_DIR && \
    rm /tmp/miniconda.sh && \
    conda config --set channel_priority strict && \
    conda clean -afy

# ------------------------------------------------------------
# Create conda environment with bioconda + conda-forge
# ------------------------------------------------------------
ENV ENV_NAME=bioinf

# Use Mamba to speed up installs
RUN conda install -y -n base -c conda-forge mamba && \
    mamba create -y -n $ENV_NAME -c conda-forge -c bioconda \
      python=3.10 \
      biopython \
      pandas \
      requests \
      sra-tools \
      sratoolkit \
      samtools \
      fastqc \
      cutadapt \
      kraken2 \
      star \
      star-fusion \
      arriba \
      fusioncatcher \
      mygene \
      pigz \
      nextflow && \
    conda clean -afy

# ------------------------------------------------------------
# Set environment activation automatically
# ------------------------------------------------------------
SHELL ["bash", "-c"]
RUN echo "source activate ${ENV_NAME}" >> ~/.bashrc
ENV PATH=$CONDA_DIR/envs/$ENV_NAME/bin:$PATH

# Default working directory
WORKDIR /data

# Test that tools are available
RUN conda run -n $ENV_NAME python -V && \
    conda run -n $ENV_NAME fastqc --version && \
    conda run -n $ENV_NAME samtools --version

# ------------------------------------------------------------
# Container metadata
# ------------------------------------------------------------
LABEL maintainer="User 1" \
      description="Docker image with bioinformatics tools installed via bioconda and conda-forge"

# ------------------------------------------------------------
# Default command
# ------------------------------------------------------------
CMD ["/bin/bash"]
