#!/usr/bin/env bash
# install.sh - install bioinformatics tools in a single conda environment (preferred)
# If conda/mamba install fails, fallback to creating a python virtualenv with python deps.
#
# Usage:
#   chmod +x install.sh
#   ./install.sh
#
# Notes:
# - This script is written to be safe and idempotent. Re-run if interrupted.
# - Installs Miniconda to $HOME/miniconda3 if conda is not present.
# - Creates conda env named "bioinf" and installs tools from bioconda + conda-forge via mamba.
# - If conda route fails, creates a venv named "pyenv-bioinf" with python packages only.
set -euo pipefail
IFS=$'\n\t'

LOG() { printf "\n[install.sh] %s\n" "$*"; }

# -------- configuration --------
CONDA_DIR="${HOME}/miniconda3"
ENV_NAME="bioinf"
PY_VENV_DIR="${HOME}/pyenv-bioinf"
# list of conda packages (channels: conda-forge, bioconda)
CONDA_PKGS=(
  "python=3.10"
  "mamba"
  "biopython"
  "pandas"
  "requests"
  "sra-tools"        # sratoolkit (contains prefetch, fasterq-dump)
  "sratoolkit"       # sometimes named sratoolkit; sra-tools is common
  "sra-tools"        # ok may be duplicate but conda resolves
  "samtools"
  "fastqc"
  "cutadapt"
  "kraken2"
  "star"
  "star-fusion"
  "arriba"
  "fusioncatcher"
  "mygene"
  "pigz"
  "nextflow"
)
# pip-only packages to ensure in env
PIP_PKGS=(
  "requests"
)

# apt packages (system-level helpers)
APT_PKGS=(curl wget git bzip2 gzip unzip build-essential default-jre-headless \
    libssl-dev libcurl4-openssl-dev procps vim locales ca-certificates \
    python3-venv)

# -------- helper functions --------
command_exists() { command -v "$1" >/dev/null 2>&1; }

install_apt_pkgs() {
  LOG "Installing apt packages (requires sudo)..."
  sudo apt-get update -y
  sudo apt-get install -y --no-install-recommends "${APT_PKGS[@]}"
  LOG "Apt packages installed."
}

install_miniconda() {
  if [[ -x "${CONDA_DIR}/bin/conda" ]]; then
    LOG "Conda already installed at ${CONDA_DIR}"
    return 0
  fi

  LOG "Downloading Miniconda installer..."
  TMPSH="$(mktemp -u)/miniconda.sh"
  TMPDIR="$(mktemp -d)"
  INSTALLER="${TMPDIR}/miniconda.sh"
  curl -fsSL "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" -o "${INSTALLER}"
  chmod +x "${INSTALLER}"
  LOG "Installing Miniconda to ${CONDA_DIR} (non-interactive)..."
  bash "${INSTALLER}" -b -p "${CONDA_DIR}"
  rm -rf "${TMPDIR}"
  LOG "Miniconda installed."

  # initialize conda for bash (affects current shell via eval)
  source "${CONDA_DIR}/etc/profile.d/conda.sh"
  conda init bash >/dev/null 2>&1 || true
  LOG "Configured conda. You may need to re-open your shell to use 'conda' normally."
}

create_conda_env() {
  LOG "Ensuring mamba is available in base..."
  # source conda
  if [[ -f "${CONDA_DIR}/etc/profile.d/conda.sh" ]]; then
    # shellcheck disable=SC1090
    source "${CONDA_DIR}/etc/profile.d/conda.sh"
  fi

  if ! command_exists conda; then
    LOG "conda not found in PATH even after installing Miniconda. Trying to source manually..."
    export PATH="${CONDA_DIR}/bin:${PATH}"
    if ! command_exists conda; then
      LOG "conda still not found. Aborting conda env creation."
      return 1
    fi
  fi

  # install mamba in base if not present
  if ! conda list -n base mamba >/dev/null 2>&1; then
    LOG "Installing mamba into base environment..."
    conda install -n base -y -c conda-forge mamba
  else
    LOG "mamba already installed in base."
  fi

  LOG "Creating conda environment '${ENV_NAME}' with bioconda/conda-forge packages..."
  # use mamba to create env
  mamba create -y -n "${ENV_NAME}" -c conda-forge -c bioconda "${CONDA_PKGS[@]}" || {
    LOG "mamba failed to create environment (non-zero exit)."
    return 2
  }

  LOG "Activating '${ENV_NAME}' and installing pip packages..."
  conda activate "${ENV_NAME}"
  python -m pip install --upgrade pip
  if [[ ${#PIP_PKGS[@]} -gt 0 ]]; then
    python -m pip install "${PIP_PKGS[@]}"
  fi

  LOG "Conda environment '${ENV_NAME}' created and configured."
  return 0
}

create_venv() {
  LOG "Creating Python virtualenv fallback at ${PY_VENV_DIR}..."
  python3 -m venv "${PY_VENV_DIR}"
  # shellcheck disable=SC1090
  source "${PY_VENV_DIR}/bin/activate"
  python -m pip install --upgrade pip
  python -m pip install biopython pandas requests
  LOG "Virtualenv created. Activate with: source ${PY_VENV_DIR}/bin/activate"
}

final_messages() {
  cat <<EOF

============================================================
INSTALLATION SUMMARY
============================================================

If conda/mamba route succeeded:
  - Activate your environment with:
      source ${CONDA_DIR}/etc/profile.d/conda.sh
      conda activate ${ENV_NAME}

  - Tools installed (via conda bioconda/conda-forge):
      ${CONDA_PKGS[*]}

  - Python packages installed:
      ${PIP_PKGS[*]}

  - sratoolkit commands (prefetch, fasterq-dump) should be available in the env.
    Example usage:
      conda activate ${ENV_NAME}
      prefetch SRRxxxxxx
      fasterq-dump SRRxxxxxx --split-files --outdir ./fastq

If conda route failed and fallback venv used:
  - Activate the venv:
      source ${PY_VENV_DIR}/bin/activate
  - Python deps available: biopython, pandas, requests
  - You still need to install heavy binary tools that are not pip-installable:
      - sratoolkit (prefetch, fasterq-dump)
      - STAR, STAR-Fusion, fusioncatcher, kraken2, etc.
    Consider installing those with conda (recommended) or using a Docker image.

Notes:
  - You may need to re-open the shell or source your ~/.bashrc to get conda initialisation.
  - To add automatic activation on login, add this to your ~/.bashrc (optional):
      source ${CONDA_DIR}/etc/profile.d/conda.sh
      conda activate ${ENV_NAME}

============================================================
EOF
}

# -------- main flow --------
LOG "Starting installation."

# 1) Install apt packages (if sudo available)
if command_exists sudo; then
  LOG "Detected sudo; installing system prerequisites via apt-get."
  install_apt_pkgs
else
  LOG "sudo not found. Skipping apt packages installation. Make sure the system has curl/wget/git/python3-venv."
fi

# 2) Try to install conda (if conda not present)
if command_exists conda; then
  LOG "Conda already installed and available in PATH."
else
  LOG "Conda not found. Installing Miniconda into ${CONDA_DIR}..."
  install_miniconda
fi

# 3) Try creating conda env
set +e
create_conda_env
RET=$?
set -e

if [[ ${RET} -eq 0 ]]; then
  LOG "Conda environment installation successful."
  final_messages
  exit 0
else
  LOG "Conda environment creation failed (exit code ${RET}). Attempting fallback: Python virtualenv."
  create_venv
  final_messages
  exit 0
fi
# End of script