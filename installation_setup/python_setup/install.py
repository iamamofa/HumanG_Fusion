#!/usr/bin/env python3
"""
install.py — unified installer for bioinformatics tools

This script installs:
  - system-level deps (via apt if available)
  - Miniconda (if missing)
  - Conda env with bioconda + conda-forge tools
  - Fallback: pure Python venv if conda fails

Usage:
  python3 install.py
"""

import os
import sys
import subprocess
import shutil
from pathlib import Path

# -----------------------------------------
# Configuration
# -----------------------------------------
CONDA_DIR = Path.home() / "miniconda3"
ENV_NAME = "bioinf"
PY_VENV_DIR = Path.home() / "pyenv-bioinf"

CONDA_PKGS = [
    "python=3.10",
    "mamba",
    "biopython",
    "pandas",
    "requests",
    "sra-tools",
    "samtools",
    "fastqc",
    "cutadapt",
    "kraken2",
    "star",
    "star-fusion",
    "arriba",
    "fusioncatcher",
    "mygene",
    "pigz",
    "nextflow",
]

PIP_PKGS = ["requests"]
APT_PKGS = [
    "curl", "wget", "git", "bzip2", "gzip", "unzip", "build-essential",
    "default-jre-headless", "libssl-dev", "libcurl4-openssl-dev", "procps",
    "vim", "locales", "ca-certificates", "python3-venv"
]

# -----------------------------------------
# Utilities
# -----------------------------------------
def log(msg: str):
    print(f"\n[install.py] {msg}", flush=True)

def run(cmd, check=True, env=None):
    log(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=check, env=env)

def command_exists(cmd):
    return shutil.which(cmd) is not None

# -----------------------------------------
# Steps
# -----------------------------------------
def install_apt_packages():
    if not command_exists("apt-get"):
        log("No apt-get found, skipping system deps.")
        return
    try:
        log("Installing system dependencies with apt-get...")
        run(["sudo", "apt-get", "update", "-y"])
        run(["sudo", "apt-get", "install", "-y", "--no-install-recommends", *APT_PKGS])
        log("Apt packages installed.")
    except subprocess.CalledProcessError:
        log("Warning: apt-get installation failed. Continuing anyway.")

def install_miniconda():
    if (CONDA_DIR / "bin" / "conda").exists():
        log(f"Conda already installed at {CONDA_DIR}")
        return
    log("Downloading and installing Miniconda...")
    url = "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    installer = Path("/tmp/miniconda.sh")
    run(["wget", "-q", url, "-O", str(installer)])
    run(["bash", str(installer), "-b", "-p", str(CONDA_DIR)])
    log("Miniconda installed.")
    (CONDA_DIR / "etc/profile.d").mkdir(parents=True, exist_ok=True)
    log("Conda setup complete.")

def create_conda_env():
    conda_path = CONDA_DIR / "bin" / "conda"
    mamba_path = CONDA_DIR / "bin" / "mamba"
    env = os.environ.copy()
    env["PATH"] = f"{CONDA_DIR}/bin:{env['PATH']}"

    if not conda_path.exists():
        raise RuntimeError("Conda not found after installation.")

    # Install mamba in base if needed
    if not mamba_path.exists():
        log("Installing mamba in base environment...")
        run([str(conda_path), "install", "-n", "base", "-y", "-c", "conda-forge", "mamba"], env=env)

    # Create env
    log(f"Creating conda env '{ENV_NAME}' with bioconda + conda-forge packages...")
    run([
        str(mamba_path), "create", "-y", "-n", ENV_NAME,
        "-c", "conda-forge", "-c", "bioconda", *CONDA_PKGS
    ], env=env)

    # Install pip packages
    log(f"Installing pip packages into conda env '{ENV_NAME}'...")
    run([
        str(conda_path), "run", "-n", ENV_NAME, "python", "-m", "pip", "install", "--upgrade", "pip"
    ], env=env)
    run([
        str(conda_path), "run", "-n", ENV_NAME, "python", "-m", "pip", "install", *PIP_PKGS
    ], env=env)

    log(f"Conda environment '{ENV_NAME}' created successfully.")

def create_venv():
    log(f"Creating fallback Python virtualenv at {PY_VENV_DIR}...")
    run(["python3", "-m", "venv", str(PY_VENV_DIR)])
    env_path = PY_VENV_DIR / "bin" / "activate"
    pip_exe = PY_VENV_DIR / "bin" / "pip"
    run([str(pip_exe), "install", "--upgrade", "pip"])
    run([str(pip_exe), "install", "biopython", "pandas", "requests"])
    log(f"Virtualenv created. Activate it with:\n  source {env_path}")

# -----------------------------------------
# Main
# -----------------------------------------
def main():
    log("Starting installation...")

    install_apt_packages()
    install_miniconda()

    try:
        create_conda_env()
        log("✅ Conda environment installation successful.")
    except Exception as e:
        log(f"⚠️ Conda installation failed: {e}")
        log("Attempting fallback to Python virtualenv.")
        create_venv()

    print("""
============================================================
INSTALLATION SUMMARY
============================================================
If conda/mamba route succeeded:
  source ~/miniconda3/etc/profile.d/conda.sh
  conda activate bioinf

If fallback venv used:
  source ~/pyenv-bioinf/bin/activate

Tools installed:
  sratoolkit, samtools, fastqc, star, star-fusion, arriba, fusioncatcher, kraken2, nextflow, etc.

============================================================
""")

if __name__ == "__main__":
    main()
