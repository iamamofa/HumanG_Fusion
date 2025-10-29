# HumanG_Fusion: RNA-Seq Fusion Detection Pipeline

![Pipeline Overview](image1.png)

**HumanG_Fusion** is a modular, Nextflow-based pipeline for **human RNA-Seq fusion detection**, featuring:
- **Preprocessing** (FastQC, trimming)
- **Optional Kraken2 decontamination**
- **STAR alignment (2-pass)**
- **Fusion calling** (STAR-Fusion, Arriba, FusionCatcher)
- **Postprocessing & merging** with gene annotation

---

## 📦 Repository Layout

```
HumanG_Fusion/
├── README.md
├── installation_setup/
│   ├── Docker_setup/
│   ├── bash_setup/
│   └── python_setup/
├── docker/
├── bin/
├── conf/
├── modules/
├── scripts/
├── bash_scripts/
├── tools/
├── sample_data/
├── main.nf
├── pipeline.config.sh
└── run_pipeline.sh
```

---

# 🚀 Getting Started

## 1️⃣ Installation Setup

Choose **one** of the following methods:

---

### 🐳 **Docker**
**Best for:** Reproducible, containerized execution.

#### Setup:
```bash
cd installation_setup/Docker_setup
# Follow README.md for build/run instructions
docker build -t humang_fusion:latest .
```
#### Or, use Docker Compose:
```bash
cd docker_compose_setup
docker compose build
```

---

### 📦 **Bash**
**Best for:** Quick setup on Linux.

#### Setup:
```bash
cd installation_setup/bash_setup
chmod +x install.sh
./install.sh
```
This will:
- Install required apt packages (if sudo)
- Install Miniconda under `~/miniconda3` (if missing)
- Create a `bioinf` environment via mamba
- Install all conda + pip packages

#### Activate environment:
```bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate bioinf
# or, if fallback venv was used:
source ~/pyenv-bioinf/bin/activate
```

---

### 🐍 **Python**
**Best for:** Python-centric environments.

#### Setup:
```bash
cd installation_setup/python_setup
chmod +x install.py
./install.py
```
This will:
- Install Miniconda (if missing)
- Create a `bioinf` environment
- Install all dependencies

#### Activate environment:
```bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate bioinf
```

---

## 2️⃣ Data Preparation

Place your **FASTQ files** in `sample_data/`:
```
sample_data/
├── sample1_R1.fastq.gz
├── sample1_R2.fastq.gz
├── sample2_R1.fastq.gz
└── sample2_R2.fastq.gz
```

<br>

---

##  2.1 Reference Preparation


**Please Visit ➡️** **[References](\references\README.md)**


<br>


## 3️⃣ Pipeline Configuration

Edit `pipeline.config.sh` with your paths:
```bash
#!/bin/bash
export READS="./sample_data/*_R1.fastq.gz"
export OUTDIR="./results"
export STAR_INDEX="/refs/STAR_index_GRCh38"
export GENOME_FASTA="/refs/GRCh38.fa"
export GTF="/refs/gencode.v41.annotation.gtf"
export STAR_FUSION_CTAT_LIB="/refs/ctat_resource_lib"
export KRAKEN2_DB="/refs/kraken2_db"
export THREADS=8
export KRAKEN_DECONTAM=true
export RUN_FUSIONCATCHER=false
```

---

## 4️⃣ Run the Pipeline

### 🏃 **Automated (Recommended)**
```bash
chmod +x run_pipeline.sh
./run_pipeline.sh
```
This script will:
- Source your config
- Run Nextflow with the correct profile (Docker, Conda, or Bash)
- Execute all steps in order

### 🛠 **Manual Step-by-Step**
If you prefer, you can run each step manually using the scripts in `bash_scripts/` or `scripts/`.

---

# 🔧 Pipeline Modules

| Module                | Description                                                                 |
|-----------------------|-----------------------------------------------------------------------------|
| `preprocess`          | FastQC, trimming (cutadapt)                                                |
| `kraken_decontam`     | Optional: Remove non-human reads                                            |
| `star_align`          | STAR 2-pass alignment                                                       |
| `fusion_callers`      | STAR-Fusion, Arriba, FusionCatcher (optional)                              |
| `postprocess`         | Merge calls, annotate genes, produce summary TSV                           |

---

# 📂 Outputs

- `results/{sample}/fastq_trimmed/`: Trimmed FASTQ
- `results/{sample}/align/`: BAM files
- `results/{sample}/fusion/`: Caller-specific outputs
- `results/reports/`: Merged, annotated fusion calls

---

# 🧪 Testing & Validation

- Test on small datasets with known fusions.
- Validate artefacts in IGV.

---

# 📝 Notes

- **STAR-Fusion CTAT resource library** must be provided.
- **Kraken2 DB** required if using decontamination.
- **mygene** is used for gene annotation.

---

# 🚀 Next Steps

- [ ] Add Singularity/Docker push automation
- [ ] Add Nextflow Tower profiles
- [ ] Create a `tests/` folder with test FASTQs

---

# 💬 Questions?

Open an issue or contact the maintainers!
```bash
Justiceoheneamofa@gmail.com / kdanquah@atu.edu.gh / kdanquah@noguchi.ug.edu.gh
```