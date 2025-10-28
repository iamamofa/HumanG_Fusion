# Nextflow Human RNA-Seq Fusion Detection Pipeline

![Pipeline Overview](image1.png)

This **Nextflow DSL2** pipeline automates **human-only RNA-Seq fusion detection** using:
- **STAR alignment (2-pass)**
- **Fusion callers**: STAR-Fusion, Arriba, FusionCatcher (optional)
- **Optional Kraken2 decontamination**
- **Postprocessing & merging** with gene annotation

---

## 📦 Repository Layout

```
nextflow-fusion-pipeline/
├── README.md
├── nextflow.config
├── main.nf
├── modules/
│   ├── preprocess/
│   ├── kraken_decontam/
│   ├── star_align/
│   ├── fusion_callers/
│   └── postprocess/
├── conf/
│   └── containers.config
├── docker/
│   ├── Dockerfile
│   └── fusion_env.Dockerfile
├── singularity/
│   └── Singularity.def
├── tools/
│   └── merge_fusions.py
├── bin/
│   └── helper_scripts/
└── sample_data/
```

---

# 🚀 Getting Started

## 1️⃣ Installation Setup

Choose **one** of the following methods to set up your environment:

---

### 🐳 **Docker**
**Best for:** Reproducible, containerized execution.

#### Build the Docker image:
```bash
docker build -t fusion-pipeline:latest docker/
```
#### Or, use Singularity:
```bash
singularity build fusion-pipeline.sif singularity/Singularity.def
```

---

### 🐍 **Python/Conda**
**Best for:** Local development or HPC with Conda.

#### Install Miniconda (if missing):
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p ~/miniconda3
source ~/miniconda3/etc/profile.d/conda.sh
conda create -n bioinf python=3.9
conda activate bioinf
```
#### Install tools:
```bash
conda install -c bioconda -c conda-forge fastqc cutadapt kraken2 star star-fusion arriba fusioncatcher mygene samtools nextflow
pip install mygene
```

---

### 📦 **Bash (Manual)**
**Best for:** Custom or existing environments.

#### Install dependencies:
```bash
sudo apt-get update && sudo apt-get install -y wget git samtools pigz python3 python3-pip
pip install mygene
```
#### Download and install each tool manually (see [SRA Toolkit](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit) for reference).

---

## 2️⃣ Data Preparation

**Place your FASTQ files in the `data/` directory:**
```
project/
└── data/
    ├── sample1_R1.fastq.gz
    ├── sample1_R2.fastq.gz
    ├── sample2_R1.fastq.gz
    └── sample2_R2.fastq.gz
```

---

## 3️⃣ Pipeline Configuration

Edit `pipeline.config.sh` with your paths:
```bash
#!/bin/bash
export READS="./data/*_R1.fastq.gz"
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

### 🐳 **Docker**
```bash
nextflow run main.nf -profile docker \
    --reads "$READS" \
    --outdir "$OUTDIR" \
    --star_index "$STAR_INDEX" \
    --genome_fasta "$GENOME_FASTA" \
    --gtf "$GTF" \
    --star_fusion_ctat_lib "$STAR_FUSION_CTAT_LIB" \
    --kraken_decontam $KRAKEN_DECONTAM \
    --kraken2_db "$KRAKEN2_DB"
```

### 🐍 **Conda/Python**
```bash
nextflow run main.nf -profile conda \
    --reads "$READS" \
    --outdir "$OUTDIR" \
    --star_index "$STAR_INDEX" \
    --genome_fasta "$GENOME_FASTA" \
    --gtf "$GTF" \
    --star_fusion_ctat_lib "$STAR_FUSION_CTAT_LIB" \
    --kraken_decontam $KRAKEN_DECONTAM \
    --kraken2_db "$KRAKEN2_DB"
```

### 📦 **Bash (Manual)**
```bash
nextflow run main.nf \
    --reads "$READS" \
    --outdir "$OUTDIR" \
    --star_index "$STAR_INDEX" \
    --genome_fasta "$GENOME_FASTA" \
    --gtf "$GTF" \
    --star_fusion_ctat_lib "$STAR_FUSION_CTAT_LIB" \
    --kraken_decontam $KRAKEN_DECONTAM \
    --kraken2_db "$KRAKEN2_DB"
```

---

# 🔧 Pipeline Modules

## **Preprocessing**
- FastQC, trimming (cutadapt)
- Output: trimmed FASTQ

## **Kraken2 Decontamination** (optional)
- Classifies reads, keeps only human/unclassified

## **STAR Alignment**
- 2-pass alignment, outputs BAM and chimeric junctions

## **Fusion Callers**
- STAR-Fusion, Arriba, FusionCatcher (optional)

## **Postprocessing**
- Merges calls, annotates genes, produces summary TSV

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
- [ ] Wrap STAR-Fusion to use STAR outputs
- [ ] Create a `tests/` folder with test FASTQs

---

# 💬 Questions?

Open an issue or contact the maintainers!
