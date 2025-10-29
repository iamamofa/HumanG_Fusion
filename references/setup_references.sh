#!/bin/bash
# ============================================================
# RNA-seq & Fusion Detection Reference Setup (GRCh38)
# Author: Justice Ohene Amofa
# Description: Downloads and prepares reference files for
# STAR, STAR-Fusion, and Kraken2.
# ============================================================

set -e  # Exit immediately if a command fails

# ----------------------------
# 1️⃣ Create Reference Directory
# ----------------------------
REF_DIR=~/refs
mkdir -p $REF_DIR
cd $REF_DIR

echo "✅ Created reference directory at: $REF_DIR"

# ----------------------------
# 2️⃣ Download Human Genome (FASTA)
# ----------------------------
echo "⬇️ Downloading GRCh38 genome (FASTA)..."
wget -c https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip -f Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38.fa
echo "✅ Genome FASTA ready at: $REF_DIR/GRCh38.fa"

# ----------------------------
# 3️⃣ Download GTF Annotation
# ----------------------------
echo "⬇️ Downloading GTF annotation (Ensembl v110)..."
wget -c https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
gunzip -f Homo_sapiens.GRCh38.110.gtf.gz
mv Homo_sapiens.GRCh38.110.gtf gencode.v41.annotation.gtf
echo "✅ GTF ready at: $REF_DIR/gencode.v41.annotation.gtf"

# ----------------------------
# 4️⃣ Build STAR Index
# ----------------------------
echo "🚀 Building STAR genome index (this may take hours)..."
mkdir -p $REF_DIR/STAR_index_GRCh38
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir $REF_DIR/STAR_index_GRCh38 \
     --genomeFastaFiles $REF_DIR/GRCh38.fa \
     --sjdbGTFfile $REF_DIR/gencode.v41.annotation.gtf \
     --sjdbOverhang 100
echo "✅ STAR index built at: $REF_DIR/STAR_index_GRCh38"

# ----------------------------
# 5️⃣ Download STAR-Fusion CTAT Resource Library
# ----------------------------
echo "⬇️ Downloading STAR-Fusion CTAT Resource Library..."
wget -c https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v41_CTAT_lib_Mar012022.plug.tar.gz
tar -xvzf GRCh38_gencode_v41_CTAT_lib_Mar012022.plug.tar.gz
mv GRCh38_gencode_v41_CTAT_lib_Mar012022.plug ctat_resource_lib
echo "✅ CTAT resource library ready at: $REF_DIR/ctat_resource_lib"

# ----------------------------
# 6️⃣ Download Kraken2 Database (standard)
# ----------------------------
echo "⬇️ Downloading Kraken2 standard database (large ~100GB)..."
mkdir -p $REF_DIR/kraken2_db
kraken2-build --standard --threads 8 --db $REF_DIR/kraken2_db
echo "✅ Kraken2 DB ready at: $REF_DIR/kraken2_db"

# ----------------------------
# 7️⃣ Export Paths (for future sessions)
# ----------------------------
echo "🧩 Exporting environment variables..."
cat <<EOF >> ~/.bashrc

# ---- RNA-seq Reference Paths ----
export STAR_INDEX="$HOME/refs/STAR_index_GRCh38"
export GENOME_FASTA="$HOME/refs/GRCh38.fa"
export GTF="$HOME/refs/gencode.v41.annotation.gtf"
export STAR_FUSION_CTAT_LIB="$HOME/refs/ctat_resource_lib"
export KRAKEN2_DB="$HOME/refs/kraken2_db"
EOF

source ~/.bashrc
echo "✅ Environment variables saved to ~/.bashrc and loaded."

echo "🎉 All references are ready! You can now run STAR or STAR-Fusion."
