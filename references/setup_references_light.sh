#!/bin/bash
# ============================================================
# Lightweight RNA-seq & Fusion Detection Reference Setup (GRCh38)
# Author: Justice Ohene Amofa
# Description: Downloads minimal references for STAR, STAR-Fusion, and Kraken2.
# Approx. total size: ~10–15 GB (instead of 150+ GB)
# ============================================================

set -e  # Stop if any command fails

# ----------------------------
# 1️⃣ Create Reference Directory
# ----------------------------
REF_DIR=~/refs
mkdir -p $REF_DIR
cd $REF_DIR
echo "✅ Reference directory created at: $REF_DIR"

# ----------------------------
# 2️⃣ Download a Small Subset of Human Genome (GRCh38)
# ----------------------------
echo "⬇️ Downloading subset of GRCh38 (chromosomes 21 & 22 for demo)..."
wget -c https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz
wget -c https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz

gunzip -f Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz
gunzip -f Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz

cat Homo_sapiens.GRCh38.dna.chromosome.*.fa > GRCh38_subset.fa
rm Homo_sapiens.GRCh38.dna.chromosome.*.fa
echo "✅ Subset genome FASTA ready at: $REF_DIR/GRCh38_subset.fa"

# ----------------------------
# 3️⃣ Download GTF Annotation (Subset)
# ----------------------------
echo "⬇️ Downloading GTF annotation (Ensembl v110)..."
wget -c https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
gunzip -f Homo_sapiens.GRCh38.110.gtf.gz
# Keep only chromosomes 21 and 22 annotations
awk '$1=="21" || $1=="22" || $1 ~ /^#/ {print}' Homo_sapiens.GRCh38.110.gtf > gencode.v41.annotation_subset.gtf
rm Homo_sapiens.GRCh38.110.gtf
echo "✅ Subset GTF ready at: $REF_DIR/gencode.v41.annotation_subset.gtf"

# ----------------------------
# 4️⃣ Build STAR Index (Subset)
# ----------------------------
echo "🚀 Building small STAR index (chromosomes 21 & 22 only)..."
mkdir -p $REF_DIR/STAR_index_GRCh38_subset
STAR --runThreadN 4 \
     --runMode genomeGenerate \
     --genomeDir $REF_DIR/STAR_index_GRCh38_subset \
     --genomeFastaFiles $REF_DIR/GRCh38_subset.fa \
     --sjdbGTFfile $REF_DIR/gencode.v41.annotation_subset.gtf \
     --sjdbOverhang 100
echo "✅ STAR index built at: $REF_DIR/STAR_index_GRCh38_subset"

# ----------------------------
# 5️⃣ Download Minimal STAR-Fusion Resource
# ----------------------------
echo "⬇️ Downloading lightweight STAR-Fusion CTAT Resource Library..."
wget -c https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v41_CTAT_lib_Mar012022.plug.tar.gz
tar -xzf GRCh38_gencode_v41_CTAT_lib_Mar012022.plug.tar.gz ctat_genome_lib_build_dir/ref_annot.gtf ctat_genome_lib_build_dir/ref_genome.fa
mkdir -p $REF_DIR/ctat_resource_lib_min
mv ctat_genome_lib_build_dir/* $REF_DIR/ctat_resource_lib_min/ || true
rm -rf ctat_genome_lib_build_dir
echo "✅ Minimal CTAT resource ready at: $REF_DIR/ctat_resource_lib_min"

# ----------------------------
# 6️⃣ Download Mini Kraken2 Database
# ----------------------------
echo "⬇️ Downloading small Kraken2 database (human-only)..."
mkdir -p $REF_DIR/kraken2_db_mini
kraken2-build --download-library human --threads 4 --db $REF_DIR/kraken2_db_mini
kraken2-build --build --threads 4 --db $REF_DIR/kraken2_db_mini
echo "✅ Mini Kraken2 DB ready at: $REF_DIR/kraken2_db_mini"

# ----------------------------
# 7️⃣ Export Paths
# ----------------------------
echo "🧩 Adding environment variables to ~/.bashrc..."
cat <<EOF >> ~/.bashrc

# ---- Lightweight RNA-seq Reference Paths ----
export STAR_INDEX="$HOME/refs/STAR_index_GRCh38_subset"
export GENOME_FASTA="$HOME/refs/GRCh38_subset.fa"
export GTF="$HOME/refs/gencode.v41.annotation_subset.gtf"
export STAR_FUSION_CTAT_LIB="$HOME/refs/ctat_resource_lib_min"
export KRAKEN2_DB="$HOME/refs/kraken2_db_mini"
EOF

source ~/.bashrc
echo "✅ Environment variables saved and loaded."

echo "🎉 Lightweight setup complete! Ready to use with STAR, STAR-Fusion, and Kraken2."
