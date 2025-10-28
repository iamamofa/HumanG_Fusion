## How to Use

###  Build the Docker image
```bash
docker compose build
```
*(or if you already built it manually, skip this)*

---

###  Run an interactive container
```bash
docker compose run bioinf
```
You’ll drop into:
```bash
root@container:/data#
```
Inside, the conda environment `bioinf` is already active.

---

### Mount data
All files in your local `./data/` directory are shared with `/data` inside the container:
```bash
mkdir data
cp some_reads.fastq data/
docker compose run bioinf fastqc some_reads.fastq
```

---

### Reuse container
To keep it running as a long session (e.g., for Nextflow):
```bash
docker compose up -d
docker exec -it bioinf bash
```

---

### Optional: Nextflow Integration
If you plan to run Nextflow pipelines directly inside, you can use the same container, e.g.:
```bash
docker compose run bioinf nextflow run my_pipeline.nf -profile docker
```
Or include in the compose file:
```yaml
command: ["nextflow", "run", "my_pipeline.nf"]
```

---

###  Directory Structure Example
```
project/
│
├── ./Docker_setup/Dockerfile/
├── docker-compose.yml
└── data/
    ├── sample_R1.fastq.gz
    ├── sample_R2.fastq.gz
    └── output/
```

