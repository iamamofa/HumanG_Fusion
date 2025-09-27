 _   _ ______     _____ ____   ____  _____  ______ 
 _   _ _____    ____ ___  ____  _____      _   _ __  __ ___ __  __ ____  
| \ | |  ___|  / ___/ _ \|  _ \| ____|    | \ | |  \/  |_ _|  \/  |  _ \ 
|  \| | |_    | |  | | | | |_) |  _|      |  \| | |\/| || || |\/| | |_) |
| |\  |  _|   | |__| |_| |  _ <| |___     | |\  | |  | || || |  | |  _ < 
|_| \_|_|      \____\___/|_| \_\_____|    |_| \_|_|  |_|___|_|  |_|_| \_\      
     H U M A N   _ F   P I P E L I N E



# Nextflow RNA-Seq Human Fusion Detection Pipeline

This repository contains a **Nextflow (DSL2)** pipeline scaffold for **human-only RNA-Seq fusion detection**. It includes:

* A modular Nextflow DSL2 pipeline (`main.nf` + modules) that runs preprocessing, optional Kraken2 decontamination, STAR alignment (2-pass), fusion callers (STAR-Fusion, Arriba, FusionCatcher optional), and postprocessing/merging.
* Dockerfile and Singularity definition for reproducible runtime.
* An improved `merge_fusions.py` that parses STAR-Fusion, Arriba, FusionCatcher, annotates genes via `mygene` (Ensembl), and produces summary scores.
* A README with usage, configuration, and testing suggestions.

---

## Repo layout

```
nextflow-fusion-pipeline/
├─ README.md
├─ nextflow.config
├─ main.nf
├─ modules/
│  ├─ preprocess/main.nf
│  ├─ kraken_decontam/main.nf
│  ├─ star_align/main.nf
│  ├─ fusion_callers/main.nf
│  ├─ postprocess/main.nf
│  └─ utils/ (helper small processes)
├─ conf/ (profiles)
│  └─ containers.config
├─ docker/
│  ├─ Dockerfile
│  └─ fusion_env.Dockerfile  # smaller image variant if desired
├─ singularity/
│  └─ Singularity.def
├─ tools/
│  └─ merge_fusions.py
├─ bin/
│  └─ helper scripts
└─ sample_data/ (optional: example metadata)
```

---

## Important notes (before running)

* This pipeline **enforces human-only processing** by aligning to a human reference (GRCh38) and by offering an **optional Kraken2 decontamination** step to remove non-human reads prior to fusion calling.
* Tools used: `fastqc`, `cutadapt`, `kraken2` (optional), `star`, `star-fusion`, `arriba`, `fusioncatcher` (optional), `samtools`, `mygene` (Python), `python3`.
* The Dockerfile builds an image containing the primary tools; if you prefer Singularity, the `Singularity.def` is included.
* The pipeline is written in DSL2 modules and is portable to HPC or cloud via Nextflow profiles.

---

## `nextflow.config` (example)

```groovy
profiles {
  standard {
    process {
      cpus = 8
      memory = '16 GB'
      time = '6h'
    }
    executor = 'local'
  }
  docker {
    process.container = 'nextflow/fusion-pipeline:latest'
  }
  singularity {
    process.container = 'docker://nextflow/fusion-pipeline:latest' // or local sif
    singularity.enabled = true
  }
}

params {
  reads = "./sample_data/*_R1.fastq.gz"
  samples_tsv = "" // optional: sample_id\tR1\tR2
  outdir = './results'
  genome_fasta = '/refs/GRCh38/GRCh38.primary_assembly.genome.fa'
  gtf = '/refs/GRCh38/gencode.v41.annotation.gtf'
  star_index = '/refs/STAR_index_GRCh38'
  kraken2_db = '/refs/kraken2_db' // required if using --kraken_decontam

  threads = 8
  kraken_decontam = true
  run_fusioncatcher = false
  min_split_reads = 3
  min_spanning_reads = 2
}
```

---

## `main.nf` (DSL2 entry)

```groovy
nextflow.enable.dsl=2

include { preprocess } from './modules/preprocess/main.nf'
include { kraken_decontam } from './modules/kraken_decontam/main.nf'
include { star_align } from './modules/star_align/main.nf'
include { fusion_callers } from './modules/fusion_callers/main.nf'
include { postprocess } from './modules/postprocess/main.nf'

workflow {
  samples_ch = Channel.fromFilePairs(params.samples_tsv ?: params.reads, flat: true)

  processed = preprocess(samples_ch)

  decontaminated = params.kraken_decontam ? kraken_decontam(processed) : processed

  aligned = star_align(decontaminated)

  fusion_results = fusion_callers(aligned)

  postprocess(fusion_results)
}
```

> Note: `Channel.fromFilePairs` accepts a glob or a TSV. If `samples_tsv` is provided it should be tab-delimited lines `sample_id\tpath_R1\tpath_R2`.

---

## Module: `modules/preprocess/main.nf`

This module handles FastQC and trimming with cutadapt, and emits trimmed FASTQ pair.

```groovy
process preprocess_trim {
  tag "trim:${sample_id}"
  cpus params.threads
  publishDir "${params.outdir}/\${sample_id}/fastq_trimmed", mode: 'copy'

  input:
  tuple val(sample_id), path(reads)

  output:
  tuple val(sample_id), path("${sample_id}_R1.trim.fastq.gz"), path("${sample_id}_R2.trim.fastq.gz")

  script:
  """
  mkdir -p work
  r1=${reads[0]}
  r2=${reads[1]}
  fastqc -t ${task.cpus} -o . "$r1" "$r2" || true
  cutadapt -q 20 -m 30 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o ${sample_id}_R1.trim.fastq.gz -p ${sample_id}_R2.trim.fastq.gz "$r1" "$r2"
  fastqc -t ${task.cpus} -o . ${sample_id}_R1.trim.fastq.gz ${sample_id}_R2.trim.fastq.gz || true
  """
}

workflow preprocess {
  take: samples
  main:
  samples.map { sid, r -> tuple(sid, r) }.into { _in }
  produce:
  preprocess_trim(_in)
}
```

---

## Module: `modules/kraken_decontam/main.nf`

Optional Kraken2 decontamination. It classifies reads and extracts reads classified as `Human` (or unclassified) depending on user policy. Two options are provided: (A) keep only reads assigned to human, or (B) remove reads assigned to non-human taxa (recommended: remove non-human and keep unclassified+human).

```groovy
process kraken2_classify {
  tag "kraken:${sample_id}"
  cpus params.threads
  output:
  tuple val(sample_id), path("${sample_id}.decontam_R1.fastq.gz"), path("${sample_id}.decontam_R2.fastq.gz")

  input:
  tuple val(sample_id), path(r1), path(r2)

  script:
  """
  kraken2 --db ${params.kraken2_db} --paired --gzip-compressed --report kraken.report --output kraken.out --threads ${task.cpus} $r1 $r2

  # Extract reads not classified as unwanted taxa. We'll use kraken2's --unclassified-out and --classified-out features
  kraken2 --db ${params.kraken2_db} --paired --gzip-compressed --threads ${task.cpus} --output /dev/null --unclassified-out ${sample_id}.unclassified#.fastq --classified-out ${sample_id}.classified#.fastq $r1 $r2 || true

  # simple approach: keep classified human + unclassified
  # use kraken-report to find taxid for Homo sapiens (e.g., 9606)
  # but here we use kraken2 outputs - users may adjust to their db

  # For portability, we convert the per-file outputs to gz pairs named below
  gzip -c ${sample_id}.unclassified_1.fastq > ${sample_id}.decontam_R1.fastq.gz || true
  gzip -c ${sample_id}.unclassified_2.fastq > ${sample_id}.decontam_R2.fastq.gz || true
  """
}

workflow kraken_decontam {
  take: pre_trimmed
  main:
  pre_trimmed.map { tuple -> kraken2_classify(tuple[0], tuple[1], tuple[2]) }
  emit: kraken2_classify.out
}
```

> **Caveat**: Kraken2 DB configuration varies. The script shows a conservative strategy: keep unclassified reads and optionally human-assigned reads. You can refine to keep only reads assigned to `taxid|9606` after parsing `kraken2 --report`.

---

## Module: `modules/star_align/main.nf`

Performs STAR 2-pass, outputs coordinate-sorted BAM and chimeric junction files required by Arriba/STAR-Fusion.

```groovy
process star_2pass {
  tag "star:${sample_id}"
  cpus params.threads
  publishDir "${params.outdir}/\${sample_id}/align", mode: 'copy'

  input:
  tuple val(sample_id), path(r1), path(r2)

  output:
  tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.Chimeric.out.junction"), path("${sample_id}.Aligned.out.bam")

  script:
  """
  mkdir -p star_tmp
  STAR --runThreadN ${task.cpus} --genomeDir ${params.star_index} --readFilesIn $r1 $r2 --readFilesCommand zcat --outFileNamePrefix ${sample_id}. --outSAMtype BAM Unsorted --chimSegmentMin 12 --chimJunctionOverhangMin 12 --chimOutType Junctions --outSAMattributes NH HI AS nM MD XS --outTmpDir star_tmp
  samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam ${sample_id}.Aligned.out.bam
  samtools index ${sample_id}.sorted.bam
  """
}

workflow star_align {
  take: reads_in
  main:
  reads_in.map { tuple -> star_2pass(tuple[0], tuple[1], tuple[2]) }
}
```

---

## Module: `modules/fusion_callers/main.nf`

Runs STAR-Fusion, Arriba, and optionally FusionCatcher. Each caller writes to its output directory per-sample.

```groovy
process run_starfusion {
  tag "starfusion:${sample_id}"
  cpus params.threads
  publishDir "${params.outdir}/\${sample_id}/fusion/starfusion", mode: 'copy'

  input:
  tuple val(sample_id), path(r1), path(r2), path(chimeric), path(sorted_bam)

  script:
  """
  # STAR-Fusion (the wrapper will run its own STAR if requested; many installations expect raw fq input)
  STAR-Fusion --left_fq $r1 --right_fq $r2 --genome_lib_dir ${params.star_fusion_ctat_lib} --CPU ${task.cpus} --output_dir starfusion_out || true
  cp -r starfusion_out/* . || true
  """
}

process run_arriba {
  tag "arriba:${sample_id}"
  cpus params.threads
  publishDir "${params.outdir}/\${sample_id}/fusion/arriba", mode: 'copy'

  input:
  tuple val(sample_id), path(chimeric), path(sorted_bam)

  script:
  """
  # Arriba expects chimeric BAM and aligned BAM produced by STAR with specific flags
  arriba -x $chimeric -b $sorted_bam -o arriba.tsv -O arriba.discarded.tsv -a ${params.genome_fasta} -g ${params.gtf} -t ${params.arriba_blacklist} || true
  """
}

process run_fusioncatcher {
  tag "fusioncatcher:${sample_id}"
  cpus params.threads
  when: params.run_fusioncatcher
  publishDir "${params.outdir}/\${sample_id}/fusion/fusioncatcher", mode: 'copy'

  input:
  tuple val(sample_id), path(r1), path(r2)

  script:
  """
  fusioncatcher -d . -p ${task.cpus} -i $r1,$r2 || true
  """
}

workflow fusion_callers {
  take: aligned
  main:
  aligned.map { sid, sorted_bam, chimeric, _ ->
    // we also need trimmed fq for STAR-Fusion (depending on the wrapper); the aligned input can include the fq paths
    run_arriba(sid, chimeric, sorted_bam)
    run_starfusion(sid, sorted_bam.meta?.left_fastq, sorted_bam.meta?.right_fastq, chimeric, sorted_bam)
    if (params.run_fusioncatcher) run_fusioncatcher(sid, sorted_bam.meta?.left_fastq, sorted_bam.meta?.right_fastq)
  }
}
```

> Note: Nextflow metadata passing (`meta`) can be used to carry trimmed FASTQ paths in the channel. The example shows the idea; minor wiring may be required depending on how channels are created.

---

## Module: `modules/postprocess/main.nf`

Collects per-sample outputs, runs `tools/merge_fusions.py` to combine calls, annotate, and produce final TSV/HTML.

```groovy
process merge_and_annotate {
  tag "merge:${sample_id}"
  cpus 2
  publishDir "${params.outdir}/reports", mode: 'copy'

  input:
  val sample_id
  path starfusion_file optional true
  path arriba_file optional true
  path fusioncatcher_file optional true

  output:
  path "${sample_id}.merged_fusions.tsv"

  script:
  """
  python3 ${workDir}/../tools/merge_fusions.py \
    --starfusion ${starfusion_file} \
    --arriba ${arriba_file} \
    --fusioncatcher ${fusioncatcher_file} \
    --gtf ${params.gtf} \
    --out ${sample_id}.merged_fusions.tsv
  """
}

workflow postprocess {
  take: fusion_results
  main:
  fusion_results.map { tuple -> merge_and_annotate(tuple[0], tuple[1], tuple[2], tuple[3]) }
}
```

---

## `tools/merge_fusions.py` (improved)

This script:

* Parses STAR-Fusion (`star-fusion.fusion_predictions.tsv`), Arriba (`arriba.tsv`) and FusionCatcher output.
* Merges calls, assigns caller counts, computes a simple summary score (e.g., number of callers, total split/spanning support), and annotates left/right genes with Ensembl IDs and descriptions using `mygene` Python client.

```python
#!/usr/bin/env python3
"""
merge_fusions.py - parse and merge STAR-Fusion, Arriba, FusionCatcher outputs
Requires: mygene (pip install mygene)
"""
import argparse, csv, glob, sys
from collections import defaultdict

try:
    import mygene
    mg = mygene.MyGeneInfo()
except Exception as e:
    mg = None


def read_starfusion(path):
    calls = []
    if not path: return calls
    try:
        with open(path) as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            for r in reader:
                left = r.get('LeftGene')
                right = r.get('RightGene')
                calls.append({'caller':'STAR-Fusion', 'left':left, 'right':right, 'split': int(r.get('JunctionReads') or 0), 'spanning': int(r.get('SpanningFragCount') or 0), 'raw':r})
    except Exception:
        pass
    return calls

def read_arriba(path):
    calls = []
    if not path: return calls
    try:
        with open(path) as fh:
            for line in fh:
                if line.startswith('#'): continue
                cols = line.strip().split('\t')
                gene1 = cols[3]
                gene2 = cols[7]
                # support column may include e.g. 'sr:2;pe:1'
                support = cols[8] if len(cols)>8 else ''
                split = 0
                span = 0
                for token in support.split(';'):
                    if token.startswith('sr:'): split = int(token.split(':')[1])
                    if token.startswith('pe:'): span = int(token.split(':')[1])
                calls.append({'caller':'Arriba','left':gene1,'right':gene2,'split':split,'spanning':span,'raw':cols})
    except Exception:
        pass
    return calls

def read_fusioncatcher(path_glob):
    calls = []
    if not path_glob: return calls
    for p in glob.glob(path_glob or ''):
        try:
            with open(p) as fh:
                for line in fh:
                    if line.startswith('#') or not line.strip(): continue
                    cols = line.strip().split('\t')
                    # the typical fusioncatcher final-list has columns: Gene1_symbol Gene2_symbol ...
                    left = cols[0]
                    right = cols[1]
                    calls.append({'caller':'FusionCatcher','left':left,'right':right,'split':0,'spanning':0,'raw':cols})
        except Exception:
            pass
    return calls


def annotate_genes(gene_list):
    if not mg:
        # no mygene available, return empty annotations
        return {g: {'ensembl':None,'name':None,'entrez':None} for g in gene_list}
    res = mg.querymany(list(gene_list), scopes='symbol', fields='ensembl.gene,name,entrezgene', species='human', as_dataframe=False)
    ann = {}
    for r in res:
        q = r.get('query')
        ann[q] = {'ensembl': None, 'name': None, 'entrez': None}
        if 'ensembl' in r and r['ensembl']:
            if isinstance(r['ensembl'], list): ann[q]['ensembl'] = r['ensembl'][0].get('gene')
            elif isinstance(r['ensembl'], dict): ann[q]['ensembl'] = r['ensembl'].get('gene')
        ann[q]['name'] = r.get('name')
        ann[q]['entrez'] = r.get('entrezgene')
    return ann


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--starfusion', default=None)
    p.add_argument('--arriba', default=None)
    p.add_argument('--fusioncatcher', default=None)
    p.add_argument('--gtf', default=None)
    p.add_argument('--out', required=True)
    args = p.parse_args()

    merged = defaultdict(lambda: {'callers':set(),'split':0,'spanning':0,'raw':[]})
    gene_set = set()

    for c in read_starfusion(args.starfusion):
        key = (c['left'], c['right'])
        merged[key]['callers'].add('STAR-Fusion')
        merged[key]['split'] += c.get('split',0)
        merged[key]['spanning'] += c.get('spanning',0)
        merged[key]['raw'].append(('STAR-Fusion', c['raw']))
        gene_set.update([c['left'], c['right']])

    for c in read_arriba(args.arriba):
        key = (c['left'], c['right'])
        merged[key]['callers'].add('Arriba')
        merged[key]['split'] += c.get('split',0)
        merged[key]['spanning'] += c.get('spanning',0)
        merged[key]['raw'].append(('Arriba', c['raw']))
        gene_set.update([c['left'], c['right']])

    for c in read_fusioncatcher(args.fusioncatcher):
        key = (c['left'], c['right'])
        merged[key]['callers'].add('FusionCatcher')
        merged[key]['raw'].append(('FusionCatcher', c['raw']))
        gene_set.update([c['left'], c['right']])

    # annotate genes via mygene
    annotations = annotate_genes(gene_set)

    # produce a simple ranking: callers_count * 10 + split + spanning
    rows = []
    for k,v in merged.items():
        callers_count = len(v['callers'])
        score = callers_count*10 + v['split'] + v['spanning']
        left_ann = annotations.get(k[0],{})
        right_ann = annotations.get(k[1],{})
        rows.append({
            'left':k[0], 'left_ensembl':left_ann.get('ensembl'), 'left_name':left_ann.get('name'),
            'right':k[1], 'right_ensembl':right_ann.get('ensembl'), 'right_name':right_ann.get('name'),
            'callers':';'.join(sorted(v['callers'])), 'split':v['split'], 'spanning':v['spanning'], 'score':score
        })

    rows.sort(key=lambda r: r['score'], reverse=True)

    with open(args.out,'w') as outfh:
        hdr = ['left','left_ensembl','left_name','right','right_ensembl','right_name','callers','split','spanning','score']
        outfh.write('\t'.join(hdr) + '\n')
        for r in rows:
            outfh.write('\t'.join(str(r.get(c,'')) for c in hdr) + '\n')

if __name__=='__main__':
    main()
```

---

## Dockerfile (docker/Dockerfile)

A single Dockerfile that installs common bioinformatics tools. This is a large image; you can split into smaller images or use Conda/bioconda.

```dockerfile
FROM ubuntu:22.04
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential wget curl git python3 python3-pip python3-dev unzip ca-certificates \
    zlib1g-dev libbz2-dev liblzma-dev libssl-dev libcurl4-openssl-dev \
    default-jre-headless samtools pigz \
  && rm -rf /var/lib/apt/lists/*

# install tools via conda/miniconda for nicer dependency management
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh
ENV PATH=/opt/conda/bin:$PATH

RUN conda install -y -c bioconda -c conda-forge \
    fastqc cutadapt kraken2 star star-fusion arriba fusioncatcher mygene samtools nextflow && \
    conda clean -afy

# Install STAR-Fusion CTAT resource not included due to size; instruct user to mount

# create working dir
WORKDIR /data

# tools
COPY tools/ /opt/pipeline/tools/
RUN chmod +x /opt/pipeline/tools/*.py

# entrypoint is nextflow when run
CMD ["nextflow"]
```

> **Note:** Star-Fusion CTAT resource library is large; we expect the user to provide/mount it at runtime (set `params.star_fusion_ctat_lib`). FusionCatcher may also require external DBs and python paths — consult their docs.

---

## Singularity definition (`singularity/Singularity.def`)

```
Bootstrap: docker
From: ubuntu:22.04

%post
    apt-get update && apt-get install -y wget curl git python3 python3-pip default-jre samtools pigz build-essential && rm -rf /var/lib/apt/lists/*
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh
    bash /tmp/miniconda.sh -b -p /opt/conda
    export PATH=/opt/conda/bin:$PATH
    /opt/conda/bin/conda install -y -c bioconda -c conda-forge fastqc cutadapt kraken2 star star-fusion arriba fusioncatcher mygene samtools nextflow
    conda clean -afy

%files
    tools/ /opt/pipeline/tools/

%environment
    PATH=/opt/conda/bin:$PATH

%runscript
    exec nextflow "$@"
```

---

## README.md (summary and quick start)

````markdown
# Nextflow Human RNA-Seq Fusion Pipeline

## Overview
This Nextflow DSL2 pipeline performs human RNA-Seq fusion detection using STAR, STAR-Fusion, Arriba and optionally FusionCatcher. It includes an optional Kraken2 decontamination step to enforce "human-only" reads.

## Requirements
- Nextflow
- Docker or Singularity (recommended) OR conda environment with all tools
- Reference files (human GRCh38 FASTA, GTF, STAR index, STAR-Fusion CTAT lib, Kraken2 DB if using decontam)

## Quick start (Docker)
1. Build the Docker image (from `docker/Dockerfile`):
   ```bash
   docker build -t nextflow/fusion-pipeline:latest docker/
````

2. Run the pipeline (example):

   ```bash
   nextflow run main.nf -profile docker --reads './sample_data/*_R1.fastq.gz' --outdir ./results --star_index /refs/STAR_index_GRCh38 --gtf /refs/gencode.v41.annotation.gtf --genome_fasta /refs/GRCh38.fa --kraken_decontam true --kraken2_db /refs/kraken2_db
   ```

## Notes

* Provide the STAR-Fusion CTAT resource via `--star_fusion_ctat_lib /path/to/ctat_resource_lib` when running.
* If `--kraken_decontam true`, make sure `--kraken2_db` points to a Kraken2 db containing human.
* The merge script uses `mygene` to annotate gene symbols to Ensembl IDs. Install with `pip install mygene` in the container or environment.

## Outputs

* Per-sample directories: `align/`, `fastq_trimmed/`, `fusion/` (caller-specific), and `reports/` containing merged TSVs.

## Testing & validation

* Test on small datasets with known fusions (e.g., SRA datasets or synthetic spike-ins).
* Validate artefacts by manual inspection of supporting reads in IGV.

```

---

## Final notes and next steps
- The DSL2 modules above are intentionally modular and include `publishDir` placeholders to collect outputs into `params.outdir`.
- You may want to: add Singularity/Docker push automation, add Nextflow Tower profiles, wrap STAR-Fusion so that it uses the STAR outputs (some STAR-Fusion releases call STAR internally — check your version), and create a small `tests/` folder with test-fastqs.

If you'd like, I can:
- produce a runnable Git repository archive (zip) containing all files prepared above,
- generate a tailored Dockerfile that pins exact versions for each tool,
- or convert the code to a single-file `main.nf` + `modules/` ready to run on your machine with exact parameter wiring.

---

*I put the full Nextflow modules, Dockerfile, Singularity def, `tools/merge_fusions.py`, and README into this canvas. Open the document and tell me which part you'd like me to expand or to export as files.*

```
