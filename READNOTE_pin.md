Notes, recommendations, and caveats

Human-only enforcement: This pipeline aligns reads to a human reference (GRCh38) and uses human gene annotation (GTF) in fusion callers. If you want to exclude reads mapping to non-human contaminants before fusion calling, add a decontamination step (e.g., Kraken2/BBMap) — I left that out for clarity but can add it.

STAR-Fusion specifics: STAR-Fusion has particular expectations for the STAR output and a CTAT resource lib; ensure you follow STAR-Fusion documentation when installing. The config variable STAR_FUSION_CTAT_LIB should point to the CTAT library for GRCh38.

Arriba expects STAR outputs created with specific flags (--chimOutType WithinBAM etc.). You may need to tweak the STAR command to produce compatible outputs.

FusionCatcher is optional; it uses its own database and can be slow but complementary.

Containerization: I strongly recommend using Docker/Singularity or conda environments to avoid dependency hell. Provide a containers/ folder with Dockerfiles if you want; I can generate them if you want.

Parallelization: The driver runs samples sequentially; to parallelize, either run multiple instances for different samples or convert to Nextflow (straightforward from this modular layout).

Testing: Test the pipeline on small sample(s) first. Use public RNA-Seq datasets with known fusions (e.g., positive controls) to validate.

Security / Reproducibility: Pin tool versions in your conda env or container. Save logs to logs/ (config already does that).
