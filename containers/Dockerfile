FROM mambaorg/micromamba:1.5.8

# Install all bioinformatics tools needed for the pipeline
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    fastqc=0.12.1 \
    fastp=0.23.4 \
    star=2.7.11b \
    subread=2.0.6 \
    samtools=1.19 \
    salmon=1.10.1 \
    multiqc=1.21 \
    r-base=4.3.1 \
    bioconductor-deseq2=1.42.0 \
    bioconductor-tximport=1.30.0 \
    r-ggplot2=3.4.4 \
    r-pheatmap=1.0.12 \
    r-rcolorbrewer=1.1_3 \
    r-dplyr=1.1.3 \
    r-readr=2.1.4 \
    && micromamba clean -a -y

# Ensure the conda environment is on the PATH
ENV PATH=/opt/conda/bin:$PATH

# Standard working directory
WORKDIR /data
