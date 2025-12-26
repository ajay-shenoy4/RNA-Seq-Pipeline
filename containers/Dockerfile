FROM mambaorg/micromamba:1.5.8

# Install bioinformatics tools (without R first)
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    fastqc=0.12.1 \
    fastp=0.23.4 \
    star=2.7.11b \
    subread=2.0.6 \
    samtools=1.19 \
    salmon=1.10.1 \
    multiqc=1.21 \
    && micromamba clean -a -y

# Install R packages separately to avoid conflicts
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    r-base=4.4.2 \
    bioconductor-deseq2 \
    bioconductor-tximport \
    r-ggplot2 \
    r-pheatmap \
    r-rcolorbrewer \
    r-dplyr \
    r-readr \
    && micromamba clean -a -y

ENV PATH=/opt/conda/bin:$PATH
WORKDIR /data
