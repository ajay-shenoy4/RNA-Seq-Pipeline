FROM mambaorg/micromamba:1.5.8

# Install all bioinformatics tools needed for the pipeline
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    fastqc=0.12.1 \
    fastp=0.23.4 \
    star=2.7.11b \
    subread=2.0.6 \
    samtools=1.19 \
    r-base=4.3.1 \
    && micromamba clean -a -y

# Ensure the conda environment is on the PATH
ENV PATH=/opt/conda/bin:$PATH

# Standard working directory
WORKDIR /data
