FROM continuumio/miniconda3:latest

LABEL maintainer="RiboLongShort"
LABEL description="Docker image for Ribo-seq preprocessing and smORF discovery"

SHELL ["/bin/bash", "-c"]

RUN apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    unzip \
    gzip \
    pigz \
    build-essential \
    default-jre \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN conda install -y -c conda-forge -c bioconda \
    sra-tools \
    fastqc \
    multiqc \
    cutadapt \
    bowtie \
    star \
    samtools \
    subread \
    seqkit \
    ribotricer \
    && conda clean -afy

WORKDIR /workspace

CMD ["/bin/bash"]
