FROM continuumio/miniconda3:latest

LABEL maintainer="RiboLongShort"
LABEL description="Container for RiboLongShort Ribo-seq pipeline"

SHELL ["/bin/bash", "-c"]

WORKDIR /opt/ribolongshort

COPY environment.yml .

RUN conda env create -f environment.yml

ENV PATH /opt/conda/envs/ribolongshort/bin:$PATH

RUN echo "source activate ribolongshort" >> ~/.bashrc

WORKDIR /workspace

CMD ["/bin/bash"]
