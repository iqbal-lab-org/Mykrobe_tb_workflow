FROM continuumio/miniconda3:4.9.2

RUN conda init bash && \
    . /root/.bashrc && \
    conda update conda && \
    conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install mykrobe=0.9.0 minimap2=2.17 samtools=1.11 nanostat=1.5 && \
    pip install pistis==0.3.3

ARG GUPPY_VERSION="4.4.2"
ARG GUPPY_URL="https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_${GUPPY_VERSION}_linux64.tar.gz"

WORKDIR /guppy

RUN wget -O - "$GUPPY_URL" | tar xzf -

ENV PATH="/guppy/ont-guppy-cpu/bin:${PATH}"

