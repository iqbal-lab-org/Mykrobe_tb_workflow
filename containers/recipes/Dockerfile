FROM continuumio/miniconda3:4.9.2-alpine

RUN conda init bash && \
    . /root/.bashrc && \
    conda update conda && \
    conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install mykrobe=0.9.0 minimap2=2.17 samtools=1.11 nanostat=1.5 && \
    pip install pistis==0.3.3

