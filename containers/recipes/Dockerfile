FROM ubuntu:20.04

RUN apt update \
    && apt install -y wget make gcc libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev python3-dev python3-pip git \
    && ln -sf $(which python3) /usr/bin/python \
    && ln -sf $(which pip3) /usr/bin/pip

ARG GUPPY_VERSION="6.2.1"
ARG GUPPY_URL="https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_${GUPPY_VERSION}_linux64.tar.gz"
ENV PATH="/guppy/ont-guppy-cpu/bin:${PATH}"
WORKDIR /guppy
RUN wget -O - "$GUPPY_URL" | tar xzf -

ARG SAMTOOLS_VERSION="1.15.1"
ARG SAMTOOLS_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2"
WORKDIR /tmp
RUN wget "$SAMTOOLS_URL" -O /tmp/samtools.tar.bz2 \
    && tar -xf samtools.tar.bz2 \
    && cd samtools-${SAMTOOLS_VERSION} \
    && ./configure --prefix /usr \
    && make \
    && make install

ARG MINIMAP_VERSION="2.24"
ARG MINIMAP_URL="https://github.com/lh3/minimap2/releases/download/v${MINIMAP_VERSION}/minimap2-${MINIMAP_VERSION}_x64-linux.tar.bz2"
WORKDIR /tmp
RUN wget "$MINIMAP_URL" -O /tmp/minimap2.tar.bz2 \
    && tar -jxf minimap2.tar.bz2 \
    && mv minimap2*/* /usr/bin/

ARG SNAKEMAKE_VERSION="7.8.5"
ARG NANOSTAT_VERSION="1.5"
ARG PISTIS_VERSION="0.3.4"
RUN python3 -m pip install --no-cache-dir -U pip \
        && python3 -m pip install --no-cache-dir requests \
            snakemake=="$SNAKEMAKE_VERSION" \
            nanostat=="$NANOSTAT_VERSION" \
            pistis=="$PISTIS_VERSION"

ARG RASUSA_VERSION="0.7.0"
ARG RASUSA_URL="https://github.com/mbhall88/rasusa/releases/download/${RASUSA_VERSION}/rasusa-${RASUSA_VERSION}-x86_64-unknown-linux-musl.tar.gz"
RUN wget "$RASUSA_URL" -O - | tar -xz && mv rasusa*/rasusa /usr/bin/rasusa

ARG MYKROBE_VERSION="0.12.1"
ARG MYKROBE_URL="https://github.com/Mykrobe-tools/mykrobe/archive/v${MYKROBE_VERSION}.tar.gz"
ARG MCCORTEX_GIT="https://github.com/phelimb/mccortex"
WORKDIR /
RUN wget "$MYKROBE_URL" -O mykrobe.tar.gz \
    && tar -xzf mykrobe.tar.gz \
    && rm mykrobe.tar.gz \
    && cd mykrobe* \
    && git clone --recursive -b geno_kmer_count "$MCCORTEX_GIT" mccortex \
    && cd mccortex \
    && make MAXK=31 \
    && cp bin/mccortex31 ../src/mykrobe/cortex/ \
    && cd ../ \
    && rm -rf mccortex \
    && python3 -m pip install --no-cache-dir . \
    && mykrobe panels update_metadata \
    && mykrobe panels update_species all

RUN apt remove -y make git wget \
    && apt autoremove -y \
    && rm -rf /var/lib/apt/lists/* \
    && rm -rf /tmp/*