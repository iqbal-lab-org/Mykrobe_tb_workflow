# Analysis pipeline for *M. tuberculosis* nanopore data

[TOC]: #

## Table of Contents
- [Install Singularity](#install-singularity)
- [Download Pipeline](#download-pipeline)
- [Install Snakemake](#install-snakemake)
- [Singularity containers](#singularity-containers)
  - [Download from Singularity Hub (recommended)](#download-from-singularity-hub-recommended)
  - [Build container locally](#build-container-locally)
- [Initial data location](#initial-data-location)
  - [Non-barcoded sample {#non_barcoded_sample}](#non-barcoded-sample-non_barcoded_sample)
  - [Barcoded sample](#barcoded-sample)
- [Configuration file](#configuration-file)
- [Cluster configuration file](#cluster-configuration-file)
- [Local](#local)
- [Cluster](#cluster)

# Overview

This pipeline is designed to analyse Oxford Nanopore Technologies sequence data. It was
developed (and will continue to be improved/maintained) as part of a project to improve
tuberculosis diagnostics using nanopore sequencing. The project is being run by a global
group of researchers and clinicians from Madagascar's National TB Program, Institute
Pasteur Madagascar (IPM), University of Oxford, European Bioinformatics Institute
(EMBL-EBI), Stony Brook University, and Hospital for Tropical Diseases - Ho Chi Minh
City. More information and a short video of the project can be found
[here](https://nanoporetech.com/about-us/news/public-health-teams-madagascar-pioneer-use-portable-real-time-dna-sequencing-fight).

The analysis run by the pipeline is:

1. Adapter trimming of the basecalled reads (and demultiplexing if required).
2. Alignment to a given reference genome (default is
   [NC\_000962.3](https://www.ncbi.nlm.nih.gov/nuccore/NC_000962.3)).
3. Plots and statistics of the data after the above.
4. A final HTML report summarising all results. [EXAMPLE
   REPORT](https://rawgit.com/iqbal-lab-org/Mykrobe_tb_workflow/master/docs/example_report.html)

![image](./docs/imgs/dag.png)

Installation
============

**Note: the following instructions assume you are working on a Linux operating system
and have Python version 3.5 or greater.**

Install Singularity
-------------------

Singularity containers can be used to package entire scientific workflows, software and
libraries, and even data. This means that you don't have to ask your cluster admin to
install anything for you - you can put it in a Singularity container and run. A
Singularity container with all the programs required to run this analysis is provided
with the pipeline, but in order to use Singularity containers you need to have
Singularity installed. If you don\'t have Singularity installed, you can find
[detailed instructions here](http://singularity.lbl.gov/install-linux).

Download Pipeline
-----------------

The first thing to do is download this repository onto the machine you want to run the
analysis on. In the spirit of making everything reproducible and tidy I would advise to
download this repository once for each nanopore experiment you want to analyse.

Let\'s create our experiment directory and clone the pipeline.

```shell
experiment=sample1
git clone https://github.com/iqbal-lab-org/Mykrobe_tb_workflow.git "$experiment"
cd "$experiment"
project_dir=$(pwd)
mkdir -p "$project_dir"/logs
```

This will download the pipeline repository into a directory named `sample1`.

Install Snakemake
-----------------

[Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) is a workflow
management system which coordinates the running of this pipeline. In order to install it
you will need to make sure you have [Python3](https://www.python.org/downloads/source/)
installed. It is best to manage all of this in a python virtual enviornment. In this
repository there is a `requirements.txt` file which can be used to set up a python
environment very easily.

```shell
cd "$project_dir"
# create a virtual environment
python3 -m venv .venv
# activate the environment
source .venv/bin/activate
# install the necessary packages
python3 -m pip install -r requirements.txt
# sometimes environment variables get cleared after activating environments
project_dir=$(pwd)
experiment=$(basename $project_dir)
```

This will install the required python package `snakemake` and activate
the virtual environment. You will need to remember to activate and deactivate the
environment each time

```shell
# activate
cd "$project_dir"
source .venv/bin/activate
# deactivate
deactivate
```

Analysis setup
==============

Singularity container
---------------------

Download the container locally

```shell
cd "$project_dir"/containers
container_name=tb.simg
singularity pull --force --name "$container_name" "docker://quay.io/iqballab/mykrobe_tb_workflow:latest"
```

If you are going to be running this pipline for many different samples on the same
machine, it is recommended to only download/build the container once, as it is about
1.5GB. Change `container_name` in the above code to a more central directory and make
sure to update the container location in `config.yaml` (see below).

Initial data location
---------------------

The pipeline expects that the data you want to analyse is placed in specific
directories. Whilst this may seem a bit rigid, it is all in the name of reproducibility.

### Non-barcoded sample
For a single sample with no barcoding (and therefore no demultiplexing required) you
just need to ensure there is a single fastq file of the basecalled reads. Generally,
when a sample has been basecalled there is multiple fastq files (the default for
instance has 4000 reads per fastq). To combine the fastq files into a single file

```shell
# change into the directory where all the fastq files are
cd /path/to/basecalled/fastq_files
cat *.fastq* | gzip > "$experiment".fastq.gz
```

Once you have this single, combined fastq file, we need to move it into the appropriate
pipeline data folder. **Note:** The combined file must have the same name as the
variable `experiment` we set earlier. It must also be `gzip`ed.

```shell
# make the directory we will move the combined file into
mkdir -p "$project_dir"/data
mv "$experiment".fastq.gz "$project_dir"/data/
cd "$project_dir"
```

### Barcoded sample

If you are working with multiplexed (barcoded) samples, then you just need to move (or
copy) the directory containing the basecalled data to the experiment directory

```shell
# make the directory we will move the reads into
mkdir -p "$project_dir"/data/basecalled
# use `cp -r` instead of `mv` if you want to copy the folder instead
mv /path/to/basecalling/dir "$project_dir"/data/basecalled/
cd "$project_dir"
```

Copying the data will double the size of the data, so moving it is recommended.

Configuration file
------------------

This is the file `config.yaml` located in the pipeline root directory.

Open this file up in a text editor and change the following fields, if necessary:

- **multiplexed** - Default is `false`. Change to `true` if sample is multiplexed. If
  set to `true` then you **MUST** enter information for `barcodes` as well (see below).
- **sample\_name** - If `multiplexed` is set to `false` then this is the name of your
  sample. **Note: this MUST be the value of** `experiment` **we defined at the start of
  the installation instructions**. If `multiplexed` is set to `true` then ignore this
  field.
- **barcodes** - If `multiplexed` is set to `true` then this needs to be a
  **space-separated** string of the expected barcodes (the ones you used in the
  experiment). An example of barcodes 01-05 is provided. These **MUST** follow the same
  format of `BC` followed by 2 digits (e.g `"BC01 BC02 BC03"`). If `multiplexed` is set
  to `false` then ignore this field.
- **reference** - The genome you would like to align the reads to. This is set by to
  default to the reference provided with the pipeline -
  [NC\_000962.3](https://www.ncbi.nlm.nih.gov/nuccore/NC_000962.3).
- **container** - If you have downloaded/built the Singularity container in a different
  location to the default (`containers/tb.simg`) then change the path for the container
  to the location you have it stored at.

Cluster configuration file
--------------------------

This is the file `cluster.yaml` located in the pipeline root directory. It holds the
settings for running the pipeline on a cluster, but also for the resource allocation for
jobs. So even if you\'re running the pipeline on a local computer the resources from
this file will be used.

The fields are pretty self-explanatory so feel free to change them you see fit. The one
section in this you **should** change is under `__defaul__`:`name` -name `JOBNAME`
something useful, such as the current value of `$experiment`.

The cluster configuration provided is also used by snakemake if it is to be run on a
cluster and is how it knows what resources to ask for for each job. This file has been
tested successfully on the LSF cluster management system. For more information on using
snakemake on other cluster management systems, [see the
documentation](https://snakemake.readthedocs.io/en/latest/snakefiles/configuration.html#cluster-configuration).
**Note:** if you change the memory parameter for a rule, ensure you also change the
value in resources in the two places with that value.

Run
===

You are all set up now. To run the pipeline simply execute the following. At the end,
all of the logs will be under `logs/`. Data will be in the appropriate subdirectories in
`data/` and the final report(s) (one for each barcode) will be under `docs/`.

Local
-----

To run the pipeline on a local computer (i.e laptop or desktop)

```{.bash}
cd "$project_dir"
snakemake --use-singularity
```

This will provide a summary of all the jobs that are to be run, and when they have been
started and finished.

Cluster
-------

This pipeline can also be run on a cluster. These instructions are for running on an LSF
cluster system. The `cluster.yaml` file *should* be general across clusters (except for
the `resources` field). The cluster submission command however is different from cluster
to cluster. We provide the command for an LSF system here. Please contact us if you use
a different cluster system and cannot figure out the command and we will see if we can
help. Additionally, if you use a different cluster management system and successfully
run it, please provide the cluster submission commands and we will add them into these
instructions for others to use.

There is script provided in the scripts directory for submitting the job to an LSF
cluster. To run this you just need to be in the pipeline root directory and provide a
name for the job (to be used by the cluster).

```{.bash}
cd "$project_dir"
JOB_NAME=snakemake_master_process
bash scripts/submit_lsf.sh "$JOB_NAME"
```

All the log files for the cluster jobs will be prefixed with `cluster_`.

