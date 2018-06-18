========================================
Analysis pipeline for *M. tuberculosis* nanopore data
========================================

.. contents::

Overview
========================================

This pipeline is designed to analyse Oxford Nanopore Technologies sequence data.

The analysis run by the pipeline is:

1. Adapter trimming of the basecalled reads (and demultiplexing if required).
2. Alignment to a given reference genome (default is `NC_000962.3`_).
3. Plots and statistics of the data after the above.
4. A final HTML report summarising all results.

.. image:: ./docs/imgs/dag.png


Installation
========================================
**Note: the following instructions assume you are working on a Linux operating system and have Python version 3.5 or greater.**

The first thing to do is download this repository onto the machine you want to
run the analysis on. In the spirit of making everything reproducible and tidy I
would advise to download this repository once for each nanopore experiment you
want to analyse.

Let's create our experiment directory and clone the pipeline.

.. code-block:: bash

    experiment=sample1
    git clone https://github.com/iqbal-lab-org/Mykrobe_tb_workflow.git "$experiment"

This will download the pipeline repository into a directory named ``sample1``.


Install Singularity
---------------------
Singularity containers can be used to package entire scientific workflows,
software and libraries, and even data. This means that you don’t have to ask
your cluster admin to install anything for you - you can put it in a Singularity
container and run. A Singularity container with all the programs required to run
this analysis is provided with the pipeline, but in order to use Singularity
containers you need to have Singularity installed. If you don't have Singularity
installed, you can find `detailed instructions here`_.


Install Snakemake
---------------------
Snakemake_ is a workflow management system which coordinates the running of this
pipeline. In order to install it you will need to make sure you have Python3_
installed. It is best to manage all of this in a python virtual enviornment. In
this repository there is a ``Pipfile`` which can be used to set up an
environment in ``pipenv``_. All you need to do is run:

.. code-block:: bash

    cd "$experiment"
    project_dir=$(pwd)
    mkdir -p ${project_dir}/logs/cluster
    pipenv install
    pipenv shell

This will install the required python packages ``snakemake`` and ``docutils``
and activate the virtual environment.

Setup
========================================
Getting Singularity containers
--------------------------------
.. image:: https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg
  :target: https://singularity-hub.org/collections/1145

I have added the paths to the appropriate Singularity containers on ``yoda`` into
the ``config.yaml`` file so there is no need to worry about this.

If you are wanting to run this somewhere else though you can get the containers
with by pulling them from Singularity Hub:

.. code-block:: bash

    singularity pull --force --name nanoporeqc.simg shub://mbhall88/Singularity_recipes:nanoporeqc

Make sure to update the ``config.yaml`` with the paths to the containers if you
pulled them down from Singularity Hub.

Moving/copying reads into correct directory
--------------------------------------------
The pipeline expects that the data is placed in specific directories. Whilst this may seem a bit rigid, it is all in the name of reproducibility.

**Non-barcoded sample**

If you have already basecalled your reads then you will only need to merge the fastq files produced by the basecaller. To do this

.. code-block:: bash

    # make the directory we will move the merged file into
    mkdir -p ${project_dir}/data/basecalled
    cd /path/to/basecalled/fastq_files
    cat *.fastq | gzip > ${project_dir}/data/basecalled/${experiment}.fastq.gz
    cd ${project_dir}

This will combine all of the fastq files into a single, compressed file named according to the experiment name and move it into our basecalled data directory.

**Barcoded sample**

If you are working with multiplexed samples (barcoded) then your directory that the basecalling was done into should contain subdirectories named after the barcode they were binned into by the basecaller. You will need to moved these directories (in exampe below) to a directory in the experiment pipeline. If you did not selected the barcoding option for basecalling, but the samples are barcoded, then do the following for the fastq files produced by the basecalling. Note: we only work with files in the "pass" directory (if there is one). Additionally, if you did not basecall the data with the demultiplexing option, then just place

.. code-block:: bash

    # make the directory we will move the reads into
    mkdir -p ${project_dir}/data/basecalled/workspace/pass
    cd ${project_dir}/data/basecalled/workspace/pass
    mv /path/to/dir/containing/barcode/dirs/* .
    cd ${project_dir}

**Basecalling required**

If basecalling is required from the pipeline then you need to do two things. First, change the ``basecall`` field to ``true`` within the config file (see below). Second, move your fast5 files into the pipeline directory.

.. code-block:: bash

    # make the directory we will move the reads into
    mkdir -p ${project_dir}/data/reads
    cd ${project_dir}/data/reads
    mv /path/to/dir/containing/fast5/files/* .
    cd ${project_dir}

If they are multiplexed then you must fill in the appropriate fields in the config file (see below).

Config file - ``config.yaml``
--------------
Open this file up in a text editor and change the following fields, if necessary:

* **multiplexed** - Default is ``false``. Change to ``true`` if sample is multiplexed. If set to ``true`` then you **MUST** enter information for ``barcodes`` as well (see below).
* **sample_name** - If ``multiplexed`` is set to ``false`` then this is the name of your sample. **Note: this MUST be the value of** ``experiment`` **we defined at the start of the installation instructions**. If ``multiplexed`` is set to ``true`` then ignore this field.
* **barcodes** - If ``multiplexed`` is set to ``true`` then this needs to be a **space-separated** string of the expected barcodes (the ones you used in the experiment). An example of barcodes 01-05 is provided. These **MUST** follow the same format of ``BC`` followed by 2 digits. If ``multiplexed`` is set to ``false`` then ignore this field.
* **basecall** - Default is ``true``. Set to ``false`` if you have already basecalled the data.
* **reference** - The genome you would like to align the reads to.
* **flowcell** - The flowcell used (if known). Default is "FLO-MIN106"
* **kit** - The sequencing kit used (if known). Default is "SQK-LSK108"
* **containers** - If you have downloaded/built the Singularity containers elsewhere as you will be using them for multiple samples then change the paths for each container to the location you have them stored at. If running this on ``yoda`` though you shouldn't need to change
this.

Cluster config file - ``cluster.yaml``
--------------------
This file holds the parameters/resources that ``snakemake`` will submit the jobs for each
rule with. The fields are pretty self-explanatory so feel free to change them as
you see fit. The one section in this you **should** change is under ``__defaul__``:``name``
you should name ``JOBNAME`` something useful, such as the current value of
``$experiment``.

**Note:** if you change the memory parameter for a rule, ensure you also change the
value in resources in the two places with that value.

Run
======
You are all set up now. To run the pipeline simply execute the following:

.. code-block:: bash

    cd ${project_dir}
    CLUSTER_CMD='"bsub -n {cluster.nCPUs} -R {cluster.resources} -M {cluster.memory} -o {cluster.output} -e {cluster.error} -J {cluster.name}"'
    bsub.py 1 logs/cluster/snakemake_master_process \
      snakemake \
        --use-singularity \
        --cluster-config cluster.yaml \
        --jobs 500 \
        --cluster "$CLUSTER_CMD"

Or if you don't have access to ``bsub.py``:

.. code-block:: bash

    bsub -R "select[mem>1000] rusage[mem=1000]" -M1000 -o logs/cluster/snakemake_master_process.o -e logs/cluster/snakemake_master_process.e -J snakemake_master_process \
      snakemake \
        --use-singularity \
        --cluster-config cluster.yaml \
        --jobs 500 \
        --cluster "$CLUSTER_CMD"

All the log files for the cluster jobs will be under ``logs/cluster`` and all
the logs for the commands themselves will be in ``logs/``. When it has all run
the data should all be in the appropriate subdirectories in ``data/``.




.. _Singularity: http://singularity.lbl.gov/
.. _`detailed instructions here`: http://singularity.lbl.gov/install-linux
.. _Snakemake: https://snakemake.readthedocs.io/en/stable/index.html
.. _Python3: https://www.python.org/downloads/source/
.. _NC_000962.3: https://www.ncbi.nlm.nih.gov/nuccore/NC_000962.3
.. _pipenv: https://docs.pipenv.org/
