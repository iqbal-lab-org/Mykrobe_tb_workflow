import os
import re
import yaml
from typing import List
from snakemake.utils import min_version

RULES_DIR = 'rules'
# need this version as minimum as it is when singularity support was added
min_version("4.2.0")


configfile: "config.yaml"

with open('cluster.yaml', 'r') as fh:
    cluster_config = yaml.load(fh)


class InvalidBarcode(Exception):
    __module__ = Exception.__module__


def barcode_parser(barcodes_string: str) -> List[str]:
    """Parses the barcodes string and ensures they follow correct format"""
    msg = "Barcode must be of the form BC01. That is, BC followed by 2 digits."
    regex = r'\bBC\d{2}\b'
    barcodes = barcodes_string.split()
    for barcode in barcodes:
        if not (len(barcode) == 4 and re.match(regex, barcode)):
            raise InvalidBarcode(barcode + '\n' + msg)
    return barcodes


# if multiplexed, add expected barcodes
MULTIPLEXED = config["multiplexed"]


if MULTIPLEXED:
    SAMPLES = barcode_parser(config["barcodes"])
else:
    SAMPLES = [config["sample_name"]]


rule all:
    input:
        expand("data/sorted/{sample}_sorted.bam.bai", sample=SAMPLES),
        expand("data/stats/{sample}_stats.txt", sample=SAMPLES),
        expand("data/plots/{sample}_plots.pdf", sample=SAMPLES)


# only run basecalling when requested
if config["basecall"]:
    include: os.path.join(RULES_DIR, 'basecall.smk')


# the snakemake files that run the different parts of the pipeline
include: os.path.join(RULES_DIR, 'porechop.smk')
include: os.path.join(RULES_DIR, 'align.smk')
include: os.path.join(RULES_DIR, 'reports.smk')
