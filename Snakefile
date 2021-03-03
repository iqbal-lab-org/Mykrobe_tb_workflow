from pathlib import Path
from snakemake.utils import min_version


# Snakemake version when Singularity support was added
min_version("6.0.0")


# ======================================================
# Config files
# ======================================================
configfile: "config.yaml"


# ======================================================
# Global variables
# ======================================================
RULES_DIR = Path("rules")
IS_MULTIPLEXED = config["multiplexed"]
GB = 1024


include: RULES_DIR / "common.smk"


if IS_MULTIPLEXED:
    SAMPLES = barcode_parser(config["barcodes"])
else:
    SAMPLES = [config["sample_name"]]


# ======================================================
# Rules
# ======================================================
rule all:
    input:
        expand("docs/report_{sample}.html", sample=SAMPLES),


# the snakemake files that run the different parts of the pipeline
if IS_MULTIPLEXED:

    include: RULES_DIR / "demux.smk"


include: RULES_DIR / "align.smk"
include: RULES_DIR / "mykrobe.smk"
include: RULES_DIR / "reports.smk"
