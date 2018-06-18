rule mykrobe:
    input:
        "data/filtered/{sample}_filtered.fastq.gz"
    output:
        "data/mykrobe/{sample}/{sample}_predict.json"
    params:
        species="tb",
        container=config["container"]
    log:
        "logs/mykrobe_{sample}.log"
    singularity:
        config["container"]
    shell:
        "scripts/run_mykrobe.sh {wildcards.sample} {params.species} {input} "
        "{output} {params.container} 2> {log}"
