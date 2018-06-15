rule plot:
    input:
        fastq="data/porechopped/{sample}.fastq.gz",
        bam="data/sorted/{sample}_sorted.bam"
    output:
        "data/plots/{sample}_plots.pdf"
    resources:
        mem_mb=cluster_config["plot"]["memory"]
    singularity:
        config["containers"]["nanoporeqc"]
    log:
        "logs/plot_{sample}.log"
    shell:
        "pistis --fastq {input.fastq} --output {output} --bam {input.bam} 2> {log} "


rule stats:
    input:
        "data/porechopped/{sample}.fastq.gz"
    output:
        "data/stats/{sample}_stats.txt"
    singularity:
        config["containers"]["nanoporeqc"]
    log:
        "logs/stats_{sample}.log"
    threads:
        cluster_config["stats"]["nCPUs"]
    resources:
        mem_mb=cluster_config["stats"]["memory"]
    shell:
        "NanoStat --fastq {input} --name {output} --threads {threads} "
        "--readtype 1D 2> {log}"
