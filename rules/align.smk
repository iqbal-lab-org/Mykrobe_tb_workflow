rule map_minimap2:
    input:
        config["reference"],
        "data/porechopped/{sample}.fastq.gz"

    output:
        temp("data/mapped/{sample}.bam")
    threads:
        cluster_config["map_minimap2"]["nCPUs"]
    resources:
        mem_mb=cluster_config["map_minimap2"]["memory"]
    log:
        "logs/minimap2_{sample}.log"
    singularity:
        config["containers"]["nanoporeqc"]
    shell:
        "(minimap2 -t {threads} -ax map-ont {input} | "
        "samtools view -b - > {output}) 2> {log}"


rule samtools_sort:
    input:
        "data/mapped/{sample}.bam"
    output:
        "data/sorted/{sample}_sorted.bam"
    threads:
        cluster_config["samtools_sort"]["nCPUs"]
    resources:
        mem_mb=cluster_config["samtools_sort"]["memory"]
    log:
        "logs/samtools_sort_{sample}.log"
    singularity:
        config["containers"]["nanoporeqc"]
    shell:
        "samtools sort -@ {threads} {input} 2> {log} > {output}"


rule samtools_index:
    input:
        "data/sorted/{sample}_sorted.bam"
    output:
        "data/sorted/{sample}_sorted.bam.bai"
    log:
        "logs/samtools_index_{sample}.log"
    singularity:
        config["containers"]["nanoporeqc"]
    shell:
        "samtools index -b {input} 2> {log}"
