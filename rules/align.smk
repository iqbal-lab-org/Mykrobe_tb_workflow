rule map_minimap2:
    input:
        ref=config["reference"],
        query="data/porechopped/{sample}.fastq.gz"

    output:
        temp("data/mapped/{sample}.bam")
    threads:
        cluster_config["map_minimap2"]["nCPUs"]
    resources:
        mem_mb=cluster_config["map_minimap2"]["memory"]
    log:
        "logs/minimap2_{sample}.log"
    shell:
        "(minimap2 -t {threads} -ax map-ont {input.ref} {input.query} | "
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
    shell:
        "samtools sort -@ {threads} {input} 2> {log} > {output}"


rule samtools_index:
    input:
        "data/sorted/{sample}_sorted.bam"
    output:
        "data/sorted/{sample}_sorted.bam.bai"
    log:
        "logs/samtools_index_{sample}.log"
    shell:
        "samtools index -b {input} 2> {log}"


rule bam_to_fastq:
    input:
        "data/sorted/{sample}_sorted.bam"
    output:
        "data/filtered/{sample}_filtered.fastq.gz"
    log:
        "logs/bam_to_fastq_{sample}.log"
    shell:
        "(samtools fastq -F 0x4 {input} | gzip > {output}) 2> {log}"
