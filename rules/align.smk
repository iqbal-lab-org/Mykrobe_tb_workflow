rule map:
    input:
        ref=config["reference"],
        query="data/{sample}.fastq.gz",
    output:
        bam=temp("data/mapped/{sample}.bam"),
    threads: 8
    resources:
        mem_mb=int(16 * GB),
    log:
        "logs/map/{sample}.log",
    params:
        extra="-aL -x map-ont --secondary=no",
    container:
        config["container"]
    shell:
        """
        (minimap2 -t {threads} {params.extra} {input.ref} {input.query} \
        | samtools view -b - > {output.bam}) 2> {log}
        """


rule sort:
    input:
        bam=rules.map.output.bam,
    output:
        bam="data/sorted/{sample}.bam",
    threads: 8
    resources:
        mem_mb=int(10 * GB),
    log:
        "logs/sort/{sample}.log",
    params:
        extra="-O BAM",
    container:
        config["container"]
    shell:
        "samtools sort {params.extra} -@ {threads} {input.bam} 2> {log} > {output.bam}"


rule index:
    input:
        bam=rules.sort.output.bam,
    output:
        idx="data/sorted/{sample}.bam.bai",
    log:
        "logs/index/{sample}.log",
    container:
        config["container"]
    shell:
        "samtools index -b {input.bam} 2> {log}"


rule bam_to_fastq:
    input:
        bam=rules.index.input.bam,
        index=rules.index.output.idx,
    output:
        fastq="data/filtered/{sample}.filtered.fastq.gz",
    log:
        "logs/bam_to_fastq/{sample}.log",
    params:
        extra="-F 2308",
    container:
        config["container"]
    shell:
        "(samtools fastq {params.extra} {input.bam} | gzip > {output.fastq}) 2> {log}"
