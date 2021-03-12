rule subsample:
    input:
        reads="data/{sample}.fastq.gz",
    output:
        reads="data/{sample}.subsampled.fq.gz",
    threads: 1
    params:
        options="",
        genome_size=4411532,
        coverage=config["max_covg"],
    log:
        "logs/subsample/{sample}.log",
    container:
        config["container"]
    shell:
        """
        rasusa {params.options} -i {input.reads} -o {output.reads} \
          -c {params.coverage} -g {params.genome_size} 2> {log}
        """
