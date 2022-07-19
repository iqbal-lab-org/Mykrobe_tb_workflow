rule mykrobe:
    input:
        reads=rules.bam_to_fastq.output.fastq,
    output:
        report="data/mykrobe/{sample}_predict.json",
    shadow:
        "shallow"
    params:
        species="tb",
        extra="--format json --ont",
    log:
        "logs/mykrobe/{sample}.log",
    container:
        config["container"]
    shell:
        """
        mykrobe predict {params.extra} --seq {input.reads} --output {output.report} \
          --sample {wildcards.sample} --species {params.species} > {log} 2>&1
        """
