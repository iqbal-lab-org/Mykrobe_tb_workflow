rule demux:
    input:
        basecalled="data/basecalled",
    output:
        fastqs=directory("data/demux"),
    threads: 4
    resources:
        mem_mb=int(8 * GB),
    params:
        extras="--recursive --trim_barcodes --compress_fastq",
        kit=config["kit"],
    log:
        "logs/demux.log",
    shell:
        """
        guppy_barcoder -i {input.basecalled} -s {output.fastqs} -t {threads} \
            --barcode_kits {params.kit} {params.extras} > {log} 2>&1
        """


rule combine_barcode_fastqs:
    input:
        demux_dir=rules.demux.output.fastqs,
    output:
        fastq=expand("data/{sample}.fastq.gz", sample=SAMPLES),
    threads: 1
    resources:
        mem_mb=int(GB),
    params:
        barcode_dir=infer_barcode_dir,
    log:
        "logs/combine_barcode_fastqs/{sample}.log",
    shell:
        "cat {input.demux_dir}/{params.barcode_dir}/*.fastq* > {output.fastq} 2> {log}"
