rule demux:
    input:
        basecalled="data/basecalled",
    output:
        fastqs=directory("data/demux"),
    threads: 8
    resources:
        mem_mb=int(8 * GB),
    params:
        extras="--recursive --trim_barcodes --compress_fastq",
        kit=config["kit"],
    log:
        "logs/demux.log",
    container:
        config["container"]
    shell:
        """
        guppy_barcoder -i {input.basecalled} -s {output.fastqs} -t {threads} \
            --barcode_kits {params.kit} {params.extras} > {log} 2>&1
        """


rule combine_barcode_fastqs:
    input:
        demux_dir=rules.demux.output.fastqs,
    output:
        fastq="data/{sample}.fastq.gz",
    threads: 1
    resources:
        mem_mb=int(GB),
    params:
        barcode_dir=infer_barcode_dir,
    container:
        config["container"]
    log:
        "logs/combine_barcode_fastqs/{sample}.log",
    shell:
        "cat {input.demux_dir}/pass/{params.barcode_dir}/*.fastq* > {output.fastq} 2> {log}"
