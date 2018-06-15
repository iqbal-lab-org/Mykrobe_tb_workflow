rule basecall_albacore:
    input:
        "data/reads"
    output:
        "data/basecalled/workspace/pass"
    params:
        kit=config["kit"],
        flowcell=config["flowcell"],
        barcoding=("--barcoding" if MULTIPLEXED else "")
    singularity:
        config["containers"]["albacore"]
    threads:
        cluster_config["basecall_albacore"]["nCPUs"]
    resources:
        mem_mb=cluster_config["basecall_albacore"]["memory"]
    shell:
        "read_fast5_basecaller.py {params.barcoding} --worker_threads {threads} "
        "--save_path data/basecalled --input {input} --flowcell {params.flowcell} "
        "--kit {params.kit} --output_format fastq --recursive"

if not MULTIPLEXED:
    rule concat_and_gzip_fastq:
        input:
            "data/basecalled/workspace/pass"
        output:
            "data/basecalled/{sample}.fastq.gz"
        shell:
            "cat {input}/*.fastq | gzip > {output}"
