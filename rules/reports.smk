rule plot_pre_filtering:
    input:
        "data/{sample}.fastq.gz",
    output:
        "data/plots/{sample}_pre_filtering.pdf",
    log:
        "logs/pistis_pre_filtering/{sample}.log",
    threads: 8
    resources:
        mem_mb=int(10 * GB),
    params:
        downsample="--downsample {}".format(config["plot_downsampling"]),
    container:
        config["container"]
    shell:
        "pistis --fastq {input} --output {output} {params.downsample} 2> {log} "


rule plot_post_filtering:
    input:
        fastq=rules.bam_to_fastq.output.fastq,
        bam=rules.sort.output.bam,
    output:
        "data/plots/{sample}_post_filtering.pdf",
    log:
        "logs/pistis_post_filtering/{sample}.log",
    threads: 8
    resources:
        mem_mb=int(10 * GB),
    params:
        downsample="--downsample {}".format(config["plot_downsampling"]),
    container:
        config["container"]
    shell:
        "pistis --fastq {input.fastq} --output {output} --bam {input.bam} "
        "{params.downsample} 2> {log} "


rule stats_pre_filtering:
    input:
        rules.plot_pre_filtering.input[0],
    output:
        "data/stats/{sample}_pre_filtering.txt",
    log:
        "logs/nanostat_pre_filtering/{sample}.log",
    threads: 4
    resources:
        mem_mb=int(10 * GB),
    container:
        config["container"]
    shell:
        "NanoStat --fastq {input} --name {output} --threads {threads} "
        "--readtype 1D 2> {log}"


rule stats_post_filtering:
    input:
        rules.bam_to_fastq.output.fastq,
    output:
        "data/stats/{sample}_post_filtering.txt",
    log:
        "logs/nanostat_post_filtering_{sample}.log",
    threads: 4
    resources:
        mem_mb=int(10 * GB),
    container:
        config["container"]
    shell:
        "NanoStat --fastq {input} --name {output} --threads {threads} "
        "--readtype 1D 2> {log}"


rule report:
    input:
        plot_pre_filter=rules.plot_pre_filtering.output[0],
        plot_post_filter=rules.plot_post_filtering.output[0],
        stats_pre_filter=rules.stats_pre_filtering.output[0],
        stats_post_filter=rules.stats_post_filtering.output[0],
        mykrobe=rules.mykrobe.output.report,
    output:
        "docs/report_{sample}.html",
    log:
        "logs/report/{sample}.log",
    container:
        config["container"]
    script:
        "../scripts/report.py"
