__author__ = "Martin Rippin, Nina Hollfelder"
__copyright__ = "Copyright 2021, Martin Rippin"
__email__ = "nina.hollfelder@scilifelab.uu.se"
__license__ = "GPL-3"


rule manta_to_tsv:
    input:
        vcf="cnv_sv/manta_run_workflow_{analysis}/{sample}.ssa.{bed}.vcf",
    output:
        tsv=temp("tsv_files/{sample}_manta_{analysis}.{bed}.tsv"),
    log:
        "tsv_files/{sample}_manta_{analysis}.{bed}.tsv.log",
    benchmark:
        repeat(
            "tsv_files/{sample}_manta_{analysis}.{bed}.tsv.benchmark.tsv",
            config.get("manta_to_tsv", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("manta_to_tsv", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("manta_to_tsv", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("manta_to_tsv", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("manta_to_tsv", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("manta_to_tsv", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("manta_to_tsv", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("manta_to_tsv", {}).get("container", config["default_container"])
    message:
        "{rule}: select relevant info from {input.vcf} and compile into tsv"
    script:
        "../scripts/manta_to_tsv.py"


rule manta_to_tsv_no_diagnostic_filtering:
    input:
        vcf="cnv_sv/manta_run_workflow_{analysis}/{sample}.ssa.vcf",
    output:
        tsv=temp("tsv_files/{sample}_manta_{analysis}.tsv"),
    log:
        "tsv_files/{sample}_manta_{analysis}.tsv.log",
    benchmark:
        repeat(
            "tsv_files/{sample}_manta_{analysis}.tsv.benchmark.tsv",
            config.get("manta_to_tsv", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("manta_to_tsv", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("manta_to_tsv", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("manta_to_tsv", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("manta_to_tsv", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("manta_to_tsv", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("manta_to_tsv", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("manta_to_tsv", {}).get("container", config["default_container"])
    message:
        "{rule}: select info from {input.vcf} without filtering for diagnostic regions and compile into tsv"
    script:
        "../scripts/manta_to_tsv.py"
