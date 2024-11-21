__author__ = "Arielle R. Munters"
__copyright__ = "Copyright 2024, Arielle R. Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"


rule fix_svdb_header:
    input:
        vcf="cnv_sv/svdb_query/{sample}_{type}.{tc_method}.svdb_query.wrong_header.vcf",
    output:
        vcf="cnv_sv/svdb_query/{sample}_{type}.{tc_method}.svdb_query.vcf",
    log:
        "cnv_sv/svdb_query/{sample}_{type}.{tc_method}.svdb_query.vcf.log",
    benchmark:
        repeat(
            "cnv_sv/svdb_query/{sample}_{type}.{tc_method}.svdb_query.vcf.tsv",
            config.get("fix_svdb_header", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("fix_svdb_header", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("fix_svdb_header", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("fix_svdb_header", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("fix_svdb_header", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("fix_svdb_header", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("fix_svdb_header", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("fix_svdb_header", {}).get("container", config["default_container"])
    message:
        "{rule}: replace SAMPLE with SAMPLES in svdb header in {input.vcf}"
    shell:
        "(sed 's/_SAMPLE,/_SAMPLES,/g' {input.vcf} > {output.vcf}) &> {log}"
