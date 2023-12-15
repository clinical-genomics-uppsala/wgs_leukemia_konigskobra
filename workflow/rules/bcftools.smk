__author__ = "Nina Hollfelder"
__copyright__ = "Copyright 2023, Nina Hollfelder"
__email__ = "nina.hollfelder@scilifelab.uu.se"
__license__ = "GPL-3"


rule bcftools_split_tn_vcf:
    input:
        vcf="parabricks/pbrun_mutectcaller_tn/{sample}.vcf",
    output:
        vcf="parabricks/pbrun_mutectcaller_tn/{sample}_N.only.vcf",
    log:
        "parabricks/pbrun_mutectcaller_tn/{sample}_N.only.vcf.log",
    benchmark:
        repeat(
            "parabricks/pbrun_mutectcaller_tn/{sample}_N.only.vcf.benchmark.tsv",
            config.get("bcftools_split_tn_vcf", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("bcftools_split_tn_vcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("bcftools_split_tn_vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcftools_split_tn_vcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bcftools_split_tn_vcf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bcftools_split_tn_vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("bcftools_split_tn_vcf", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("bcftools_split_tn_vcf", {}).get("container", config["default_container"])
    message:
        "{rule}: Create fam file for all samples for peddy input"
    shell:
        "bcftools view {input.vcf} -s {wildcards.sample}_N > {output.vcf}"

