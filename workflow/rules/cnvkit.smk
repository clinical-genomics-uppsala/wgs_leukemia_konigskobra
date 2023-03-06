__author__ = "Nina Hollfelder"
__copyright__ = "Copyright 2021, Nina Hollfelder"
__email__ = "nina.hollfelder@scilifelab.uu.se"
__license__ = "GPL-3"


rule cnvkit_call_no_custom_purity:
    input:
        segment="cnv_sv/cnvkit_batch/{sample}/{sample}_{type}.cns",
        vcf="parabricks/pbrun_mutectcaller_t/{sample}_{type}.vep.filter.germline.vcf",
    output:
        segment=temp("cnv_sv/cnvkit_call/{sample}_{type}.loh.cns"),
    log:
        "cnv_sv/cnvkit_call/{sample}_{type}.loh.cns.log",
    benchmark:
        repeat(
            "cnv_sv/cnvkit_call/{sample}_{type}.loh.cns.benchmark.tsv",
            config.get("cnvkit_call", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnvkit_call", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("cnvkit_call", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvkit_call", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("cnvkit_call", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvkit_call", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvkit_call", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("cnvkit_call", {}).get("container", config["default_container"])
    message:
        "{rule}: Call cnvs with loh info into cnv_sv/cnvkit_call/{wildcards.sample}_{wildcards.type}.loh.cns"
    shell:
        "(cnvkit.py call {input.segment} "
        "-v {input.vcf} "
        "-o {output.segment}) &> {log}"
