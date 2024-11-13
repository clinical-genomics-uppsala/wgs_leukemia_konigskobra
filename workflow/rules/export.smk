__author__ = "Arielle R. Munters"
__copyright__ = "Copyright 2024, Arielle R. Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"


rule export_to_xlsx_snvs:
    input:
        vcfs=expand(
            "parabricks/pbrun_mutectcaller_{{analysis}}/{{sample_type}}.normalized.vep.ratio.filter.somatic.include.{bed}.vcf.gz",
            bed=["all", "aml", "tm"],
        ),
        vcf_pindel="cnv_sv/pindel_vcf/{sample_type}_T.no_tc.vep_annotated.vcf",
        all_bed=config["bcftools_SNV"]["all"],
        aml_bed=config["bcftools_SNV"]["aml"],
        tm_bed=config["bcftools_SNV"]["tm"],
        pindel_bed=config["pindel_call"]["include_bed"],
    output:
        xlsx=temp("export_to_xlsx/{analysis}/{sample_type}.snvs.xlsx"),
    params:
        filterfile = config["filter_vcf"]["somatic"],
        extra=config.get("export_to_xlsx_snvs", {}).get("extra", ""),
    log:
        "export_to_xlsx/{analysis}/{sample_type}.snvs.xslx.log",
    benchmark:
        repeat(
            "export_to_xlsx/{analysis}/{sample_type}.snvs.xslx.benchmark.tsv",
            config.get("export_to_xlsx_snvs", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("export_to_xlsx_snvs", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("export_to_xlsx_snvs", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("export_to_xlsx_snvs", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("export_to_xlsx_snvs", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("export_to_xlsx_snvs", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("export_to_xlsx_snvs", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("export_to_xlsx_snvs", {}).get("container", config["default_container"])
    message:
        "{rule}: merge {input.vcfs} into {output.xlsx}"
    script:
        "../scripts/export_to_xlsx_snvs.py"


rule export_to_xlsx_manta:
    input:
        vcf="cnv_sv/manta_run_workflow_{analysis}/{sample}.ssa.vcf",
        vcfs_bed=expand("cnv_sv/manta_run_workflow_{{analysis}}/{{sample}}.ssa.include.{bed}.vcf.gz", bed=["all", "aml"]), #ska tm med?
        all_bed=config["bcftools_SV"]["all"],
        aml_bed=config["bcftools_SV"]["aml"],
    output:
        temp("export_to_xlsx/{analysis}/{sample}.manta.xlsx"),
    params:
        extra=config.get("export_to_xlsx_manta", {}).get("extra", ""),
    log:
        "export_to_xlsx/{analysis}/{sample}.manta.xlsx.log",
    benchmark:
        repeat(
            "export_to_xlsx/{analysis}/{sample}.manta.xlsx.benchmark.tsv",
            config.get("export_to_xlsx_manta", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("export_to_xlsx_manta", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("export_to_xlsx_manta", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("export_to_xlsx_manta", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("export_to_xlsx_manta", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("export_to_xlsx_manta", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("export_to_xlsx_manta", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("export_to_xlsx_manta", {}).get("container", config["default_container"])
    message:
        "{rule}: merge {input.vcfs_bed} and {input.vcf} into {output.xlsx}"
    script:
        "../scripts/export_to_xlsx_manta.py"
