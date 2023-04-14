__author__ = "Nina Hollfelder"
__copyright__ = "Copyright 2023, Nina Hollfelder"
__email__ = "nina.hollfelder@scilifelab.uu.se"
__license__ = "GPL-3"


rule peddy_create_ped:
    input:
        config["samples"],
    output:
        ["qc/peddy/" + sample + ".peddy.fam" for sample in get_samples(samples)],
    log:
        "qc/peddy/peddy_create_ped.fam.log",
    benchmark:
        repeat(
            "qc/peddy/peddy_create_ped.fam.benchmark.tsv",
            config.get("peddy_create_ped", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("peddy_create_ped", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("peddy_create_ped", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("peddy_create_ped", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("peddy_create_ped", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("peddy_create_ped", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("peddy_create_ped", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("peddy_create_ped", {}).get("container", config["default_container"])
    message:
        "{rule}: Create fam file for all samples for peddy input"
    script:
        "../scripts/peddy_create_ped.py"
