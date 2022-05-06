__author__ = "Martin Rippin"
__copyright__ = "Copyright 2022, Martin Rippin"
__email__ = "martin.rippin@igp.uu.se"
__license__ = "GPL-3"

import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *

min_version("6.10.0")

### Set and validate config file


configfile: "config.yaml"


validate(config, schema="../schemas/config.schema.yaml")
config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")


### Read and validate samples file

samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

### Read and validate units file

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "type", "flowcell", "lane", "barcode"], drop=False).sort_index()
validate(units, schema="../schemas/units.schema.yaml")

### Set wildcard constraints


wildcard_constraints:
    sample="|".join(samples.index),
    unit="N|T|R",


def compile_output_list(wildcards):
    chromosomes = ["X", "Y"]
    chromosomes.extend(range(1, 23))
    output_list = ["qc/multiqc/multiqc.html"]
    output_list.append(["cnv_sv/cnvkit_vcf/%s_T.vcf" % (sample) for sample in get_samples(samples)])
    output_list.append(["cnv_sv/pindel/%s.vcf" % (sample) for sample in get_samples(samples)])
    output_list.append(["cnv_sv/cnvkit_diagram/%s_T.png" % (sample) for sample in get_samples(samples)])
    output_list.append(
        [
            "cnv_sv/cnvkit_scatter/%s_T_chr%s.png" % (sample, chromosome)
            for sample in get_samples(samples)
            for chromosome in chromosomes
        ]
    )
    output_list.append(
        [
            "cnv_sv/manta_run_workflow_tn/%s.ssa.%s.vcf" % (sample, diagnosis)
            for sample in get_samples(samples)
            for diagnosis in ["aml", "all"]
        ]
    )
    output_list.append(
        [
            "parabricks/pbrun_mutectcaller_tn/%s.vep.%s.vcf" % (sample, diagnosis)
            for sample in get_samples(samples)
            for diagnosis in ["aml", "all"]
        ]
    )
    output_list.append(
        [
            "compression/spring/%s_%s_%s_%s_%s.spring" % (sample, flowcell, lane, barcode, t)
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
            for flowcell in set([u.flowcell for u in units.loc[(sample,t,)].dropna().itertuples()])
            for barcode in set([u.barcode for u in units.loc[(sample,t,)].dropna().itertuples()])
            for lane in set([u.lane for u in units.loc[(sample,t,)].dropna().itertuples()])
        ]
    )
    output_list.append(
        [
            "compression/crumble/%s_%s.crumble.cram%s" % (sample, t, ext)
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
            for ext in ["", ".crai"]
        ]
    )
    return output_list
