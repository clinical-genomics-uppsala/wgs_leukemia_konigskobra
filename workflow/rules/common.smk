__author__ = "Martin Rippin, Arielle R. Munters, Nina Hollfelder"
__copyright__ = "Copyright 2022, Martin Rippin"
__email__ = "arielle.munters@scilifelab.uu.se, nina.hollfelder@scilifelab.uu.se"
__license__ = "GPL-3"

import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *

min_version("7.8.0")


### Set and validate config file
configfile: "config.yaml"


validate(config, schema="../schemas/config.schema.yaml")
config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")


### Read and validate samples file
samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

### Read and validate units file
units = (
    pd.read_table(config["units"], dtype=str)
    .set_index(["sample", "type", "flowcell", "lane", "barcode"], drop=False)
    .sort_index()
)
validate(units, schema="../schemas/units.schema.yaml")


### Read output_files for cp rules
with open(config["output"]) as output:
    output_json = json.load(output)


### Set wildcard constraints
wildcard_constraints:
    sample="|".join(samples.index),
    unit="N|T|R",
    bed="aml|all",


def type_generator(types):
    if "N" in types and "T" in types:
        types.add("TN")
        return types
    else:
        return types

def get_bam_input(wildcards, t_n=None, use_sample_wildcard=True):
    if use_sample_wildcard is True and t_n is None:
        sample_str = "{}_{}".format(wildcards.sample, wildcards.type)
    elif use_sample_wildcard is True and t_n:
        sample_str = "{}_{}".format(wildcards.sample, t_n)
    else:
        sample_str = wildcards.file


    aligner = config.get("aligner", None)
    if aligner is None:
        sys.exit("aligner missing from config, valid options: bwa_gpu or bwa_sentieon")
    elif aligner == "bwa_gpu":
        bam_input = "parabricks/pbrun_fq2bam/{}.bam".format(sample_str)
    elif aligner == "bwa_sentieon":
        bam_input = "sentieon/realign/{}_REALIGNED.bam".format(sample_str)
    else:
        sys.exit("valid options for aligner are: bwa_gpu or bwa_sentieon")

    bai_input = "{}.bai".format(bam_input)

    return (bam_input, bai_input)

def get_vcf_input(wildcards):
    aligner = config.get("aligner", None)
    if aligner is None:
        sys.exit("aligner missing from config, valid options: bwa_gpu or bwa_sentieon")
    elif aligner == "bwa_gpu":
        vcf_input = "parabricks/pbrun_mutectcaller_t/{}_{}.vep.filter.germline.vcf".format(wildcards.sample, wildcards.type)
    elif aligner == "bwa_sentieon":
        vcf_input = "sentieon/tnscope/{}_TNscope_tn_ML.vcf".format(wildcards.sample)
    else:
        sys.exit("valid options for aligner are: bwa_gpu or bwa_sentieon")

    return vcf_input

def compile_output_list(wildcards):
    output_files = []
    types = type_generator(set([unit.type for unit in units.itertuples()]))
    chromosome_numbers = ["X", "Y"]
    chromosome_numbers.extend(range(1, 23))
    for output in output_json:
        output_files += set(
            [
                output.format(
                    sample=sample,
                    type=unit_type,
                    project=samples.loc[(sample)]["project"],
                    chr=chromosome_number,
                    flowcell=flowcell,
                    barcode=barcode,
                    lane=lane,
                )
                for chromosome_number in chromosome_numbers
                for sample in get_samples(samples)
                for unit_type in get_unit_types(units, sample)
                if unit_type in set(output_json[output]["types"])
                for flowcell in set([u.flowcell for u in units.loc[(sample, unit_type)].dropna().itertuples()])
                for barcode in set([u.barcode for u in units.loc[(sample, unit_type)].dropna().itertuples()])
                for lane in set([u.lane for u in units.loc[(sample, unit_type)].dropna().itertuples()])
            ]
        )
    for output in output_json:
        output_files += set(
            [
                output.format(
                    sample=sample,
                    type=unit_type,
                    project=samples.loc[(sample)]["project"],
                    chr=chromosome_number,
                )
                for chromosome_number in chromosome_numbers
                for sample in get_samples(samples)
                for unit_type in type_generator(get_unit_types(units, sample))
                if unit_type in set(output_json[output]["types"]) and unit_type == "TN"
            ]
        )
    return list(set(output_files))


def generate_copy_code(workflow, output_json):
    code = ""
    for result, values in output_json.items():
        if values["file"] is not None:
            input_file = values["file"]
            output_file = result
            rule_name = values["name"]
            mem_mb = config.get("_copy", {}).get("mem_mb", config["default_resources"]["mem_mb"])
            mem_per_cpu = config.get("_copy", {}).get("mem_mb", config["default_resources"]["mem_mb"])
            partition = config.get("_copy", {}).get("partition", config["default_resources"]["partition"])
            threads = config.get("_copy", {}).get("threads", config["default_resources"]["threads"])
            time = config.get("_copy", {}).get("time", config["default_resources"]["time"])
            copy_container = config.get("_copy", {}).get("container", config["default_container"])
            result_file = os.path.basename(output_file)
            code += f'@workflow.rule(name="{rule_name}")\n'
            code += f'@workflow.input("{input_file}")\n'
            code += f'@workflow.output("{output_file}")\n'
            if "{project}" in output_file:
                if "{chr}" in output_file:
                    code += f'@workflow.log("logs/{rule_name}_{{project}}_{result_file}_chr{{chr}}.log")\n'
                else:
                    code += f'@workflow.log("logs/{rule_name}_{{project}}_{result_file}.log")\n'
            else:
                code += f'@workflow.log("logs/{rule_name}_{result_file}.log")\n'
            code += f'@workflow.container("{copy_container}")\n'
            code += f'@workflow.conda("../envs/copy_result.yaml")\n'
            code += f'@workflow.resources(time = "{time}", threads = {threads}, mem_mb = {mem_mb}, mem_per_cpu = {mem_per_cpu}, partition = "{partition}")\n'
            code += '@workflow.shellcmd("cp {input} {output}")\n\n'
            code += "@workflow.run\n"
            code += (
                f"def __rule_{rule_name}(input, output, params, wildcards, threads, resources, log, version, rule, "
                "conda_env, container_img, singularity_args, use_singularity, env_modules, bench_record, jobid, is_shell, "
                "bench_iteration, cleanup_scripts, shadow_dir, edit_notebook, conda_base_path, basedir, runtime_sourcecache_path, "
                "__is_snakemake_rule_func=True):\n"
                '\tshell ( "(cp {input[0]} {output[0]}) &> {log}" , bench_record=bench_record, bench_iteration=bench_iteration)\n\n'
            )
    exec(compile(code, "result_to_copy", "exec"), workflow.globals)


generate_copy_code(workflow, output_json)
