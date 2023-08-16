#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Martin Rippin"
__copyright__ = "Copyright 2021, Martin Rippin"
__email__ = "martin.rippin@igp.uu.se"
__license__ = "GPL-3"


import logging
import pysam
import re
import sys
import gzip


# identify caller software from input path
def getCaller(path: str):
    pathParts = path.split("/")
    if len(pathParts) == 3:
        return pathParts[1]
    else:
        raise ValueError(
            "{} is not a valid input for this script. Required:"
            "'snv_indels/caller/sample_type.vcf'.".format(path)
        )


# modify vcf header if necessary
def modifyHeader(caller: str, header: pysam.libcbcf.VariantHeader):
    if (caller == "pisces" or caller == "gatk_mutect2" or
       caller == "pbrun_deepvariant" or caller == "pbrun_mutectcaller_t"):
        header.info.add("AF", "A", "Float", "DescriptionDescription")
    elif caller == "varscan":
        header.info.add(
            "AF",
            "A",
            "Float",
            (
                "Allel count divided on depth"
                "(Quality of bases: Phred score >= 15)"
            )
        )
    return header


# fix af field in freebayes vcf entries
def fixFreebayes(
        header: pysam.libcbcf.VariantHeader,
        row: pysam.libcbcf.VariantRecord
):
    sample = header.samples[0]
    ads = row.samples[sample].get("AD")
    af = []
    for ad in ads:
        af.append(ad/sum(ads))
    return tuple(af[1:])


# loop through input vcf and write modified entries to new vcf
def writeNewVcf(
        path: str,
        header: pysam.libcbcf.VariantHeader,
        vcf: pysam.libcbcf.VariantFile,
        caller: str
):
    new_vcf = pysam.VariantFile(path, "w", header=header)
    for row in vcf.fetch():
        if caller == "freebayes":
            row.info["AF"] = fixFreebayes(header, row)
        elif caller == "haplotypecaller":
            row.info["AF"] = row.samples[0].get("AF")
        elif caller == "gatk_mutect2":
            row.info["AF"] = row.samples[0].get("AF")
        elif caller == "pisces":
            row.info["AF"] = row.samples[0].get("VF")
        elif caller == "vardict":
            row.info["AF"] = row.samples[0].get("AF")
        elif caller == "pbrun_mutectcaller_t":
            row.info["AF"] = row.samples[0].get("AF")
        elif caller == "pbrun_deepvariant":
            row.info["AF"] = row.samples[0].get("VAF")
        elif caller == "gatk_select_variants_final":
            row.info["AF"] = row.samples[0].get("AF")
        elif caller == "varscan":
            row.info["AF"] = row.samples[0].get("AD")/row.samples[0].get("DP")
        else:
            raise ValueError(
                "{} is not a valid caller for this script. Choose between: "
                "freebayes, haplotypecaller, gatk_mutect2, gatk_select_variants_final, "
                "pbrun_deepvariant, pisces, vardict, varscan.".format(caller)
            )
        new_vcf.write(row)
    return


# call function
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, filename=snakemake.log[0])
    logging.info("Read %s", snakemake.input.vcf)
    vcf = pysam.VariantFile(snakemake.input.vcf)
    logging.info("Determine caller...")
    caller = getCaller(snakemake.input.vcf)
    logging.info("Caller is %s", caller)
    logging.info("Add info to header if necessary")
    header = modifyHeader(caller, vcf.header)
    logging.info("Start writing to %s", snakemake.output.vcf)
    writeNewVcf(snakemake.output.vcf, header, vcf, caller)
    logging.info(
        "Successfully written vcf file %s",
        snakemake.output.vcf,
    )
