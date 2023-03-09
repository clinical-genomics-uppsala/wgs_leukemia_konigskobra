#!/bin/python3

import sys
import csv


header_row = "Sample_ID,Sample_Name,Description,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project\n"
header_row_split = header_row.strip().split(",")

samples = {}
header = False
with open(snakemake.input.sample_sheet, 'r') as samplesheet:
    for lline in samplesheet:
        if len(lline.split()) == 0:  # skip blank lines
            continue
        if header:
            line = lline.strip().split(",")
            if line[7] == "WGSWP2":
                samples[line[0]] = {
                                        header_row_split[1]: line[1],
                                        "Type": line[2].split("_")[0],
                                        "Sex": line[2].split("_")[1],
                                        "Pedegree_id": line[2].split("_")[2],
                                        "Seq_run": line[2].split("_")[3],
                                        "Project": line[2].split("_")[4],
                                        header_row_split[3]: line[3],
                                        header_row_split[4]: line[4],
                                        header_row_split[5]: line[5],
                                        header_row_split[6]: line[6],
                                        header_row_split[7]: line[7],
                                    }
        if lline == header_row:
            header = True


if len(samples) == 0:
    raise Exception("No samples found, has the header in SampleSheet changed?")

# replacement files for rna and dna samples seperate to get sample order correct
with open(snakemake.output.replacement_dna, "w+") as replacement_tsv_dna, \
     open(snakemake.output.replacement_rna, "w+") as replacement_tsv_rna, \
     open(snakemake.output.order_dna, "w+") as order_tsv_dna, \
     open(snakemake.output.order_rna, "w+") as order_tsv_rna:
    order_tsv_dna.write("\t".join(["Sample Order", "Pedegree ID", "DNA number"])+"\n")
    order_tsv_rna.write("\t".join(["Sample Order", "Pedegree ID", "RNA number"])+"\n")
    i = 1
    j = 1
    for sample in samples.values():
        sample["Type"] = "R" if sample["Type"] == "Heltranskriptom" else sample["Type"]
        replacement_line = [sample["Pedegree_id"] + "_" + sample["Type"]]
        order_line = [sample["Pedegree_id"] + "_" + sample["Type"], sample["Sample_Name"]]
        if sample["Type"] == "T" or sample["Type"] == "N":
            replacement_tsv_dna.write("\t".join(replacement_line + ["sample_"+str(f"{i:03}")]) + "\n")
            order_tsv_dna.write("\t".join(["sample_"+str(f"{i:03}")] + order_line) + "\n")
            i += 1
        elif sample["Type"] == "R":
            replacement_tsv_rna.write("\t".join(replacement_line + ["sample_"+str(f"{j:03}")]) + "\n")
            order_tsv_rna.write("\t".join(["sample_"+str(f"{j:03}")] + order_line) + "\n")
            j += 1
