#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import sys
import vcf

input_file = snakemake.input.vcf
analysis = snakemake.params.analysis
sample = snakemake.params.sample
output_file = snakemake.output.tsv

vcf = vcf.Reader(open(input_file, "r"))

with open(output_file, "wt") as tsv:
    tsv_writer = csv.writer(tsv, delimiter='\t')
    tsv_writer.writerow(["#SAMPLE", "CHROMOSOME", "POSITION", "REFERENCE", "ALTERNATIVE", "ALLELEFREQUENCY", "DEPTH", "GENE", "TRANSCRIPT", "NUCLEOTIDESEQ", "PROTEIN", "PROTEINSEQ", "CONSEQUENCE"])
    for row in vcf:
        genes = row.INFO["CSQ"][0].split("|")
        transcript_field = genes[10].split(":")
        if len(transcript_field) == 1:
            transcript_field = ["", ""]
        protein_field = genes[11].split(":")
        if len(protein_field) == 1:
            protein_field = [genes[24], ""]
        i = 0
        if analysis == "tn":
            i = 1
        af = row.samples[i]["AF"]
        tsv_writer.writerow([sample, row.CHROM, row.POS, row.REF, row.ALT, af, row.INFO["DP"], genes[3], transcript_field[0], transcript_field[1], protein_field[0], protein_field[1], genes[1]])
