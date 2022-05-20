#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import sys
import vcf

input_file = snakemake.input.vcf
output_file = snakemake.output.tsv

vcf = vcf.Reader(open(input_file, "r"))

with open(output_file, "wt") as tsv:
    tsv_writer = csv.writer(tsv, delimiter='\t')
    tsv_writer.writerow(["#POSITION1", "POSITION2", "GENES", "DEPTH"])
    for row in vcf:
        if "MantaBND" in row.ID:
            genes = row.INFO["ANN"][0].split("|")
            tsv_writer.writerow([row.CHROM + ":" + str(row.POS), row.ALT, genes[3] + "(" + genes[4] + ")", row.INFO["BND_DEPTH"]])
