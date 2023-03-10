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
    tsv_writer.writerow(["#POSITION1", "MANTAID","BREAKEND", "GENES", "DEPTH", "ANNOTATIONINFO"])
    for row in vcf:
        if "MantaBND" in row.ID and not any(x in ["MinQUAL", "MinGQ", "MinSomaticScore",
                                                  "Ploidy", "MaxDepth", "MaxMQ0Frac", "NoPairSupport",
                                                  "SampleFT", "HomRef"] for x in row.FILTER):
            try:
                genes = row.INFO["ANN"][0].split("|")
            except KeyError:
                genes = ["NA", "NA", "NA", "NA", "NA"]
            manta_id = ":".join(row.ID.split(":")[0:2])
            tsv_writer.writerow([row.CHROM + ":" + str(row.POS), manta_id, str(row.ALT)[1:-1],
                                genes[3] + "(" + genes[4] + ")", row.INFO["BND_DEPTH"], row.FILTER])
