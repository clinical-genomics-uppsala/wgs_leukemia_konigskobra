#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import sys
import vcf

input_file = snakemake.input.vcf
output_file_dels = snakemake.output.dels
output_file_dup = snakemake.output.dups
output_file_ins = snakemake.output.ins

datei = vcf.Reader(open(input_file, "r"))

with open(output_file_dels, "wt") as tsv:
    tsv_writer = csv.writer(tsv, delimiter='\t')
    tsv_writer.writerow(["#POSITION1", "POSITION2", "LENGTH", "MANTAID", "GENES"])
    for row in datei:
        if "MantaDEL" in row.ID and not bool(row.FILTER):
            genes = row.INFO["ANN"][0].split("|")
            dellength = row.INFO["SVLEN"]
            pos2 = row.INFO["END"]
            manta_id = ":".join(row.ID.split(":")[0:2])
            tsv_writer.writerow([row.CHROM + ":" + str(row.POS), row.CHROM + ":" + str(pos2),
                                str(dellength)[1:-1], manta_id, genes[3] + "(" + genes[4] + ")"])


datei = vcf.Reader(open(input_file, "r"))

with open(output_file_dup, "wt") as tsv:
    tsv_writer = csv.writer(tsv, delimiter='\t')
    tsv_writer.writerow(["#POSITION1", "POSITION2", "LENGTH", "MANTAID", "GENES", "HOMLENGTH", "HOMSEQ"])
    for row in datei:
        if "MantaDUP" in row.ID and not bool(row.FILTER):
            genes = row.INFO["ANN"][0].split("|")
            dellength = row.INFO["SVLEN"]
            pos2 = row.INFO["END"]
            manta_id = ":".join(row.ID.split(":")[0:3])
            try:
                homlen = row.INFO["HOMLEN"]
                homseq = row.INFO["HOMSEQ"]
            except KeyError:
                homlen = "NA"
                homseq = "NA"
            tsv_writer.writerow([row.CHROM + ":" + str(row.POS), row.CHROM + ":" + str(pos2), str(dellength)[1:-1],
                                manta_id, genes[3] + "(" + genes[4] + ")", str(homlen)[1:-1], str(homseq)[1:-1]])


datei = vcf.Reader(open(input_file, "r"))

with open(output_file_ins, "wt") as tsv:
    tsv_writer = csv.writer(tsv, delimiter='\t')
    tsv_writer.writerow(["#POSITION", "REFERENCE" "ALTERNATIVE", "LENGTH", "MANTAID", "GENES", "HOMLENGTH", "HOMSEQ"])
    for row in datei:
        if "MantaINS" in row.ID and not bool(row.FILTER):
            genes = row.INFO["ANN"][0].split("|")
            dellength = row.INFO["SVLEN"]
            pos2 = row.INFO["END"]
            try:
                homlen = row.INFO["HOMLEN"]
                homseq = row.INFO["HOMSEQ"]
            except KeyError:
                homlen = "NA"
                homseq = "NA"
            manta_id = ":".join(row.ID.split(":")[0:2])
            tsv_writer.writerow([row.CHROM + ":" + str(row.POS), row.REF, str(row.ALT)[1:-1], str(dellength)[1:-1],
                                 manta_id, genes[3] + "(" + genes[4] + ")", str(homlen)[1:-1], str(homseq)[1:-1]])
