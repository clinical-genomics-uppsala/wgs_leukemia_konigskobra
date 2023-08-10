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
    tsv_writer.writerow(["#POSITION1", "POSITION2", "LENGTH", "MANTAID", "GENES", "ANNOTATIONINFO", "PR_ALT_FREQ", "SR_ALT_FREQ"])
    for row in datei:
        if "MantaDEL" in row.ID and not any(xxx in ["MinQUAL", "MinGQ", "MinSomaticScore",
                                                    "Ploidy", "MaxDepth", "MaxMQ0Frac", "NoPairSupport", "SampleFT", "HomRef"]
                                            for xxx in row.FILTER):
            genes = row.INFO["ANN"][0].split("|")
            dellength = row.INFO["SVLEN"]
            pos2 = row.INFO["END"]
            manta_id = ":".join(row.ID.split(":")[0:2])
            # get frequency of paired and split alternate reads
            last_sample_index = len(row.samples) - 1
            last_sample = row.samples[last_sample_index]
            pr_values = last_sample.data.PR if hasattr(last_sample.data, 'PR') else None
            sr_values = last_sample.data.SR if hasattr(last_sample.data, 'SR') else None
            if pr_values:
                pr_denominator, pr_numerator = pr_values
                pr_frequency = pr_numerator / (pr_denominator + pr_numerator) if pr_denominator + pr_numerator != 0 else None
            else:
                pr_frequency = None

            if sr_values:
                sr_denominator, sr_numerator = sr_values
                sr_frequency = sr_numerator / (sr_denominator + sr_numerator) if sr_denominator + sr_numerator != 0 else None
            else:
                sr_frequency = None

            # Format frequencies only if they are not None
            pr_formatted = f"{pr_frequency:.4f}" if pr_frequency is not None else 'NA'
            sr_formatted = f"{sr_frequency:.4f}" if sr_frequency is not None else 'NA'
            if dellength[0] <= -100:
                tsv_writer.writerow([row.CHROM + ":" + str(row.POS), row.CHROM + ":" + str(pos2),
                                    str(dellength)[1:-1], manta_id, genes[3] + "(" + genes[4] + ")",
                                    row.FILTER, str(pr_formatted), str(sr_formatted)])


datei = vcf.Reader(open(input_file, "r"))

with open(output_file_dup, "wt") as tsv:
    tsv_writer = csv.writer(tsv, delimiter='\t')
    tsv_writer.writerow(["#POSITION1", "POSITION2", "LENGTH", "MANTAID", "GENES", "HOMLENGTH", "HOMSEQ", "ANNOTATIONINFO",
                         "PR_ALT_FREQ", "SR_ALT_FREQ"])
    for row in datei:
        if "MantaDUP" in row.ID and not any(xxx in ["MinQUAL", "MinGQ", "MinSomaticScore", "Ploidy",
                                                    "MaxDepth", "MaxMQ0Frac", "NoPairSupport", "SampleFT",
                                                    "HomRef"] for xxx in row.FILTER):
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
            # get frequency of paired and split alternate reads
            last_sample_index = len(row.samples) - 1
            last_sample = row.samples[last_sample_index]
            pr_values = last_sample.data.PR if hasattr(last_sample.data, 'PR') else None
            sr_values = last_sample.data.SR if hasattr(last_sample.data, 'SR') else None
            if pr_values:
                pr_denominator, pr_numerator = pr_values
                pr_frequency = pr_numerator / (pr_denominator + pr_numerator) if pr_denominator + pr_numerator != 0 else None
            else:
                pr_frequency = None

            if sr_values:
                sr_denominator, sr_numerator = sr_values
                sr_frequency = sr_numerator / (sr_denominator + sr_numerator) if sr_denominator + sr_numerator != 0 else None
            else:
                sr_frequency = None

            # Format frequencies only if they are not None
            pr_formatted = f"{pr_frequency:.4f}" if pr_frequency is not None else 'NA'
            sr_formatted = f"{sr_frequency:.4f}" if sr_frequency is not None else 'NA'
            if dellength[0] >= 100:
                tsv_writer.writerow([row.CHROM + ":" + str(row.POS), row.CHROM + ":" + str(pos2), str(dellength)[1:-1],
                                    manta_id, genes[3] + "(" + genes[4] + ")", str(homlen)[1:-1], str(homseq)[1:-1],
                                    row.FILTER, str(pr_formatted), str(sr_formatted)])


datei = vcf.Reader(open(input_file, "r"))

with open(output_file_ins, "wt") as tsv:
    tsv_writer = csv.writer(tsv, delimiter='\t')
    tsv_writer.writerow(["#POSITION", "REFERENCE", "ALTERNATIVE", "LENGTH", "MANTAID", "GENES",
                         "HOMLENGTH", "HOMSEQ", "ANNOTATIONINFO", "PR_ALT_FREQ", "SR_ALT_FREQ"])
    for row in datei:
        if "MantaINS" in row.ID and not any(xxx in ["MinQUAL", "MinGQ", "MinSomaticScore", "Ploidy",
                                                    "MaxDepth", "MaxMQ0Frac", "NoPairSupport", "SampleFT", "HomRef"]
                                            for xxx in row.FILTER):
            try:
                genes = row.INFO["ANN"][0].split("|")
            except KeyError:
                genes = ["NA", "NA", "NA", "NA", "NA"]
            try:
                dellength = row.INFO["SVLEN"]
            except KeyError:
                dellength = "NA"
            try:
                homlen = row.INFO["HOMLEN"]
                homseq = row.INFO["HOMSEQ"]
            except KeyError:
                homlen = "NA"
                homseq = "NA"
            manta_id = ":".join(row.ID.split(":")[0:2])
            # get frequency of paired and split alternate reads
            last_sample_index = len(row.samples) - 1
            last_sample = row.samples[last_sample_index]
            pr_values = last_sample.data.PR if hasattr(last_sample.data, 'PR') else None
            sr_values = last_sample.data.SR if hasattr(last_sample.data, 'SR') else None
            if pr_values:
                pr_denominator, pr_numerator = pr_values
                pr_frequency = pr_numerator / (pr_denominator + pr_numerator) if pr_denominator + pr_numerator != 0 else None
            else:
                pr_frequency = None

            if sr_values:
                sr_denominator, sr_numerator = sr_values
                sr_frequency = sr_numerator / (sr_denominator + sr_numerator) if sr_denominator + sr_numerator != 0 else None
            else:
                sr_frequency = None

            # Format frequencies only if they are not None
            pr_formatted = f"{pr_frequency:.4f}" if pr_frequency is not None else 'NA'
            sr_formatted = f"{sr_frequency:.4f}" if sr_frequency is not None else 'NA'
            if dellength == "NA" or dellength[0] >= 100:
                tsv_writer.writerow([row.CHROM + ":" + str(row.POS), row.REF, str(row.ALT)[1:-1], str(dellength)[1:-1],
                                     manta_id, genes[3] + "(" + genes[4] + ")", str(homlen)[1:-1], str(homseq)[1:-1],
                                     row.FILTER, str(pr_formatted), str(sr_formatted)])
