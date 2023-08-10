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
    tsv_writer.writerow(["#POSITION1", "MANTAID", "BREAKEND", "GENES", "DEPTH", "ANNOTATIONINFO", "PR_ALT_FREQ", "SR_ALT_FREQ"])
    for row in vcf:
        if "MantaBND" in row.ID and not any(x in ["MinQUAL", "MinGQ", "MinSomaticScore",
                                                  "Ploidy", "MaxDepth", "MaxMQ0Frac", "NoPairSupport",
                                                  "SampleFT", "HomRef"] for x in row.FILTER):
            try:
                genes = row.INFO["ANN"][0].split("|")
            except KeyError:
                genes = ["NA", "NA", "NA", "NA", "NA"]
            manta_id = ":".join(row.ID.split(":")[0:2])
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
            tsv_writer.writerow([row.CHROM + ":" + str(row.POS), manta_id, str(row.ALT)[1:-1],
                                genes[3] + "(" + genes[4] + ")", row.INFO["BND_DEPTH"],
                                row.FILTER, str(pr_formatted), str(sr_formatted)])
