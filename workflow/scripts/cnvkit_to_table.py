#!/bin/python3.6

import sys
import xlsxwriter

# import bedfile for CNA genes from bedfile.
bedtable = []
bed_header = ['chr', 'start', 'end', 'annotation']
with open(snakemake.input.gene_interest) as bedfile:
    #next(bedfile) if header
    for line in bedfile:
        bedtable.append(line.strip().split('\t'))

relevant_cnvs = []
relevant_cnvs_header = ['Chromosome', 'Start', 'End', 'Log2', 'BAF', 'Depth', 'Probes', 'Weight']
with open(snakemake.input.cns, 'r+') as cnsfile:
    cns_header = next(cnsfile).rstrip().split("\t")
    for cnv in cnsfile:
        cnv_chr = cnv[cns_header.index('chromosome')]
        cnv_start = int(cnv[cns_header.index('start')])
        cnv_end = int(cnv[cns_header.index('end')])
        if (cnv_end - cnv_start) >= 100000:
            outline = [cnv_chr, cnv_start, cnv_end, cnv[cns_header.index('log2')], cnv[cns_header.index('baf')],
                       cnv[cns_header.index('depth')], cnv[cns_header.index('probes')], cnv[cns_header.index('weight')]]
            relevant_cnvs.append(outline)
        else:
            for bedline in bedtable:
                if (cnv_chr == bedline[bed_header.index('chr')]):
                    bed_start = int(bedline[bed_header.index('start')])
                    bed_end = int(bedline[bed_header.index('end')])
                    if (cnv_start >= bed_start and cnv_end <= bed_end) or
                       (cnv_start >= bed_start and cnv_start <= bed_end ) or
                       (cnv_end >= bed_start and cnv_end <= bed_end) or
                       (cnv_start <= bed_start and cnv_end >= bed_end):
                            outline = [cnv_chr, cnv_start, cnv_end, cnv[cns_header.index('log2')], cnv[cns_header.index('baf')],
                                       cnv[cns_header.index('depth')], cnv[cns_header.index('probes')], cnv[cns_header.index('weight')]]
                            relevant_cnvs.append(outline)

''' Creating xlsx file '''
sample = str(snakemake.input.cns).split("/")[-1].split("_")[0]
workbook = xlsxwriter.Workbook(snakemake.output[0])
# add sheet with png?
worksheet_calls = workbook.add_worksheet('Calls')
heading_format = workbook.add_format({'bold': True, 'font_size': 18})
tablehead_format = workbook.add_format({'bold': True, 'text_wrap': True})

worksheet_calls.write('A1', 'CNVkit calls', heading_format)
worksheet_calls.write('A3', 'Sample: '+str(sample))
worksheet_calls.write('A5', 'Calls larger than 100 kb or in CNA bedfile included')
# Add link to bedfile

worksheet_calls.write('A7', relevant_cnvs_header, tablehead_format)
row = 7
col = 0
for line in relevant_cnvs:
    worksheet_calls.write_row(row, col, line)
    row += 1
