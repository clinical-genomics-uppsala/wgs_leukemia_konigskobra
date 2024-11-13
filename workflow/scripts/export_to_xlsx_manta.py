#!/bin/python3

import xlsxwriter
import datetime
from pysam import VariantFile

""" Prepping data """
sample_name = snakemake.output.xlsx.split("/")[-1].split(".manta.xlsx")[0]

bedfiles = [snakemake.input.all_bed, snakemake.input.aml_bed]

filter_flags = ["MinQUAL", "MinGQ", "MinSomaticScore", "Ploidy", "MaxDepth", "MaxMQ0Frac", "NoPairSupport", "SampleFT", "HomRef"]
manta_tables_full = create_manta_tables(snakemake.input.vcf, filter_flags)


""" Creating xlsx file """
workbook = xlsxwriter.Workbook(snakemake.output.xlsx)
worksheet_overview = workbook.add_worksheet("Overview")
worksheet_deletions = workbook.add_worksheet("Deletions")
worksheet_insertions = workbook.add_worksheet("Insertions")
worksheet_duplications = workbook.add_worksheet("Duplications")
worksheet_bnd = workbook.add_worksheet("Translocations")
worksheet_bnd_all = workbook.add_worksheet("Translocations ALL")
worksheet_bnd_aml = workbook.add_worksheet("Translocations AML")
worksheet_bnd_tm = workbook.add_worksheet("Translocations TM")

format_heading = workbook.add_format({"bold": True, "font_size": 18})
format_bold = workbook.add_format({"bold": True, "text_wrap": True})

#overview
worksheet_overview.write(0, 0, sample_name, format_heading)
worksheet_overview.write(1, 0, "Processing date: " + datetime.datetime.now().strftime("%d %B, %Y"))

worksheet_overview.write(4, 0, "Created by: ")
worksheet_overview.write(4, 4, "Valid from: ")
worksheet_overview.write(5, 0, "Signed by: ")
worksheet_overview.write(5, 4, "Document nr: ")

worksheet_overview.write(7, 0, "Sheets:", format_bold)
worksheet_overview.write_url(8, 0, "internal:'Deletions'!A1", string="Manta Deletions")
worksheet_overview.write_url(9, 0, "internal:'Insertions'!A1", string="Manta Insertions")
worksheet_overview.write_url(10, 0, "internal:'Duplications'!A1", string="Manta Duplications")
worksheet_overview.write_url(11,0, "internal: 'Translocations'!A1", string="Manta Translocations or breakpoints")
i_overview = 12

worksheet_overview.write_url(13, 0, "internal:'Translocations AML'!A1", string="Manta Translocations in AML genes")

worksheet_overview.write(17, 0, "ALL bedfile: " + snakemake.input.all_bed)
worksheet_overview.write(18, 0, "AML bedfile: " + snakemake.input.aml_bed)

worksheet_overview(22, 0, "Only calls NOT containing the following annotation are included: " + ", ".join(filter_flags))

# del
worksheet_deletions.set_column("B:B", 12)
worksheet_deletions.write("A1", "Deletions found by Manta", format_heading)
worksheet_deletions.write("A3", "Sample: " + str(sample_name))
worksheet_deletions.write("A5", "Only variants longer than 100 bp included.")

i = 7
if len(manta_tables_full["del"]["data"]) > 0:
    table_area = "A" + str(i) + ":K" + str(len(manta_tables_full["del"]["data"]) + i) # for TN borde det vara separata?
else:
    table_area = "A" + str(i) + ":K" + str(i + 1) # for TN borde det vara separata?

worksheet_deletions.add_table(table_area, {"columns": manta_tables_full["del"]["headers"], "data": manta_tables_full["del"]["data"], "style": "Table Style Light 1"})

# ins
worksheet_insertions.set_column("B:B", 12)
worksheet_insertions.write("A1", "Insertions found by Manta", format_heading)
worksheet_insertions.write("A3", "Sample: " + str(sample_name))

i = 5
if len(manta_tables_full["ins"]["data"]) > 0:
    table_area = "A" + str(i) + ":N" + str(len(manta_tables_full["ins"]["data"]) + i) # for TN borde det vara separata?
else:
    table_area = "A" + str(i) + ":N" + str(i + 1) # for TN borde det vara separata?

worksheet_deletions.add_table(table_area, {"columns": manta_tables_full["ins"]["headers"], "data": manta_tables_full["ins"]["data"], "style": "Table Style Light 1"})

# dup
worksheet_duplications.set_column("B:C", 12)
worksheet_duplications.write("A1", "Duplications found by Manta", format_heading)
worksheet_duplications.write("A3", "Sample: " + str(sample_name))

i = 5
if len(manta_tables_full["dup"]["data"]) > 0:
    table_area = "A" + str(i) + ":M" + str(len(manta_tables_full["dup"]["data"]) + i) # for TN borde det vara separata?
else:
    table_area = "A" + str(i) + ":M" + str(i + 1) # for TN borde det vara separata?

worksheet_deletions.add_table(table_area, {"columns": manta_tables_full["dup"]["headers"], "data": manta_tables_full["dup"]["data"], "style": "Table Style Light 1"})

# bnd
worksheet_bnd.set_column("B:B", 12)
worksheet_bnd.write("A1", "Translocations found by Manta", format_heading)
worksheet_bnd.write("A3", "Sample: " + str(sample_name))

i = 5
if len(manta_tables_full["bnd"]["data"]) > 0:
    table_area = "A" + str(i) + ":K" + str(len(manta_tables_full["bnd"]["data"]) + i) # for TN borde det vara separata?
else:
    table_area = "A" + str(i) + ":K" + str(i + 1) # for TN borde det vara separata?

worksheet_deletions.add_table(table_area, {"columns": manta_tables_full["bnd"]["headers"], "data": manta_tables_full["bnd"]["data"], "style": "Table Style Light 1"})

# Each panel
for vcf in snakemake.input.vcfs_bed:
    panel = vcf.split(".")[-2]
    panel_tables = create_manta_tables(vcf, filter_flags)
    worksheet = workbook.add_worksheet("Translocations " + panel.upper())
    worksheet_overview.write_url(i_overview, 0, "internal:'Translocations " + panel.upper() +"'!A1", string="Manta Translocations in " + panel.upper() + " genes")
    i_overview += 1
    worksheet.set_column("B:B", 12)
    worksheet.write("A1", "Translocations in " + panel.upper() + " genes", format_heading)
    worksheet.write("A3", "Sample: " + str(sample_name))

    i = 5
    if len(panel_tables["bnd"]["data"]) > 0:
        table_area = "A" + str(i) + ":K" + str(len(panel_tables["bnd"]["data"]) + i) # for TN borde det vara separata?
    else:
        table_area = "A" + str(i) + ":K" + str(i + 1) # for TN borde det vara separata?    

    worksheet_deletions.add_table(table_area, {"columns": panel_tables["bnd"]["headers"], "data": panel_tables["bnd"]["data"], "style": "Table Style Light 1"})

workbook.close()