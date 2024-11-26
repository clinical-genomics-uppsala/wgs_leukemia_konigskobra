#!/bin/python3

from export_to_xlsx_create_tables import *
import xlsxwriter
import datetime
from pysam import VariantFile
import sys
import logging

logging.basicConfig(
    format="{asctime} - {levelname} - {message}",
    style="{",
    datefmt="%Y-%m-%d %H:%M",
    level=logging.INFO,
)


def convert_columns_to_letter(nr_columns):
    # Function to convert number of columns to alphabetical coordinates for xlsx-sheets
    if nr_columns < 27:
        letter = chr(nr_columns + 64)
    elif nr_columns < 703:
        i = int((nr_columns - 1) / 26)
        letter = chr(i + 64) + chr(nr_columns - (i * 26) + 64)
    else:
        logging.error(f"Nr columns has to be less than 703, does not support three letter column-index for tables {nr_columns=}")
        sys.exit()
    return letter


""" Prepping data """
logging.info(f"Prepping data, such as loading {snakemake.input.vcf}=")
sample_name = snakemake.output.xlsx.split("/")[-1].split(".manta.xlsx")[0]

bedfiles = [snakemake.input.all_bed, snakemake.input.aml_bed]

filter_flags = ["MinQUAL", "MinGQ", "MinSomaticScore", "Ploidy", "MaxDepth", "MaxMQ0Frac", "NoPairSupport", "SampleFT", "HomRef"]
manta_tables_full = create_manta_tables(snakemake.input.vcf, filter_flags)

""" Creating xlsx file """
workbook = xlsxwriter.Workbook(snakemake.output.xlsx)
logging.info(f"Creating xlsx workbook {snakemake.output.xlsx}=")
worksheet_overview = workbook.add_worksheet("Overview")
worksheet_deletions = workbook.add_worksheet("Deletions")
worksheet_insertions = workbook.add_worksheet("Insertions")
worksheet_duplications = workbook.add_worksheet("Duplications")
worksheet_bnd = workbook.add_worksheet("Translocations")

format_heading = workbook.add_format({"bold": True, "font_size": 18})
format_bold = workbook.add_format({"bold": True, "text_wrap": True})

# overview
logging.debug(f"Creating Overview sheet")
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
worksheet_overview.write_url(11, 0, "internal: 'Translocations'!A1", string="Manta Translocations or breakpoints")
i_overview = 12

worksheet_overview.write_url(13, 0, "internal:'Translocations AML'!A1", string="Manta Translocations in AML genes")

worksheet_overview.write(17, 0, "ALL bedfile: " + snakemake.input.all_bed)
worksheet_overview.write(18, 0, "AML bedfile: " + snakemake.input.aml_bed)

worksheet_overview.write(22, 0, "Only calls NOT containing the following annotation are included: " + ", ".join(filter_flags))

# del
logging.debug(f"Creating Deletions sheet")
worksheet_deletions.set_column("B:C", 12)
worksheet_deletions.set_column("E:E", 12)
worksheet_deletions.write("A1", "Deletions found by Manta", format_heading)
worksheet_deletions.write("A3", "Sample: " + str(sample_name))
worksheet_deletions.write("A5", "Only calls NOT containing the following annotation are included: " + ", ".join(filter_flags))
worksheet_deletions.write("A6", "Calls have to be longer than 100 bp to be included.")

i = 8
column_end = ":" + convert_columns_to_letter(len(manta_tables_full["del"]["headers"]))
if len(manta_tables_full["del"]["data"]) > 0:
    table_area = "A" + str(i) + column_end + str(len(manta_tables_full["del"]["data"]) + i)
else:
    table_area = "A" + str(i) + column_end + str(i + 1)

worksheet_deletions.add_table(
    table_area,
    {"columns": manta_tables_full["del"]["headers"], "data": manta_tables_full["del"]["data"], "style": "Table Style Light 1"},
)

# ins
logging.debug(f"Creating Insertions sheet")
worksheet_insertions.set_column("B:B", 12)
worksheet_insertions.set_column("F:F", 12)
worksheet_insertions.write("A1", "Insertions found by Manta", format_heading)
worksheet_insertions.write("A3", "Sample: " + str(sample_name))
worksheet_insertions.write("A5", "Only calls NOT containing the following annotation are included: " + ", ".join(filter_flags))

i = 7
column_end = ":" + convert_columns_to_letter(len(manta_tables_full["ins"]["headers"]))
if len(manta_tables_full["ins"]["data"]) > 0:
    table_area = "A" + str(i) + column_end + str(len(manta_tables_full["ins"]["data"]) + i)
else:
    table_area = "A" + str(i) + column_end + str(i + 1)


worksheet_insertions.add_table(
    table_area,
    {"columns": manta_tables_full["ins"]["headers"], "data": manta_tables_full["ins"]["data"], "style": "Table Style Light 1"},
)

# dup
logging.debug(f"Creating Duplications sheet")
worksheet_duplications.set_column("B:C", 12)
worksheet_duplications.set_column("E:E", 12)
worksheet_duplications.write("A1", "Duplications found by Manta", format_heading)
worksheet_duplications.write("A3", "Sample: " + str(sample_name))
worksheet_duplications.write("A5", "Only calls NOT containing the following annotation are included: " + ", ".join(filter_flags))

column_end = ":" + convert_columns_to_letter(len(manta_tables_full["dup"]["headers"]))
i = 7
if len(manta_tables_full["dup"]["data"]) > 0:
    table_area = "A" + str(i) + column_end + str(len(manta_tables_full["dup"]["data"]) + i)
else:
    table_area = "A" + str(i) + column_end + str(i + 1)

worksheet_duplications.add_table(
    table_area,
    {"columns": manta_tables_full["dup"]["headers"], "data": manta_tables_full["dup"]["data"], "style": "Table Style Light 1"},
)

# bnd
logging.debug(f"Creating bnd sheet")
worksheet_bnd.set_column("B:B", 12)
worksheet_bnd.set_column("C:D", 15)
worksheet_bnd.write("A1", "Translocations found by Manta", format_heading)
worksheet_bnd.write("A3", "Sample: " + str(sample_name))
worksheet_bnd.write("A5", "Only calls NOT containing the following annotation are included: " + ", ".join(filter_flags))

i = 7
column_end = ":" + convert_columns_to_letter(len(manta_tables_full["bnd"]["headers"]))
if len(manta_tables_full["bnd"]["data"]) > 0:
    table_area = "A" + str(i) + column_end + str(len(manta_tables_full["bnd"]["data"]) + i)
else:
    table_area = "A" + str(i) + column_end + str(i + 1)

worksheet_bnd.add_table(
    table_area,
    {"columns": manta_tables_full["bnd"]["headers"], "data": manta_tables_full["bnd"]["data"], "style": "Table Style Light 1"},
)

# Each panel
for vcf in snakemake.input.vcfs_bed:
    panel = vcf.split(".")[-3]
    logging.debug(f"Creating {panel} sheet")
    panel_tables = create_manta_tables(vcf, filter_flags)
    worksheet = workbook.add_worksheet("Translocations " + panel.upper())
    worksheet_overview.write_url(
        i_overview,
        0,
        "internal:'Translocations " + panel.upper() + "'!A1",
        string="Manta Translocations in " + panel.upper() + " genes",
    )
    i_overview += 1
    worksheet.set_column("B:B", 12)
    worksheet.set_column("C:D", 15)
    worksheet.write("A1", "Translocations in " + panel.upper() + " genes", format_heading)
    worksheet.write("A3", "Sample: " + str(sample_name))
    worksheet.write("A5", "Only calls NOT containing the following annotation are included: " + ", ".join(filter_flags))

    i = 7
    column_end = ":" + convert_columns_to_letter(len(panel_tables["bnd"]["headers"]))
    if len(panel_tables["bnd"]["data"]) > 0:
        table_area = "A" + str(i) + column_end + str(len(panel_tables["bnd"]["data"]) + i)
    else:
        table_area = "A" + str(i) + column_end + str(i + 1)

    worksheet.add_table(
        table_area,
        {"columns": panel_tables["bnd"]["headers"], "data": panel_tables["bnd"]["data"], "style": "Table Style Light 1"},
    )

workbook.close()
logging.info("All done")
