#!/bin/python3

from export_to_xlsx_create_tables import *
import xlsxwriter
import datetime
from pysam import VariantFile
import yaml


def index_vep(variantfile):
    csq_index = []
    for x in variantfile.header.records:
        if "CSQ" in str(x):
            csq_index = str(x).split("Format: ")[1].strip().strip('">').split("|")
    return csq_index


""" Prepping input data """
bedfiles = {}
bedfiles["all"] = snakemake.input.all_bed
bedfiles["aml"] = snakemake.input.aml_bed
bedfiles["tm"] = snakemake.input.tm_bed

vcfs = {}
vcfs["all"] = [x for x in snakemake.input.vcfs if "include.all" in x][0]
vcfs["aml"] = [x for x in snakemake.input.vcfs if "include.aml" in x][0]
vcfs["tm"] = [x for x in snakemake.input.vcfs if "include.tm" in x][0]
subsections = ["all", "aml", "tm"]

sample_name = snakemake.output.xlsx.split("/")[-1].split(".snvs.xlsx")[0]

snv_tables = {}
for subsection in subsections:
    snv_tables[subsection] = create_snv_table(vcfs[subsection])
    # pindel_table = create_pindel_table(pindel)


# Adding bedfiles
bed_tables = {}
for subsection in subsections:
    bed_table = []
    with open(bedfiles[subsection], "r") as file:
        for line in file:
            bed_table.append(line.strip().split("\t"))
    bed_tables[subsection] = bed_table


# Add filters info
filters = []
with open(snakemake.params.filterfile, "r") as filter_file:
    filters_dict = yaml.safe_load(filter_file)

for key, items in filters_dict["filters"].items():
    filters.append(items["soft_filter_flag"] + ": " + items["description"])


""" Creating xlsx file """
workbook = xlsxwriter.Workbook(snakemake.output.xlsx)
worksheet_overview = workbook.add_worksheet("Overview")

format_heading = workbook.add_format({"bold": True, "font_size": 18})
format_bold = workbook.add_format({"bold": True, "text_wrap": True})
format_orange = workbook.add_format({"bg_color": "#ffd280"})

# Overview sheet
worksheet_overview.write(0, 0, sample_name, format_heading)
worksheet_overview.write(1, 0, "Processing date: " + datetime.datetime.now().strftime("%d %B, %Y"))

worksheet_overview.write(4, 0, "Created by: ")
worksheet_overview.write(4, 4, "Valid from: ")
worksheet_overview.write(5, 0, "Signed by: ")
worksheet_overview.write(5, 4, "Document nr: ")


worksheet_overview.write(7, 0, "Sheets:", format_bold)
worksheet_overview.write_url(8, 0, "internal:'ALL'!A1", string="Variants in ALL genes")
worksheet_overview.write_url(9, 0, "internal:'AML'!A1", string="Variants in AML genes")
worksheet_overview.write_url(10, 0, "internal:'TM'!A1", string="Variants in TM exons")
worksheet_overview.write_url(11, 0, "internal:'ALL bedfile'!A1", string="Gene regions included in ALL bedfile")
worksheet_overview.write_url(12, 0, "internal:'AML bedfile'!A1", string="Gene regions included in AML bedfile")
worksheet_overview.write_url(13, 0, "internal:'TM bedfile'!A1", string="Gene regions included in TM bedfile")

worksheet_overview.write(16, 0, "ALL bedfile: " + bedfiles["all"])
worksheet_overview.write(17, 0, "AML bedfile: " + bedfiles["aml"])
worksheet_overview.write(18, 0, "TM exons bedfile: " + bedfiles["tm"])

# Add snv variants sheets
for sheet in subsections:
    data_table = snv_tables[sheet]
    worksheet = workbook.add_worksheet(sheet.upper())
    worksheet.set_column("B:B", 12)
    worksheet.set_column("E:E", 10)
    worksheet.set_column("K:K", 15)
    worksheet.write("A1", "Variants in " + sheet.upper() + " regions", format_heading)
    worksheet.write("A3", "Sample: " + str(sample_name))
    worksheet.write("A4", "Databases used: " + data_table["vep_line"])

    worksheet.write("A6", "Filters: ", format_orange)
    for i, filter_txt in enumerate(filters):
        i += 7
        worksheet.write("B" + str(i), filter_txt, format_orange)

    i = len(filters) + 7 + 1
    worksheet.write_rich_string("A" + str(i), "Only variants with filter-flag ", format_bold, "PASS", " shown by default.")
    i += 1
    worksheet.write(
        "A" + str(i),
        "To see all variants; put marker on header row, then click on 'Standard Filter' and remove any values. "
        + "You can then use the drop-downs in the header row to filter to your liking.",
    )

    i += 2
    if len(data_table["data"]) > 0:
        table_area = "A" + str(i) + ":T" + str(len(data_table["data"]) + i)
        table_area_data = "A" + str(i + 1) + ":T" + str(len(data_table["data"]) + i)
    else:
        table_area = "A" + str(i) + ":T" + str(i + 1)
        table_area_data = "A" + str(i + 1) + ":T" + str(i + 1)

    worksheet.add_table(table_area, {"columns": data_table["headers"], "style": "Table Style Light 1"})

    cond_formula = "=LEFT($A" + str(i + 1) + ', 4)<>"PASS"'
    worksheet.conditional_format(table_area_data, {"type": "formula", "criteria": cond_formula, "format": format_orange})

    worksheet.autofilter(table_area)
    worksheet.filter_column("A", "Filter != PASS")
    # worksheet.filter_column("I", "AF >= 0.05")
    for row_data in data_table["data"]:
        if row_data[0] == "PASS":
            pass
        else:
            worksheet.set_row(i, options={"hidden": True})
        worksheet.write_row(i, 0, row_data)
        i += 1


# Add bedfile sheets
for sheet in subsections:
    bed_data = bed_tables[sheet]
    worksheet = workbook.add_worksheet(sheet.upper() + " bedfile")
    worksheet.set_column("B:C", 10)
    worksheet.write("A1", "Bedfile used to filter out " + sheet.upper() + " regions", format_heading)
    tableheading = ["Chr", "Start", "End", "Region"]
    heading_list = [{"header": a} for a in tableheading]
    worksheet.add_table(
        "A4:D" + str(len(bed_data) + 5), {"data": bed_data, "columns": heading_list, "style": "Table Style Light 1"}
    )

workbook.close()
