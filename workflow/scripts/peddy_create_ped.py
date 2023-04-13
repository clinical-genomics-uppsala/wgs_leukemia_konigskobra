#!/usr/bin/env python
# coding: utf-8


import sys


input_file = snakemake.input[0]


with open(input_file, 'r') as samplesheet :
    next(samplesheet)
    for lline in samplesheet :
        line=lline.strip().split("\t")
        if line[3] == "M" :
        	sex = "1"
        elif line[3] == "K" :
        	sex = "2"
        else : 
        	sex="0"
        with open("qc/peddy/" + line[0] + ".peddy.fam", 'w+') as pedfile: 
        	pedfile.write("\t".join([line[0], line[0] + "_T", "0", "0", sex, "-9"])+"\n")



