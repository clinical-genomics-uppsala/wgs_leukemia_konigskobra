
## [export_to_xlsx](url_to_tool)
Introduction to export_to_xlsx

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__export__export_to_xlsx#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__export__export_to_xlsx#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__export_to_xlsx#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__export_to_xlsx#

## [annotate_normal_ratio]
Add normal ratio to variants for TN samples 

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__annotate__annotate_normal_ratio#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__annotate__annotate_normal_ratio#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__annotate_normal_ratio#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__annotate_normal_ratio#

## [fix_svdb_header]
With svdb 2.6.0 there is a spelling mistake in the vcf header, using sed to replace _SAMPLE, with _SAMPLES, to match format used in variants

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__fix_svdb_header__fix_svdb_header#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__fix_svdb_header__fix_svdb_header#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__fix_svdb_header#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__fix_svdb_header#
