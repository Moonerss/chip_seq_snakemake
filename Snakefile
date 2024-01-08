import pandas as pd
from modules.basic import CHIPipeLine
from scripts.check_samplesheet import check_samplesheet

## set config file
CONFIG_FILE = "config/config.yaml"
configfile: CONFIG_FILE

## check sample sheet
check_samplesheet(config['sample_sheet'], 'config/samplesheet_valid.csv')

## create CHIPipeLine object
SAMPLES_DF = pd.read_csv('config/samplesheet_valid.csv')
BASE = config['base']
PIPELINE = CHIPipeLine(CONFIG_FILE)

## get final outputs
def get_final_outputs(wildcards):
    files = []
    for _, row in SAMPLES_DF.iterrows():
        sample_id = row["sample"]
        antibody_id = row["antibody"]
        paired = row["paired_end"]
        files.extend(PIPELINE.output_files(sample_id, antibody_id, paired))
    return files


## get final outputs
rule all:
    input: get_final_outputs
    threads: config["threads"]


## Include all other rules files
## This order is also important as some of the functions are reused in other files.
include: 'rules/common.smk'
include: 'rules/link.smk'
include: 'rules/fastqc.smk'
include: 'rules/trimmed.smk'
include: 'rules/build_index.smk'
include: 'rules/align.smk'
include: 'rules/markdup.smk'
include: 'rules/filter.smk'
include: 'rules/bam_qc.smk'
include: 'rules/bam_to_bigwig.smk'
include: 'rules/deeptools_qc.smk'
include: 'rules/call_peak.smk'
include: 'rules/homer_anno.smk'

