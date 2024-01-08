from snakemake.utils import validate
import pandas as pd
import os
from smart_open import open
import yaml

samples = pd.read_csv('config/samplesheet_valid.csv', sep=",", dtype = str).set_index("sample", drop=False)

antibody = set([x for x in samples.antibody.tolist() if str(x) != 'nan'])

##### wildcard constraints #####
wildcard_constraints:
    sample_id = "|".join(samples.index),
    antibody_id = "|".join(antibody)

def get_individual_fastq(wildcards):
    """Get individual raw FASTQ files from unit sheet, based on a read (end) wildcard"""
    if samples.loc[wildcards.sample_id, "paired_end"]:
        return {"fastq_1": samples.loc[wildcards.sample_id, "fastq_1"], "fastq_2":samples.loc[wildcards.sample_id, "fastq_2"]}
    else:
        return {"fastq_1": samples.loc[wildcards.sample_id, "fastq_1"]}

def get_paired(wildcards):
    if samples.loc[wildcards.sample_id, "paired_end"]:
        return True
    else:
        return False

def get_antibody_samples(wildcards):
    """Get all bam files for multiBamSummary"""
    sub_samples = samples[samples['antibody'] == wildcards.antibody_id]
    control_samples = sub_samples.control.tolist()
    sub_sample_id = sub_samples.index.tolist()
    all_samples = list(set(sub_sample_id + control_samples))
    return all_samples

def only_sample_not_input(wildcards):
    """Get sample not input"""
    if str(samples.loc[wildcards.sample_id, "control"]) != 'nan':
        return True
    else:
        return False

def get_paired_samples(wildcards):
    """Get each paired sample"""
    sample = wildcards.sample_id
    sample_input = samples.loc[wildcards.sample_id, "control"]
    return [sample, sample_input]

