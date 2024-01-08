'''
Author: 
TODO: what
LastEditors: zhaoej
LastEditTime: 2023-12-25 10:18:20
Date: 2023-12-25 10:10:50
'''
"""
CHIPipeLine (c) 2024

Author: Jeason Zhao, Qian Lab, NJMU

This deals with PipeLine class which handles all the output file name
generations.
"""

import yaml
import os
import pandas as pd
from modules.constants import *

class CHIPipeLine:
    def __init__(self, config_file: str):
        self.filename = config_file
        self._data = None
        self.tasks = []
        self.de_methods = []
        ## the last task number, this need to change when add new task
        self._last_task = 12

    @property
    def data(self) -> dict:
        if self._data is None:
            with open(self.filename) as f:
                self._data = yaml.load(f, Loader=yaml.SafeLoader)
        return self._data

    @property
    def base(self) -> str:
        return self.data["base"]

    def _generate_tasks(self):
        self.tasks = []
        if isinstance(self.data['tasks'], int):
            if self.data['tasks'] == TASK_ALL:
                self.tasks = list(range(1, self._last_task + 1))
            else:
                self.tasks.append(self.data['tasks'])
        else:
            tsk = self.data['tasks'].split(',')
            correct_tasks = []
            for t in tsk:
                try:
                    correct_tasks.append(int(t))
                except ValueError as e:
                    r_range = str(t).strip().split("-")
                    if len(r_range) <= 1:
                        raise ValueError(f"'tasks' should be integer, list of "
                                         f"integers or range. e.g. tasks: 0 or "
                                         f"tasks: 1, 2, 3 ... or 2-5") from e
                    elif len(r_range) > 2:
                        raise ValueError(f"Range should be in "
                                         f"START_TASK-END_TASK format."
                                         f" e.g. 1-5, 11-15") from e
                    else:
                        try:
                            tmp = sorted([int(x) for x in r_range])
                            correct_tasks.extend(
                                list(range(tmp[0], tmp[1] + 1)))
                        except ValueError:
                            raise ValueError(
                                f"Range should have only integers, "
                                f"you have provided '{r_range}'")
            self.tasks = sorted(list(set(correct_tasks)))
    
    def _get_macs2_file(self, sample_name = 'a', peak_type = 'broad'):
        if peak_type == 'broad':
            file_end = ['_peaks.xls', '_peaks.broadPeak', '_summits.bed']
        elif peak_type == 'narrow':
            file_end = ['_peaks.xls', '_peaks.narrowPeak', '_summits.bed']
        else:
            raise ValueError("The `peak_type` in sample.csv file must be `broad` or `narrow`")
        return [sample_name + y for y in file_end]

    def get_peak_files(self, sample_id:str, antibody_id:str):
        sample_sheet_file = self.data['sample_sheet']
        sample_sheet_data = pd.read_csv(sample_sheet_file, sep=",", dtype = str).set_index("sample", drop=False)
        sample_data = sample_sheet_data[sample_sheet_data['antibody'] == antibody_id]
        sample_data = sample_data[sample_data['sample'] == sample_id]
        if sample_data.shape[0] > 1:
            raise ValueError("There are two samples with same info in sample.csv, please check ...")
        peak_type = list(sample_data['peak_type'])[0]
        if peak_type == 'broad':
            bf_extend = 'bowtie2/macs2/broad_peak/'
        else:
            bf_extend = 'bowtie2/macs2/narrow_peak/'
        return_file = [bf_extend + x for x in self._get_macs2_file(os.path.join(antibody_id, sample_id), peak_type = peak_type)]
        return return_file

    def get_anno_peak(self, sample_id:str, antibody_id:str):
        sample_sheet_file = self.data['sample_sheet']
        sample_sheet_data = pd.read_csv(sample_sheet_file, sep=",", dtype = str).set_index("sample", drop=False)
        sample_data = sample_sheet_data[sample_sheet_data['antibody'] == antibody_id]
        sample_data = sample_data[sample_data['sample'] == sample_id]
        if sample_data.shape[0] > 1:
            raise ValueError("There are two samples with same info in sample.csv, please check ...")
        peak_type = list(sample_data['peak_type'])[0]
        if peak_type == 'broad':
            bf_extend = 'bowtie2/macs2/broad_peak'
        else:
            bf_extend = 'bowtie2/macs2/narrow_peak'
        return_file = os.path.join(bf_extend, antibody_id, sample_id + '_peaks.annotatePeaks.txt')
        return return_file

    def output_files(self, sample_id: str, antibody_id: str, is_paired: bool):
        self._generate_tasks()
        out = []
        ## start with linked file
        if is_paired:
            out.append(f"originFASTQ/{sample_id}_R1.fastq.gz")
            out.append(f"originFASTQ/{sample_id}_R2.fastq.gz")
        else:
            out.append(f"originFASTQ/{sample_id}.fastq.gz")
        for t in self.tasks:
            if t == TASK_RAW_QC:
                if is_paired:
                    out.append(f"fastqc/zips/{sample_id}_R1_fastqc.zip")
                    out.append(f"fastqc/zips/{sample_id}_R2_fastqc.zip")
                else:
                    out.append(f"fastqc/zips/{sample_id}_fastqc.zip")
            elif t == TASK_TRIMMED:
                if is_paired:
                    out.append(f"trimmed/{sample_id}_R1.fastq.gz")
                    out.append(f"trimmed/{sample_id}_R2.fastq.gz")
                    out.append(f"trimmed/{sample_id}_R1_trimming_report.txt")
                    out.append(f"trimmed/{sample_id}_R2_trimming_report.txt")
                else:
                    out.append(f"trimmed/{sample_id}.fastq.gz")
                    out.append(f"trimmed/{sample_id}_trimming_report.txt")
            elif t == TASK_TRIMMED_QC:
                if is_paired:
                    out.append(f"trimmed_fastqc/zips/{sample_id}_R1_fastqc.zip")
                    out.append(f"trimmed_fastqc/zips/{sample_id}_R2_fastqc.zip")
                else:
                    out.append(f"trimmed_fastqc/zips/{sample_id}_fastqc.zip")
            elif t == TASK_BUILD_INDEX:
                if self.data['bowtie2_index'] is None or self.data['bowtie2_index'] == '':
                    ## bowtie2 index
                    out.append(f"genome/bowtie2_index/genome.1.bt2")
                    out.append(f"genome/bowtie2_index/genome.2.bt2")
                    out.append(f"genome/bowtie2_index/genome.3.bt2")
                    out.append(f"genome/bowtie2_index/genome.4.bt2")
                    out.append(f"genome/bowtie2_index/genome.4.bt2")
                    out.append(f"genome/bowtie2_index/genome.rev.1.bt2")
                    out.append(f"genome/bowtie2_index/genome.rev.2.bt2")
                    ## static
                    out.append(f"genome/genome.fa")
                    out.append(f"genome/genome.fa.fai")
                    out.append(f"genome/genome.fa.sizes")
                    out.append(f"genome/genome.fa.effectiveSize")
                    out.append(f"genome/genome.bed")
                    ## black list
                    out.append(f"genome/genome.include_regions.bed")
                else:
                    print('You supplied the `bowtie2_index` in config file, we use it ...')
            elif t == TASK_ALIGN:
                if is_paired:
                    out.append(f"bowtie2/align/unmapped/{sample_id}.unmapped_1.fastq.gz")
                    out.append(f"bowtie2/align/unmapped/{sample_id}.unmapped_2.fastq.gz")
                else:
                    out.append(f"bowtie2/align/unmapped/{sample_id}.unmapped.fastq.gz")
                out.append(f"bowtie2/align/{sample_id}.bam")
                out.append(f"bowtie2/align/{sample_id}.sorted.bam")
                out.append(f"bowtie2/align/{sample_id}.sorted.bam.bai")
                out.append(f"bowtie2/align/samtools_stat/{sample_id}.sorted.bam.stats")
                out.append(f"bowtie2/align/samtools_stat/{sample_id}.sorted.bam.flagstat")
                out.append(f"bowtie2/align/samtools_stat/{sample_id}.sorted.bam.idxstats")
            elif t == TASK_MARKDUP:
                out.append(f"bowtie2/markdup/{sample_id}.mkD.sorted.bam")
                out.append(f"bowtie2/markdup/picard_metrics/{sample_id}.mkD.sorted.metrics.txt")
                out.append(f"bowtie2/markdup/{sample_id}.mkD.sorted.bam.bai")
                out.append(f"bowtie2/markdup/samtools_stat/{sample_id}.mkD.sorted.bam.stats")
                out.append(f"bowtie2/markdup/samtools_stat/{sample_id}.mkD.sorted.bam.flagstat")
                out.append(f"bowtie2/markdup/samtools_stat/{sample_id}.mkD.sorted.bam.idxstats")
            elif t == TASK_FILTER:
                out.append(f"bowtie2/filter/{sample_id}.flT.sorted.bam")
                out.append(f"bowtie2/filter/{sample_id}.clN.sorted.bam")
                out.append(f"bowtie2/filter/{sample_id}.clN.sorted.bam.bai")
                out.append(f"bowtie2/filter/samtools_stat/{sample_id}.clN.sorted.bam.stats")
                out.append(f"bowtie2/filter/samtools_stat/{sample_id}.clN.sorted.bam.flagstat")
                out.append(f"bowtie2/filter/samtools_stat/{sample_id}.clN.sorted.bam.idxstats")
                out.append(f"bowtie2/filter/picard_metrics/{sample_id}.CollectMultipleMetrics.alignment_summary_metrics")
                out.append(f"bowtie2/filter/picard_metrics/{sample_id}.CollectMultipleMetrics.base_distribution_by_cycle_metrics")
                out.append(f"bowtie2/filter/picard_metrics/{sample_id}.CollectMultipleMetrics.insert_size_metrics")
                out.append(f"bowtie2/filter/picard_metrics/{sample_id}.CollectMultipleMetrics.quality_by_cycle_metrics")
                out.append(f"bowtie2/filter/picard_metrics/{sample_id}.CollectMultipleMetrics.quality_distribution_metrics")
                out.append(f"bowtie2/filter/picard_metrics/pdf/{sample_id}.CollectMultipleMetrics.insert_size_histogram.pdf")
                out.append(f"bowtie2/filter/picard_metrics/pdf/{sample_id}.CollectMultipleMetrics.quality_by_cycle.pdf")
                out.append(f"bowtie2/filter/picard_metrics/pdf/{sample_id}.CollectMultipleMetrics.read_length_histogram.pdf")
                out.append(f"bowtie2/filter/picard_metrics/pdf/{sample_id}.CollectMultipleMetrics.quality_distribution.pdf")
            elif t == TASK_QC:
                out.append(f"bowtie2/phantompeakqualtools/{sample_id}.spp.pdf")
                # out.append(f"bowtie2/phantompeakqualtools/{sample_id}.spp.Rdata")
                out.append(f"bowtie2/phantompeakqualtools/{sample_id}.spp.out")
                out.append(f"bowtie2/phantompeakqualtools/{sample_id}.spp_correlation_mqc.tsv")
                out.append(f"bowtie2/phantompeakqualtools/{sample_id}.spp_nsc_mqc.tsv")
                out.append(f"bowtie2/phantompeakqualtools/{sample_id}.spp_rsc_mqc.tsv")
            elif t == TASK_BAM2BIGWIG:
                out.append(f"bowtie2/bam_coverage/{sample_id}.bigwig")
            elif t == TASK_DEEPTOOLS_QC:
                if str(antibody_id) != 'nan':
                    out.append(f"bowtie2/deeptools_qc/multiBamSummary/{antibody_id}.read_coverage.bins.npz")
                    out.append(f"bowtie2/deeptools_qc/multiBamSummary/{antibody_id}.scaling_factors.txt")
                    out.append(f"bowtie2/deeptools_qc/multiBamSummary/{antibody_id}.readCounts.tab")
                    out.append(f"bowtie2/deeptools_qc/plotPCA/{antibody_id}.PCA.read_coverage.pdf")
                    out.append(f"bowtie2/deeptools_qc/plotPCA/{antibody_id}.PCA.read_coverage.tsv")
                    out.append(f"bowtie2/deeptools_qc/plotCoverage/{antibody_id}.read_coverage.pdf")
                    out.append(f"bowtie2/deeptools_qc/plotCoverage/{antibody_id}.read_coverage.tsv")
                    out.append(f"bowtie2/deeptools_qc/plotCorrelation/{antibody_id}.pearson.heatmap.pdf")
                    out.append(f"bowtie2/deeptools_qc/plotCorrelation/{antibody_id}.spearman.heatmap.pdf")
                    out.append(f"bowtie2/deeptools_qc/plotCorrelation/{antibody_id}.pearson.heatmap.tsv")
                    out.append(f"bowtie2/deeptools_qc/plotCorrelation/{antibody_id}.spearman.heatmap.tsv")
                    out.append(f"bowtie2/deeptools_qc/plotFingerprint/{antibody_id}.Fingerprint.pdf")
                    out.append(f"bowtie2/deeptools_qc/plotFingerprint/{antibody_id}.Fingerprint.tsv")
                    out.append(f"bowtie2/deeptools_qc/bamPEFragmentSize/{antibody_id}.bamPEFragmentSize.pdf")
                    out.append(f"bowtie2/deeptools_qc/bamPEFragmentSize/{antibody_id}.bamPEFragmentSize.tsv")
                    out.append(f"bowtie2/enrichedHeatmap/computeMatrix/reference_point/{antibody_id}/{sample_id}.matrix.gz")
                    out.append(f"bowtie2/enrichedHeatmap/computeMatrix/reference_point/{antibody_id}/{sample_id}.region.bed")
                    out.append(f"bowtie2/enrichedHeatmap/computeMatrix/reference_point/{antibody_id}/{sample_id}.matrix.tsv")
                    out.append(f"bowtie2/enrichedHeatmap/computeMatrix/scale_regions/{antibody_id}/{sample_id}.matrix.gz")
                    out.append(f"bowtie2/enrichedHeatmap/computeMatrix/scale_regions/{antibody_id}/{sample_id}.region.bed")
                    out.append(f"bowtie2/enrichedHeatmap/computeMatrix/scale_regions/{antibody_id}/{sample_id}.matrix.tsv")
                    out.append(f"bowtie2/enrichedHeatmap/plotHeatmap/reference_point/{antibody_id}/{sample_id}.heatmap.pdf")
                    out.append(f"bowtie2/enrichedHeatmap/plotHeatmap/scale_regions/{antibody_id}/{sample_id}.heatmap.pdf")
            elif t == TASK_CALL_PEAK:
                if str(antibody_id) != 'nan':
                    out.extend(self.get_peak_files(sample_id, antibody_id))
            elif t == TASK_PEAK_ANNO:
                if str(antibody_id) != 'nan':
                    out.append(self.get_anno_peak(sample_id, antibody_id))
            elif t < self._last_task + 1:
                pass
            else:
                raise KeyError(f"Invalid task ID: '{t}'. Check README.md "
                               f"file to know what are the task IDs available")
        return [f"{self.base}/{x}" for x in list(set(out))]