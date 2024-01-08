rule picard_markdup:
    input:
        bam_file = '{BASE_FOLDER}/bowtie2/align/{sample_id}.sorted.bam'
    output:
        out_bam = '{BASE_FOLDER}/bowtie2/markdup/{sample_id}.mkD.sorted.bam',
        out_metrics = '{BASE_FOLDER}/bowtie2/markdup/picard_metrics/{sample_id}.mkD.sorted.metrics.txt'
    params:
        picard = config['tools']['picard']
    shell:
        """
        {params.picard} MarkDuplicates \
        --ASSUME_SORTED true --REMOVE_DUPLICATES false \
        --VALIDATION_STRINGENCY LENIENT \
        --INPUT {input.bam_file} \
        --OUTPUT {output.out_bam} \
        --METRICS_FILE {output.out_metrics}
        """

rule bam_index_mkD:
    input:
        bam_file = '{BASE_FOLDER}/bowtie2/markdup/{sample_id}.mkD.sorted.bam'
    output:
        bam_index_file = '{BASE_FOLDER}/bowtie2/markdup/{sample_id}.mkD.sorted.bam.bai'
    params:
        samtools = config['tools']['samtools']
    threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
    shell:
        """
        {params.samtools} index --threads {threads} {input.bam_file}
        """

rule bam_idxstats_mkD:
    input:
        bam_file = '{BASE_FOLDER}/bowtie2/markdup/{sample_id}.mkD.sorted.bam'
    output:
        idxstats_file = '{BASE_FOLDER}/bowtie2/markdup/samtools_stat/{sample_id}.mkD.sorted.bam.idxstats'
    params:
        samtools = config['tools']['samtools']
    threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
    shell:
        """
        {params.samtools} idxstats --threads {threads} {input.bam_file} > {output.idxstats_file}
        """

rule bam_flagstat_mkD:
    input:
        bam_file = '{BASE_FOLDER}/bowtie2/markdup/{sample_id}.mkD.sorted.bam'
    output:
        flagstat_file = '{BASE_FOLDER}/bowtie2/markdup/samtools_stat/{sample_id}.mkD.sorted.bam.flagstat'
    params:
        samtools = config['tools']['samtools']
    threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
    shell:
        """
        {params.samtools} flagstat --threads {threads} {input.bam_file} > {output.flagstat_file}
        """

rule bam_stats_mkD:
    input:
        bam_file = '{BASE_FOLDER}/bowtie2/markdup/{sample_id}.mkD.sorted.bam'
    output:
        stats_file = '{BASE_FOLDER}/bowtie2/markdup/samtools_stat/{sample_id}.mkD.sorted.bam.stats'
    params:
        samtools = config['tools']['samtools']
    threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
    shell:
        """
        {params.samtools} stats --threads {threads} {input.bam_file} > {output.stats_file}
        """