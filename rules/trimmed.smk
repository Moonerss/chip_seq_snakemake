rule trim_paired:
    input:
        fq1 = '{BASE_FOLDER}/originFASTQ/{sample_id}_R1.fastq.gz',
        fq2 = '{BASE_FOLDER}/originFASTQ/{sample_id}_R2.fastq.gz'
    output:
        out_fq1 = '{BASE_FOLDER}/trimmed/{sample_id}_R1.fastq.gz',
        out_fq2 = '{BASE_FOLDER}/trimmed/{sample_id}_R2.fastq.gz',
        trim_file1 = '{BASE_FOLDER}/trimmed/{sample_id}_R1_trimming_report.txt',
        trim_file2 = '{BASE_FOLDER}/trimmed/{sample_id}_R2_trimming_report.txt'
    params:
        trim_galore = config['tools']['trim_galore'],
        cutadapt = config['tools']['cutadapt'],
        trim_galore_extra = config['extra']['trim_galore'] if config['extra']['trim_galore'] is not None or config['extra']['trim_galore'] != '' else config['extra']['trim_galore']
    threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
    log: '{BASE_FOLDER}/logs/trimmed/raw/{sample_id}.log'
    shell:
        """
        {params.trim_galore} {params.trim_galore_extra} --cores {threads} \
        --gzip --paired --path_to_cutadapt {params.cutadapt} \
        --output_dir {wildcards.BASE_FOLDER}/trimmed {input.fq1} {input.fq2} &> {log}
        mv {wildcards.BASE_FOLDER}/trimmed/{wildcards.sample_id}_R1_val_1.fq.gz {output.out_fq1}
        mv {wildcards.BASE_FOLDER}/trimmed/{wildcards.sample_id}_R2_val_2.fq.gz {output.out_fq2}
        mv {wildcards.BASE_FOLDER}/trimmed/{wildcards.sample_id}_R1.fastq.gz_trimming_report.txt {output.trim_file1}
        mv {wildcards.BASE_FOLDER}/trimmed/{wildcards.sample_id}_R2.fastq.gz_trimming_report.txt {output.trim_file2}
        """

rule trim_singled:
    input:
        fq = '{BASE_FOLDER}/originFASTQ/{sample_id}.fastq.gz'
    output:
        out_fq = '{BASE_FOLDER}/trimmed/{sample_id}.fastq.gz',
        trim_file = '{BASE_FOLDER}/trimmed/{sample_id}_trimming_report.txt'
    params:
        trim_galore = config['tools']['trim_galore'],
        cutadapt = config['tools']['cutadapt'],
        trim_galore_extra = config['extra']['trim_galore'] if config['extra']['trim_galore'] is not None or config['extra']['trim_galore'] != '' else config['extra']['trim_galore']
    threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
    log: '{BASE_FOLDER}/logs/trimmed/raw/{sample_id}.log'
    shell:
        """
        {params.trim_galore} {params.trim_galore_extra} --cores {threads} \
        --gzip --paired --path_to_cutadapt {params.cutadapt} \
        --output_dir {wildcards.BASE_FOLDER}/trimmed &> {log}
        mv {wildcards.BASE_FOLDER}/trimmed/{wildcards.sample_id}_trimmed.fq.gz {output.out_fq}
        mv {wildcards.BASE_FOLDER}/trimmed/{wildcards.sample_id}.fastq.gz_trimming_report.txt {output.trim_file}
        """


rule trimmed_fastqc_raw_paired:
    input:
        fq1 = '{BASE_FOLDER}/trimmed/{sample_id}_R1.fastq.gz',
        fq2 = '{BASE_FOLDER}/trimmed/{sample_id}_R2.fastq.gz'
    output:
        fqc_1 = '{BASE_FOLDER}/trimmed_fastqc/zips/{sample_id}_R1_fastqc.zip',
        fqc_2 = '{BASE_FOLDER}/trimmed_fastqc/zips/{sample_id}_R2_fastqc.zip'
    params:
        fastqc = config['tools']['fastqc'],
        fastqc_extra = config['extra']['fastqc'] if config['extra']['fastqc'] is not None or config['extra']['fastqc'] != '' else config['extra']['fastqc']
    threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
    log: '{BASE_FOLDER}/logs/fastqc/trimmed/{sample_id}.log'
    shell:
        """
        {params.fastqc} {params.fastqc_extra} --threads {threads} {input.fq1} {input.fq2} \
        --outdir {wildcards.BASE_FOLDER}/trimmed_fastqc &> {log}
        mv {wildcards.BASE_FOLDER}/trimmed_fastqc/{wildcards.sample_id}_R1_fastqc.zip {wildcards.BASE_FOLDER}/trimmed_fastqc/zips
        mv {wildcards.BASE_FOLDER}/trimmed_fastqc/{wildcards.sample_id}_R2_fastqc.zip {wildcards.BASE_FOLDER}/trimmed_fastqc/zips
        """

rule trimmed_fastqc_raw_singled:
    input:
        fq1 = '{BASE_FOLDER}/trimmed/{sample_id}.fastq.gz'
    output:
        fqc_1 = '{BASE_FOLDER}/trimmed_fastqc/zips/{sample_id}_fastqc.zip'
    params:
        fastqc = config['tools']['fastqc'],
        fastqc_extra = config['extra']['fastqc'] if config['extra']['fastqc'] is not None or config['extra']['fastqc'] != '' else config['extra']['fastqc']
    threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
    log: '{BASE_FOLDER}/logs/fastqc/trimmed/{sample_id}.log'
    shell:
        """
        {params.fastqc} {params.fastqc_extra} --threads {threads} {input.fq1} \
        --outdir {wildcards.BASE_FOLDER}/trimmed_fastqc &> {log}
        mv {wildcards.BASE_FOLDER}/trimmed_fastqc/{wildcards.sample_id}_fastqc.zip {wildcards.BASE_FOLDER}/trimmed_fastqc/zips
        """