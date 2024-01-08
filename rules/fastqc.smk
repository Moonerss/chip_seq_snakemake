rule fastqc_raw_paired:
    input:
        fq1 = '{BASE_FOLDER}/originFASTQ/{sample_id}_R1.fastq.gz',
        fq2 = '{BASE_FOLDER}/originFASTQ/{sample_id}_R2.fastq.gz'
    output:
        fqc_1 = '{BASE_FOLDER}/fastqc/zips/{sample_id}_R1_fastqc.zip',
        fqc_2 = '{BASE_FOLDER}/fastqc/zips/{sample_id}_R2_fastqc.zip'
    params:
        fastqc = config['tools']['fastqc'],
        fastqc_extra = config['extra']['fastqc'] if config['extra']['fastqc'] is not None or config['extra']['fastqc'] != '' else config['extra']['fastqc']
    threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
    log: '{BASE_FOLDER}/logs/fastqc/raw/{sample_id}.log'
    shell:
        """
        {params.fastqc} {params.fastqc_extra} --threads {threads} {input.fq1} {input.fq2} \
        --outdir {wildcards.BASE_FOLDER}/fastqc &> {log}
        mv {wildcards.BASE_FOLDER}/fastqc/{wildcards.sample_id}_R1_fastqc.zip {wildcards.BASE_FOLDER}/fastqc/zips
        mv {wildcards.BASE_FOLDER}/fastqc/{wildcards.sample_id}_R2_fastqc.zip {wildcards.BASE_FOLDER}/fastqc/zips
        """

rule fastqc_raw_singled:
    input:
        fq1 = '{BASE_FOLDER}/originFASTQ/{sample_id}.fastq.gz'
    output:
        fqc_1 = '{BASE_FOLDER}/fastqc/zips/{sample_id}_fastqc.zip'
    params:
        fastqc = config['tools']['fastqc'],
        fastqc_extra = config['extra']['fastqc'] if config['extra']['fastqc'] is not None or config['extra']['fastqc'] != '' else config['extra']['fastqc']
    threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
    log: '{BASE_FOLDER}/logs/fastqc/raw/{sample_id}.log'
    shell:
        """
        {params.fastqc} {params.fastqc_extra} --threads {threads} {input.fq1} \
        --outdir {wildcards.BASE_FOLDER}/fastqc &> {log}
        mv {wildcards.BASE_FOLDER}/fastqc/{wildcards.sample_id}_fastqc.zip {wildcards.BASE_FOLDER}/fastqc/zips
        """