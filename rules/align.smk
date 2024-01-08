
if get_paired:
    rule bowtie2_align_paired:
        input:
           fq1 = '{BASE_FOLDER}/trimmed/{sample_id}_R1.fastq.gz',
           fq2 = '{BASE_FOLDER}/trimmed/{sample_id}_R2.fastq.gz'
        output:
            unmap_1 = '{BASE_FOLDER}/bowtie2/align/unmapped/{sample_id}.unmapped_1.fastq.gz',
            unmap_2 = '{BASE_FOLDER}/bowtie2/align/unmapped/{sample_id}.unmapped_2.fastq.gz',
            bam_file = '{BASE_FOLDER}/bowtie2/align/{sample_id}.bam'
        params:
            bowtie2 = config['tools']['bowtie2'],
            samtools = config['tools']['samtools']
        threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
        log: '{BASE_FOLDER}/logs/bowtie2/align/{sample_id}.log'
        shell:
            """
            {params.bowtie2} -x {wildcards.BASE_FOLDER}/genome/bowtie2_index/genome \
            -1 {input.fq1} -2 {input.fq2} --threads {threads} \
            --un-conc-gz {wildcards.BASE_FOLDER}/bowtie2/align/{wildcards.sample_id}.unmapped.fastq.gz 2> {log} | \
            {params.samtools} view --threads {threads} -o {output.bam_file}
            mv {wildcards.BASE_FOLDER}/bowtie2/align/{wildcards.sample_id}.unmapped.fastq.1.gz {output.unmap_1}
            mv {wildcards.BASE_FOLDER}/bowtie2/align/{wildcards.sample_id}.unmapped.fastq.2.gz {output.unmap_2}
            """
else:
    rule bowtie2_align_singled:
        input:
           fq = '{BASE_FOLDER}/trimmed/{sample_id}.fastq.gz'
        output:
            unmap = '{BASE_FOLDER}/bowtie2/align/unmapped/{sample_id}.unmapped.fastq.gz',
            bam_file = '{BASE_FOLDER}/bowtie2/align/{sample_id}.bam'
        params:
            bowtie2 = config['tools']['bowtie2'],
            samtools = config['tools']['samtools']
        threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
        log: '{BASE_FOLDER}/logs/bowtie2/align/{sample_id}.log'
        shell:
            """
            {params.bowtie2} -x {wildcards.BASE_FOLDER}/genome/bowtie2_index/genome \
            -U {input.fq} --threads {threads} \
            --un-conc-gz {wildcards.BASE_FOLDER}/bowtie2/align/{wildcards.sample_id}.unmapped.fastq.gz 2> {log} | \
            {params.samtools} view --threads {threads} -o {output.bam_file}
            mv {wildcards.BASE_FOLDER}/bowtie2/align/{wildcards.sample_id}.unmapped.fastq.gz {output.unmap}
            """

rule sort_bam:
    input:
        bam_file = '{BASE_FOLDER}/bowtie2/align/{sample_id}.bam'
    output:
        sort_bam_file = '{BASE_FOLDER}/bowtie2/align/{sample_id}.sorted.bam'
    params:
        samtools = config['tools']['samtools']
    threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
    shell:
        """
        {params.samtools} sort --threads {threads} \
        -o {output.sort_bam_file} -T {wildcards.sample_id}.sorted {input.bam_file}
        """

rule bam_index:
    input:
        bam_file = '{BASE_FOLDER}/bowtie2/align/{sample_id}.sorted.bam'
    output:
        bam_index_file = '{BASE_FOLDER}/bowtie2/align/{sample_id}.sorted.bam.bai'
    params:
        samtools = config['tools']['samtools']
    threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
    shell:
        """
        {params.samtools} index --threads {threads} {input.bam_file}
        """

rule bam_stats:
    input:
        bam_file = '{BASE_FOLDER}/bowtie2/align/{sample_id}.sorted.bam'
    output:
        stats_file = '{BASE_FOLDER}/bowtie2/align/samtools_stat/{sample_id}.sorted.bam.stats'
    params:
        samtools = config['tools']['samtools']
    threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
    shell:
        """
        {params.samtools} stats --threads {threads} {input.bam_file} > {output.stats_file}
        """

rule bam_flagstat:
    input:
        bam_file = '{BASE_FOLDER}/bowtie2/align/{sample_id}.sorted.bam'
    output:
        flagstat_file = '{BASE_FOLDER}/bowtie2/align/samtools_stat/{sample_id}.sorted.bam.flagstat'
    params:
        samtools = config['tools']['samtools']
    threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
    shell:
        """
        {params.samtools} flagstat --threads {threads} {input.bam_file} > {output.flagstat_file}
        """

rule bam_idxstats:
    input:
        bam_file = '{BASE_FOLDER}/bowtie2/align/{sample_id}.sorted.bam'
    output:
        idxstats_file = '{BASE_FOLDER}/bowtie2/align/samtools_stat/{sample_id}.sorted.bam.idxstats'
    params:
        samtools = config['tools']['samtools']
    threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
    shell:
        """
        {params.samtools} idxstats --threads {threads} {input.bam_file} > {output.idxstats_file}
        """
