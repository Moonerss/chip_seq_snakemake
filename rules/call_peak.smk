rule macs2_broad_peak:
    input:
        bam1 = lambda wildcards: expand('{{BASE_FOLDER}}/bowtie2/filter/{sss}.clN.sorted.bam', sss = get_antibody_samples(wildcards)[0]),
        bam2 = lambda wildcards: expand('{{BASE_FOLDER}}/bowtie2/filter/{sss}.clN.sorted.bam', sss = get_antibody_samples(wildcards)[1]),
        gsizes = '{BASE_FOLDER}/genome/genome.fa.effectiveSize'
    output:
        expand('{{BASE_FOLDER}}/bowtie2/macs2/broad_peak/{{antibody_id}}/{{sample_id}}{end}', end = ['_peaks.xls', '_peaks.broadPeak', '_summits.bed'])
    params:
        macs2 = config['tools']['macs2'],
        broad_peak = config['extra']['macs2']['broad']
    log: '{BASE_FOLDER}/logs/macs2/broad/{antibody_id}/{sample_id}.log'
    shell:
        """
        gsize=`cat {input.gsizes}`
        {params.macs2} callpeak \
        --treatment {input.bam1} --control {input.bam2} \
        --gsize $gsize {params.broad_peak} \
        --name {wildcards.sample_id} \
        --outdir {wildcards.BASE_FOLDER}/bowtie2/macs2/broad_peak/{wildcards.antibody_id} \
        &> {log}
        """

rule macs2_narrow_peak:
    input:
        bam1 = lambda wildcards: expand('{{BASE_FOLDER}}/bowtie2/filter/{sss}.clN.sorted.bam', sss = get_antibody_samples(wildcards)[0]),
        bam2 = lambda wildcards: expand('{{BASE_FOLDER}}/bowtie2/filter/{sss}.clN.sorted.bam', sss = get_antibody_samples(wildcards)[1]),
        gsizes = '{BASE_FOLDER}/genome/genome.fa.effectiveSize'
    output:
        expand('{{BASE_FOLDER}}/bowtie2/macs2/narrow_peak/{{antibody_id}}/{{sample_id}}{end}', end = ['_peaks.xls', '_peaks.narrowPeak', '_summits.bed'])
    params:
        macs2 = config['tools']['macs2'],
        narrow_peak = config['extra']['macs2']['narrow']
    log: '{BASE_FOLDER}/logs/macs2/narrow/{antibody_id}/{sample_id}.log'
    shell:
        """
        gsize=`cat {input.gsizes}`
        {params.macs2} callpeak \
        --treatment {input.bam1} --control {input.bam2} \
        --gsize $gsize {params.narrow_peak} \
        --name {wildcards.sample_id} \
        --outdir {wildcards.BASE_FOLDER}/bowtie2/macs2/narrow_peak/{wildcards.antibody_id} \
        &> {log}
        """

