rule bam_coverage:
    input:
        bam_file = '{BASE_FOLDER}/bowtie2/filter/{sample_id}.clN.sorted.bam',
        effect_size = '{BASE_FOLDER}/genome/genome.fa.effectiveSize'
    output:
        out_bigwig = '{BASE_FOLDER}/bowtie2/bam_coverage/{sample_id}.bigwig'
    params:
        deeptools = config['tools']['deeptools'],
        extra_bam_coverage = config['extra']['bam_coverage']
    threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
    log: '{BASE_FOLDER}/logs/deeptools/bamCoverage/{sample_id}.log'
    shell:
        """
        bam_coverage=`dirname {params.deeptools}`/bamCoverage
        effect_size_num=`cat {input.effect_size}`
        $bam_coverage --bam {input.bam_file} \
        --outFileName {output.out_bigwig} \
        --outFileFormat bigwig \
        {params.extra_bam_coverage} \
        --numberOfProcessors {threads} \
        --effectiveGenomeSize $effect_size_num \
        --verbose &> {log}
        """