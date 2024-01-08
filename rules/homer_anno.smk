rule homer_anno_narrow:
    input:
        peak = '{BASE_FOLDER}/bowtie2/macs2/narrow_peak/{antibody_id}/{sample_id}_peaks.narrowPeak'
    output:
        anno_peak = '{BASE_FOLDER}/bowtie2/macs2/narrow_peak/{antibody_id}/{sample_id}_peaks.annotatePeaks.txt'
    params:
        homer = config['tools']['homer'],
        fa = '{BASE_FOLDER}/genome/genome.fa',
        gtf = config['reference']['gtf']
    threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
    log: '{BASE_FOLDER}/logs/homer_anno/narrow_peak/{antibody_id}/{sample_id}.log'
    shell:
        """
        cmd_annotatePeaks=`dirname {params.homer}`/annotatePeaks.pl
        $cmd_annotatePeaks {input.peak} \
        {params.fa} -gid -gtf {params.gtf} \
        -cpu {threads} > {output.anno_peak} 2> {log}
        """

rule homer_anno_broad:
    input:
        peak = '{BASE_FOLDER}/bowtie2/macs2/broad_peak/{antibody_id}/{sample_id}_peaks.broadPeak'
    output:
        anno_peak = '{BASE_FOLDER}/bowtie2/macs2/broad_peak/{antibody_id}/{sample_id}_peaks.annotatePeaks.txt'
    params:
        homer = config['tools']['homer'],
        fa = '{BASE_FOLDER}/genome/genome.fa',
        gtf = config['reference']['gtf']
    threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
    log: '{BASE_FOLDER}/logs/homer_anno/broad_peak/{antibody_id}/{sample_id}.log'
    shell:
        """
        cmd_annotatePeaks=`dirname {params.homer}`/annotatePeaks.pl
        $cmd_annotatePeaks {input.peak} \
        {params.fa} -gid -gtf {params.gtf} \
        -cpu {threads} > {output.anno_peak} 2> {log}
        """

