rule link_fastq_paired:
    input:
        unpack(get_individual_fastq)
    output:
        fq1 = '{BASE_FOLDER}/originFASTQ/{sample_id}_R1.fastq.gz',
        fq2 = '{BASE_FOLDER}/originFASTQ/{sample_id}_R2.fastq.gz'
    shell:
        """
        ln -s {input.fastq_1} {output.fq1}
        ln -s {input.fastq_2} {output.fq2}
        """

rule link_fastq_singled:
    input:
        unpack(get_individual_fastq)
    output:
        fq = '{BASE_FOLDER}/originFASTQ/{sample_id}.fastq.gz'
    shell:
        """
        ln -s {input.fastq_1} {output.fq}
        """