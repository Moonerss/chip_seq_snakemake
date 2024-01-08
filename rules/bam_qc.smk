rule phantompeakqualtools:
    input:
        bam_file = '{BASE_FOLDER}/bowtie2/filter/{sample_id}.clN.sorted.bam'
    output:
        pdf_file = '{BASE_FOLDER}/bowtie2/phantompeakqualtools/{sample_id}.spp.pdf',
        data_file = temp('{BASE_FOLDER}/bowtie2/phantompeakqualtools/{sample_id}.spp.Rdata'),
        out_file = '{BASE_FOLDER}/bowtie2/phantompeakqualtools/{sample_id}.spp.out'
    params:
        Rscript = config['tools']['Rscript'],
        samtools = config['tools']['samtools'],
        run_spp = 'scripts/run_spp_modified.R'
    threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
    log: '{BASE_FOLDER}/logs/phantompeakqualtools/{sample_id}.log'
    shell:
        """
        {params.Rscript} {params.run_spp} \
        -c={input.bam_file} \
        -savp={output.pdf_file} -savd={output.data_file} -out={output.out_file} \
        -p={threads} -samtoolsPath={params.samtools} &> {log}
        """

rule write_spp_multiqc:
    input:
        data_file = '{BASE_FOLDER}/bowtie2/phantompeakqualtools/{sample_id}.spp.Rdata',
        out_file = '{BASE_FOLDER}/bowtie2/phantompeakqualtools/{sample_id}.spp.out',
        spp_correlation_header = 'assests/spp_correlation_header.txt',
        spp_nsc_header = 'assests/spp_nsc_header.txt',
        spp_rsc_header = 'assests/spp_rsc_header.txt'
    output:
        spp_correlation_mqc = '{BASE_FOLDER}/bowtie2/phantompeakqualtools/{sample_id}.spp_correlation_mqc.tsv',
        spp_nsc_mqc = '{BASE_FOLDER}/bowtie2/phantompeakqualtools/{sample_id}.spp_nsc_mqc.tsv',
        spp_rsc_mqc = '{BASE_FOLDER}/bowtie2/phantompeakqualtools/{sample_id}.spp_rsc_mqc.tsv'
    params:
        Rscript = config['tools']['Rscript']
    shell:
        """
        cp {input.spp_correlation_header} {output.spp_correlation_mqc}
        
        {params.Rscript} --max-ppsize=500000 \
        -e "load('{input.data_file}'); write.table(crosscorr\$cross.correlation, file='{output.spp_correlation_mqc}', sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE,append=TRUE)"

        awk -v OFS='\\t' '{{print "{wildcards.sample_id}", $9}}'  {input.out_file} | cat {input.spp_nsc_header} - > {output.spp_nsc_mqc}
        awk -v OFS='\\t' '{{print "{wildcards.sample_id}", $10}}'  {input.out_file} | cat {input.spp_rsc_header} - > {output.spp_rsc_mqc}
        """


# cp spp_correlation_header.txt H33_DMSO_1.spp_correlation_mqc.tsv
# Rscript --max-ppsize=500000 -e "load('H33_DMSO_1.spp.Rdata'); write.table(crosscorr\$cross.correlation, file=\"H33_DMSO_1.spp_correlation_mqc.tsv\", sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE,append=TRUE)"
# 
# awk -v OFS='	' '{print "H33_DMSO_1", $9}'  H33_DMSO_1.spp.out | cat spp_nsc_header.txt - > H33_DMSO_1.spp_nsc_mqc.tsv
# awk -v OFS='	' '{print "H33_DMSO_1", $10}' H33_DMSO_1.spp.out | cat spp_rsc_header.txt - > H33_DMSO_1.spp_rsc_mqc.tsv
