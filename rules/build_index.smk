rule gtf2bed:
    input:
        gtf = config['reference']['gtf']
    output:
        bed =  '{BASE_FOLDER}/genome/genome.bed'
    params:
        gtf2bed = 'scripts/gtf2bed'
    shell:
        """
        {params.gtf2bed} {input.gtf} > {output.bed}
        """

rule bowtie2_index:
    input:
        fa = config['reference']['fasta']
    output:
        multiext("{BASE_FOLDER}/genome/bowtie2_index/genome",
                    ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2",
                    ".rev.1.bt2", ".rev.2.bt2")
    params:
        bowtie2 = config['tools']['bowtie2']
    threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
    log: '{BASE_FOLDER}/logs/bowtie2/index/bowtie2_index.log'
    shell:
        """
        {params.bowtie2}-build --threads {threads} {input.fa} \
        {wildcards.BASE_FOLDER}/genome/bowtie2_index/genome &> {log}
        """

rule fasta_fai:
    input:
        fa = config['reference']['fasta']
    output:
        fa = '{BASE_FOLDER}/genome/genome.fa',
        fa_fai = '{BASE_FOLDER}/genome/genome.fa.fai',
        fa_size = '{BASE_FOLDER}/genome/genome.fa.sizes'
    params:
        samtools = config['tools']['samtools']
    shell:
        """
        if [ ! -f {wildcards.BASE_FOLDER}/genome/genome.fa ]; then
            cp {input.fa} {output.fa}
        fi
        {params.samtools} faidx {output.fa}
        cut -f 1,2 {output.fa_fai} > {output.fa_size}
        """

if config['reference']['blacklist'] is None or config['reference']['blacklist'] != '':
    rule include_regions:
        input:
            blacklist_file = config['reference']['blacklist'],
            fa_size = '{BASE_FOLDER}/genome/genome.fa.sizes'
        output:
            include = '{BASE_FOLDER}/genome/genome.include_regions.bed'
        params:
            bedtools = config['tools']['bedtools']
        shell:
            """
            {params.bedtools} sort -i {input.blacklist_file} -g {input.fa_size} | \
            {params.bedtools} complement -i stdin -g {input.fa_size} > {output.include}
            """
else:
    rule include_regions:
        input:
            fa_size = '{BASE_FOLDER}/genome/genome.fa.sizes'
        output:
            include = '{BASE_FOLDER}/genome/genome.include_regions.bed'
        shell:
            """
            awk '{{print $1, '0' , $2}}' OFS='\\t' {input.fa_size} > {output.include}
            """

rule get_effect_size:
    input:
        fa = '{BASE_FOLDER}/genome/genome.fa'
    output:
        effect_size = '{BASE_FOLDER}/genome/genome.fa.effectiveSize'
    params:
        seqtk = config['tools']['seqtk']
    log: '{BASE_FOLDER}/logs/seqtk/effectiveSize.log'
    shell:
        """
        {params.seqtk} comp {input.fa} | \
        awk '{{tot += $3 + $4 + $5 + $6}}END{{print tot}}' > {output.effect_size} 2> {log}
        """

