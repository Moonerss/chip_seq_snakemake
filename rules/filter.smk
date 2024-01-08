if get_paired:
    rule filter_bam_paired:
        input:
            bam_file = '{BASE_FOLDER}/bowtie2/markdup/{sample_id}.mkD.sorted.bam',
            include = '{BASE_FOLDER}/genome/genome.include_regions.bed',
            bamtools_filter_pe = 'assests/bamtools_filter_pe.json'
        output:
            filter_bam = '{BASE_FOLDER}/bowtie2/filter/{sample_id}.flT.sorted.bam'
        params:
            samtools = config['tools']['samtools'],
            bamtools = config['tools']['bamtools']
        threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
        shell:
            """
            {params.samtools} view \
            -F 0x004 -F 0x0008 -f 0x001 -F 0x0400 -q 1 \
            -L {input.include} -b {input.bam_file} \
            | {params.bamtools} filter -out {output.filter_bam} \
            -script {input.bamtools_filter_pe}
            """
else:
    rule filter_bam_paired:
        input:
            bam_file = '{BASE_FOLDER}/bowtie2/markdup/{sample_id}.mkD.sorted.bam',
            include = '{BASE_FOLDER}/genome/genome.include_regions.bed',
            bamtools_filter_pe = 'assests/bamtools_filter_se.json'
        output:
            filter_bam = '{BASE_FOLDER}/bowtie2/filter/{sample_id}.flT.sorted.bam'
        params:
            samtools = config['tools']['samtools'],
            bamtools = config['tools']['bamtools']
        threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
        shell:
            """
            {params.samtools} view --threads {threads} \
            -F 0x004 -F 0x0008 -f 0x001 -F 0x0400 -q 1 \
            -L {input.include} -b {input.bam_file} \
            | {params.bamtools} filter -out {output.filter_bam} \
            -script {input.bamtools_filter_se}
            """
if get_paired:
    rule name_sort_bam_paired:
        input:
            bam_file = '{BASE_FOLDER}/bowtie2/filter/{sample_id}.flT.sorted.bam',
        output:
            sorted_bam = temp('{BASE_FOLDER}/bowtie2/filter/{sample_id}.clN.name.sorted.bam'),
            out_bam = temp('{BASE_FOLDER}/bowtie2/filter/{sample_id}.clN.bam'),
            out_sort_bam = '{BASE_FOLDER}/bowtie2/filter/{sample_id}.clN.sorted.bam'
        params:
            samtools = config['tools']['samtools'],
            bampe_rm_orphan = 'scripts/bampe_rm_orphan.py',
            python = config['tools']['python']
        threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
        shell:
            """
            {params.samtools} sort -n --threads {threads} \
            -o {output.sorted_bam} \
            -T {wildcards.BASE_FOLDER}/bowtie2/filter/{wildcards.sample_id}.mLb.clN.name.sorted \
            {input.bam_file}

            {params.python} {params.bampe_rm_orphan} \
            {output.sorted_bam} {output.out_bam} --only_fr_pairs

            {params.samtools} sort --threads {threads} \
            -o {output.out_sort_bam} \
            -T {wildcards.BASE_FOLDER}/bowtie2/filter/{wildcards.sample_id}.clN.sorted \
            {output.out_bam}
            """
else:
    rule name_sort_bam_singled:
        input:
            bam_file = '{BASE_FOLDER}/bowtie2/filter/{sample_id}.flT.sorted.bam',
        output:
            out_bam = temp('{BASE_FOLDER}/bowtie2/filter/{sample_id}.clN.bam'),
            out_sort_bam = '{BASE_FOLDER}/bowtie2/filter/{sample_id}.clN.sorted.bam'
        params:
            samtools = config['tools']['samtools']
        threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
        shell:
            """
            cp {input.bam_file} {output.out_bam}
            
            {params.samtools} sort --threads {threads} \
            -o {output.sorted_bam} \
            -T {wildcards.BASE_FOLDER}/bowtie2/filter/{wildcards.sample_id}.clN.sorted \
            {input.out_bam}
            """

rule bam_index_clN:
    input:
        bam_file = '{BASE_FOLDER}/bowtie2/filter/{sample_id}.clN.sorted.bam'
    output:
        index_file = '{BASE_FOLDER}/bowtie2/filter/{sample_id}.clN.sorted.bam.bai'
    params:
        samtools = config['tools']['samtools']
    threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
    shell:
        """
        {params.samtools} index --threads {threads} {input.bam_file}
        """

rule bam_idxstats_clN:
    input:
        bam_file = '{BASE_FOLDER}/bowtie2/filter/{sample_id}.clN.sorted.bam'
    output:
        idxstats_file = '{BASE_FOLDER}/bowtie2/filter/samtools_stat/{sample_id}.clN.sorted.bam.idxstats'
    params:
        samtools = config['tools']['samtools']
    threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
    shell:
        """
        {params.samtools} idxstats --threads {threads} {input.bam_file} > {output.idxstats_file}
        """

rule bam_flagstat_clN:
    input:
        bam_file = '{BASE_FOLDER}/bowtie2/filter/{sample_id}.clN.sorted.bam'
    output:
        flagstat_file = '{BASE_FOLDER}/bowtie2/filter/samtools_stat/{sample_id}.clN.sorted.bam.flagstat'
    params:
        samtools = config['tools']['samtools']
    threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
    shell:
        """
        {params.samtools} flagstat --threads {threads} {input.bam_file} > {output.flagstat_file}
        """

rule bam_stats_clN:
    input:
        bam_file = '{BASE_FOLDER}/bowtie2/filter/{sample_id}.clN.sorted.bam'
    output:
        stats_file = '{BASE_FOLDER}/bowtie2/filter/samtools_stat/{sample_id}.clN.sorted.bam.stats'
    params:
        samtools = config['tools']['samtools']
    threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
    shell:
        """
        {params.samtools} stats --threads {threads} {input.bam_file} > {output.stats_file}
        """

rule picard_CollectMultipleMetrics:
    input:
        bam_file = '{BASE_FOLDER}/bowtie2/filter/{sample_id}.clN.sorted.bam',
        fa = '{BASE_FOLDER}/genome/genome.fa'
    output:
        multiext("{BASE_FOLDER}/bowtie2/filter/picard_metrics/{sample_id}.CollectMultipleMetrics",
                    ".alignment_summary_metrics", ".base_distribution_by_cycle_metrics",
                    ".insert_size_metrics", ".quality_by_cycle_metrics",
                    ".quality_distribution_metrics"),
        multiext("{BASE_FOLDER}/bowtie2/filter/picard_metrics/pdf/{sample_id}.CollectMultipleMetrics",
                    ".insert_size_histogram.pdf", ".quality_by_cycle.pdf",
                    ".read_length_histogram.pdf", ".quality_distribution.pdf")
    params:
        picard = config['tools']['picard']
    shell:
        """
        {params.picard} CollectMultipleMetrics \
        --VALIDATION_STRINGENCY LENIENT \
        --INPUT {input.bam_file} \
        --TMP_DIR {wildcards.BASE_FOLDER}/bowtie2/filter/picard_metrics \
        --OUTPUT {wildcards.BASE_FOLDER}/bowtie2/filter/picard_metrics/{wildcards.sample_id}.CollectMultipleMetrics \
        --REFERENCE_SEQUENCE {input.fa}

        mv {wildcards.BASE_FOLDER}/bowtie2/filter/picard_metrics/{wildcards.sample_id}.CollectMultipleMetrics*pdf \
        {wildcards.BASE_FOLDER}/bowtie2/filter/picard_metrics/pdf
        """

