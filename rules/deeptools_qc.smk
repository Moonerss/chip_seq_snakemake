rule multiBamSummary:
    input:
        bams = lambda wildcards: expand('{{BASE_FOLDER}}/bowtie2/filter/{sample_id}.clN.sorted.bam', sample_id = get_antibody_samples(wildcards))
    output:
        npz = '{BASE_FOLDER}/bowtie2/deeptools_qc/multiBamSummary/{antibody_id}.read_coverage.bins.npz',
        scalingFactors = '{BASE_FOLDER}/bowtie2/deeptools_qc/multiBamSummary/{antibody_id}.scaling_factors.txt',
        outRawCounts = '{BASE_FOLDER}/bowtie2/deeptools_qc/multiBamSummary/{antibody_id}.readCounts.tab'
    params:
        deeptools = config['tools']['deeptools'],
        extra_multiBamSummary = config['extra']['multiBamSummary'],
        labels = lambda wildcards: get_antibody_samples(wildcards)
    threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
    log: '{BASE_FOLDER}/logs/deeptools/multiBamSummary/{antibody_id}.log'
    shell:
        """
        cmd_multiBamSummary=`dirname {params.deeptools}`/multiBamSummary
        $cmd_multiBamSummary bins --bamfiles {input.bams} \
        --outFileName {output.npz} \
        --labels {params.labels} \
        --numberOfProcessors {threads} \
        --outRawCounts {output.outRawCounts} \
        --scalingFactors {output.scalingFactors} \
        {params.extra_multiBamSummary} \
        --verbose &> {log} 
        """

rule plotPCA:
    input:
        npz = '{BASE_FOLDER}/bowtie2/deeptools_qc/multiBamSummary/{antibody_id}.read_coverage.bins.npz'
    output:
        plotFile = '{BASE_FOLDER}/bowtie2/deeptools_qc/plotPCA/{antibody_id}.PCA.read_coverage.pdf',
        outFile = '{BASE_FOLDER}/bowtie2/deeptools_qc/plotPCA/{antibody_id}.PCA.read_coverage.tsv'
    params:
        deeptools = config['tools']['deeptools']
    shell:
        """
        cmd_plotPCA=`dirname {params.deeptools}`/plotPCA
        $cmd_plotPCA -in {input.npz} \
        --plotFile {output.plotFile} \
        --outFileNameData {output.outFile}
        """

rule plotCoverage:
    input:
        bams = lambda wildcards: expand('{{BASE_FOLDER}}/bowtie2/filter/{sample_id}.clN.sorted.bam', sample_id = get_antibody_samples(wildcards))
    output:
        plotFile = '{BASE_FOLDER}/bowtie2/deeptools_qc/plotCoverage/{antibody_id}.read_coverage.pdf',
        outFile = '{BASE_FOLDER}/bowtie2/deeptools_qc/plotCoverage/{antibody_id}.read_coverage.tsv'
    params:
        deeptools = config['tools']['deeptools'],
        labels = lambda wildcards: get_antibody_samples(wildcards)
    threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
    log: '{BASE_FOLDER}/logs/deeptools/plotCoverage/{antibody_id}.log'
    shell:
        """
        cmd_plotCoverage=`dirname {params.deeptools}`/plotCoverage
        $cmd_plotCoverage --bamfiles {input.bams} \
        --labels {params.labels} --plotFile {output.plotFile} \
        --skipZeros --outRawCounts {output.outFile} \
        --numberOfProcessors {threads} &> {log}
        """

rule plotCorrelation:
    input:
        npz = '{BASE_FOLDER}/bowtie2/deeptools_qc/multiBamSummary/{antibody_id}.read_coverage.bins.npz'
    output:
        plotFile1 = '{BASE_FOLDER}/bowtie2/deeptools_qc/plotCorrelation/{antibody_id}.pearson.heatmap.pdf',
        plotFile2 = '{BASE_FOLDER}/bowtie2/deeptools_qc/plotCorrelation/{antibody_id}.spearman.heatmap.pdf',
        outFile1 = '{BASE_FOLDER}/bowtie2/deeptools_qc/plotCorrelation/{antibody_id}.pearson.heatmap.tsv',
        outFile2 = '{BASE_FOLDER}/bowtie2/deeptools_qc/plotCorrelation/{antibody_id}.spearman.heatmap.tsv'
    params:
        deeptools = config['tools']['deeptools'],
        labels = lambda wildcards: get_antibody_samples(wildcards)
    shell:
        """
        cmd_plotCorrelation=`dirname {params.deeptools}`/plotCorrelation

        $cmd_plotCorrelation --corData {input.npz} \
        --corMethod pearson --whatToPlot heatmap --removeOutliers \
        --plotFile {output.plotFile1} --skipZeros \
        --labels {params.labels} --outFileCorMatrix {output.outFile1}

        $cmd_plotCorrelation --corData {input.npz} \
        --corMethod spearman --whatToPlot heatmap \
        --plotFile {output.plotFile2} --skipZeros \
        --labels {params.labels} --outFileCorMatrix {output.outFile2}
        """

rule plotFingerprint:
    input:
        bams = lambda wildcards: expand('{{BASE_FOLDER}}/bowtie2/filter/{sample_id}.clN.sorted.bam', sample_id = get_antibody_samples(wildcards))
    output:
        plotFile = '{BASE_FOLDER}/bowtie2/deeptools_qc/plotFingerprint/{antibody_id}.Fingerprint.pdf',
        outFile = '{BASE_FOLDER}/bowtie2/deeptools_qc/plotFingerprint/{antibody_id}.Fingerprint.tsv'
    params:
        deeptools = config['tools']['deeptools'],
        labels = lambda wildcards: get_antibody_samples(wildcards)
    threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
    shell:
        """
        cmd_plotFingerprint=`dirname {params.deeptools}`/plotFingerprint
        $cmd_plotFingerprint --bamfiles {input.bams} \
        --plotFile {output.plotFile} --outRawCounts {output.outFile} \
        --labels {params.labels} --skipZeros --numberOfProcessors {threads}
        """

rule bamPEFragmentSize:
    input:
        bams = lambda wildcards: expand('{{BASE_FOLDER}}/bowtie2/filter/{sample_id}.clN.sorted.bam', sample_id = get_antibody_samples(wildcards))
    output:
        plotFile = '{BASE_FOLDER}/bowtie2/deeptools_qc/bamPEFragmentSize/{antibody_id}.bamPEFragmentSize.pdf',
        outFile = '{BASE_FOLDER}/bowtie2/deeptools_qc/bamPEFragmentSize/{antibody_id}.bamPEFragmentSize.tsv'
    params:
        deeptools = config['tools']['deeptools'],
        labels = lambda wildcards: get_antibody_samples(wildcards)
    threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
    shell:
        """
        cmd_bamPEFragmentSize=`dirname {params.deeptools}`/bamPEFragmentSize
        $cmd_bamPEFragmentSize --bamfiles {input.bams} --histogram {output.plotFile} \
        --numberOfProcessors {threads} --samplesLabel {params.labels} \
        --binSize 100 --table {output.outFile}
        """

rule computeMatrix_point:
    input:
        bigwigs = lambda wildcards: expand('{{BASE_FOLDER}}/bowtie2/bam_coverage/{all_damples}.bigwig', all_damples = get_paired_samples(wildcards)),
        regions = '{BASE_FOLDER}/genome/genome.bed'
    output:
        gz_file = '{BASE_FOLDER}/bowtie2/enrichedHeatmap/computeMatrix/reference_point/{antibody_id}/{sample_id}.matrix.gz',
        bed_file = '{BASE_FOLDER}/bowtie2/enrichedHeatmap/computeMatrix/reference_point/{antibody_id}/{sample_id}.region.bed',
        tab_file = '{BASE_FOLDER}/bowtie2/enrichedHeatmap/computeMatrix/reference_point/{antibody_id}/{sample_id}.matrix.tsv'
    params:
        deeptools = config['tools']['deeptools'],
        labels = lambda wildcards: get_paired_samples(wildcards),
        extra_computeMatrix = config['extra']['computeMatrix_reference_point']
    threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
    log: '{BASE_FOLDER}/logs/deeptools/computeMatrix/reference_point/{antibody_id}_{sample_id}.log'
    shell:
        """
        cmd_computeMatrix=`dirname {params.deeptools}`/computeMatrix
        $cmd_computeMatrix reference-point \
        --regionsFileName {input.regions} \
        --scoreFileName {input.bigwigs} \
        --outFileName {output.gz_file} \
        --outFileNameMatrix {output.tab_file} \
        --outFileSortedRegions {output.bed_file} \
        {params.extra_computeMatrix} \
        --missingDataAsZero --skipZeros \
        --samplesLabel {params.labels} \
        --numberOfProcessors {threads} \
        --verbose &> {log}
        """

rule computeMatrix_region:
    input:
        bigwigs = lambda wildcards: expand('{{BASE_FOLDER}}/bowtie2/bam_coverage/{all_damples}.bigwig', all_damples = get_paired_samples(wildcards)),
        regions = '{BASE_FOLDER}/genome/genome.bed'
    output:
        gz_file = '{BASE_FOLDER}/bowtie2/enrichedHeatmap/computeMatrix/scale_regions/{antibody_id}/{sample_id}.matrix.gz',
        bed_file = '{BASE_FOLDER}/bowtie2/enrichedHeatmap/computeMatrix/scale_regions/{antibody_id}/{sample_id}.region.bed',
        tab_file = '{BASE_FOLDER}/bowtie2/enrichedHeatmap/computeMatrix/scale_regions/{antibody_id}/{sample_id}.matrix.tsv'
    params:
        deeptools = config['tools']['deeptools'],
        labels = lambda wildcards: get_paired_samples(wildcards),
        extra_computeMatrix = config['extra']['computeMatrix_scale_regions']
    threads: lambda wildcards: 10 if 10 < config['threads'] else config['threads']
    log: '{BASE_FOLDER}/logs/deeptools/computeMatrix/scale_regions/{antibody_id}_{sample_id}.log'
    shell:
        """
        cmd_computeMatrix=`dirname {params.deeptools}`/computeMatrix
        $cmd_computeMatrix scale-regions \
        --regionsFileName {input.regions} \
        --scoreFileName {input.bigwigs} \
        --outFileName {output.gz_file} \
        --outFileNameMatrix {output.tab_file} \
        --outFileSortedRegions {output.bed_file} \
        {params.extra_computeMatrix} \
        --missingDataAsZero --skipZeros \
        --samplesLabel {params.labels} \
        --numberOfProcessors {threads} \
        --verbose &> {log}
        """

rule plotHeatmap_point:
    input:
        gz_file = '{BASE_FOLDER}/bowtie2/enrichedHeatmap/computeMatrix/reference_point/{antibody_id}/{sample_id}.matrix.gz'
    output:
        plot = '{BASE_FOLDER}/bowtie2/enrichedHeatmap/plotHeatmap/reference_point/{antibody_id}/{sample_id}.heatmap.pdf'
    params:
        deeptools = config['tools']['deeptools'],
        extra_plotHeatmap = config['extra']['plotHeatmap']
    shell:
        """
        cmd_plotHeatmap=`dirname {params.deeptools}`/plotHeatmap
        $cmd_plotHeatmap --matrixFile {input.gz_file} \
        --outFileName {output.plot} \
        {params.extra_plotHeatmap}
        """

rule plotHeatmap_region:
    input:
        gz_file = '{BASE_FOLDER}/bowtie2/enrichedHeatmap/computeMatrix/scale_regions/{antibody_id}/{sample_id}.matrix.gz'
    output:
        plot = '{BASE_FOLDER}/bowtie2/enrichedHeatmap/plotHeatmap/scale_regions/{antibody_id}/{sample_id}.heatmap.pdf'
    params:
        deeptools = config['tools']['deeptools'],
        extra_plotHeatmap = config['extra']['plotHeatmap']
    shell:
        """
        cmd_plotHeatmap=`dirname {params.deeptools}`/plotHeatmap
        $cmd_plotHeatmap --matrixFile {input.gz_file} \
        --outFileName {output.plot} \
        {params.extra_plotHeatmap}
        """