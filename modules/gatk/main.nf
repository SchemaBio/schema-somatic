// GATK 变异检测模块
// 用途：标记重复、合并比对、质控统计、变异检测
// 包含的process：
//   MARKDUPLICATES - 标记PCR重复
//   MARKDUPLICATESUMI - 基于UMI标记重复
//   MERGEBAMALIGNMENT - 合并比对和未比对BAM
//   COLLECTQCMETRICS - 收集质控指标
//   LEFTALIGNANDTRIMVARIANTS - 左对齐和标准化变异
//   MUTECT2 - 体细胞突变检测

process MARKDUPLICATES {
    tag "MARKDUPLICATES on $sample_id"
    label 'markduplicates'
    publishDir "${params.output}/02.Alignment", mode: 'copy'

    input:
        tuple val(sample_id), path(cram), path(crai)
        val   fasta

    output:
        path("${sample_id}.markdup.cram"), emit: cram
        path("${sample_id}.markdup.metrics.txt"), emit: metrics

    script:
    """
    gatk MarkDuplicates \\
        -I ${cram} \\
        -O ${sample_id}.markdup.cram \\
        -M ${sample_id}.markdup.metrics.txt \\
        --CREATE_INDEX false \\
        --REFERENCE_SEQUENCE ${fasta}
    """
}

process MARKDUPLICATESUMI {
    tag "MARKDUPLICATESUMI on $sample_id"
    label 'markduplicatesumi'
    publishDir "${params.output}/02.Alignment", mode: 'copy'

    input:
        tuple val(sample_id), path(cram), path(crai)
        val   fasta

    output:
        path("${sample_id}.markdup.cram"), emit: cram
        path("${sample_id}.markdup.metrics.txt"), emit: metrics
        path("${sample_id}.markdup.metrics.umi.txt"), emit: umi_metrics

    script:
    """
    gatk UmiAwareMarkDuplicatesWithMateCigar \\
        -I ${cram} \\
        -O ${sample_id}.markdup.cram \\
        -M ${sample_id}.markdup.metrics.txt \\
        --CREATE_INDEX false \\
        --REFERENCE_SEQUENCE ${fasta} \\
        --UMI_METRICS_FILE ${sample_id}.markdup.umi.metrics.txt
    """
}

process MERGEBAMALIGNMENT {
    tag "MERGEBAMALIGNMENT on $sample_id"
    label 'mergebamalignment'
    publishDir "${params.output}/02.Alignment", mode: 'copy'

    input:
        tuple val(sample_id), path(mapped_bam), path(unmapped_bam)
        val   fasta

    output:
        path("${sample_id}.bam"), emit: merged_bam

    script:
    """
    gatk MergeBamAlignment \\
        -R ${fasta} \\
        -ALIGNED ${mapped_bam} -UNMAPPED ${unmapped_bam} \\
        -O ${sample_id}.merged.bam \\
        --ALIGNER_PROPER_PAIR_FLAGS true \\
        --ATTRIBUTES_TO_RETAIN XS
    """
}

process COLLECTQCMETRICS {
    tag "COLLECTQCMETRICS on $sample_id"
    label 'collectqcmeterics'
    publishDir "${params.output}/01.QC", mode: 'copy'

    input:
        tuple val(sample_id), path(cram), path(crai)
        val   fasta
        val   fasta_dict
        val   bed

    output:
        path("${sample_id}.hs_metrics.txt"), emit: hs_metrics
        path("${sample_id}.insert_size_metrics.txt"), emit: insert_size_metrics
        path("${sample_id}.insert_size_histogram.pdf"), emit: insert_size_histogram
        path("${sample_id}.est_lib_complex_metrics.txt"), emit: lib_metrics

    script:
    """
    gatk BedToIntervalList \\
        --VALIDATION_STRINGENCY SILENT \\
        --TMP_DIR tmp \\
        --SEQUENCE_DICTIONARY ${fasta_dict} \\
        --INPUT ${bed} \\
        --OUTPUT Bait.interval_list

    gatk CollectHsMetrics \\
        --VALIDATION_STRINGENCY SILENT \\
        --TMP_DIR tmp \\
        --BAIT_INTERVALS Bait.interval_list \\
        --TARGET_INTERVALS Bait.interval_list \\
        --INPUT ${cram} \\
        --OUTPUT ${sample_id}.hs_metrics.txt

    gatk CollectInsertSizeMetrics \\
        --VALIDATION_STRINGENCY SILENT \\
        --TMP_DIR tmp \\
        --INPUT ${cram} \\
        --OUTPUT ${sample_id}.insert_size_metrics.txt \\
        --Histogram_FILE ${sample_id}.insert_size_histogram.pdf

    gatk EstimateLibraryComplexity \\
        --MAX_RECORDS_IN_RAM 303942330 \\
        --VALIDATION_STRINGENCY SILENT \\
        --TMP_DIR tmp \\
        --INPUT ${cram} \\
        --OUTPUT ${sample_id}.est_lib_complex_metrics.txt
    """
}

process LEFTALIGNANDTRIMVARIANTS {
    tag "LEFTALIGNANDTRIMVARIANTS on $sample_id"
    label 'leftalignandtrimvariants'
    publishDir "${params.vcf_output_dir}", mode: 'copy'

    input:
        tuple val(sample_id), path(vcf)
        val   fasta

    output:
        path("${sample_id}.normalized.vcf.gz"), emit: left_vcf

    script:
    """
    gatk LeftAlignAndTrimVariants \\
        --variant ${vcf} \\
        --reference ${fasta} \\
        --output ${sample_id}.normalized.vcf.gz \\
        --split-multi-allelics
    """
}

process MUTECT2 {
    tag "MUTECT2 on $sample_id"
    label 'mutect2'
    publishDir "${params.output}/03.Mutations/SNV_InDel", mode: 'copy'

    input:
        tuple val(sample_id), path(cram), path(crai)
        path fasta
        path fasta_index
        val normal_cram   // 可选，正常样本CRAM
        val normal_crai  // 可选，正常样本索引
        val gnomad       // 可选，germline resource
        val pon          // 可选，panel of normals
        val bed          // 可选，捕获区域bed文件

    output:
        path("${sample_id}.mutect2.vcf.gz"), emit: vcf
        path("${sample_id}.mutect2.vcf.gz.tbi"), emit: tbi
        path("${sample_id}.mutect2.filtered.normalized.vcf.gz"), emit: left_vcf
        path("${sample_id}.mutect2.filtered.normalized.vcf.gz.tbi"), emit: left_tbi

    script:
    def normal_input = normal_cram ? "-I ${normal_cram}" : ""
    def germline_resource = gnomad ? "--germline-resource ${gnomad}" : ""
    def panel_of_normals = pon ? "-pon ${pon}" : ""
    def capture_bed = bed ? "-L ${bed}" : ""
    def threads = task.cpus
    def downsample = params.maxReadsPerAlignmentStart ? "--max-reads-per-alignment-start ${params.maxReadsPerAlignmentStart}" : "0"

    """
    gatk Mutect2 \\
        -R ${fasta} \\
        -I ${cram} \\
        ${normal_input} \\
        -O ${sample_id}.mutect2.vcf.gz \\
        -tumor ${sample_id} \\
        ${germline_resource} \\
        ${panel_of_normals} \\
        --native-pair-hmm-threads ${threads} \\
        ${capture_bed} \\
        -A Coverage \\
        -A GenotypeSummaries \\
        -mbq 15 \\
        --force-active true \\
        --callable-depth 50 \\
        ${downsample} \\
        --f1r2-tar-gz ${sample_id}.f1r2.tar.gz

    gatk LearnReadOrientationModel \\
        -I ${sample_id}.f1r2.tar.gz \\
        -O ${sample_id}.f1r2.model.tar.gz
    
    gatk FilterMutectCalls \\
        --min-slippage-length 5 \\
        -O ${sample_id}.mutect2.filtered.vcf.gz \\
        -R ${fasta} \\
        --orientation-bias-artifact-priors ${sample_id}.f1r2.model.tar.gz \\
        -V ${sample_id}.mutect2.vcf.gz

    gatk LeftAlignAndTrimVariants \\
        --variant ${sample_id}.mutect2.filtered.vcf.gz \\
        --reference ${fasta} \\
        --output ${sample_id}.mutect2.filtered.normalized.vcf.gz \\
        --split-multi-allelics
    """
}
