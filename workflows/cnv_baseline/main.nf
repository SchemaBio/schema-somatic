// ============== CNV基线建立流程 ==============
// 导入模块
include { CNVKITBASELINE; CNVKITNORMAL } from '../modules/cnvkit/main'

// 参考基因组
ref_fasta = file(params.reference)
ref_index = Channel.fromPath("${params.reference}*").collect()
genome_version = params.genome_version ?: 'GRCh38'

workflow CNV_BASELINE {
    main:
    // 输入多个正常样本的BAM文件
    Channel.fromFilePairs(params.normal_bams + "/*.bam", flat: true)
        .ifEmpty { exit 1, "No BAM files found in ${params.normal_bams}" }
        .map { sample_id, bam, bai -> tuple(sample_id, bam, bai) }
        .set { normal_samples }

    // 建立CNV基线
    if (params.cnv_baseline_type == 'multiple') {
        // 多个样本建立基线
        CNVKITBASELINE(
            normal_samples.map { it[1] }.collect(),
            normal_samples.map { it[2] }.collect(),
            ref_fasta,
            ref_index,
            genome_version,
            params.user_bed ?: null,
            params.basename ?: 'cnv_baseline'
        )

        emit:
        baseline = CNVKITBASELINE.out.baseline
    } else {
        // 单样本建立基线
        CNVKITNORMAL(
            normal_samples.map { it[1] }.first(),
            normal_samples.map { it[2] }.first(),
            ref_fasta,
            ref_index,
            genome_version,
            params.user_bed ?: null,
            params.basename ?: 'cnv_baseline'
        )

        emit:
        baseline = CNVKITNORMAL.out.baseline
    }
}