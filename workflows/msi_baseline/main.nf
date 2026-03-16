// ============== MSI基线建立流程 ==============
// 导入模块
include { MSISENSORPRO_SCAN; MSISENSORPRO_BASELINE } from '../modules/msisensorpro/main'

// 参考基因组
ref_fasta = file(params.reference)
ref_index = Channel.fromPath("${params.reference}*").collect()

workflow MSI_BASELINE {
    main:
    // 扫描参考基因组获取MSI位点
    MSISENSORPRO_SCAN(
        ref_fasta,
        ref_index,
        params.msi_scan_name ?: 'msi_scan'
    )

    // 输入多个正常样本的BAM文件
    Channel.fromFilePairs(params.normal_bams + "/*.bam", flat: true)
        .ifEmpty { exit 1, "No BAM files found in ${params.normal_bams}" }
        .map { sample_id, bam, bai -> tuple(sample_id, bam, bai) }
        .set { normal_samples }

    // 建立MSI基线
    MSISENSORPRO_BASELINE(
        normal_samples.map { it[1] }.collect(),
        normal_samples.map { it[2] }.collect(),
        ref_fasta,
        ref_index,
        MSISENSORPRO_SCAN.out.loci,
        params.basename ?: 'msi_baseline'
    )

    emit:
    baseline = MSISENSORPRO_BASELINE.out.baseline
    loci = MSISENSORPRO_SCAN.out.loci
}