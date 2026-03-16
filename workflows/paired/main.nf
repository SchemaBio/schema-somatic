// ============== 配对样本分析流程 ==============
// 导入模块
include { FASTP } from '../modules/fastp/main'
include { BWAMEM } from '../modules/bwamem/main'
include { MUTECT2 } from '../modules/gatk/main'
include { MANTA } from '../modules/manta/main'
include { GRIDSS } from '../modules/gridss/main'
include { CNVKIT } from '../modules/cnvkit/main'
include { MSISENSORPRO } from '../modules/msisensorpro/main'
include { LOHHLA } from '../modules/lohhla/main'

// 参考基因组
ref_fasta = file(params.reference)
ref_index = Channel.fromPath("${params.reference}*").collect()
genome_version = params.genome_version ?: 'GRCh38'

workflow PAIRED {
    main:
    // 输入肿瘤和正常样本FASTQ
    tumor_reads = Channel.fromFilePairs(params.tumor_input + "/*_{1,2}.fastq.gz", flat: true)
        .ifEmpty { exit 1, "No tumor input files found" }
    normal_reads = Channel.fromFilePairs(params.normal_input + "/*_{1,2}.fastq.gz", flat: true)
        .ifEmpty { exit 1, "No normal input files found" }

    // 1. 质控
    tumor_qc = FASTP(tumor_reads)
    normal_qc = FASTP(normal_reads)

    // 2. 比对
    tumor_aligned = BWAMEM(tumor_qc.clean_reads, ref_fasta, ref_index, 'cram', 'tumor')
    normal_aligned = BWAMEM(normal_qc.clean_reads, ref_fasta, ref_index, 'cram', 'normal')

    // 3. 体细胞突变检测 - SNP/Indel
    if (params.run_mutect2) {
        MUTECT2(
            tumor_aligned.cram.map { file -> tuple(file.baseName, file, file + '.crai') },
            ref_fasta,
            ref_index,
            normal_aligned.cram,
            normal_aligned.crai,
            params.gnomad ?: null,
            params.pon ?: null,
            params.bed ?: null
        )
    }

    // 4. 结构变异
    if (params.run_manta) {
        MANTA(
            tumor_aligned.cram.map { file -> tuple(file.baseName, file, file + '.crai') },
            ref_fasta,
            ref_index
        )
    }

    if (params.run_gridss) {
        GRIDSS(
            tumor_aligned.cram.map { file -> tuple(file.baseName, file, file + '.crai') },
            ref_fasta,
            ref_index
        )
    }

    // 5. CNV分析
    if (params.run_cnvkit) {
        if (params.cnv_baseline) {
            CNVKIT(
                tumor_aligned.cram.map { file -> tuple(file.baseName, file, file + '.crai') },
                ref_fasta,
                ref_index,
                genome_version,
                params.user_bed ?: null,
                params.cnv_baseline
            )
        }
    }

    // 6. MSI分析
    if (params.run_msisensorpro) {
        if (params.msi_baseline) {
            MSISENSORPRO(
                tumor_aligned.cram.map { file -> tuple(file.baseName, file, file + '.crai') },
                ref_fasta,
                ref_index,
                params.msi_baseline
            )
        }
    }

    // 7. LOHHLA分析
    if (params.run_lohhla && params.hla_results) {
        tumor_bams = tumor_aligned.cram.map { file -> tuple(file.baseName, file, file + '.crai') }
        normal_bams = normal_aligned.cram.map { file -> tuple(file.baseName, file, file + '.crai') }

        tumor_bams.cross(normal_bams) { it[0] }.map { tumor, normal ->
            tuple(tumor[0], tumor[1], tumor[2], normal[1], normal[2])
        }.set { paired_bams }

        if (params.cnv_baseline) {
            LOHHLA(
                paired_bams.map { sample_id, tumor_bam, tumor_bai, normal_bam, normal_bai ->
                    tuple(sample_id, tumor_bam, tumor_bai)
                },
                paired_bams.map { sample_id, tumor_bam, tumor_bai, normal_bam, normal_bai ->
                    tuple(sample_id, normal_bam, normal_bai)
                },
                params.cnv_baseline,
                params.hla_results,
                genome_version
            )
        }
    }

    emit:
    tumor_aligned = tumor_aligned
    normal_aligned = normal_aligned
    tumor_qc = tumor_qc
    normal_qc = normal_qc
}