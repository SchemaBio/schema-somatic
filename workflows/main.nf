// Schema-Somatic 主流程
// 支持单样本和配对样本的体细胞突变分析

// 导入模块
include { FASTP } from '../modules/fastp/main'
include { BWAMEM; BWAMEM2 } from '../modules/bwamem/main'
include { INDEXALIGNMENT } from '../modules/samtools/main'
include { MARKDUPLICATES; MUTECT2 } from '../modules/gatk/main'
include { WHATSHAP; EGFRHAP } from '../modules/whatshap/main'
include { DEEPVARIANT } from '../modules/deepvariant/main'
include { VEP } from '../modules/vep/main'
include { MANTA } from '../modules/manta/main'
include { GRIDSS } from '../modules/gridss/main'
include { CNVKIT } from '../modules/cnvkit/main'

// 参考基因组
ref_fasta = file(params.reference)
ref_index = Channel.fromPath("${params.reference}*").collect()

// 检查输入
if (!params.input) {
    exit 1, "No input specified. Use --input to provide input files."
}

// ============== 单样本分析流程 ==============
workflow SINGLE {
    main:
    // 输入FASTQ文件
    Channel.fromFilePairs(params.input + "/*_{1,2}.fastq.gz", flat: true)
        .ifEmpty { exit 1, "No input files found in ${params.input}" }
        .set { reads }

    // 1. 质控
    FASTP(reads)

    // 2. 比对
    aligned = BWAMEM(
        FASTP.out.clean_reads,
        ref_fasta,
        ref_index,
        params.output_format ?: 'cram',
        params.rgid ?: 'RG1'
    )

    // 3. 索引
    INDEXALIGNMENT(
        aligned.cram.map { file -> tuple(file.baseName, file) },
        params.output + "/02.Alignment"
    )

    // 4. 变异检测
    if (params.run_mutect2) {
        MUTECT2(
            aligned.cram.map { file -> tuple(file.baseName, file, file + '.crai') },
            ref_fasta,
            ref_index,
            null, null, null, null, null
        )
    }

    if (params.run_deepvariant) {
        DEEPVARIANT(
            aligned.cram.map { file -> tuple(file.baseName, file, file + '.crai') },
            ref_fasta,
            ref_index,
            params.model_type ?: 'WGS'
        )
    }

    // 5. 结构变异
    if (params.run_manta) {
        MANTA(
            aligned.cram.map { file -> tuple(file.baseName, file, file + '.crai') },
            ref_fasta,
            ref_index
        )
    }

    emit:
    aligned = aligned
}

// ============== 配对样本分析流程 ==============
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

    // 3. 体细胞突变检测
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

    emit:
    tumor_aligned = tumor_aligned
    normal_aligned = normal_aligned
}

// ============== 主入口 ==============
workflow {
    if (!params.analysis_mode) {
        exit 1, "No analysis mode specified. Use --analysis_mode single or paired"
    }

    switch (params.analysis_mode) {
        case 'single':
            SINGLE()
            break
        case 'paired':
            PAIRED()
            break
        default:
            exit 1, "Invalid analysis mode: ${params.analysis_mode}. Use 'single' or 'paired'"
    }
}