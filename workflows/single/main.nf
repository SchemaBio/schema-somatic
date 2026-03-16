// ============== 单样本分析流程 ==============
// 导入模块
include { FASTP } from '../modules/fastp/main'
include { BWAMEM } from '../modules/bwamem/main'
include { INDEXALIGNMENT } from '../modules/samtools/main'
include { MUTECT2 } from '../modules/gatk/main'
include { DEEPVARIANT } from '../modules/deepvariant/main'
include { MANTA } from '../modules/manta/main'
include { GRIDSS } from '../modules/gridss/main'
include { CNVKIT } from '../modules/cnvkit/main'
include { MSISENSORPRO } from '../modules/msisensorpro/main'

// 参考基因组
ref_fasta = file(params.reference)
ref_index = Channel.fromPath("${params.reference}*").collect()
genome_version = params.genome_version ?: 'GRCh38'

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

    // 4. 变异检测 - SNP/Indel
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

    if (params.run_gridss) {
        GRIDSS(
            aligned.cram.map { file -> tuple(file.baseName, file, file + '.crai') },
            ref_fasta,
            ref_index
        )
    }

    // 6. CNV分析
    if (params.run_cnvkit) {
        if (params.cnv_baseline) {
            CNVKIT(
                aligned.cram.map { file -> tuple(file.baseName, file, file + '.crai') },
                ref_fasta,
                ref_index,
                genome_version,
                params.user_bed ?: null,
                params.cnv_baseline
            )
        }
    }

    // 7. MSI分析
    if (params.run_msisensorpro) {
        if (params.msi_baseline) {
            MSISENSORPRO(
                aligned.cram.map { file -> tuple(file.baseName, file, file + '.crai') },
                ref_fasta,
                ref_index,
                params.msi_baseline
            )
        }
    }

    emit:
    aligned = aligned
    qc = FASTP.out
}