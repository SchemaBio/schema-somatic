// ============== 单样本分析流程 ==============
// 导入模块
include { FASTP } from '../modules/fastp/main'
include { BWAMEM } from '../modules/bwamem/main'
include { INDEXALIGNMENT } from '../modules/samtools/main'
include { MUTECT2 } from '../modules/gatk/main'
include { MUTECT2 as HOTSPOT } from '../modules/gatk/main'
include { DEEPVARIANT } from '../modules/deepvariant/main'
include { MANTA } from '../modules/manta/main'
include { GRIDSS } from '../modules/gridss/main'
include { CNVKIT } from '../modules/cnvkit/main'
include { MSISENSORPRO } from '../modules/msisensorpro/main'
include { WHATSHAP } from '../modules/whatshap/main'
include { EGFRHAP } from '../modules/whatshap/main'
include { VEP } from '../modules/vep/main'
include { VEP as VEP_HOTSPOT } from '../modules/vep/main'
include { VEP as VEP_GERMLINE } from '../modules/vep/main'

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
            null, null, null, null, null, null, null
        )

        // 4.2 单倍型分型
        WHATSHAP(
            MUTECT2.out.left_vcf.join(aligned.cram.map { file -> tuple(file.baseName, file, file + '.crai') }),
            ref_fasta,
            ref_index
        )

        // 4.3 EGFR区域变异合并
        EGFRHAP(
            WHATSHAP.out.phased_vcf,
            ref_fasta,
            ref_index,
            genome_version
        )

        // 4.4 EGFR区域VEP注释
        if (params.vep_database) {
            VEP(
                EGFRHAP.out.merged_vcf,
                ref_fasta,
                ref_index,
                params.vep_database,
                genome_version,
                'egfr'
            )
        }
    }

    // 4.1 热点区域变异检测
    if (params.run_hotspot && params.hotspot) {
        HOTSPOT(
            aligned.cram.map { file -> tuple(file.baseName, file, file + '.crai') },
            ref_fasta,
            ref_index,
            null, null, null, null, null,
            params.hotspot,
            'hotspot'
        )

        // 4.5 热点区域VEP注释
        if (params.vep_database) {
            VEP_HOTSPOT(
                HOTSPOT.out.left_vcf,
                ref_fasta,
                ref_index,
                params.vep_database,
                genome_version,
                'hotspot'
            )
        }
    }

    if (params.run_deepvariant) {
        DEEPVARIANT(
            aligned.cram.map { file -> tuple(file.baseName, file, file + '.crai') },
            ref_fasta,
            ref_index,
            params.model_type ?: 'WGS'
        )

        // 4.6 种系变异VEP注释
        if (params.vep_database) {
            VEP_GERMLINE(
                DEEPVARIANT.out.vcf,
                ref_fasta,
                ref_index,
                params.vep_database,
                genome_version,
                'germline'
            )
        }
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