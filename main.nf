// schema-somatic Nextflow Main Pipeline
// A Nextflow pipeline for somatic mutation analysis

nextflow.enable.dsl = 2

// ============== 导入子流程 ==============
include { SINGLE } from './workflows/single/main'
include { PAIRED } from './workflows/paired/main'
include { CNV_BASELINE } from './workflows/cnv_baseline/main'
include { MSI_BASELINE } from './workflows/msi_baseline/main'

// ============== 参数设置 ==============
// Print help message
def helpMessage() {
    log.info """
    ============================================
    Schema-Somatic: Somatic Mutation Analysis Pipeline
    ============================================

    Usage:
        nextflow run . [options]

    Options:
        --input         Input directory containing FASTQ files (single mode)
        --tumor_input   Tumor FASTQ directory (paired mode)
        --normal_input  Normal FASTQ directory (paired mode)
        --output        Output directory (default: results)
        --reference     Reference genome fasta file
        --genome_version Reference genome version (GRCh38 or GRCh37)

    Analysis Modes:
        --analysis_mode single      Single sample analysis
        --analysis_mode paired     Paired sample analysis (tumor-normal)

    Baseline Building:
        --build_baseline cnv       Build CNV baseline
        --build_baseline msi       Build MSI baseline

    Optional Steps:
        --run_mutect2              Run Mutect2 for SNP/Indel detection
        --run_deepvariant          Run DeepVariant for SNP/Indel detection
        --run_manta                Run Manta for structural variants
        --run_gridss               Run GRIDSS for structural variants
        --run_cnvkit               Run CNVkit for CNV analysis
        --run_msisensorpro         Run MSI sensor pro for MSI analysis
        --run_lohhla               Run LOHHLA for HLA loss analysis

    Additional Parameters:
        --cnv_baseline             CNV baseline file
        --msi_baseline             MSI baseline file
        --hla_results              HLA typing results file
        --user_bed                 User-provided BED file for CNVkit
        --normal_bams              Normal BAM files directory for baseline building

    Profiles:
        -profile standard          Standard execution
        -profile docker            Docker container execution
        -profile singularity       Singularity container execution
    """
}

// Show help message if --help is specified
if (params.help) {
    helpMessage()
    exit 0
}

// ============== 参考基因组 ==============
ref_fasta = file(params.reference)
ref_index = Channel.fromPath("${params.reference}*").collect()
genome_version = params.genome_version ?: 'GRCh38'

// 检查输入
if (!params.input && !params.tumor_input && !params.build_baseline) {
    exit 1, "No input specified. Use --input (single), --tumor_input (paired), or --build_baseline"
}

// ============== 主入口 ==============
workflow {
    // 建立基线模式
    if (params.build_baseline) {
        switch (params.build_baseline) {
            case 'cnv':
                CNV_BASELINE()
                break
            case 'msi':
                MSI_BASELINE()
                break
            default:
                exit 1, "Invalid baseline type: ${params.build_baseline}. Use 'cnv' or 'msi'"
        }
    } else {
        // 分析模式
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
}