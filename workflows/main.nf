// Main pipeline workflow

// Reference files
ref_fasta = file(params.reference)
ref_index = Channel.fromPath("${params.reference}*").collect()

// Input channel
if (!params.input) {
    exit 1, "No input specified. Use --input to provide input files."
}

// ============== 单样本分析流程 ==============
workflow SINGLE {
    main:
    // Create input channel from FASTQ files (single sample)
    Channel.fromFilePairs(params.input + "/*_{1,2}.fastq.gz", flat: true)
        .ifEmpty { exit 1, "No input files found in ${params.input}" }
        .set { reads }

    // Quality control
    FASTQC(reads)

    // Trimming
    trimmed = TRIMGALORE(reads)

    // Alignment
    aligned = BWAMEM(trimmed, ref_fasta, ref_index, 'cram')

    // Index alignment
    INDEXALIGNMENT(aligned.cram, ref_fasta, ref_index, null)

    // Variant calling
    GATKHAPLOTYPECALLER(aligned.cram, ref_fasta)

    emit:
    aligned = aligned
}

// ============== 配对样本分析流程 ==============
workflow PAIRED {
    main:
    // Create input channel from FASTQ files (paired: tumor + normal)
    // Expected format: samples/Tumor_SampleID/*_{1,2}.fastq.gz and samples/Normal_SampleID/*_{1,2}.fastq.gz
    tumor_reads = Channel.fromFilePairs(params.tumor_input + "/*_{1,2}.fastq.gz", flat: true)
        .ifEmpty { exit 1, "No tumor input files found in ${params.tumor_input}" }
    normal_reads = Channel.fromFilePairs(params.normal_input + "/*_{1,2}.fastq.gz", flat: true)
        .ifEmpty { exit 1, "No normal input files found in ${params.normal_input}" }

    // Quality control
    FASTQC(tumor_reads)
    FASTQC(normal_reads)

    // Trimming
    tumor_trimmed = TRIMGALORE(tumor_reads)
    normal_trimmed = TRIMGALORE(normal_reads)

    // Alignment
    tumor_aligned = BWAMEM(tumor_trimmed, ref_fasta, ref_index, 'cram')
    normal_aligned = BWAMEM(normal_trimmed, ref_fasta, ref_index, 'cram')

    // Index alignment
    INDEXALIGNMENT(tumor_aligned.cram, ref_fasta, ref_index, null)
    INDEXALIGNMENT(normal_aligned.cram, ref_fasta, ref_index, null)

    // Variant calling (somatic)
    GATKMUTECT2(
        tumor_aligned.cram,
        normal_aligned.cram,
        ref_fasta
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