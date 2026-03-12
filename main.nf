// schema-somatic Nextflow Main Pipeline
// A Nextflow pipeline for somatic mutation analysis

nextflow.enable.dsl = 2

// Print help message
def helpMessage() {
    log.info """
    ============================================
    Schema-Somatic: Somatic Mutation Analysis Pipeline
    ============================================

    Usage:
        nextflow run . [options]

    Options:
        --input         Input directory containing FASTQ/BAM files
        --output        Output directory (default: results)
        --reference     Reference genome fasta file
        --sample_sheet  CSV file with sample information

    Profiles:
        -profile standard      Standard execution
        -profile docker        Docker container execution
        -profile singularity   Singularity container execution
        -profile conda         Conda environment execution
    """
}

// Show help message if --help is specified
if (params.help) {
    helpMessage()
    exit 0
}

// Include modules
include { fastqc } from './modules/fastqc/main'
include { trimGalore } from './modules/trimgalore/main'
include { bwaMem } from './modules/bwamem/main'
include { gatkHaplotypeCaller } from './modules/gatk/main'

// Include main workflow
include { pipeline } from './workflows/main'

// Main entry point
workflow {
    pipeline()
}