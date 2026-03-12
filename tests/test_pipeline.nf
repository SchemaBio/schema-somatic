#!/usr/bin/env nextflow

// Test workflow for schema-somatic pipeline

nextflow.enable.dsl = 2

include { fastqc } from '../modules/fastqc/main'

workflow test_fastqc_module {
    // Create test input channel
    Channel.fromPath("tests/data/test_*.fastq.gz")
        .map { file -> tuple(file.baseName.replaceAll(/_\d+$/, ''), file) }
        .groupTuple()
        .set { test_reads }

    fastqc(test_reads)
}

workflow {
    test_fastqc_module()
}