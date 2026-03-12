// FGBIO module for UMI processing (Complete Workflow)

process FGBIOUMI2FQ {
    tag "FGBIOUMI2FQ on $sample_id"
    label 'fgbioumi2fq'
    publishDir "${params.output}/01.QC/Fgbio", mode: 'copy'

    input:
        tuple val(sample_id), path(read1), path(read2)
        val read_structure  // e.g. "5M5S+T +T"
        path fasta
        path fasta_index

    output:
        path("${sample_id}.unmapped.umi.bam"), emit: unmapped_bam

    script:
    """
    fgbio FastqToBam \\
        --input ${read1} ${read2} \\
        --read-structures ${read_structure} \\
        --output ${sample_id}.unmapped.bam

    fgbio ExtractUmisFromBam \\
        -i ${sample_id}.unmapped.bam \\
        -o ${sample_id}.unmapped.umi.bam \\
        -r ${read_structure} \\
        -s RX \\
        -t ZA ZB
    """
}