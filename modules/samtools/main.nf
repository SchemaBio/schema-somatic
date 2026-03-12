process INDEXALIGNMENT {
    tag "INDEXALIGNMENT on $sample_id"
    label 'indexalignment'
    if (publish_dir) {
        publishDir "${publish_dir}", mode: 'copy'
    }

    input:
        tuple val(sample_id), path(alignment_file)
        val(publish_dir)

    output:
        path("*.bai"), emit: bai, optional: true
        path("*.crai"), emit: crai, optional: true

    script:
    def threads = task.cpus
    def ext = alignment_file.extension
    if (ext == 'cram') {
        """
        samtools index -@ ${threads} -C ${alignment_file}
        """
    } else {
        """
        samtools index -@ ${threads} ${alignment_file}
        """
    }
}

process BAM2FASTQ {
    tag "BAM2FASTQ on $sample_id"
    label 'bam2fastq'
    if (publish_dir) {
        publishDir "${publish_dir}", mode: 'copy'
    }

    input:
        tuple val(sample_id), path(bam)
        val(publish_dir)

    output:
        path("${sample_id}.clean.umi_*.fq.gz"), emit: reads

    script:
    """
    samtools fastq ${bam} \\
        -1 ${sample_id}_1.fq.gz \\
        -2 ${sample_id}_2.fq.gz
    """
}