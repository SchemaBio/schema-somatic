// BWA MEM module for read alignment

process BWAMEM {
    tag "BWAMEM on $sample_id"
    label 'bwamem'
    publishDir "${params.output}/02.Alignment", mode: 'copy'

    input:
        tuple val(sample_id), path(reads)
        path fasta
        path fasta_index
        val output_format
        val rgid

    output:
        path("${sample_id}.cram"), emit: cram, optional: true
        path("${sample_id}.cram.crai"), emit: crai, optional: true
        path("${sample_id}.bam"), emit: bam, optional: true
        path("${sample_id}.bam.bai"), emit: bai, optional: true

    script:
    def threads = task.cpus
    if (output_format == 'bam') {
        """
        bwa mem -t ${threads} -M -R "@RG\\tID:${rgid}\\tSM:${sample_id}\\tPL:SCHEMABIO\\tPU:Somatic" $fasta $reads \\
            | samtools sort -@ ${threads} -O bam -o ${sample_id}.bam -
        samtools index -@ ${threads} ${sample_id}.bam
        """
    } else {
        """
        bwa mem -t ${threads} -M -R "@RG\\tID:${rgid}\\tSM:${sample_id}\\tPL:SCHEMABIO\\tPU:Somatic" $fasta $reads \\
            | samtools sort -@ ${threads} \\
            --reference ${fasta} -O cram --output-fmt-option version=3.0 -o ${sample_id}.cram -
        samtools index -@ ${threads} ${sample_id}.cram
        """
    }
}

process BWAMEM2 {
    tag "BWAMEM2 on $sample_id"
    label 'bwamem2'
    publishDir "${params.output}/02.Alignment", mode: 'copy'

    input:
        tuple val(sample_id), path(reads)
        path fasta
        path fasta_index
        val output_format
        val rgid

    output:
        path("${sample_id}.cram"), emit: cram, optional: true
        path("${sample_id}.cram.crai"), emit: crai, optional: true
        path("${sample_id}.bam"), emit: bam, optional: true
        path("${sample_id}.bam.bai"), emit: bai, optional: true

    script:
    def threads = task.cpus
    if (output_format == 'bam') {
        """
        bwa-mem2 mem -t ${threads} -M -R "@RG\\tID:${rgid}\\tSM:${sample_id}\\tPL:SCHEMABIO\\tPU:Somatic" $fasta $reads \\
            | samtools sort -@ ${threads} -O bam -o ${sample_id}.bam -
        samtools index -@ ${threads} ${sample_id}.bam
        """
    } else {
        """
        bwa-mem2 mem -t ${threads} -M -R "@RG\\tID:${rgid}\\tSM:${sample_id}\\tPL:SCHEMABIO\\tPU:Somatic" $fasta $reads \\
            | samtools sort -@ ${threads} \\
            --reference ${fasta} -O cram --output-fmt-option version=3.0 -o ${sample_id}.cram -
        samtools index -@ ${threads} ${sample_id}.cram
        """
    }
}
