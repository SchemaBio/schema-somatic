// DeepVariant module for variant calling

process DEEPVARIANT {
    tag "DEEPVARIANT on $sample_id"
    label 'deepvariant'
    publishDir "${params.output}/03.Mutations/Germline", mode: 'copy'

    input:
        tuple val(sample_id), path(bam), path(bai)
        path fasta
        path fasta_index
        val model_type  // 'WGS', 'WES', or 'PACBIO'

    output:
        tuple val(sample_id), path("${sample_id}.deepvariant.vcf.gz"), path("${sample_id}.deepvariant.vcf.gz.tbi"), emit: vcf
        path("${sample_id}.deepvariant.gvcf.gz"), emit: gvcf, optional: true

    script:
    def threads = task.cpus
    """
    /opt/deepvariant/bin/run_deepvariant \\
        --model_type=${model_type} \\
        --ref=${fasta} \\
        --reads=${bam} \\
        --output_vcf=${sample_id}.deepvariant.vcf.gz \\
        --output_gvcf=${sample_id}.deepvariant.gvcf.gz \\
        --num_shards=${threads}

    tabix -p vcf ${sample_id}.deepvariant.vcf.gz
    """
}
