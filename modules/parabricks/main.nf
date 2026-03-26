// Parabricks GPU加速工具模块
// 用途：使用NVIDIA Clara Parabricks进行GPU加速的生物信息学分析
// 用法：
//   PB_BAM2FQ - 将BAM/CRAM文件转换为FASTQ格式
//   PB_DEEPVARIANT - GPU加速的DeepVariant变异检测
//   PB_FQ2BAM - GPU加速的BWA MEM比对

process PB_BAM2FQ {
    tag "PB_BAM2FQ on $sample_id"
    label 'pb_bam2fq'
    if (publish_dir) {
        publishDir "${publish_dir}", mode: 'copy'
    }

    input:
        tuple val(sample_id), path(bam), path(bai)
        val(publish_dir)

    output:
        tuple val(sample_id), path("${sample_id}_1.fq.gz"), path("${sample_id}_2.fq.gz"), emit: reads

    script:
    def threads = task.cpus
    """
    pbrun bam2fq \\
        --input ${bam} \\
        --output1 ${sample_id}_1.fq.gz \\
        --output2 ${sample_id}_2.fq.gz \\
        --threads ${threads}
    """
}

process PB_DEEPVARIANT {
    tag "PB_DEEPVARIANT on $sample_id"
    label 'pb_deepvariant'
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
    def gpu_devices = task.gpu ? task.gpu : 'all'
    """
    pbrun deepvariant \\
        --ref ${fasta} \\
        --in-bam ${bam} \\
        --out-variants ${sample_id}.deepvariant.vcf.gz \\
        --out-gvcf ${sample_id}.deepvariant.gvcf.gz \\
        --model-type ${model_type} \\
        --threads ${threads} \\
        --gpu-devices ${gpu_devices}

    tabix -p vcf ${sample_id}.deepvariant.vcf.gz
    """
}

process PB_FQ2BAM {
    tag "PB_FQ2BAM on $sample_id"
    label 'pb_fq2bam'
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
    def gpu_devices = task.gpu ? task.gpu : 'all'
    if (output_format == 'bam') {
        """
        pbrun fq2bam \\
            --ref ${fasta} \\
            --in-fq ${reads} \\
            --out-bam ${sample_id}.bam \\
            --rg-id "${rgid}" \\
            --rg-sm "${sample_id}" \\
            --rg-pl "SCHEMABIO" \\
            --rg-pu "Somatic" \\
            --threads ${threads} \\
            --gpu-devices ${gpu_devices} \\
            --gpusort

        samtools index -@ ${threads} ${sample_id}.bam
        """
    } else {
        """
        pbrun fq2bam \\
            --ref ${fasta} \\
            --in-fq ${reads} \\
            --out-bam ${sample_id}.cram \\
            --rg-id "${rgid}" \\
            --rg-sm "${sample_id}" \\
            --rg-pl "SCHEMABIO" \\
            --rg-pu "Somatic" \\
            --threads ${threads} \\
            --gpu-devices ${gpu_devices} \\
            --gpusort

        samtools index -@ ${threads} ${sample_id}.cram
        """
    }
}