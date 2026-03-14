// DeepVariant 变异检测模块
// 用途：使用深度学习进行高精度的种系变异检测
// 用法：
//   - 输入：样本ID、CRAM文件及索引、参考基因组、模型类型
//   - 模型类型：WGS(全基因组)、WES(外显子)、PACBIO(长读长)
//   - 输出：VCF文件和gVCF文件
//   - 自动使用多线程加速

process DEEPVARIANT {
    tag "DEEPVARIANT on $sample_id"
    label 'deepvariant'
    publishDir "${params.output}/03.Mutations/Germline", mode: 'copy'

    input:
        tuple val(sample_id), path(cram), path(crai)
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
        --reads=${cram} \\
        --output_vcf=${sample_id}.deepvariant.vcf.gz \\
        --output_gvcf=${sample_id}.deepvariant.gvcf.gz \\
        --num_shards=${threads}

    tabix -p vcf ${sample_id}.deepvariant.vcf.gz
    """
}
