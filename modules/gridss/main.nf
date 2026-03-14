// GRIDSS 结构变异检测模块
// 用途：高精度结构变异和断点检测
// 用法：
//   - 输入：样本ID、BAM/CRAM文件及索引、参考基因组
//   - 输出：结构变异VCF文件和assembly BAM
//   - 单样本模式，支持复杂结构变异检测
//   - 支持BAM和CRAM格式

process GRIDSS {
    tag "GRIDSS on $sample_id"
    label 'gridss'
    publishDir "${params.output}/03.Mutations/SV", mode: 'copy'

    input:
        tuple val(sample_id), path(alignment), path(index)
        path fasta
        path fasta_index

    output:
        tuple val(sample_id), path("${sample_id}.gridss.vcf.gz"), path("${sample_id}.gridss.vcf.gz.tbi"), emit: vcf
        path("${sample_id}.gridss.assembly.bam"), emit: assembly_bam

    script:
    def threads = task.cpus
    def memory = task.memory.toGiga()
    """
    gridss \\
        --reference ${fasta} \\
        --output ${sample_id}.gridss.vcf.gz \\
        --assembly ${sample_id}.gridss.assembly.bam \\
        --threads ${threads} \\
        --jvmheap ${memory}g \\
        ${alignment}

    tabix -p vcf ${sample_id}.gridss.vcf.gz
    """
}
