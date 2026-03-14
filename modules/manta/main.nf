// Manta 结构变异检测模块
// 用途：检测结构变异(SV)，包括缺失、插入、倒位、易位等
// 用法：
//   - 输入：样本ID、BAM/CRAM文件及索引、参考基因组
//   - 输出：候选SV和高置信度SV的VCF文件
//   - 单样本模式，适用于种系变异检测
//   - 支持BAM和CRAM格式

process MANTA {
    tag "MANTA on $sample_id"
    label 'manta'
    publishDir "${params.output}/03.Mutations/SV", mode: 'copy'

    input:
        tuple val(sample_id), path(alignment), path(index)
        path fasta
        path fasta_index

    output:
        tuple val(sample_id), path("${sample_id}.manta.vcf.gz"), path("${sample_id}.manta.vcf.gz.tbi"), emit: vcf
        path("${sample_id}.manta.candidateSV.vcf.gz"), emit: candidate_vcf

    script:
    def threads = task.cpus
    """
    configManta.py \\
        --bam ${alignment} \\
        --referenceFasta ${fasta} \\
        --runDir manta_work

    manta_work/runWorkflow.py \\
        -m local \\
        -j ${threads}

    mv manta_work/results/variants/diploidSV.vcf.gz ${sample_id}.manta.vcf.gz
    mv manta_work/results/variants/diploidSV.vcf.gz.tbi ${sample_id}.manta.vcf.gz.tbi
    mv manta_work/results/variants/candidateSV.vcf.gz ${sample_id}.manta.candidateSV.vcf.gz
    """
}
