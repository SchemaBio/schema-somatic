// MSIsensor-pro 微卫星不稳定性检测模块
// 用途：检测肿瘤样本的微卫星不稳定性(MSI)状态
// 用法：
//   - 输入：样本ID、肿瘤BAM/CRAM、正常BAM/CRAM（可选）、微卫星位点列表
//   - 输出：MSI评分和分类结果
//   - 支持单样本和配对样本模式

process MSISENSORPRO_SCAN {
    tag "MSISENSORPRO_SCAN"
    label 'msisensorpro'
    publishDir "${params.output}/MSI_Reference", mode: 'copy'

    input:
        path fasta
        path fasta_index

    output:
        path("microsatellites.list"), emit: msi_list

    script:
    """
    msisensor-pro scan -d ${fasta} -o microsatellites.list
    """
}

process MSISENSORPRO {
    tag "MSISENSORPRO on $sample_id"
    label 'msisensorpro'
    publishDir "${params.output}/03.Mutations/MSI", mode: 'copy'

    input:
        tuple val(sample_id), path(tumor_bam), path(tumor_bai)
        path normal_bam
        path normal_bai
        path msi_list

    output:
        path("${sample_id}_msisensor"), emit: msi_result
        path("${sample_id}_msisensor_dis"), emit: msi_distribution

    script:
    def normal_arg = normal_bam.name != 'NO_FILE' ? "-n ${normal_bam}" : ""
    def threads = task.cpus
    """
    msisensor-pro msi \\
        -d ${msi_list} \\
        -t ${tumor_bam} \\
        ${normal_arg} \\
        -o ${sample_id}_msisensor \\
        -b ${threads}
    """
}

process MSISENSORPRO_BASELINE {
    tag "MSISENSORPRO_BASELINE"
    label 'msisensorpro'
    publishDir "${params.output}/MSI_Reference", mode: 'copy'

    input:
        path normal_bams
        path normal_bais
        path msi_list

    output:
        path("baseline"), emit: baseline

    script:
    def threads = task.cpus
    def bam_list = normal_bams.collect { it.toString() }.join(' ')
    """
    msisensor-pro baseline \\
        -d ${msi_list} \\
        -i ${bam_list} \\
        -o baseline \\
        -b ${threads}
    """
}
