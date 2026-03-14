// BCFtools 过滤模块
// 用途：根据等位基因深度(AD)和频率(AF)过滤VCF变异
// 用法：
//   - 输入：样本ID、VCF文件、最小AD值、最小AF值
//   - 输出：过滤后的VCF文件
//   - 过滤条件：AD[0:1] >= minAD 且 AF >= minAF

process BCFTOOLS {
    tag "BCFTOOLS on $sample_id"
    label 'bcftools'
    if (publish_dir) {
        publishDir "${publish_dir}", mode: 'copy'
    }

    input:
        tuple val(sample_id), path(vcf)
        val(publish_dir)
        val(minAD)
        val(minAF)

    output:
        path("${sample}.filtered.vcf"), emit: filtered_vcf

    script:
    """
    bcftools view -e "FORMAT/AD[0:1]<${minAD} || FORMAT/AF<${minAF}" ${vcf} > ${sample}.filtered.vcf
    """
}