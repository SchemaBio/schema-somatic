// WhatsHap 单倍型分析模块
// 用途：计算单倍型并合并EGFR关键区域的相邻变异
// 包含的process：
//   WHATSHAP - 基于reads进行单倍型分型
//   EGFRHAP - 合并EGFR 19del和20ins区域的MNP（跳过T790M和C797S）
// 用法：
//   - 支持GRCh37和GRCh38两种参考基因组
//   - 只对EGFR 19del和20ins区域进行变异合并
//   - 自动跳过T790M和C797S位点的合并，保留为独立变异

process WHATSHAP {
    tag "WHATSHAP on $sample_id"
    label 'whatshap'
    publishDir "${params.output}/03.Mutations/SNV_InDel", mode: 'copy'

    input:
        tuple val(sample_id), path(vcf), path(vcf_index), path(bam), path(bai)
        path fasta
        path fasta_index

    output:
        tuple val(sample_id), path("${sample_id}.phased.vcf.gz"), path("${sample_id}.phased.vcf.gz.tbi"), emit: phased_vcf

    script:
    """
    whatshap phase \\
        --reference ${fasta} \\
        --output ${sample_id}.phased.vcf.gz \\
        ${vcf} \\
        ${bam}

    tabix -p vcf ${sample_id}.phased.vcf.gz
    """
}

process EGFRHAP {
    tag "EGFRHAP on $sample_id"
    label 'egfrhap'
    publishDir "${params.output}/03.Mutations/SNV_InDel", mode: 'copy'

    input:
        tuple val(sample_id), path(vcf), path(vcf_index)
        path fasta
        path fasta_index
        val genome_build  // 'GRCh37' or 'GRCh38'

    output:
        tuple val(sample_id), path("${sample_id}.egfr.merged.vcf.gz"), path("${sample_id}.egfr.merged.vcf.gz.tbi"), emit: merged_vcf

    script:
    """
    # 合并EGFR 19del和20ins区域的MNP，跳过T790M和C797S
    python /usr/local/bin/merge_egfr_mnp.py \\
        ${vcf} \\
        ${fasta} \\
        --out_file ${sample_id}.egfr.merged.vcf.gz \\
        --genome ${genome_build}
    """
}
