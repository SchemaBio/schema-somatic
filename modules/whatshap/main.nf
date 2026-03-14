// WhatsHap 单倍型分析模块
// 用途：计算单倍型并合并EGFR关键区域的相邻变异
// 包含的process：
//   WHATSHAP - 基于reads进行单倍型分型
//   EGFRHAP - 合并EGFR 19del和20ins区域的MNP
// 用法：
//   - 支持GRCh37和GRCh38两种参考基因组
//   - 只对EGFR特定区域进行变异合并，其他区域保持不变

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
    def coords = genome_build == 'GRCh38' ?
        [chr: 'chr7', del: '55174771-55174791', ins: '55181320-55181380'] :
        [chr: 'chr7', del: '55242464-55242484', ins: '55249071-55249131']
    """
    # 提取EGFR 19del区域
    bcftools view -r ${coords.chr}:${coords.del} ${vcf} -O z -o ${sample_id}.19del.vcf.gz
    tabix -p vcf ${sample_id}.19del.vcf.gz

    # 提取EGFR 20ins区域
    bcftools view -r ${coords.chr}:${coords.ins} ${vcf} -O z -o ${sample_id}.20ins.vcf.gz
    tabix -p vcf ${sample_id}.20ins.vcf.gz

    # 合并19del区域
    python /usr/local/bin/merge_mnp.py \\
        ${sample_id}.19del.vcf.gz \\
        ${fasta} \\
        --out_file ${sample_id}.19del.merged.vcf

    # 合并20ins区域
    python /usr/local/bin/merge_mnp.py \\
        ${sample_id}.20ins.vcf.gz \\
        ${fasta} \\
        --out_file ${sample_id}.20ins.merged.vcf

    # 提取其他区域
    bcftools view -T ^<(echo -e "${coords.chr}\\t${coords.del.split('-')[0]}\\t${coords.del.split('-')[1]}\\n${coords.chr}\\t${coords.ins.split('-')[0]}\\t${coords.ins.split('-')[1]}") ${vcf} -O z -o ${sample_id}.other.vcf.gz

    # 合并所有区域
    bcftools concat -a \\
        ${sample_id}.19del.merged.vcf \\
        ${sample_id}.20ins.merged.vcf \\
        ${sample_id}.other.vcf.gz \\
        -O z -o ${sample_id}.egfr.merged.vcf.gz

    tabix -p vcf ${sample_id}.egfr.merged.vcf.gz
    """
}
