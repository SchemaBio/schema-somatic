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

// BCFtools MPileup 模块
// 用途：快速统计给定VCF坐标的深度和基因型信息（胚系位点）
// 应用场景：
//   1. 基因组范围内的BAF绘图
//   2. 化疗位点的基因型确认
//   3. 样本识别位点的基因型确认
// 参数说明：
//   - max_depth: 最大深度限制，默认1000X，避免高深度区域计算过慢
// 输入：
//   - BAM/CRAM文件及其索引
//   - 参考基因组及其索引
//   - 目标位点VCF文件（包含待检测的坐标）
// 输出：
//   - 包含深度和基因型信息的TSV文件

process MPILEUP {
    tag "MPILEUP on $sample_id"
    label 'bcftools'

    input:
        tuple val(sample_id), path(bam), path(bai)
        path fasta
        path fasta_index  // .fai 和 .dict 文件
        path target_vcf   // 目标位点VCF文件

    output:
        tuple val(sample_id), path("${sample_id}.pileup.tsv"), emit: pileup_tsv

    script:
    def max_depth = params.max_depth ?: 1000
    """
    # 使用bcftools mpileup获取目标位点的深度信息
    # -a DP,AD: 输出总深度(DP)和等位基因深度(AD)
    # -d ${max_depth}: 最大深度限制为1000X，适用于胚系位点
    # -Q 0: 最低碱基质量设为0，不过滤低质量碱基
    # -q 0: 最低映射质量设为0
    # -f: 参考基因组
    # -T: 目标位点文件
    bcftools mpileup \\
        -a DP,AD \\
        -d ${max_depth} \\
        -Q 0 \\
        -q 0 \\
        -f ${fasta} \\
        -T ${target_vcf} \\
        ${bam} \\
        -O u \\
    | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%DP\\t%AD]\\n' \\
    > ${sample_id}.pileup.tsv

    # 添加表头
    echo -e "CHROM\\tPOS\\tREF\\tALT\\tDP\\tAD" > ${sample_id}.pileup.tsv.header
    cat ${sample_id}.pileup.tsv >> ${sample_id}.pileup.tsv.header
    mv ${sample_id}.pileup.tsv.header ${sample_id}.pileup.tsv
    """
}

// MPILEUP 批量处理模块（带基因型calling）
// 用途：批量处理多个样本的MPileup分析并calling基因型
// 适用于需要比较多样本的场景（如样本识别位点确认）
// 参数说明：
//   - max_depth: 最大深度限制，默认1000X

process MPILEUP_BATCH {
    tag "MPILEUP_BATCH on ${sample_id}"
    label 'bcftools'

    input:
        tuple val(sample_id), path(bam), path(bai)
        path fasta
        path fasta_index
        path target_vcf

    output:
        tuple val(sample_id), path("${sample_id}.genotype.tsv"), emit: genotype_tsv

    script:
    def max_depth = params.max_depth ?: 1000
    """
    # 生成mpileup并计算基因型
    # -d ${max_depth}: 最大深度限制为1000X
    bcftools mpileup \\
        -a DP,AD,SP \\
        -d ${max_depth} \\
        -Q 0 \\
        -q 0 \\
        -f ${fasta} \\
        -T ${target_vcf} \\
        ${bam} \\
        -O u \\
    | bcftools call -mv -O u \\
    | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT\\t%DP\\t%AD]\\n' \\
    > ${sample_id}.genotype.tmp

    # 添加表头并保存
    echo -e "CHROM\\tPOS\\tREF\\tALT\\tGT\\tDP\\tAD" > ${sample_id}.genotype.tsv
    cat ${sample_id}.genotype.tmp >> ${sample_id}.genotype.tsv
    """
}

// MPILEUP 简化版（仅统计深度，不calling）
// 用途：快速统计目标位点的覆盖深度
// 适用于BAF绘图等只需深度信息的场景

process MPILEUP_DEPTH {
    tag "MPILEUP_DEPTH on $sample_id"
    label 'bcftools'

    input:
        tuple val(sample_id), path(bam), path(bai)
        path fasta
        path fasta_index
        path target_bed   // 目标区域BED文件

    output:
        tuple val(sample_id), path("${sample_id}.depth.tsv"), emit: depth_tsv

    script:
    def max_depth = params.max_depth ?: 1000
    """
    # 使用bcftools mpileup获取目标区域的深度信息
    # 输出格式: CHROM POS REF ALT DP
    bcftools mpileup \\
        -a DP \\
        -d ${max_depth} \\
        -Q 0 \\
        -q 0 \\
        -f ${fasta} \\
        -R ${target_bed} \\
        ${bam} \\
        -O u \\
    | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%DP\\n' \\
    > ${sample_id}.depth.tmp

    # 添加表头
    echo -e "CHROM\\tPOS\\tREF\\tDP" > ${sample_id}.depth.tsv
    cat ${sample_id}.depth.tmp >> ${sample_id}.depth.tsv
    """
}