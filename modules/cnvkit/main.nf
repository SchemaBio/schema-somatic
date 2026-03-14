// CNVkit 拷贝数变异检测模块
// 用途：检测拷贝数变异(CNV)，支持WGS和靶向测序
// 用法：
//   - 输入：样本ID、BAM/CRAM文件及索引、参考基因组、目标区域BED（可选）
//   - 输出：CNV结果文件(.cns)、分段文件(.cnr)、可视化图表
//   - 单样本模式，可使用参考样本或平坦参考

process CNVKIT {
    tag "CNVKIT on $sample_id"
    label 'cnvkit'
    publishDir "${params.output}/03.Mutations/CNV", mode: 'copy'

    input:
        tuple val(sample_id), path(alignment), path(index)
        path fasta
        path fasta_index
        path target_bed  // optional, use 'NO_FILE' for WGS
        path antitarget_bed  // 已调整好的antitarget格式bed
        path baseline

    output:
        path("${sample_id}.call.cns"), emit: call_result
        path("${sample_id}.seg.call.cns"), emit: segment_result

    script:
    def threads = task.cpus
    def target_arg = target_bed.name != 'NO_FILE' ? "--targets ${target_bed}" : ""
    """
    cnvkit.py coverage \\
        ${alignment} ${target} -o ${sample_id}.targetcoverage.cnn -p ${threads}
    cnvkit.py coverage \\
        ${alignment} ${antitarget_bed} -o ${sample_id}.antitargetcoverage.cnn -p ${threads}
    cnvkit.py fix ${sample_id}.targetcoverage.cnn ${sample_id}.antitargetcoverage.cnn \\
        ${baseline} -o ${sample_id}.cnr
    cnvkit.py segment ${sample_id}.cnr -o ${sample_id}.seg.cns -p ${threads}
    cnvkit.py call ${sample_id}.cnr -o ${sample_id}.call.cns
    cnvkit.py call ${sample_id}.seg.cns -o ${sample_id}.seg.call.cns.tmp
    head -1 ${sample_id}.seg.call.cns.tmp > ${sample_id}.seg.call.cns
    tail -n+2 ${sample_id}.seg.call.cns.tmp \\
        | sort -t $'\t' -k 6 -r -n \\
        | awk '{if($4 != "-") print $0}' >> ${sample_id}.seg.call.cns
    """
}

process CNVKITBASELINE {
    tag "CNVKITBASELINE"
    label 'cnvkitbaseline'
    publishDir "${params.output}/CNVkit_Reference", mode: 'copy'

    input:
        path normal_bams  // 多个正常样本的BAM/CRAM文件
        path normal_bais  // 对应的索引文件
        path fasta
        path fasta_index
        path target_bed      // 已调整好的target格式bed
        path antitarget_bed  // 已调整好的antitarget格式bed
        val  prefix

    output:
        path("${prefix}.cnv.reference.cnn"), emit: baseline

    script:
    def threads = task.cpus
    """
    # 对每个正常样本计算coverage
    for bam in ${normal_bams}; do
        cnvkit.py coverage \$bam ${target_bed} -o \${bam%.bam}.targetcoverage.cnn -p ${threads}
        cnvkit.py coverage \$bam ${antitarget_bed} -o \${bam%.bam}.antitargetcoverage.cnn -p ${threads}
    done

    # 创建基线参考
    cnvkit.py reference *.targetcoverage.cnn *.antitargetcoverage.cnn \\
        --fasta ${fasta} \\
        -o ${prefix}.cnv.reference.cnn
    """
}

process CNVKITNORMAL {
    tag "CNVKITBASELINE"
    label 'cnvkitnormal'
    publishDir "${params.output}/03.Mutations/CNV", mode: 'copy'

    input:
        path normal_bam  // 多个正常样本的BAM/CRAM文件
        path normal_bai  // 对应的索引文件
        path fasta
        path fasta_index
        path target_bed      // 已调整好的target格式bed
        path antitarget_bed  // 已调整好的antitarget格式bed
        val  prefix

    output:
        path("${prefix}.cnv.reference.cnn"), emit: baseline

    script:
    def threads = task.cpus
    """
    cnvkit.py coverage \\
        ${normal_bam} ${target_bed} -o ${prefix}.targetcoverage.cnn -p ${threads}
    cnvkit.py coverage \\
        ${normal_bam} ${antitarget_bed} -o ${prefix}.antitargetcoverage.cnn -p ${threads}
    cnvkit.py reference \\
        ${prefix}.targetcoverage.cnn ${prefix}.antitargetcoverage.cnn \\
        --fasta ${fasta} -o ${prefix}.cnv.reference.cnn
    """
}

process CNVKITSEX {
    tag "CNVKITSEX"
    label 'cnvkitsex'
    publishDir "${params.output}/01.QC", mode: 'copy'

    input:
        path cnr
        val  sample_id

    output:
        path("${sample_id}.gender.txt"), emit: cnvkit_sex

    script:
    def threads = task.cpus
    """
    cnvkit.py sex ${cnr} -o ${sample_id}.gender.txt
    """
}
