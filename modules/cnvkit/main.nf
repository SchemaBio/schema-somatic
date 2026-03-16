// CNVkit 拷贝数变异检测模块
// 用途：检测拷贝数变异(CNV)，支持WGS和靶向测序
// 用法：
//   - 输入：样本ID、BAM/CRAM文件及索引、参考基因组版本参数
//   - 输出：CNV结果文件(.cns)、分段文件(.cnr)、可视化图表
//   - 单样本模式，可使用参考样本或平坦参考
//   - 预置BED文件位于 assets/cnvkit/ 目录
//   - 使用方式：params.cnvkit_target_bed 和 params.cnvkit_antitarget_bed 指定预置文件

process CNVKIT {
    tag "CNVKIT on $sample_id"
    label 'cnvkit'
    publishDir "${params.output}/03.Mutations/CNV", mode: 'copy'

    input:
        tuple val(sample_id), path(alignment), path(index)
        path fasta
        path fasta_index
        val genome_version  // 参考基因组版本：GRCh38 或 GRCh37
        path user_bed  // optional, 用户自定义的target bed（只和target取交集）
        path baseline

    output:
        path("${sample_id}.call.cns"), emit: call_result
        path("${sample_id}.seg.call.cns"), emit: segment_result

    script:
    def threads = task.cpus
    def genome_ver = genome_version ?: 'GRCh38'

    // 根据genome_version确定预置bed文件路径
    def preset_target_beds = params.cnvkit_target_beds ?: ['GRCh38': "${projectDir}/assets/cnvkit/GRCh38_target.bed", 'GRCh37': "${projectDir}/assets/cnvkit/GRCh37_target.bed"]
    def preset_antitarget_beds = params.cnvkit_antitarget_beds ?: ['GRCh38': "${projectDir}/assets/cnvkit/GRCh38_antitarget.bed", 'GRCh37': "${projectDir}/assets/cnvkit/GRCh37_antitarget.bed"]

    def preset_target = preset_target_beds[genome_ver] ?: preset_target_beds['GRCh38']
    def preset_antitarget = preset_antitarget_beds[genome_ver] ?: preset_antitarget_beds['GRCh38']

    // 确定最终使用的target bed
    def target_bed_file = "${sample_id}.target.bed"
    def antitarget_bed_file = preset_antitarget  // antitarget直接使用预置

    // 用户提供了自定义bed时，与预置target取交集；否则直接使用预置
    def target_bed_setup = ""
    if (user_bed) {
        target_bed_setup = """
        # 用户自定义bed与预置target取交集，保留预置bed的完整格式(-wb)
        COLS=\$(head -n1 ${preset_target} | awk '{print NF}')
        START_COL=\$((\$COLS + 1))
        bedtools intersect -a ${user_bed} -b ${preset_target} -wb | cut -f\$START_COL- > ${target_bed_file}
        """
    } else {
        target_bed_setup = """
        cp ${preset_target} ${target_bed_file}
        """
    }
    """
    # 设置target bed
    ${target_bed_setup}

    cnvkit.py coverage \\
        ${alignment} ${target_bed_file} -o ${sample_id}.targetcoverage.cnn -p ${threads}
    cnvkit.py coverage \\
        ${alignment} ${antitarget_bed_file} -o ${sample_id}.antitargetcoverage.cnn -p ${threads}
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
        val genome_version  // 参考基因组版本：GRCh38 或 GRCh37
        path user_bed  // optional, 用户自定义的bed
        val  prefix

    output:
        path("${prefix}.cnv.reference.cnn"), emit: baseline

    script:
    def threads = task.cpus
    def genome_ver = genome_version ?: 'GRCh38'

    // 根据genome_version确定预置bed文件路径
    def preset_target_beds = params.cnvkit_target_beds ?: ['GRCh38': "${projectDir}/assets/cnvkit/GRCh38_target.bed", 'GRCh37': "${projectDir}/assets/cnvkit/GRCh37_target.bed"]
    def preset_antitarget_beds = params.cnvkit_antitarget_beds ?: ['GRCh38': "${projectDir}/assets/cnvkit/GRCh38_antitarget.bed", 'GRCh37': "${projectDir}/assets/cnvkit/GRCh37_antitarget.bed"]

    def preset_target = preset_target_beds[genome_ver] ?: preset_target_beds['GRCh38']
    def preset_antitarget = preset_antitarget_beds[genome_ver] ?: preset_antitarget_beds['GRCh38']

    // 确定最终使用的target bed
    def target_bed_file = "target.bed"
    def antitarget_bed_file = preset_antitarget  // antitarget直接使用预置

    // 用户提供了自定义bed时，与预置target取交集；否则直接使用预置
    def target_bed_setup = ""
    if (user_bed) {
        target_bed_setup = """
        # 用户自定义bed与预置target取交集，保留预置bed的完整格式(-wb)
        COLS=\$(head -n1 ${preset_target} | awk '{print NF}')
        START_COL=\$((\$COLS + 1))
        bedtools intersect -a ${user_bed} -b ${preset_target} -wb | cut -f\$START_COL- > ${target_bed_file}
        """
    } else {
        target_bed_setup = """
        cp ${preset_target} ${target_bed_file}
        """
    }
    """
    # 设置target bed
    ${target_bed_setup}

    # 对每个正常样本计算coverage
    for bam in ${normal_bams}; do
        cnvkit.py coverage \$bam ${target_bed_file} -o \${bam%.bam}.targetcoverage.cnn -p ${threads}
        cnvkit.py coverage \$bam ${antitarget_bed_file} -o \${bam%.bam}.antitargetcoverage.cnn -p ${threads}
    done

    # 创建基线参考
    cnvkit.py reference *.targetcoverage.cnn *.antitargetcoverage.cnn \\
        --fasta ${fasta} \\
        -o ${prefix}.cnv.reference.cnn
    """
}

process CNVKITNORMAL {
    tag "CNVKITNORMAL"
    label 'cnvkitnormal'
    publishDir "${params.output}/03.Mutations/CNV", mode: 'copy'

    input:
        path normal_bam  // 正常样本的BAM/CRAM文件
        path normal_bai  // 对应的索引文件
        path fasta
        path fasta_index
        val genome_version  // 参考基因组版本：GRCh38 或 GRCh37
        path user_bed  // optional, 用户自定义的bed
        val  prefix

    output:
        path("${prefix}.cnv.reference.cnn"), emit: baseline

    script:
    def threads = task.cpus
    def genome_ver = genome_version ?: 'GRCh38'

    // 根据genome_version确定预置bed文件路径
    def preset_target_beds = params.cnvkit_target_beds ?: ['GRCh38': "${projectDir}/assets/cnvkit/GRCh38_target.bed", 'GRCh37': "${projectDir}/assets/cnvkit/GRCh37_target.bed"]
    def preset_antitarget_beds = params.cnvkit_antitarget_beds ?: ['GRCh38': "${projectDir}/assets/cnvkit/GRCh38_antitarget.bed", 'GRCh37': "${projectDir}/assets/cnvkit/GRCh37_antitarget.bed"]

    def preset_target = preset_target_beds[genome_ver] ?: preset_target_beds['GRCh38']
    def preset_antitarget = preset_antitarget_beds[genome_ver] ?: preset_antitarget_beds['GRCh38']

    // 确定最终使用的target bed
    def target_bed_file = "target.bed"
    def antitarget_bed_file = preset_antitarget  // antitarget直接使用预置

    // 用户提供了自定义bed时，与预置target取交集；否则直接使用预置
    def target_bed_setup = ""
    if (user_bed) {
        target_bed_setup = """
        # 用户自定义bed与预置target取交集，保留预置bed的完整格式(-wb)
        COLS=\$(head -n1 ${preset_target} | awk '{print NF}')
        START_COL=\$((\$COLS + 1))
        bedtools intersect -a ${user_bed} -b ${preset_target} -wb | cut -f\$START_COL- > ${target_bed_file}
        """
    } else {
        target_bed_setup = """
        cp ${preset_target} ${target_bed_file}
        """
    }
    """
    # 设置target bed
    ${target_bed_setup}

    cnvkit.py coverage \\
        ${normal_bam} ${target_bed_file} -o ${prefix}.targetcoverage.cnn -p ${threads}
    cnvkit.py coverage \\
        ${normal_bam} ${antitarget_bed_file} -o ${prefix}.antitargetcoverage.cnn -p ${threads}
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
        path("${sample_id}.cnvkit.gender.txt"), emit: cnvkit_sex

    script:
    def threads = task.cpus
    """
    cnvkit.py sex ${cnr} -o ${sample_id}.cnvkit.gender.txt
    """
}
