// OptiType HLA分型模块
// 用途：基于RNA-seq或WES数据进行HLA分型
// 用法：
//   - 输入：样本ID、HLA区域的FASTQ文件
//   - 输出：HLA分型结果
//   - 自动调整线程数配置

process OPTITYPE {
    tag "OPTITYPE on $sample_id"
    label 'optitype'
    publishDir "${params.output}/03.Mutations/HLA", mode: 'copy'

    input:
        tuple val(sample_id), path(read1), path(read2)

    output:
        path("${sample_id}.optitype.txt"), emit: hla_result
        path("${sample_id}_coverage_plot.pdf"), emit: coverage_plot, optional: true

    script:
    def threads = task.cpus
    """
    mkdir ${sample_id}_tmp
    sed -i 's/threads=1/threads=${threads}/g' /usr/local/bin/OptiType/config.ini
    OptiTypePipeline.py \\
        -i ${read1} ${read2} \\
        -d -o ${sample_id}_tmp -p ${sample_id} -v
    cp ${sample_id}_tmp/${sample_id}_result.tsv ${sample_id}.optitype.txt
    [ -f ${sample_id}_tmp/${sample_id}_coverage_plot.pdf ] && cp ${sample_id}_tmp/${sample_id}_coverage_plot.pdf .
    """
}
