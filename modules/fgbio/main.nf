// FGBIO UMI处理模块
// 用途：提取和处理分子标签(UMI)用于去重
// 用法：
//   - 输入：样本ID、FASTQ文件、UMI结构（如"5M5S+T +T"）、参考基因组
//   - 输出：包含UMI信息的未比对BAM文件
//   - UMI结构说明：M=分子标签，S=跳过碱基，T=模板序列

process EXTRACTUMI {
    tag "EXTRACTUMI on $sample_id"
    label 'extractumi'
    publishDir "${params.output}/01.QC/Fgbio", mode: 'copy'

    input:
        tuple val(sample_id), path(read1), path(read2)
        val read_structure  // e.g. "5M5S+T +T"
        path fasta
        path fasta_index

    output:
        path("${sample_id}.unmapped.umi.bam"), emit: unmapped_bam

    script:
    """
    fgbio FastqToBam \\
        --input ${read1} ${read2} \\
        --read-structures ${read_structure} \\
        --output ${sample_id}.unmapped.bam

    fgbio ExtractUmisFromBam \\
        -i ${sample_id}.unmapped.bam \\
        -o ${sample_id}.unmapped.umi.bam \\
        -r ${read_structure} \\
        -s RX \\
        -t ZA ZB
    """
}