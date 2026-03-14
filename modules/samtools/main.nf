// Samtools 工具模块
// 用途：BAM/CRAM文件索引和格式转换
// 用法：
//   INDEXALIGNMENT - 为比对文件创建索引（.bai或.crai）
//   BAM2FASTQ - 将BAM文件转换回FASTQ格式
//   SRYSEX - 基于SRY基因判断性别
//   EXTRACTHLA - 提取HLA区域reads

process INDEXALIGNMENT {
    tag "INDEXALIGNMENT on $sample_id"
    label 'indexalignment'
    if (publish_dir) {
        publishDir "${publish_dir}", mode: 'copy'
    }

    input:
        tuple val(sample_id), path(alignment_file)
        val(publish_dir)

    output:
        path("*.bai"), emit: bai, optional: true
        path("*.crai"), emit: crai, optional: true

    script:
    def threads = task.cpus
    def ext = alignment_file.extension
    if (ext == 'cram') {
        """
        samtools index -@ ${threads} -C ${alignment_file}
        """
    } else {
        """
        samtools index -@ ${threads} ${alignment_file}
        """
    }
}

process BAM2FASTQ {
    tag "BAM2FASTQ on $sample_id"
    label 'bam2fastq'
    if (publish_dir) {
        publishDir "${publish_dir}", mode: 'copy'
    }

    input:
        tuple val(sample_id), path(bam)
        val(publish_dir)

    output:
        path("${sample_id}.clean.umi_*.fq.gz"), emit: reads

    script:
    """
    samtools fastq ${bam} \\
        -1 ${sample_id}_1.fq.gz \\
        -2 ${sample_id}_2.fq.gz
    """
}

process SRYSEX {
    tag "SRYSEX on $sample_id"
    label 'srysex'
    publishDir "${params.output}/01.QC", mode: 'copy'

    input:
        tuple val(sample_id), path(alignment), path(index)
        path fasta
        path fasta_index
        val genome_build  // 'GRCh37' or 'GRCh38'

    output:
        path("${sample_id}.sry.sex.txt"), emit: sry_sex

    script:
    def coords = genome_build == 'GRCh38' ? 'chrY:2786855-2787699' : 'chrY:2654896-2655740'
    """
    # 统计SRY区域reads数
    sry_reads=\$(samtools view -c ${alignment} ${coords})

    # 统计总reads数
    total_reads=\$(samtools view -c ${alignment})

    # 计算比例并判断性别
    if [ \$total_reads -gt 0 ]; then
        ratio=\$(echo "scale=6; \$sry_reads / \$total_reads" | bc)
        if (( \$(echo "\$ratio > 0.00001" | bc -l) )); then
            sex="Male"
        else
            sex="Female"
        fi
    else
        sex="Unknown"
        ratio=0
    fi

    echo -e "Sample\tSRY_reads\tTotal_reads\tRatio\tSex" > ${sample_id}.sry.sex.txt
    echo -e "${sample_id}\t\$sry_reads\t\$total_reads\t\$ratio\t\$sex" >> ${sample_id}.sry.sex.txt
    """
}

process EXTRACTHLA {
    tag "EXTRACTHLA on $sample_id"
    label 'extracthla'
    publishDir "${params.output}/HLA_Reads", mode: 'copy'

    input:
        tuple val(sample_id), path(alignment), path(index)
        path fasta
        path fasta_index
        val genome_build  // 'GRCh37' or 'GRCh38'

    output:
        tuple val(sample_id), path("${sample_id}.hla_1.fq.gz"), path("${sample_id}.hla_2.fq.gz"), emit: hla_reads

    script:
    def region = genome_build == 'GRCh38' ? 'chr6:28510120-33480577' : 'chr6:28477797-33448354'
    def threads = task.cpus
    """
    samtools view -@ ${threads} -b ${alignment} ${region} | \\
    samtools sort -@ ${threads} -n -o ${sample_id}.hla.sorted.bam - && \\
    samtools fastq -@ ${threads} \\
        -1 ${sample_id}.hla_1.fq.gz \\
        -2 ${sample_id}.hla_2.fq.gz \\
        -0 /dev/null -s /dev/null -n ${sample_id}.hla.sorted.bam
    """
}
