version 1.2

# fgbio UMI 处理模块
# 用途：将 FASTQ 转为未比对 BAM，并提取分子标签（UMI）供后续比对和去重使用
# UMI 结构示例："5M5S+T +T"，其中 M=分子标签，S=跳过碱基，T=模板序列

task FastqToBam {
    input {
        String prefix
        File read_1
        File read_2
        String read_structure_1
        String read_structure_2
    }

    command <<<
        set -euo pipefail

        fgbio FastqToBam \
            --input ~{read_1} ~{read_2} \
            --read-group-id "~{prefix}" --sample "~{prefix}" --library "~{prefix}" --sort true \
            --read-structures ~{read_structure_1} ~{read_structure_2} \
            --output ~{prefix}.unmapped.bam
    >>>

    output {
        File unmapped_bam = "~{prefix}.unmapped.bam"
    }

    runtime {
        cpu: 4
        memory: "8G"
        docker: "schemabio/fgbio:3.1.1"
    }
}
