version 1.2

task Gridss {
    input {
        String prefix
        String fasta
        Directory ref_dir
        File tumor_bam
        File tumor_bai
        File? normal_bam
        File? normal_bai
        Int threads
    }

    command <<<
        set -euo pipefail
        ln -s "~{tumor_bam}" tumor.bam
        ln -s "~{tumor_bai}" tumor.bam.bai
        if [[ -n "~{default="" normal_bam}" ]]; then
            ln -s "~{default="" normal_bam}" normal.bam
            ln -s "~{default="" normal_bai}" normal.bam.bai
        fi
        inputs=(tumor.bam)
        if [[ -e normal.bam ]]; then inputs=(normal.bam tumor.bam); fi
        gridss --reference "~{ref_dir}/~{fasta}" --output "~{prefix}.gridss.vcf" \
            --assembly "~{prefix}.gridss.assembly.bam" --threads ~{threads} "${inputs[@]}"
    >>>

    output {
        File vcf = "~{prefix}.gridss.vcf"
        File assembly_bam = "~{prefix}.gridss.assembly.bam"
    }

    requirements {
        container: "schemabio/gridss:2.13.1"
        cpu: threads
        memory: "32 GiB"
    }
}
