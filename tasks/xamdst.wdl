version 1.2

task Xamdst {
    input {
        String prefix
        File bam
        File bai
        File bed
        Int threads
    }

    Int memory_gb = threads * 2

    command <<<
        set -euo pipefail
        ln -s "~{bam}" input.bam
        ln -s "~{bai}" input.bam.bai
        xamdst -p "~{bed}" -o result --threads ~{threads} --cutoffdepth 1000 input.bam
        cp result/coverage.report.json ~{prefix}.xamdst.report.json
    >>>

    output {
        File xamdst_report = "~{prefix}.xamdst.report.json"
    }

    runtime {
        cpu: threads
        memory: "~{memory_gb}G"
        docker: "schemabio/mapping:v1.1.0"
    }
}
