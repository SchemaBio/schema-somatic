version 1.2

task MsiSingle {
    input {
        String prefix
        File tumor_bam
        File tumor_bai
        File msi_baseline
        Int threads
    }

    command <<<
        set -euo pipefail
        ln -s "~{tumor_bam}" tumor.bam
        ln -s "~{tumor_bai}" tumor.bam.bai
        msisensor-pro pro -t tumor.bam -d "~{msi_baseline}" -b ~{threads} -o "~{prefix}.msi"
    >>>

    output {
        File score = "~{prefix}.msi"
    }

    requirements {
        container: "quay.io/biocontainers/msisensor-pro:1.3.0--hd979922_1"
        cpu: threads
        memory: "8 GiB"
    }
}

task MsiPair {
    input {
        String prefix
        File tumor_bam
        File tumor_bai
        File normal_bam
        File normal_bai
        File msi_sites
        Int threads
    }

    command <<<
        set -euo pipefail
        ln -s "~{tumor_bam}" tumor.bam
        ln -s "~{tumor_bai}" tumor.bam.bai
        ln -s "~{normal_bam}" normal.bam
        ln -s "~{normal_bai}" normal.bam.bai
        msisensor-pro msi -n normal.bam -t tumor.bam -d "~{msi_sites}" -b ~{threads} -o "~{prefix}.msi"
    >>>

    output {
        File score = "~{prefix}.msi"
    }

    requirements {
        container: "quay.io/biocontainers/msisensor-pro:1.3.0--hd979922_1"
        cpu: threads
        memory: "8 GiB"
    }
}
