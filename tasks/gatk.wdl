version 1.2

task MergeBAMAlignment {
    input {
        String prefix
        File unmapped_bam
        File aligned_bam
        File ref_fasta
        Directory ref_dir
    }

    command <<<
        set -euo pipefail

        gatk MergeBamAlignment \
            -R ~{ref_fasta} \
            -UNMAPPED ~{unmapped_bam} \
            -ALIGNED ~{aligned_bam} \
            -O ~{prefix}.merged.bam \
            --ALIGNER_PROPER_PAIR_FLAGS true \
            --ATTRIBUTES_TO_RETAIN XS
    >>>

    output {
        File merged_bam = "~{prefix}.merged.bam"
        File merged_bai = "~{prefix}.merged.bai"
    }

    runtime {
        cpu: 4
        memory: "16G"
        docker: "schemabio/gatk:4.6.2.0"
    }
}

task MarkDuplicates {
    input {
        String prefix
        File bam
        File bai
    }

    command <<<
        set -e

        mkdir -p ./tmp

        gatk MarkDuplicates \
            -I ~{bam} \
            -O ~{prefix}.markdup.bam \
            -M ~{prefix}.markdup.metrics.txt \
            --TMP_DIR ./tmp \
            --CREATE_INDEX true \
            --BARCODE_TAG RX
        
        mv ~{prefix}.markdup.bai ~{prefix}.markdup.bam.bai

    >>>

    output {
        File markdup_bam = "~{prefix}.markdup.bam"
        File markdup_bai = "~{prefix}.markdup.bam.bai"
        File metrics = "~{prefix}.markdup.metrics.txt"
    }

    runtime {
        docker: "schemabio/gatk:4.6.2.0"
    }
}


