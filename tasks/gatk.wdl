version 1.2

task MergeBAMAlignment {
    input {
        String prefix
        File unmapped_bam
        File aligned_bam
        String ref_fasta
        Directory ref_dir
    }

    command <<<
        set -euo pipefail

        gatk MergeBamAlignment \
            -R "~{ref_dir}/~{ref_fasta}" \
            -UNMAPPED ~{unmapped_bam} \
            -ALIGNED ~{aligned_bam} \
            -O ~{prefix}.merged.bam \
            --ALIGNER_PROPER_PAIR_FLAGS true \
            --ATTRIBUTES_TO_RETAIN XS \
            --CREATE_INDEX true
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

task CollectQCMetrics {
    input {
        String prefix
        File bam
        File bai
        File bed
        String fasta
        Int threads
        Directory ref_dir
    }

    Int memory_gb = threads * 2

    command <<<
        set -euo pipefail

        gatk CollectAlignmentSummaryMetrics \
            -I ~{bam} \
            -R ~{ref_dir}/~{fasta} \
            -O ~{prefix}.metrics.txt

        gatk BedToIntervalList \
            -I ~{bed} \
            -O ~{prefix}.interval_list \
            -SD "~{ref_dir}/~{fasta}"

        gatk CollectHsMetrics \
            -BI ~{prefix}.interval_list \
            -TI ~{prefix}.interval_list \
            -I ~{bam} \
            -O ~{prefix}.hs.txt
    >>>

    output {
        File summary = "~{prefix}.metrics.txt"
        File hs_metric = "~{prefix}.hs.txt"
    }

    runtime {
        cpu: threads
        memory: "~{memory_gb}G"
        docker: "schemabio/gatk:4.6.2.0"
    }
}


task Mutect2 {
    input {
        String prefix
        String fasta
        Directory ref_dir
        File tumor_bam
        File tumor_bai
        File? normal_bam
        File? normal_bai
        String? normal_sample_name
        File bed
        File germline_resource_vcf
        File germline_resource_index
        File pon_vcf
        File pon_index
    }

    command <<<
        set -euo pipefail
        ln -s "~{tumor_bam}" tumor.bam
        ln -s "~{tumor_bai}" tumor.bam.bai
        if [[ -n "~{default="" normal_bam}" ]]; then
            ln -s "~{default="" normal_bam}" normal.bam
            ln -s "~{default="" normal_bai}" normal.bam.bai
        fi
        ln -s "~{germline_resource_vcf}" germline.vcf.gz
        ln -s "~{germline_resource_index}" germline.vcf.gz.tbi
        ln -s "~{pon_vcf}" pon.vcf.gz
        ln -s "~{pon_index}" pon.vcf.gz.tbi
        normal_args=()
        if [[ -e normal.bam ]]; then
            normal_args=(-I normal.bam -normal "~{default="" normal_sample_name}")
        fi
        gatk Mutect2 -R "~{ref_dir}/~{fasta}" -L "~{bed}" \
            -I tumor.bam "${normal_args[@]}" --germline-resource germline.vcf.gz \
            --panel-of-normals pon.vcf.gz -O "~{prefix}.unfiltered.vcf.gz"
    >>>

    output {
        File vcf = "~{prefix}.unfiltered.vcf.gz"
        File index = "~{prefix}.unfiltered.vcf.gz.tbi"
        File stats = "~{prefix}.unfiltered.vcf.gz.stats"
    }

    requirements {
        container: "schemabio/gatk:4.6.2.0"
        cpu: 2
        memory: "8 GiB"
    }
}

task FilterMutectCalls {
    input {
        String prefix
        String fasta
        Directory ref_dir
        File vcf
        File index
        File stats
    }

    command <<<
        set -euo pipefail
        ln -s "~{vcf}" calls.vcf.gz
        ln -s "~{index}" calls.vcf.gz.tbi
        gatk FilterMutectCalls -R "~{ref_dir}/~{fasta}" -V calls.vcf.gz \
            --stats "~{stats}" -O "~{prefix}.filtered.vcf.gz"
    >>>

    output {
        File filtered_vcf = "~{prefix}.filtered.vcf.gz"
        File filtered_index = "~{prefix}.filtered.vcf.gz.tbi"
    }

    requirements {
        container: "schemabio/gatk:4.6.2.0"
        memory: "8 GiB"
    }
}
