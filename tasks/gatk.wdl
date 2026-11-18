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
        # gatk CollectInsertSizeMetrics \
        #     -I ~{bam} \
        #     -O ~{prefix}.insertsize.txt \
        #     -H ~{prefix}.histogram.pdf
        
        # gatk EstimateLibraryComplexity \
        #     -I ~{bam} \
        #     -O ~{prefix}.complexity_metrics.txt

        gatk CollectAlignmentSummaryMetrics \
            -I ~{bam} \
            -R ~{ref_dir}/~{fasta} \
            -O ~{prefix}.metrics.txt

        gatk BedToIntervalList \
            -I ~{bed} \
            -O ~{prefix}.interval_list \
            -SD ~{ref_dir}/~{fasta}

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

Mutect2Pair {
    input {
        String prefix
        File tumor_bam
        File tumor_bai
        String tumor_sample_name
        File normal_bam
        File normal_bai
        String normal_sample_name
        File ref_fasta
        File germline_resource_vcf
        File pon_vcf
    }

    command <<<
        set -euo pipefail

        gatk Mutect2 \
            -R ~{ref_fasta} \
            -I ~{tumor_bam} \
            -I ~{normal_bam} \
            -tumor ~{tumor_sample_name} \
            -normal ~{normal_sample_name} \
            --germline-resource ~{germline_resource_vcf} \
            --panel-of-normals ~{pon_vcf} \
            -O ~{prefix}.unfiltered.vcf
    >>>

    output {
        File unfiltered_vcf = "~{prefix}.unfiltered.vcf"
        File unfiltered_vcf_index = "~{prefix}.unfiltered.vcf.idx"
    }
}


Mutect2 {

}

