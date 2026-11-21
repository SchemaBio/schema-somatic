version 1.2

task Manta {
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
        normal_args=()
        result=tumorSV
        if [[ -e normal.bam ]]; then
            normal_args=(--normalBam normal.bam)
            result=somaticSV
        fi
        configManta.py --exome --referenceFasta "~{ref_dir}/~{fasta}" \
            --tumorBam tumor.bam "${normal_args[@]}" --runDir manta
        manta/runWorkflow.py -m local -j ~{threads}
        cp "manta/results/variants/$result.vcf.gz" "~{prefix}.manta.vcf.gz"
        cp "manta/results/variants/$result.vcf.gz.tbi" "~{prefix}.manta.vcf.gz.tbi"
    >>>

    output {
        File vcf = "~{prefix}.manta.vcf.gz"
        File index = "~{prefix}.manta.vcf.gz.tbi"
    }

    requirements {
        container: "schemabio/manta:1.6.0"
        cpu: threads
        memory: "32 GiB"
    }
}
