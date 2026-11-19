version 1.2

task Manta {
    input {
        String prefix
        File tumor_bam
        File tumor_bai
        File normal_bam
        File normal_bai
        File fasta
        Directory ref_dir
        Int threads
    }

    command {
        set -euo pipefail

        configManta.py \
            --normalBam ${normal_bam} \
            --tumorBam ${tumor_bam} \
            --referenceFasta ${fasta} \
            --runDir ${prefix}.manta_run

        ${prefix}.manta_run/runWorkflow.py -m local -j ${threads}

        zcat ${prefix}.manta_run/results/variants/diploidSV.vcf.gz | bgzip -c > ${prefix}.manta.vcf.gz
    }

    output {
        File vcf = "${prefix}.manta.vcf.gz"
    }

    runtime {
        docker: "schemabio/manta:1.6.0"
        memory: "32G"
        cpu: threads
    }
}