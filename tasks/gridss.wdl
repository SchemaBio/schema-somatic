version 1.2

task Gridss {
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

        gridss \
            --reference ${fasta} \
            --output ${prefix}.gridss.vcf.gz \
            --assembly ${prefix}.gridss.assembly.bam \
            --threads ${threads} \
            --tumor ${tumor_bam} \
            --normal ${normal_bam}
    }

    output {
        File vcf = "${prefix}.gridss.vcf.gz"
        File assembly_bam = "${prefix}.gridss.assembly.bam"
    }

    runtime {
        docker: "schemabio/gridss:2.13.1"
        memory: "32G"
        cpu: threads
    }
}

