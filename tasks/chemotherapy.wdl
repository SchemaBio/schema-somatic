version 1.2

task ChemotherapySites {
    input {
        String prefix
        String role
        File bam
        File bai
        Directory ref_dir
        String fasta
        String assembly
        File script
        File catalog
    }

    command <<<
        set -euo pipefail
        ln -s "~{bam}" input.bam
        ln -s "~{bai}" input.bam.bai
        python3 "~{script}" --sample "~{prefix}" --role "~{role}" --bam input.bam \
            --assembly "~{assembly}" --fasta "~{ref_dir}/~{fasta}" --catalog "~{catalog}" \
            --output "~{prefix}.chemotherapy.tsv"
    >>>

    output {
        File report = "~{prefix}.chemotherapy.tsv"
    }

    requirements {
        container: "schemabio/gatk:4.6.2.0"
    }
}
