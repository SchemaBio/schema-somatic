version 1.2

task SamtoolsSexCheck {
    input {
        String prefix
        File bam
        File bai
        String assembly
    }

    command <<<
        if [ "~{assembly}" == "GRCh38" ]; then
            samtools view ~{bam} Y:2786855-2787682 | wc -l > ~{prefix}.SRY.count.txt
        elif [ "~{assembly}" == "GRCh37" ]; then
            samtools view ~{bam} Y:2654896-2655723 | wc -l > ~{prefix}.SRY.count.txt
        else
            echo "Unsupported assembly: ~{assembly}" >&2
            exit 1
        fi
    >>>

    output {
        File SRY_count = "~{prefix}.SRY.count.txt"
    }

    runtime {
        docker: "schemabio/mapping:v1.1.0"
    }
}

task SamToFastq {
    input {
        String prefix
        File bam
    }

    command <<<
        samtools fastq -1 ~{prefix}_R1.fastq.gz -2 ~{prefix}_R2.fastq.gz -0 /dev/null -s /dev/null -n ~{bam}
    >>>

    output {
        File fastq_1 = "~{prefix}_R1.fastq.gz"
        File fastq_2 = "~{prefix}_R2.fastq.gz"
    }

    runtime {
        docker: "schemabio/mapping:v1.1.0"
    }
}

