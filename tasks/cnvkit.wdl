version 1.2


task AntiTarget {
    input {
        File target_bed
        File access_bed
    }

    command <<<
        cnvkit.py antitarget ~{target_bed} -g ~{access_bed} -o antitarget.bed
    >>>

    output {
        File antitarget_bed = "antitarget.bed"
    }

    runtime {
        docker: "schemabio/cnvkit:0.9.13.4"
    }
}

task Coverage {
    input {
        String prefix
        File bam
        File bai
        File target_bed
        File antitarget_bed
    }

    command <<<
        set -euo pipefail
        ln -s "~{bam}" input.bam
        ln -s "~{bai}" input.bam.bai
        cnvkit.py coverage input.bam ~{target_bed} -o ~{prefix}.targetcoverage.cnn
        cnvkit.py coverage input.bam ~{antitarget_bed} -o ~{prefix}.antitargetcoverage.cnn
    >>>

    output {
        File target_coverage = "~{prefix}.targetcoverage.cnn"
        File antitarget_coverage = "~{prefix}.antitargetcoverage.cnn"
    }

    runtime {
        docker: "schemabio/cnvkit:0.9.13.4"
    }
}

task Reference {
    input {
        String prefix
        File target_coverage
        File antitarget_coverage
        String fasta
        Directory ref_dir
    }

    command <<<
        cnvkit.py reference ~{target_coverage} ~{antitarget_coverage} -f "~{ref_dir}/~{fasta}" -o ~{prefix}.reference.cnn
    >>>

    output {
        File reference_cnn = "~{prefix}.reference.cnn"
    }

    runtime {
        docker: "schemabio/cnvkit:0.9.13.4"
    }
}

task CnvkitCore {
    input {
        File target_coverage
        File antitarget_coverage
        File reference_cnn
        String prefix
        String seg_method = "cbs"
    }

    command <<<
        set -euo pipefail
        cnvkit.py fix ~{target_coverage} ~{antitarget_coverage} ~{reference_cnn} -o ~{prefix}.cnr
        cnvkit.py segment ~{prefix}.cnr -o ~{prefix}.cns -m ~{seg_method}
        cnvkit.py call ~{prefix}.cns -o ~{prefix}.call.cns

    >>>

    output {
        File call_cns = "~{prefix}.call.cns"
        File cnr = "~{prefix}.cnr"
        File segments = "~{prefix}.cns"
    }

    runtime {
        docker: "schemabio/cnvkit:0.9.13.4"
    }
}


task CnvkitSingle {
    input {
        String prefix
        File bam
        File bai
        File reference_cnn
        Int threads
    }

    command <<<
        set -euo pipefail
        ln -s "~{bam}" "~{prefix}.bam"
        ln -s "~{bai}" "~{prefix}.bam.bai"
        cnvkit.py batch "~{prefix}.bam" -r "~{reference_cnn}" -p ~{threads} -d result
        cnvkit.py call "result/~{prefix}.cns" -o "result/~{prefix}.call.cns"
    >>>

    output {
        File cnr = "result/~{prefix}.cnr"
        File segments = "result/~{prefix}.cns"
        File call_cns = "result/~{prefix}.call.cns"
    }

    requirements {
        container: "schemabio/cnvkit:0.9.13.4"
        cpu: threads
        memory: "8 GiB"
    }
}
