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
        cnvkit.py coverage ~{bam} -b ~{target_bed} -o ~{prefix}.targetcoverage.cnn
        cnvkit.py coverage ~{bam} -b ~{antitarget_bed} -o ~{prefix}.antitargetcoverage.cnn
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
        Array[File] target_coverage
        Array[File] antitarget_coverage
        File fasta
    }

    command <<<
        cnvkit.py reference ~{target_coverage} ~{antitarget_coverage} -f ~{fasta} -o ~{prefix}.reference.cnn
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
        cnvkit.py fix ~{target_coverage} ~{antitarget_coverage} -r ~{reference_cnn} -o ~{prefix}.cnr
        cnvkit.py segment ~{prefix}.cnr -o ~{prefix}.cns -m ~{seg_method}
        cnvkit.py call ~{prefix}.cns -o ~{prefix}.call.cns
        cnvkit.py call ~{prefix}.cnr -o ~{prefix}.call.cnr
    >>>

    output {
        File call_cns = "~{prefix}.call.cns"
        File call_cnr = "~{prefix}.call.cnr"
    }

    runtime {
        docker: "schemabio/cnvkit:0.9.13.4"
    }
}


