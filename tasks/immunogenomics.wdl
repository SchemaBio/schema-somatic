version 1.2

task OptiType {
    input {
        String prefix
        File read_1
        File read_2
        Int threads
    }

    command <<<
        set -euo pipefail
        optitype run -i "~{read_1}" -i "~{read_2}" --dna --threads ~{threads} -p "~{prefix}" -o typing
        cp typing/*_result.tsv "~{prefix}.hla.tsv"
        python - "~{prefix}.hla.tsv" <<'PY' > alleles.txt
        import csv
        import sys
        with open(sys.argv[1]) as table:
            row = next(csv.DictReader(table, delimiter="\t"))
        if float(row["Reads"]) <= 0:
            sys.exit("No HLA-supporting reads")
        alleles = [row[key] for key in ("A1", "A2", "B1", "B2", "C1", "C2")]
        if not all("*" in allele for allele in alleles):
            sys.exit("Incomplete HLA-I genotype")
        print(",".join("HLA-" + allele for allele in dict.fromkeys(alleles)))
        PY
    >>>

    output {
        File genotype = "~{prefix}.hla.tsv"
        String alleles = read_string("alleles.txt")
    }

    requirements {
        container: "quay.io/biocontainers/optitype:1.5.0--pyhdfd78af_1"
        cpu: threads
        memory: "16 GiB"
    }
}

task PvacPlugins {
    command <<<
        mkdir plugins && pvacseq install_vep_plugin plugins
    >>>
    output {
        Directory plugins = "plugins"
    }
    requirements {
        container: "griffithlab/pvactools:7.1.5"
    }
}

task Neoantigens {
    input {
        String prefix
        File annotated_vcf
        String alleles
        String? normal_sample_name
        Int threads
    }

    command <<<
        set -euo pipefail
        normal_args=()
        if [[ -n "~{normal_sample_name}" ]]; then normal_args=(--normal-sample-name "~{normal_sample_name}"); fi
        pvacseq run "~{annotated_vcf}" "~{prefix}" "~{alleles}" MHCflurry pvacseq \
            --pass-only -e1 8,9,10,11 -t ~{threads} "${normal_args[@]}"
    >>>

    output {
        Directory results = "pvacseq"
        Array[File] filtered = glob("pvacseq/MHC_Class_I/*.filtered.tsv")
    }

    requirements {
        container: "griffithlab/pvactools:7.1.5"
        cpu: threads
        memory: "16 GiB"
    }
}
