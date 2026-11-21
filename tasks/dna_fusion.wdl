version 1.2

# One row per annotated breakpoint and transcript pair; the original VCF ID
# links each prediction to the caller's genotype and read-support fields.
task DnaFusionCandidates {
    input {
        String prefix
        String caller
        File sv_vcf
        # SnpEff database directory named for the build, containing snpEffectPredictor.bin.
        Directory snpeff_database
    }

    command <<<
        set -euo pipefail
        db=$(basename "~{snpeff_database}")
        test -s "~{snpeff_database}/snpEffectPredictor.bin"
        mkdir data
        ln -s "~{snpeff_database}" "data/$db"
        printf '%s.genome : %s\n' "$db" "$db" > snpEff.config

        snpEff -Xmx6g ann -c snpEff.config -dataDir "$PWD/data" -nodownload -noStats \
            "$db" "~{sv_vcf}" |
        awk -v caller="~{caller}" 'BEGIN {
            FS=OFS="\t"
            print "caller","chrom","pos","id","ref","alt","end","filter","genes","effect","feature_id"
        }
        !/^#/ {
            n=split($8, info, ";"); end="."; ann=""
            for (i=1; i<=n; i++) {
                if (info[i] ~ /^END=/) end=substr(info[i],5)
                if (info[i] ~ /^ANN=/) ann=substr(info[i],5)
            }
            n=split(ann, annotations, ",")
            for (i=1; i<=n; i++) {
                split(annotations[i], a, /[|]/)
                if (a[2] ~ /(^|&)(gene_fusion|bidirectional_gene_fusion)(&|$)/)
                    print caller,$1,$2,$3,$4,$5,end,$7,a[4],a[2],a[7]
            }
        }' > "~{prefix}.~{caller}.dna_fusions.tsv"
    >>>

    output {
        File candidates = "~{prefix}.~{caller}.dna_fusions.tsv"
    }

    requirements {
        container: "quay.io/biocontainers/snpeff:5.1--hdfd78af_2"
        memory: "8 GiB"
    }
}
