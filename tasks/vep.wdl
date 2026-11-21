version 1.2

task VEP {
    input {
        String prefix
        Directory pvac_plugins
        File vcf
        Directory cache_dir
        Directory schema_bundle
        Directory ref_dir
        String fasta
        String assembly
        String clinvar_version = "20260415"
        Int threads
    }

    command <<<
        set -euo pipefail
        case "~{assembly}" in
            GRCh37) build=hg19 ;;
            GRCh38) build=hg38 ;;
            *) exit 1 ;;
        esac
        vep --offline --cache --merged --dir_cache "~{cache_dir}" --assembly "~{assembly}" \
            --fasta "~{ref_dir}/~{fasta}" --fork ~{threads} --symbol --hgvs --terms SO --tsl --mane_select --canonical --biotype \
            --dir_plugins "~{pvac_plugins}" --plugin Frameshift --plugin Wildtype --safe \
            --custom "file=~{schema_bundle}/${build}_clinvar_~{clinvar_version}.vcf.gz,short_name=ClinVar,format=vcf,type=exact,fields=CLNSIG%CLNDN%CLNREVSTAT" \
            --format vcf --vcf -i "~{vcf}" -o "~{prefix}.vep.vcf"
    >>>

    output {
        File annotated_vcf = "~{prefix}.vep.vcf"
    }

    requirements {
        container: "schemabio/vep:116.0"
        cpu: threads
        memory: "8 GiB"
    }
}
