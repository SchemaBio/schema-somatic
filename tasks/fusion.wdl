version 1.2

task StarAlign {
    input {
        File read_1
        File read_2
        Directory genome_lib_dir
        Int threads
    }

    command <<<
        STAR --runThreadN ~{threads} --genomeDir "~{genome_lib_dir}/ref_genome.fa.star.idx" \
            --readFilesIn "~{read_1}" "~{read_2}" --readFilesCommand zcat \
            --outSAMtype BAM Unsorted --outSAMunmapped Within \
            --chimOutType Junctions WithinBAM HardClip --chimOutJunctionFormat 1 \
            --chimSegmentMin 10 --chimJunctionOverhangMin 10 --chimMultimapNmax 50 \
            --outFilterMultimapNmax 50 --alignSJstitchMismatchNmax 5 -1 5 5 \
            --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 \
            --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3
    >>>

    output {
        File bam = "Aligned.out.bam"
        File junctions = "Chimeric.out.junction"
    }

    requirements {
        container: "trinityctat/starfusion:1.14.0"
        cpu: threads
        memory: "48 GiB"
    }
}

task StarFusion {
    input {
        String prefix
        File junctions
        Directory genome_lib_dir
        Int threads
    }

    command <<<
        set -euo pipefail
        STAR-Fusion --genome_lib_dir "~{genome_lib_dir}" -J "~{junctions}" \
            --CPU ~{threads} --output_dir fusion
        cp fusion/star-fusion.fusion_predictions.tsv "~{prefix}.star-fusion.tsv"
    >>>

    output {
        File fusions = "~{prefix}.star-fusion.tsv"
    }

    requirements {
        container: "trinityctat/starfusion:1.14.0"
        cpu: threads
        memory: "8 GiB"
    }
}

task Arriba {
    input {
        String prefix
        File bam
        Directory genome_lib_dir
        File blacklist
    }

    command <<<
        /arriba_v2.5.1/arriba -x "~{bam}" -a "~{genome_lib_dir}/ref_genome.fa" \
            -g "~{genome_lib_dir}/ref_annot.gtf" -b "~{blacklist}" \
            -o "~{prefix}.arriba.tsv" -O "~{prefix}.arriba.discarded.tsv"
    >>>

    output {
        File fusions = "~{prefix}.arriba.tsv"
        File discarded = "~{prefix}.arriba.discarded.tsv"
    }

    requirements {
        container: "uhrigs/arriba:2.5.1"
        memory: "16 GiB"
    }
}
