version 1.2

import "tasks/fastp.wdl" as FASTP
import "tasks/fusion.wdl" as FUSION

workflow RnaFusion {
    input {
        String prefix
        File read_1
        File read_2
        Directory genome_lib_dir
        File arriba_blacklist
        Int threads
    }

    call FASTP.Fastp {
        input:
            prefix = prefix,
            read_1 = read_1,
            read_2 = read_2,
            threads = threads
    }

    call FUSION.StarAlign {
        input:
            threads = threads,
            genome_lib_dir = genome_lib_dir,
            read_1 = Fastp.clean_read_1,
            read_2 = Fastp.clean_read_2
    }

    call FUSION.StarFusion {
        input:
            prefix = prefix,
            threads = threads,
            genome_lib_dir = genome_lib_dir,
            junctions = StarAlign.junctions
    }

    call FUSION.Arriba {
        input:
            prefix = prefix,
            genome_lib_dir = genome_lib_dir,
            bam = StarAlign.bam,
            blacklist = arriba_blacklist
    }

    output {
        File star_fusion_tsv = StarFusion.fusions
        File arriba_tsv = Arriba.fusions
        File arriba_discarded_tsv = Arriba.discarded
        File fastp_report = Fastp.json_report
    }
}
