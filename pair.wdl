version 1.2

import "tasks/somatic.wdl" as SOMATIC
import "tasks/samtools.wdl" as SAMTOOLS
import "preprocess_bam.wdl" as PREPROCESSBAM

workflow PairCancer {
    input {
        String prefix
        File read_1
        File read_2
        String? umi_struct_1
        String? umi_struct_2
        String normal
        File normal_1
        File normal_2
        String? normal_umi_struct_1
        String? normal_umi_struct_2
        File bed
        String fasta
        String assembly
        Directory ref_dir
        Directory cache_dir
        Directory schema_bundle
        Int threads
    }

    # 参数调整
    String ref_fasta_name = basename(fasta)
    String clinvar_version = '20260415'

    call PREPROCESSBAM.PreBam as PreprocessTumor {
        input:
            prefix = prefix,
            read_1 = read_1,
            read_2 = read_2,
            umi_struct_1 = umi_struct_1,
            umi_struct_2 = umi_struct_2,
            fasta = fasta,
            ref_dir = ref_dir,
            threads = threads
    }

    call PREPROCESSBAM.PreBam as PreprocessNormal {
        input:
            prefix = normal,
            read_1 = normal_1,
            read_2 = normal_2,
            umi_struct_1 = normal_umi_struct_1,
            umi_struct_2 = normal_umi_struct_2,
            fasta = fasta,
            ref_dir = ref_dir,
            threads = threads
    }



}




