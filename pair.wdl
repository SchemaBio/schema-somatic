version 1.2

import "tasks/somatic.wdl" as GERMLINE
import "tasks/fastp.wdl" as FASTP

workflow PairCancer {
    input {
        String prefix
        File read_1
        File read_2
        String umi_struct
        String normal
        File normal_1
        File normal_2
        String normal_umi_struct
        File bed
        String fasta
        String assembly
        Directory ref_dir
        Directory cache_dir
        Directory schema_bundle
    }

    # 参数调整
    String ref_fasta_name = basename(fasta)
    String clinvar_version = '20260415'

    # bed 文件调整
    call SOMATIC.FixBed as FixBed {
        input:
            bed = bed
    }

    # fastp
    call FASTP.Fastp as Fastp {
        input:
            prefix = prefix,
            read_1 = read_1,
            read_2 = read_2,
            threads = 16
    }

    call FASTP.Fastp as NormalFastp {
        input:
            prefix = normal,
            read_1 = normal_1,
            read_2 = normal_2,
            threads = 16
    }

    # 比对
    call BWAMEM.BwaAlign as BwaAlign {
        input:
            prefix = prefix,
            read_1 = Fastp.clean_read_1,
            read_2 = Fastp.clean_read_2,
            ref_dir = ref_dir,
            ref_fasta_name = ref_fasta_name,
            threads = 16
    }
    call SAMBAMBA.SambambaMarkdup as Markdup {
        input:
            prefix = prefix,
            bam = BwaAlign.out_bam,
            bai = BwaAlign.out_bai,
            threads = 16
    }
    call BWAMEM.BwaAlign as NormalBwaAlign {
        input:
            prefix = normal,
            read_1 = NormalFastp.clean_read_1,
            read_2 = NormalFastp.clean_read_2,
            ref_dir = ref_dir,
            ref_fasta_name = ref_fasta_name,
            threads = 16
    }
    call SAMBAMBA.SambambaMarkdup as NormalMarkdup {
        input:
            prefix = normal,
            bam = NormalBwaAlign.out_bam,
            bai = NormalBwaAlign.out_bai,
            threads = 16
    }

    # 质控线路
    # 包含性别检测，但须确认bed是否有SRY
    # 包含配对样本的来源一致性检测


    # 化疗


    # 热点检测


    # mutect2


    # 注释


    # cnv


    # fusion



}




