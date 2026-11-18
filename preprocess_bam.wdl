version 1.2

import "tasks/fastp.wdl" as FASTP
import "tasks/bwamem.wdl" as BWAMEM
import "tasks/sambamba.wdl" as SAMBAMBA
import "tasks/fgbio.wdl" as FGBIO
import "tasks/samtools.wdl" as SAMTOOLS
import "tasks/gatk.wdl" as GATK

workflow PreBam {
    input {
        String prefix
        File read_1
        File read_2
        String? umi_struct_1
        String? umi_struct_2
        String fasta
        Directory ref_dir
        Int threads
    }

    String ref_fasta_name = basename(fasta)

    # fgbio 要求每个 read 都有结构，两个都非空才走 UMI 流程
    String umi_1 = select_first([umi_struct_1, ""])
    String umi_2 = select_first([umi_struct_2, ""])
    Boolean has_umi = umi_1 != "" && umi_2 != ""

    call FASTP.Fastp as Fastp {
        input:
            prefix = prefix,
            read_1 = read_1,
            read_2 = read_2,
            threads = threads
    }

    # WDL 没有 if 块，用空数组 scatter 实现条件执行。
    # 两个数组互补，恰好一个长度为 1，另一个为 0：未配置 UMI 时 umi_shard 为空，
    # 配置 UMI 时 plain_shard 为空，对应的 task 一个都不会启动。
    Array[String] umi_shard = if has_umi then [prefix] else []
    Array[String] plain_shard = if has_umi then [] else [prefix]

    scatter (shard in umi_shard) {
        call FGBIO.FastqToBam as FastqToBam {
            input:
                prefix = shard,
                read_1 = Fastp.clean_read_1,
                read_2 = Fastp.clean_read_2,
                read_structure_1 = umi_1,
                read_structure_2 = umi_2
        }

        call FGBIO.ExtractUmisFromBam as ExtractUmisFromBam {
            input:
                prefix = shard,
                bam = FastqToBam.unmapped_bam,
                read_structure_1 = umi_1,
                read_structure_2 = umi_2
        }

        call SAMTOOLS.SamToFastq as SamToFastq {
            input:
                prefix = shard,
                bam = ExtractUmisFromBam.umi_bam
        }

        call BWAMEM.BwaAlign as BwaAlign {
            input:
                prefix = shard,
                read_1 = SamToFastq.fastq_1,
                read_2 = SamToFastq.fastq_2,
                ref_dir = ref_dir,
                ref_fasta_name = ref_fasta_name,
                threads = threads
        }

        call GATK.MergeBAMAlignment as MergeBAMAlignment {
            input:
                prefix = shard,
                unmapped_bam = FastqToBam.unmapped_bam,
                aligned_bam = BwaAlign.out_bam,
                ref_fasta = fasta,
                ref_dir = ref_dir
        }

        call GATK.MarkDuplicates as MarkDuplicates {
            input:
                prefix = shard,
                bam = MergeBAMAlignment.merged_bam,
                bai = MergeBAMAlignment.merged_bai
        }

    }

    scatter (shard in plain_shard) {
        call BWAMEM.BwaAlign as BwaAlign {
            input:
                prefix = shard,
                read_1 = Fastp.clean_read_1,
                read_2 = Fastp.clean_read_2,
                ref_dir = ref_dir,
                ref_fasta_name = ref_fasta_name,
                threads = threads
        }

        call SAMBAMBA.SambambaMarkdup as Markdup {
            input:
                prefix = shard,
                bam = BwaAlign.out_bam,
                bai = BwaAlign.out_bai,
                threads = threads
        }
    }

    # 两条分支的输出都是长度 0 或 1 的数组，合并后恰好一个元素
    Array[File] result_bam = flatten([MarkDuplicates.markdup_bam, Markdup.markdup_bam])
    Array[File] result_bai = flatten([MarkDuplicates.markdup_bai, Markdup.markdup_bai])
    Array[File] stat_file = flatten([MarkDuplicates.metrics, Markdup.sambamba_stats])

    output {
        File fastq_report = Fastp.json_report
        File bam = result_bam[0]
        File bai = result_bai[0]
        File stats = stat_file[0]
    }
}
