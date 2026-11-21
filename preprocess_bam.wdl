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

    Boolean has_umi = defined(umi_struct_1) || defined(umi_struct_2)

    call FASTP.Fastp as Fastp {
        input:
            prefix = prefix,
            read_1 = read_1,
            read_2 = read_2,
            threads = threads
    }

    if (has_umi) {
        call FGBIO.FastqToBam as FastqToBam {
            input:
                prefix = prefix,
                read_1 = Fastp.clean_read_1,
                read_2 = Fastp.clean_read_2,
                read_structure_1 = select_first([umi_struct_1]),
                read_structure_2 = select_first([umi_struct_2])
        }

        call SAMTOOLS.SamToFastq as SamToFastq {
            input:
                prefix = prefix,
                bam = FastqToBam.unmapped_bam
        }

        call BWAMEM.BwaAlign as AlignUmi {
            input:
                prefix = prefix,
                read_1 = SamToFastq.fastq_1,
                read_2 = SamToFastq.fastq_2,
                ref_dir = ref_dir,
                ref_fasta_name = ref_fasta_name,
                threads = threads
        }

        call GATK.MergeBAMAlignment as MergeBAMAlignment {
            input:
                prefix = prefix,
                unmapped_bam = FastqToBam.unmapped_bam,
                aligned_bam = AlignUmi.out_bam,
                ref_fasta = ref_fasta_name,
                ref_dir = ref_dir
        }

        call GATK.MarkDuplicates as MarkDuplicates {
            input:
                prefix = prefix,
                bam = MergeBAMAlignment.merged_bam,
                bai = MergeBAMAlignment.merged_bai
        }

    }

    if (!has_umi) {
        call BWAMEM.BwaAlign as BwaAlign {
            input:
                prefix = prefix,
                read_1 = Fastp.clean_read_1,
                read_2 = Fastp.clean_read_2,
                ref_dir = ref_dir,
                ref_fasta_name = ref_fasta_name,
                threads = threads
        }

        call SAMBAMBA.SambambaMarkdup as Markdup {
            input:
                prefix = prefix,
                bam = BwaAlign.out_bam,
                bai = BwaAlign.out_bai,
                threads = threads
        }
    }

    output {
        File template_read_1 = select_first([SamToFastq.fastq_1, Fastp.clean_read_1])
        File template_read_2 = select_first([SamToFastq.fastq_2, Fastp.clean_read_2])
        File fastq_report = Fastp.json_report
        File bam = select_first([MarkDuplicates.markdup_bam, Markdup.markdup_bam])
        File bai = select_first([MarkDuplicates.markdup_bai, Markdup.markdup_bai])
        File stats = select_first([MarkDuplicates.metrics, Markdup.sambamba_stats])
    }
}
