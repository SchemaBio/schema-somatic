version 1.2

import "tasks/samtools.wdl" as SAMTOOLS
import "tasks/xamdst.wdl" as XAMDST
import "tasks/somatic.wdl" as SOMATIC
import "tasks/gatk.wdl" as GATK

workflow PreQC {
    input {
        String prefix
        File bam
        File bai
        File bed
        File json_report
        File markdup_stats
        Int sry_sex_cutoff
        Int threads
        String assembly
        String fasta
        Directory ref_dir
    }


    call SAMTOOLS.SamtoolsSexCheck as SamtoolsSexCheck {
        input:
            prefix = prefix,
            bam = bam,
            bai = bai,
            assembly = assembly
    }
    call XAMDST.Xamdst as Xamdst {
        input:
            prefix = prefix,
            bam = bam,
            bai = bai,
            bed = bed,
            threads = threads
    }

    call GATK.CollectQCMetrics as CollectQCMetrics {
        input:
            prefix = prefix,
            bam = bam,
            bai = bai,
            bed = bed,
            fasta = fasta,
            threads = threads,
            ref_dir = ref_dir
    }

    call SOMATIC.QCReport as QCReport {
        input:
            prefix = prefix,
            fastp_stats = json_report,
            xamdst_json = Xamdst.xamdst_report,
            gatk_metric = CollectQCMetrics.summary,
            hs_metric = CollectQCMetrics.hs_metric,
            sry_result = SamtoolsSexCheck.SRY_count,
            sry_cutoff = sry_sex_cutoff,
            stats = markdup_stats
    }

    output {
        File qc_report = QCReport.qc_report
        File xamdst_report = Xamdst.xamdst_report
    }

}
