version 1.2

import "preprocess_bam.wdl" as PREPROCESSBAM
import "preprocess_qc.wdl" as PREPROCESSQC
import "somatic_calls.wdl" as CALLS
import "tasks/somatic.wdl" as SOMATIC
import "tasks/cnvkit.wdl" as CNV

import "tasks/msisensor_pro.wdl" as MSI
import "tasks/immunogenomics.wdl" as IMMUNE

import "tasks/chemotherapy.wdl" as CHEMO

workflow SingleCancer {
    input {
        String prefix
        File read_1
        File read_2
        String? umi_struct_1
        String? umi_struct_2
        File msi_baseline
        File chemo_script = "scripts/chemotherapy_sites.py"
        File chemo_catalog = "data/chemotherapy_sites.tsv"
        File bed
        String fasta
        String assembly
        Directory ref_dir
        Directory cache_dir
        Directory schema_bundle
        Directory snpeff_database
        String clinvar_version = "20260415"
        File germline_resource_vcf
        File germline_resource_index
        File pon_vcf
        File pon_index
        Int threads
        Int sry_sex_cutoff = 40
        File reference_cnn
    }

    String ref_name = basename(fasta)

    call SOMATIC.FixBed {
        input:
            bed = bed,
            ref_dir = ref_dir,
            fasta = ref_name
    }

    call PREPROCESSBAM.PreBam as PreprocessBamTumor {
        input:
            prefix = prefix,
            read_1 = read_1,
            read_2 = read_2,
            umi_struct_1 = umi_struct_1,
            umi_struct_2 = umi_struct_2,
            ref_dir = ref_dir,
            fasta = ref_name,
            threads = threads
    }

    call PREPROCESSQC.PreQC as PreprocessQCTumor {
        input:
            prefix = prefix,
            bam = PreprocessBamTumor.bam,
            bai = PreprocessBamTumor.bai,
            bed = FixBed.fixed_bed,
            json_report = PreprocessBamTumor.fastq_report,
            markdup_stats = PreprocessBamTumor.stats,
            assembly = assembly,
            sry_sex_cutoff = sry_sex_cutoff,
            threads = threads,
            ref_dir = ref_dir,
            fasta = ref_name
    }

    call CALLS.SomaticCalls {
        input:
            prefix = prefix,
            tumor_bam = PreprocessBamTumor.bam,
            tumor_bai = PreprocessBamTumor.bai,
            bed = FixBed.fixed_bed,
            ref_dir = ref_dir,
            fasta = ref_name,
            assembly = assembly,
            cache_dir = cache_dir,
            schema_bundle = schema_bundle,
            snpeff_database = snpeff_database,
            clinvar_version = clinvar_version,
            threads = threads,
            germline_resource_vcf = germline_resource_vcf,
            germline_resource_index = germline_resource_index,
            pon_vcf = pon_vcf,
            pon_index = pon_index
    }

    call CNV.CnvkitSingle as CopyNumber {
        input:
            prefix = prefix,
            reference_cnn = reference_cnn,
            threads = threads,
            bam = PreprocessBamTumor.bam,
            bai = PreprocessBamTumor.bai
    }

    call MSI.MsiSingle as Msi {
        input:
            prefix = prefix,
            tumor_bam = PreprocessBamTumor.bam,
            tumor_bai = PreprocessBamTumor.bai,
            msi_baseline = msi_baseline,
            threads = threads
    }

    call IMMUNE.OptiType as Hla {
        input:
            prefix = prefix,
            read_1 = PreprocessBamTumor.template_read_1,
            read_2 = PreprocessBamTumor.template_read_2,
            threads = threads
    }

    call IMMUNE.Neoantigens {
        input:
            prefix = prefix,
            annotated_vcf = SomaticCalls.annotated_vcf,
            alleles = Hla.alleles,
            threads = threads
    }

    call CHEMO.ChemotherapySites {
        input:
            prefix = prefix,
            role = "tumor",
            bam = PreprocessBamTumor.bam,
            bai = PreprocessBamTumor.bai,
            assembly = assembly,
            ref_dir = ref_dir,
            fasta = ref_name,
            script = chemo_script,
            catalog = chemo_catalog
    }

    output {
        File chemotherapy_report = ChemotherapySites.report
        File msi_score = Msi.score
        File hla_genotype = Hla.genotype
        Directory neoantigen_results = Neoantigens.results
        Array[File] neoantigen_filtered = Neoantigens.filtered
        File tumor_bam = PreprocessBamTumor.bam
        File tumor_bai = PreprocessBamTumor.bai
        File tumor_qc = PreprocessQCTumor.qc_report
        File tumor_xamdst = PreprocessQCTumor.xamdst_report
        File filtered_vcf = SomaticCalls.filtered_vcf
        File filtered_index = SomaticCalls.filtered_index
        File annotated_vcf = SomaticCalls.annotated_vcf
        File manta_vcf = SomaticCalls.manta_vcf
        File manta_index = SomaticCalls.manta_index
        File gridss_candidates = SomaticCalls.gridss_candidates
        File gridss_assembly = SomaticCalls.gridss_assembly
        File manta_fusions = SomaticCalls.manta_fusions
        File gridss_fusions = SomaticCalls.gridss_fusions
        File cnv_ratios = CopyNumber.cnr
        File cnv_segments = CopyNumber.segments
        File cnv_calls = CopyNumber.call_cns
    }
}
