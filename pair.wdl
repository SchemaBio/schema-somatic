version 1.2

import "preprocess_bam.wdl" as PREPROCESSBAM
import "preprocess_qc.wdl" as PREPROCESSQC
import "somatic_calls.wdl" as CALLS
import "tasks/somatic.wdl" as SOMATIC
import "tasks/cnvkit.wdl" as CNV

import "tasks/msisensor_pro.wdl" as MSI
import "tasks/immunogenomics.wdl" as IMMUNE

import "tasks/chemotherapy.wdl" as CHEMO

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
        File msi_sites
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
        File access_bed
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

    call PREPROCESSBAM.PreBam as PreprocessBamNormal {
        input:
            prefix = normal,
            read_1 = normal_1,
            read_2 = normal_2,
            umi_struct_1 = normal_umi_struct_1,
            umi_struct_2 = normal_umi_struct_2,
            ref_dir = ref_dir,
            fasta = ref_name,
            threads = threads
    }

    call PREPROCESSQC.PreQC as PreprocessQCNormal {
        input:
            prefix = normal,
            bam = PreprocessBamNormal.bam,
            bai = PreprocessBamNormal.bai,
            bed = FixBed.fixed_bed,
            json_report = PreprocessBamNormal.fastq_report,
            markdup_stats = PreprocessBamNormal.stats,
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
            normal_bam = PreprocessBamNormal.bam,
            normal_bai = PreprocessBamNormal.bai,
            normal_sample_name = normal,
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

    call CNV.AntiTarget {
        input:
            target_bed = FixBed.fixed_bed,
            access_bed = access_bed
    }

    call CNV.Coverage as TumorCoverage {
        input:
            prefix = prefix,
            bam = PreprocessBamTumor.bam,
            bai = PreprocessBamTumor.bai,
            target_bed = FixBed.fixed_bed,
            antitarget_bed = AntiTarget.antitarget_bed
    }

    call CNV.Coverage as NormalCoverage {
        input:
            prefix = normal,
            bam = PreprocessBamNormal.bam,
            bai = PreprocessBamNormal.bai,
            target_bed = FixBed.fixed_bed,
            antitarget_bed = AntiTarget.antitarget_bed
    }

    call CNV.Reference {
        input:
            prefix = normal,
            target_coverage = NormalCoverage.target_coverage,
            antitarget_coverage = NormalCoverage.antitarget_coverage,
            ref_dir = ref_dir,
            fasta = ref_name
    }

    call CNV.CnvkitCore as CopyNumber {
        input:
            prefix = prefix,
            target_coverage = TumorCoverage.target_coverage,
            antitarget_coverage = TumorCoverage.antitarget_coverage,
            reference_cnn = Reference.reference_cnn
    }

    call MSI.MsiPair as Msi {
        input:
            prefix = prefix,
            tumor_bam = PreprocessBamTumor.bam,
            tumor_bai = PreprocessBamTumor.bai,
            normal_bam = PreprocessBamNormal.bam,
            normal_bai = PreprocessBamNormal.bai,
            msi_sites = msi_sites,
            threads = threads
    }

    call IMMUNE.OptiType as Hla {
        input:
            prefix = normal,
            read_1 = PreprocessBamNormal.template_read_1,
            read_2 = PreprocessBamNormal.template_read_2,
            threads = threads
    }

    call IMMUNE.Neoantigens {
        input:
            prefix = prefix,
            annotated_vcf = SomaticCalls.annotated_vcf,
            alleles = Hla.alleles,
            normal_sample_name = normal,
            threads = threads
    }

    call CHEMO.ChemotherapySites {
        input:
            prefix = normal,
            role = "normal",
            bam = PreprocessBamNormal.bam,
            bai = PreprocessBamNormal.bai,
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
        File normal_bam = PreprocessBamNormal.bam
        File normal_bai = PreprocessBamNormal.bai
        File normal_qc = PreprocessQCNormal.qc_report
        File normal_xamdst = PreprocessQCNormal.xamdst_report
        File pair_qc = write_json(object {tumor: read_json(tumor_qc), normal: read_json(normal_qc)})
        File cnv_reference = Reference.reference_cnn
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
