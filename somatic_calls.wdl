version 1.2

import "tasks/gatk.wdl" as GATK
import "tasks/vep.wdl" as ANNOTATION
import "tasks/manta.wdl" as MANTA
import "tasks/gridss.wdl" as GRIDSS

import "tasks/dna_fusion.wdl" as DNAFUSION

import "tasks/immunogenomics.wdl" as IMMUNE

workflow SomaticCalls {
    input {
        String prefix
        File tumor_bam
        File tumor_bai
        File? normal_bam
        File? normal_bai
        String? normal_sample_name
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
    }

    call GATK.Mutect2 {
        input:
            prefix = prefix,
            ref_dir = ref_dir,
            fasta = fasta,
            tumor_bam = tumor_bam,
            tumor_bai = tumor_bai,
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            normal_sample_name = normal_sample_name,
            bed = bed,
            germline_resource_vcf = germline_resource_vcf,
            germline_resource_index = germline_resource_index,
            pon_vcf = pon_vcf,
            pon_index = pon_index
    }

    call GATK.FilterMutectCalls {
        input:
            prefix = prefix,
            ref_dir = ref_dir,
            fasta = fasta,
            vcf = Mutect2.vcf,
            index = Mutect2.index,
            stats = Mutect2.stats
    }

    call IMMUNE.PvacPlugins {}

    call ANNOTATION.VEP {
        input:
            prefix = prefix,
            ref_dir = ref_dir,
            fasta = fasta,
            assembly = assembly,
            cache_dir = cache_dir,
            schema_bundle = schema_bundle,
            clinvar_version = clinvar_version,
            threads = threads,
            pvac_plugins = PvacPlugins.plugins,
            vcf = FilterMutectCalls.filtered_vcf
    }

    call MANTA.Manta {
        input:
            prefix = prefix,
            ref_dir = ref_dir,
            fasta = fasta,
            tumor_bam = tumor_bam,
            tumor_bai = tumor_bai,
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            threads = threads
    }

    call GRIDSS.Gridss {
        input:
            prefix = prefix,
            ref_dir = ref_dir,
            fasta = fasta,
            tumor_bam = tumor_bam,
            tumor_bai = tumor_bai,
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            threads = threads
    }

    call DNAFUSION.DnaFusionCandidates as MantaFusions {
        input:
            prefix = prefix,
            caller = "Manta",
            sv_vcf = Manta.vcf,
            snpeff_database = snpeff_database
    }

    call DNAFUSION.DnaFusionCandidates as GridssFusions {
        input:
            prefix = prefix,
            caller = "Gridss",
            sv_vcf = Gridss.vcf,
            snpeff_database = snpeff_database
    }

    output {
        File filtered_vcf = FilterMutectCalls.filtered_vcf
        File filtered_index = FilterMutectCalls.filtered_index
        File annotated_vcf = VEP.annotated_vcf
        File manta_vcf = Manta.vcf
        File manta_index = Manta.index
        File gridss_candidates = Gridss.vcf
        File gridss_assembly = Gridss.assembly_bam
        File manta_fusions = MantaFusions.candidates
        File gridss_fusions = GridssFusions.candidates
    }
}
