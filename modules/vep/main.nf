// VEP 变异注释模块
// 用途：使用Ensembl VEP对变异进行功能注释
// 用法：
//   - 输入：样本ID、VCF文件、参考基因组、VEP数据库路径、基因组版本
//   - 输出：注释后的VCF文件
//   - 支持GRCh37和GRCh38
//   - 包含ClinVar、InterVar、dbNSFP、SpliceAI等数据库

process VEP {
    tag "VEP on $sample_id"
    label 'vep'
    publishDir "${params.output}/04.Annotation", mode: 'copy'

    input:
        tuple val(sample_id), path(vcf), path(vcf_index)
        path fasta
        path fasta_index
        val vep_database
        val genome_build  // 'GRCh37' or 'GRCh38'

    output:
        tuple val(sample_id), path("${sample_id}.vep.vcf"), emit: annotated_vcf
        path("${sample_id}.vep.html"), emit: html_report

    script:
    def threads = task.cpus
    def assembly = genome_build == 'GRCh38' ? 'GRCh38' : 'GRCh37'
    def cache_str = "Uploaded_variation,Location,REF_ALLELE,Allele,Consequence,IMPACT,DOMAINS,Feature,DISTANCE,EXON,INTRON,SYMBOL,STRAND,HGNC_ID,HGVSc,HGVSp,HGVSg,MAX_AF,Protein_position,Amino_acids,Codons,PUBMED,Existing_variation"
    def custom_str = "cytoBand,ClinVar_CLNSIG,ClinVar_CLNREVSTAT,ClinVar_CLNDN,ClinVar_CLNHGVS,CIViC"
    """
    vep --offline --cache --dir_cache ${vep_database} --refseq \\
        --dir_plugins ${vep_database}/Plugins \\
        --force_overwrite --fork ${threads} \\
        -i ${vcf} -o ${sample_id}.vep.vcf \\
        --format vcf --vcf --fa ${fasta} \\
        --shift_3prime 1 --assembly ${assembly} --no_escape --check_existing --exclude_predicted \\
        --uploaded_allele --show_ref_allele --numbers --domains \\
        --pick --pick_order mane_plus_clinical,mane_select,refseq,canonical,ensembl --pick_allele \\
        --transcript_filter "stable_id match N[MR]_" \\
        --total_length --hgvs --hgvsg --symbol --ccds --uniprot --max_af --pubmed \\
        --custom file=${vep_database}/clinvar/clinvar.vcf.gz,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN%CLNHGVS \\
        --custom file=${vep_database}/cytoband/cytoBand.bed.gz,short_name=cytoBand,format=bed,type=overlap,coords=0 \\
        --plugin Pangolin,${vep_database}/pangolin/pangolin.vcf.gz \\
        --fields ${cache_str},${custom_str},${dbnsfp_str},${spliceai_str} \\
        --stats_file ${sample_id}.vep.html
    """
}
