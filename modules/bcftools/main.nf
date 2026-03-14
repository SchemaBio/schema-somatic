process BCFTOOLS {
    tag "BCFTOOLS on $sample_id"
    label 'bcftools'
    if (publish_dir) {
        publishDir "${publish_dir}", mode: 'copy'
    }

    input:
        tuple val(sample_id), path(vcf)
        val(publish_dir)
        val(minAD)
        val(minAF)

    output:
        path("${sample}.filtered.vcf"), emit: filtered_vcf

    script:
    """
    bcftools view -e "FORMAT/AD[0:1]<${minAD} || FORMAT/AF<${minAF}" ${vcf} > ${sample}.filtered.vcf
    """
}