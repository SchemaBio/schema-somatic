process BGZIPANDINDEX {
    tag "BGZIPANDINDEX"
    label 'somatic'
    if (publish_dir) {
        publishDir "${publish_dir}", mode: 'copy'
    }

    input:
        val(file)
        val(index_tag)

    output:
        path("*.gz"), emit: zip_file
        path("*.tbi"), emit: tbi

    script:
    """
    bgzip ${file}
    tabix -p ${index_tag} ${file}.gz
    """
}