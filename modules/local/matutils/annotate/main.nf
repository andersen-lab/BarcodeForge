process MATUTILS_ANNOTATE {
    tag "matutils_annotate"
    label "process_low"

    conda "${moduleDir}/environment.yml"

    input:
    path protobuf_tree_file
    path clades
    val matUtils_overlap

    output:
    path 'annotated_tree.pb', emit: annotated_tree_file
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    matUtils \\
        annotate \\
        ${args} \\
        --set-overlap ${matUtils_overlap} \\
        -i ${protobuf_tree_file} \\
        -c ${clades} \\
        -o annotated_tree.pb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        usher: \$(usher --version 2>&1)
    END_VERSIONS
    """

    stub:
    """
    touch tree.pb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        usher: \$(usher --version 2>&1)
    END_VERSIONS
    """
}
