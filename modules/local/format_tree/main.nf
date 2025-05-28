process FORMAT_TREE {
    tag "format_tree"
    label "process_medium"

    conda "${moduleDir}/environment.yml"

    input:
    path tree_file
    val tree_file_format

    output:
    path "formatted_tree.nwk", emit: formatted_tree
    path "versions.yml", emit: versions

    script:
    """
        format_tree.py \\
            --input ${tree_file} \\
            --format ${tree_file_format} \\
            --reformat \\
            --output formatted_tree.nwk
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1)
        END_VERSIONS
        """

    stub:
    """
        touch formatted_tree.nwk

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1)
        END_VERSIONS
        """
}
