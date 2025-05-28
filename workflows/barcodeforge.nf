/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { FATOVCF                } from '../modules/local/faToVcf/main.nf'
include { USHER                  } from '../modules/local/usher/main.nf'
include { MATUTILS_ANNOTATE      } from '../modules/local/matutils/annotate/main.nf'
include { MATUTILS_EXTRACT       } from '../modules/local/matutils/extract/main.nf'
include { ADD_REF_MUTS           } from '../modules/local/add_ref_muts/main.nf'
include { GENERATE_BARCODES      } from '../modules/local/generate_barcodes/main.nf'
include { FORMAT_TREE            } from '../modules/local/format_tree/main.nf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BARCODEFORGE {
    main:

    ch_versions = Channel.empty()

    FATOVCF(Channel.fromPath(params.alignment, checkIfExists: true))
        .set { ch_faToVcf_versions }

    FORMAT_TREE(
        Channel.fromPath(params.tree_file, checkIfExists: true),
        params.tree_file_format,
    )

    USHER(FORMAT_TREE.out.formatted_tree, FATOVCF.out.vcf)

    MATUTILS_ANNOTATE(
        USHER.out.protobuf_tree,
        Channel.fromPath(params.lineages, checkIfExists: true),
        params.matUtils_overlap,
    )

    MATUTILS_EXTRACT(
        MATUTILS_ANNOTATE.out.annotated_tree_file
    )

    ADD_REF_MUTS(
        Channel.fromPath(params.reference_genome, checkIfExists: true),
        MATUTILS_EXTRACT.out.sample_paths_file,
        MATUTILS_EXTRACT.out.lineage_definition_file,
        Channel.fromPath(params.alignment, checkIfExists: true),
    )

    GENERATE_BARCODES(ADD_REF_MUTS.out.modified_lineage_paths, params.barcode_prefix)

    ch_versions = ch_versions.mix(
        FORMAT_TREE.out.versions,
        USHER.out.versions,
        MATUTILS_ANNOTATE.out.versions,
        MATUTILS_EXTRACT.out.versions,
    )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'barcodeforge_software_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }

    emit:
    versions = ch_versions // channel: [ path(versions.yml) ]
}
