// Caller-independent variant tables and the compact JSON consumed by the LLM report.
process SUMMARIZE_VARIANTS {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.11 conda-forge::pandas=2.2.2"

    input:
    tuple val(meta), path(references), path(fasta), path(coverage), path(callable), path(calls)
    val   caller

    output:
    tuple val(meta), path("${prefix}_variants.tsv")       , emit: variants
    tuple val(meta), path("${prefix}_variant_summary.tsv"), emit: summary
    tuple val(meta), path("${prefix}_variants.json")      , emit: json
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    summarize_variants.py \\
        --sample ${meta.id} \\
        --caller ${caller} \\
        --platform ${params.sequencing_platform} \\
        --calls ${calls} \\
        --references ${references} \\
        --ref-fasta ${fasta} \\
        --coverage ${coverage} \\
        --callable ${callable} \\
        --min-depth ${params.variant_min_depth} \\
        --min-af ${params.variant_min_af} \\
        --min-alt-reads ${params.variant_min_alt_reads} \\
        --consensus-af ${params.variant_consensus_af} \\
        --iontorrent-max-hp ${params.variant_iontorrent_max_hp} \\
        --max-report-per-ref ${params.variant_max_report_per_ref} \\
        --prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
