// iVar: frequency/depth-threshold caller on samtools mpileup, the de facto standard for
// amplicon (ARTIC-like) viral data. Base quality is filtered by iVar itself (-q),
// so mpileup runs with -Q 0 and without BAQ.
process IVAR_VARIANTS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::ivar=1.4.4 bioconda::samtools=1.21"

    input:
    tuple val(meta), path(bam), path(bai), path(fasta)

    output:
    tuple val(meta), path("${prefix}.ivar.tsv"), emit: tsv
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -euo pipefail
    samtools faidx ${fasta}

    samtools mpileup -aa -A -d 0 -B -Q 0 -q ${params.variant_min_mapq} \\
        --reference ${fasta} ${bam} \\
        | ivar variants \\
            -p ${prefix}.ivar \\
            -q ${params.variant_min_baseq} \\
            -t ${params.variant_min_af} \\
            -m ${params.variant_min_depth} \\
            -r ${fasta}

    # header-only table when nothing piled up
    if [[ ! -s ${prefix}.ivar.tsv ]]; then
        printf 'REGION\\tPOS\\tREF\\tALT\\tREF_DP\\tREF_RV\\tREF_QUAL\\tALT_DP\\tALT_RV\\tALT_QUAL\\tALT_FREQ\\tTOTAL_DP\\tPVAL\\tPASS\\n' > ${prefix}.ivar.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ivar: \$(ivar version | head -1 | sed 's/iVar version //')
    END_VERSIONS
    """
}
