// Keep only reads usable for variant calling (primary, mapped, MAPQ >= cutoff) and
// measure per-reference coverage and the number of callable bases
// (depth >= variant_min_depth), the denominator for SNV density / identity.
process VARIANT_BAM_QC {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.21"

    input:
    tuple val(meta), path(bam), path(fasta)

    output:
    tuple val(meta), path("${prefix}.filtered.bam"), path("${prefix}.filtered.bam.bai"), emit: bam
    tuple val(meta), path("${prefix}.ref_coverage.txt")                                , emit: coverage
    tuple val(meta), path("${prefix}.callable.tsv")                                    , emit: callable
    path "versions.yml"                                                                , emit: versions

    when:
    task.ext.when == nll || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -euo pipefail

    # 0x904 = unmapped | secondary | supplementary
    samtools view -b -q ${params.variant_min_mapq} -F 0x904 \\
        --threads ${task.cpus} -o ${prefix}.filtered.bam ${bam}
    samtools index ${prefix}.filtered.bam

    samtools coverage -Q ${params.variant_min_baseq} ${prefix}.filtered.bam > ${prefix}.ref_coverage.txt

    # -aa: every position of every reference, including references without reads
    samtools depth -aa -q ${params.variant_min_baseq} ${prefix}.filtered.bam \\
        | awk -v m=${params.variant_min_depth} 'BEGIN{OFS="\\t"} {n[\$1] += 0; if (\$3 >= m) n[\$1]++}
                                                END{for (r in n) print r, n[r]}' \\
        > ${prefix}.callable.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
