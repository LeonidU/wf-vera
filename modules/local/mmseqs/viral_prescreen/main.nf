// Fast pre-screen of translated ORFs against a small, viral-only protein DB.
//
// Purpose: MMSEQS_PROTEINS (easy-taxonomy against the full UniRef90) spends its time
// in the prefilter, and that time grows with the number of query ORFs. Only contigs
// carrying at least one ORF with a viral hit can ever be selected downstream, so we
// find those contigs here (seconds to minutes on a viral-only DB) and send ONLY their
// ORFs to the expensive full-DB search.
//
// All ORFs of a candidate contig are kept (not just the ORFs that hit), because
// PARSE_MMSEQS_PROTEIN_RESULTS looks at every ORF of a contig to decide between
// "clean" / "unstable" / "mixed" lineage.
process MMSEQS_VIRAL_PRESCREEN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/../easytaxonomy/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2:15.6f452--pl5321h6a68c12_3':
        'biocontainers/mmseqs2:15.6f452--pl5321h6a68c12_3' }"

    input:
    tuple val(meta) , path(orfs)
    tuple val(meta2), path(db_target)

    output:
    tuple val(meta), path("${prefix}_candidate_orfs.fa.gz") , emit: orfs
    tuple val(meta), path("${prefix}_candidate_contigs.txt"), emit: contigs
    tuple val(meta), path("${prefix}_viral_hits.m8")        , emit: hits
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: "*.dbtype"
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Same DB prefix detection as in MMSEQS_EASYTAXONOMY
    DB_TARGET_PATH_NAME=\$(find -L "$db_target/" -maxdepth 1 -name "$args2" | sed 's/\\.[^.]*\$//' | sed -e 'N;s/^\\(.*\\).*\\n\\1.*\$/\\1\\n\\1/;D' )

    mmseqs \\
        easy-search \\
        $orfs \\
        \$DB_TARGET_PATH_NAME \\
        ${prefix}_viral_hits.m8 \\
        tmp_prescreen \\
        $args \\
        --threads ${task.cpus} \\
        --db-load-mode 2

    # ORF ids look like <contig>_frame=<f>_begin=<b>_end=<e>  ->  contig id
    cut -f1 ${prefix}_viral_hits.m8 | sed 's/_frame=.*//' | sort -u > ${prefix}_candidate_contigs.txt

    # Keep every ORF that belongs to a candidate contig
    gzip -cdf $orfs \\
        | awk 'NR==FNR { keep[\$1] = 1; next }
               /^>/    { id = substr(\$1, 2); sub(/_frame=.*/, "", id); p = (id in keep) }
               p' ${prefix}_candidate_contigs.txt - \\
        | gzip -c > ${prefix}_candidate_orfs.fa.gz

    echo "pre-screen: \$(wc -l < ${prefix}_candidate_contigs.txt) candidate contigs," \\
         "\$(gzip -cd ${prefix}_candidate_orfs.fa.gz | grep -c '>' || true) of" \\
         "\$(gzip -cdf $orfs | grep -c '>' || true) ORFs go to the full protein search"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_candidate_orfs.fa.gz
    touch ${prefix}_candidate_contigs.txt
    touch ${prefix}_viral_hits.m8

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """
}
