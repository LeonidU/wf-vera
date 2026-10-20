// Build the per-sample reference set for read-based variant calling.
//
//   mode 'mmseqs'    : references = best MMseqs2 nucleotide hits of the sample's viral
//                      contigs, sequences pulled from the same MMseqs2 nucleotide DB
//                      (no separate genome collection needed)
//   mode 'fasta'     : same selection, sequences pulled from a multi-FASTA
//                      (e.g. locally downloaded NCBI genomes)
//   mode 'fasta_all' : every sequence of the multi-FASTA (a fixed panel, e.g. bacterial
//                      genomes); tophit/taxa inputs are ignored
//
// Near-identical candidates are then collapsed (MMseqs2 easy-linclust: made for high
// identity, and unlike easy-cluster it needs no 4^15 nucleotide k-mer table in RAM) so
// that reads are not split between two strains of the same virus during competitive mapping.
// No output is emitted for a sample without usable references.
process PREPARE_VARIANT_REFERENCES {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::mmseqs2=18.8cc5c bioconda::seqkit=2.8.2 conda-forge::pandas=2.2.2"

    input:
    tuple val(meta), path(tophit), path(taxa)
    path  ref_source
    val   mode

    output:
    tuple val(meta), path("${prefix}_references.fa") , emit: fasta, optional: true
    tuple val(meta), path("${prefix}_references.tsv"), emit: table, optional: true
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def taxa_arg = taxa ? "--taxa ${taxa}" : ""
    """
    set -euo pipefail

    # ── 1. candidate references ──────────────────────────────────────────────
    if [[ "${mode}" == "fasta_all" ]]; then
        printf "reference\\tref_length\\n" > cand_ref_candidates.tsv
        seqkit fx2tab -n -i -l ${ref_source} | awk -F'\\t' 'BEGIN{OFS="\\t"} {print \$1, \$2}' >> cand_ref_candidates.tsv
        tail -n +2 cand_ref_candidates.tsv | cut -f1 > cand_ref_candidates.txt
    else
        select_variant_references.py \\
            --tophit ${tophit} \\
            ${taxa_arg} \\
            --min-aln-len ${params.variant_ref_min_aln_len} \\
            --min-fident ${params.variant_ref_min_fident} \\
            --max-refs ${params.variant_ref_max} \\
            --prefix cand
    fi

    if [[ ! -s cand_ref_candidates.txt ]]; then
        echo "No candidate references for ${meta.id}; skipping variant calling" >&2
    else
        # ── 2. extract sequences ─────────────────────────────────────────────
        if [[ "${mode}" == "mmseqs" ]]; then
            # same DB prefix detection as MMSEQS_EASYTAXONOMY
            DB=\$(find -L "${ref_source}/" -maxdepth 1 -name "*.dbtype" | sed 's/\\.[^.]*\$//' | sed -e 'N;s/^\\(.*\\).*\\n\\1.*\$/\\1\\n\\1/;D')
            for f in "\${DB}.lookup" "\${DB}_h.dbtype"; do
                [[ -e "\$f" ]] || { echo "ERROR: \$f not found - cannot extract references from the MMseqs2 DB" >&2; exit 1; }
            done
            # accession -> internal DB key
            awk 'NR==FNR {want[\$1]=1; next} (\$2 in want) {print \$1}' cand_ref_candidates.txt "\${DB}.lookup" > keys.txt
            mkdir -p sub
            mmseqs createsubdb keys.txt "\${DB}"   sub/refs   -v 1
            mmseqs createsubdb keys.txt "\${DB}_h" sub/refs_h -v 1
            mmseqs convert2fasta sub/refs cand.fa -v 1
        else
            seqkit grep -f cand_ref_candidates.txt ${ref_source} > cand.fa
        fi
        seqkit seq -n -i cand.fa > present.txt

        missing=\$(grep -vxFf present.txt cand_ref_candidates.txt | wc -l || true)
        [[ "\$missing" -gt 0 ]] && echo "WARNING: \$missing candidate reference(s) not found in ${ref_source}" >&2

        # ── 3. collapse near-identical references ────────────────────────────
        n=\$(grep -c '^>' cand.fa || true)
        if [[ "\$n" -gt 1 ]] && awk 'BEGIN{exit !(${params.variant_ref_derep_id} > 0)}'; then
            mmseqs easy-linclust cand.fa clu tmp_clu \\
                --min-seq-id ${params.variant_ref_derep_id} \\
                -c 0.8 --cov-mode 1 \\
                --max-seq-len ${params.mmseqs_max_seq_len} \\
                --threads ${task.cpus} -v 1
            dereplicate_references.py --candidates cand_ref_candidates.tsv --clusters clu_cluster.tsv \\
                --present present.txt --prefix ${prefix}
        elif [[ "\$n" -gt 0 ]]; then
            dereplicate_references.py --candidates cand_ref_candidates.tsv \\
                --present present.txt --prefix ${prefix}
        fi

        if [[ -s ${prefix}_references.txt ]]; then
            # upper-case bases: soft-masked references confuse some callers
            seqkit grep -f ${prefix}_references.txt cand.fa | seqkit seq -u > ${prefix}_references.fa
        else
            rm -f ${prefix}_references.tsv
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs version)
        seqkit: \$(seqkit version | sed 's/seqkit v//')
    END_VERSIONS
    """
}
