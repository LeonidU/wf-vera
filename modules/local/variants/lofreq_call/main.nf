// LoFreq: quality-aware caller, sensitive down to ~1% allele frequency.
//   viterbi    - realign around indels (fewer false SNVs next to indels)
//   indelqual  - indel qualities: Dindel model for Illumina, uniform for Ion Torrent
//   call       - SNVs (+ indels) with LoFreq's default dynamic strand-bias filter
process LOFREQ_CALL {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::lofreq=2.1.5 bioconda::samtools=1.21"

    input:
    tuple val(meta), path(bam), path(bai), path(fasta)

    output:
    tuple val(meta), path("${prefix}.lofreq.vcf"), emit: vcf
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def indelqual = params.sequencing_platform == 'iontorrent' ?
        "--uniform ${params.variant_indel_uniform_q}" : "--dindel -f ${fasta}"
    def indels = params.variant_call_indels ? "--call-indels" : ""
    def caller = task.cpus > 1 ? "call-parallel --pp-threads ${task.cpus}" : "call"
    """
    set -euo pipefail
    samtools faidx ${fasta}

    if [[ \$(samtools view -c ${bam}) -eq 0 ]]; then
        printf '##fileformat=VCFv4.0\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n' > ${prefix}.lofreq.vcf
    else
        lofreq viterbi -f ${fasta} ${bam} | samtools sort -@ ${task.cpus} -o realn.bam -
        lofreq indelqual ${indelqual} -o iq.bam realn.bam
        samtools index iq.bam

        lofreq ${caller} \\
            -f ${fasta} \\
            ${indels} \\
            --min-bq ${params.variant_min_baseq} \\
            --min-alt-bq ${params.variant_min_baseq} \\
            --min-mq ${params.variant_min_mapq} \\
            --min-cov ${params.variant_min_depth} \\
            -o ${prefix}.lofreq.vcf \\
            iq.bam
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq: \$(lofreq version | grep version | sed 's/version: //')
    END_VERSIONS
    """
}
