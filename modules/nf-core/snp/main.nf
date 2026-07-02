process SNP_CALLING {
    input:
    tuple val(meta), path(input), path(input_index)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    
    output:
    tuple val(meta), path("*.vcf"), emit: coverage

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    lofreq call -f ${fasta} -o ${prefix}.vcf ${input}
    """
}
