ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
// 
// MODULE: Local to the pipeline
// 
include { INPUT_CHECK                            } from '../subworkflows/local/input_check'

//
// MODULE: Installed directly from nf-core/modules
//

// input quality control
include { FASTQC as FASTQC_RAW                                  } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIMMED                              } from '../modules/nf-core/fastqc/main'
// input taxonomy control
include { KRAKEN2_KRAKEN2 as KRAKEN2_RAW                        } from '../modules/nf-core/kraken2/kraken2/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_CLEAN                      } from '../modules/nf-core/kraken2/kraken2/main'
// cleaning from host and contaminants
include { BOWTIE2_BUILD as BOWTIE2_HOST_BUILD  } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN as BOWTIE2_HOST_ALIGN  } from '../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_BUILD as BOWTIE2_CONTAM_BUILD  } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN as BOWTIE2_CONTAM_ALIGN  } from '../modules/nf-core/bowtie2/align/main'

include { FASTP                                                 } from '../modules/nf-core/fastp/main'
include { CAT_FASTQ                                             } from '../modules/nf-core/cat/fastq/main'

include { SPADES } from '../modules/nf-core/spades/main'
include { GUNZIP as GUNZIP_CONTIGS } from '../modules/nf-core/gunzip/main'

include { SEQKIT_SEQ } from '../modules/nf-core/seqkit/seq/main'
include { SEQKIT_TRANSLATE } from '../modules/nf-core/seqkit/translate/main'
include { SEQKIT_GREP as SEQKIT_GREP_FOR_NUCL } from '../modules/nf-core/seqkit/grep/main'
include { SEQKIT_GREP as SEQKIT_GREP_FOR_EXTENSION } from '../modules/nf-core/seqkit/grep/main'

include { MMSEQS_EASYTAXONOMY as MMSEQS_PROTEINS } from '../modules/local/mmseqs/easytaxonomy/main'
include { MMSEQS_EASYTAXONOMY as MMSEQS_NUCLEOTIDES } from '../modules/local/mmseqs/easytaxonomy/main'

include { PARSE_MMSEQS_PROTEIN_RESULTS } from '../modules/local/parsers/mmseqs2table_protein/main'
include { PARSE_MMSEQS_NUCLEOTIDE_RESULTS } from '../modules/local/parsers/mmseqs2table/main'

include { BOWTIE2_BUILD as BOWTIE2_CONTIGS_BUILD  } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN as BOWTIE2_CONTIGS_ALIGN  } from '../modules/nf-core/bowtie2/align/main'
include { SAMTOOLS_COVERAGE } from '../modules/nf-core/samtools/coverage/main'
include { SAMTOOLS_FAIDX } from '../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index/main'

include { CREATE_EXTENDED_TABLE } from '../modules/local/parsers/create_extended_table/main'

include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { MULTIQC                                } from '../modules/nf-core/multiqc/main'

workflow VERA {

    // 
    // Объявляем разные каналы
    // 
    
    // каналы для версий и multiqc
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    INPUT_CHECK ()
    
    // входной канал для сырых ридов
    ch_raw_short_reads = INPUT_CHECK.out.raw_short_reads

    // fastqc для сырых ридов
    FASTQC_RAW (
        ch_raw_short_reads
    )
    ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect{it[1]}.ifEmpty([]))

    if ( params.do_kraken2 ) {
            // kraken2 для сырых ридов
    KRAKEN2_RAW (
        ch_raw_short_reads,
        params.kraken2_db,
        params.kraken2_save_output_fastqc,
        params.kraken2_save_reads_assignment
    )
    ch_versions = ch_versions.mix(KRAKEN2_RAW.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2_RAW.out.report.collect{it[1]}.ifEmpty([]))
    }


    // очистка от хоста
    if ( params.host_index && params.host_fasta ) {
        // host_index is a directory of prebuilt *.bt2 files; BOWTIE2_ALIGN
        // globs the staged dir for *.rev.1.bt2 and derives the prefix itself.
        ch_host_index = Channel.fromPath(params.host_index, type: 'dir').map{ [ [], it ] }
        ch_host_fasta = Channel.value(file(params.host_fasta)).map { [ [], it ] }
        BOWTIE2_HOST_ALIGN ( 
            ch_raw_short_reads,
            // возможно, именно тут он начинает бухтеть на то что не нужно first(),
            // но без этого аргумента он выполняет только одно выравнивание ридов!
            ch_host_index.first(),
            ch_host_fasta,
            true,
            false
        )
        ch_versions = ch_versions.mix(BOWTIE2_HOST_ALIGN.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_HOST_ALIGN.out.log.collect{it[1]}.ifEmpty([]))

        ch_no_host_short_reads = BOWTIE2_HOST_ALIGN.out.fastq
    } else {
        ch_no_host_short_reads = ch_raw_short_reads
    }

    // очистка от контаминации
    if ( params.contam_fasta ) {
        // канал для очистки от контаминации
        ch_contam_fasta = Channel.value(file(params.contam_fasta)).map { [ [], it ] }

        BOWTIE2_CONTAM_BUILD ( ch_contam_fasta )
        ch_versions = ch_versions.mix(BOWTIE2_CONTAM_BUILD.out.versions.first())

        // вот тут если не добавлять collect(), то он выполнит только для первого аргумента
        // почитал документацию, если мы возвращаем канал с одним элементом, то он выполняет только для первого
        // элемента, а вот если мы возвращаем список, то он воспринимает список как значение,
        // а значит выполняет операцию с таким значением (инедксом) с каждым элементом из канала ридов 
        BOWTIE2_CONTAM_ALIGN ( 
            ch_no_host_short_reads,
            BOWTIE2_CONTAM_BUILD.out.index,
            ch_contam_fasta,
            true,
            false
        )
        ch_versions = ch_versions.mix(BOWTIE2_CONTAM_ALIGN.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_CONTAM_ALIGN.out.log.collect{it[1]}.ifEmpty([]))

        ch_clean_short_reads = BOWTIE2_CONTAM_ALIGN.out.fastq
    } else {
        ch_clean_short_reads = ch_no_host_short_reads
    }

    // чистим и режем риды
    FASTP (
        ch_clean_short_reads,
        [],
        params.fastp_discard_trimmed_fail,
        params.fastp_save_trimmed_fail,
        params.fastp_save_merged
    )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]}.ifEmpty([]))

    // получаем очищенные риды
    ch_trimmed_short_reads = FASTP.out.reads

    FASTQC_TRIMMED (
        ch_trimmed_short_reads
    )
    ch_versions = ch_versions.mix(FASTQC_TRIMMED.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMMED.out.zip.collect{it[1]}.ifEmpty([]))

    // объединяю риды с разных запусков
    CAT_FASTQ (
        // тут нужно объединить риды по id образца
        ch_trimmed_short_reads
        // для этого выкидываю информацию про run
        .map { 
            meta, reads -> 
                def new_meta = meta - meta.subMap('run') 
            [new_meta, reads]
        }
        // группирую файлы по имени последовательности по id образца
        .groupTuple()
        // объединяю все пути до ридов в один "плоский" список
        .map {
            meta, reads -> [ meta, reads.flatten() ]
        }
    )
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())

    ch_final_reads = CAT_FASTQ.out.reads

    if ( params.do_kraken2 ) {
        // и запускаю кракен для слитых "чистых" ридов
        KRAKEN2_CLEAN (
            ch_final_reads,
            params.kraken2_db,
            params.kraken2_save_output_fastqc,
            params.kraken2_save_reads_assignment
        )
        ch_versions = ch_versions.mix(KRAKEN2_CLEAN.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2_CLEAN.out.report.collect{it[1]}.ifEmpty([]))
    }


    // запускаем сборки
    // чтобы заработал спадес, взял с гита кусок кода, хз почему так работает
    if (params.spades_hmm_path) { ch_spades_hmm = file(params.spades_hmm_path) } else { ch_spades_hmm = [] }

    SPADES (
        // они на вход хотят сразу с нанопорой или пакбио, поэтому модифицируем канал
        ch_final_reads
        .map {
            meta, reads -> [meta, reads, [], []]
        },
        [],
        ch_spades_hmm
    )
    ch_versions = ch_versions.mix(SPADES.out.versions.first())

    ch_assembled_contigs = SPADES.out.contigs.filter { 
        meta, contig -> 
            contig.countFasta() > 0
        }

    // фильтруем контиги по длине
    SEQKIT_SEQ (
        ch_assembled_contigs
    )
    ch_versions = ch_versions.mix(SEQKIT_SEQ.out.versions.first())

    ch_length_filtered_contigs = SEQKIT_SEQ.out.fastx.filter { 
        meta, contig -> 
            contig.countFasta() > 0
        }
    // ch_length_filtered_contigs.view()

    // транслируем контиги в orf рамки
    SEQKIT_TRANSLATE (
        ch_length_filtered_contigs
    )
    ch_versions = ch_versions.mix(SEQKIT_TRANSLATE.out.versions.first())

    ch_generated_orfs = SEQKIT_TRANSLATE.out.fastx.filter { 
        meta, contig -> 
            contig.countFasta() > 0
        }

    // читаем белковую базу данных mmseqs
    ch_mmseqs_protein_db = Channel
        .fromPath(params.mmseqs_easy_taxonomy_protein_db, checkIfExists: true)
        .map { 
            db -> 
                def meta = [:]
                meta.id = "mmseqs_easy_taxonomy_protein_db"
                [ meta, db ]
        }

    // запускаем поиск по белковым последовательностям
    MMSEQS_PROTEINS (
        ch_generated_orfs, 
        ch_mmseqs_protein_db.collect()
    )
    ch_versions = ch_versions.mix(MMSEQS_PROTEINS.out.versions.first())

    // запускаем заключение из белковых последовательностей
    PARSE_MMSEQS_PROTEIN_RESULTS (
        MMSEQS_PROTEINS.out.lca
    )
    ch_versions = ch_versions.mix(PARSE_MMSEQS_PROTEIN_RESULTS.out.versions.first())

    ch_mixed_for_grep = SEQKIT_SEQ.out.fastx.combine(
        PARSE_MMSEQS_PROTEIN_RESULTS.out.as_list, by: 0
    ).filter {
        it -> it[2].countLines() > 0
    }

    ch_fastx = ch_mixed_for_grep.multiMap {
        it -> 
            fastx: [it[0], it[1]]
            list: [it[2]]
    }

    // извлекаем последовательности, которые удалось характеризовать как вирусные
    SEQKIT_GREP_FOR_NUCL (
        ch_fastx.fastx,
        ch_fastx.list
    )
    ch_versions = ch_versions.mix(SEQKIT_GREP_FOR_NUCL.out.versions.first())

    ch_potential_viruses = SEQKIT_GREP_FOR_NUCL.out.filter.filter { 
        meta, contig -> 
            contig.countFasta() > 0
        }

    // читаем нуклеотидную базу данных mmseqs
    ch_mmseqs_viral_nt_db = Channel
        .fromPath(params.mmseqs_easy_taxonomy_refseq_viral_db, checkIfExists: true)
        .map { 
            db -> 
                def meta = [:]
                meta.id = "mmseqs_easy_taxonomy_refseq_viral_db"
                [ meta, db ]
        }

    // для избранных по белковому поиску контигов запускаем поиск нуклеотидный 
    MMSEQS_NUCLEOTIDES (
        ch_potential_viruses, 
        ch_mmseqs_viral_nt_db.collect()
    )
    ch_versions = ch_versions.mix(MMSEQS_NUCLEOTIDES.out.versions.first())

    ch_parse_nucl = MMSEQS_NUCLEOTIDES.out.lca.combine (
        MMSEQS_NUCLEOTIDES.out.tophit_aln, by: 0
    )

    // парсим результаты нуклеотидного поиска
    PARSE_MMSEQS_NUCLEOTIDE_RESULTS (
        ch_parse_nucl.map { it -> [it[0], it[1]]},
        ch_parse_nucl.map { it -> [it[0], it[2]]}
    )
    ch_versions = ch_versions.mix(PARSE_MMSEQS_NUCLEOTIDE_RESULTS.out.versions.first())

    GUNZIP_CONTIGS (
        ch_assembled_contigs
    )

    // вычисление реального покрытия у всех собранных контигов
    // сперва индексируем контиги для bt2
    BOWTIE2_CONTIGS_BUILD ( GUNZIP_CONTIGS.out.gunzip )
    ch_versions = ch_versions.mix(BOWTIE2_CONTIGS_BUILD.out.versions.first())

    ch_for_coverage = ch_final_reads.combine(
        BOWTIE2_CONTIGS_BUILD.out.index, by: 0
    ).combine(
        GUNZIP_CONTIGS.out.gunzip, by:0
    )

    ch_reads = ch_for_coverage.map {
        meta, reads, index, contigs -> [meta, reads]
    }
    ch_index = ch_for_coverage.map {
        meta, reads, index, contigs -> [meta, index]
    }
    ch_coverage_contigs = ch_for_coverage.map {
        meta, reads, index, contigs -> [meta, contigs]
    }

    // затем запускаем картирование ридов, которые использовали для сборки, на избранные контиги
    BOWTIE2_CONTIGS_ALIGN (
        ch_reads,
        ch_index,
        ch_coverage_contigs,
        false,
        true
    )
    ch_versions = ch_versions.mix(BOWTIE2_CONTIGS_ALIGN.out.versions.first())
    // индексируем сами контиги для fai
    SAMTOOLS_FAIDX (
        GUNZIP_CONTIGS.out.gunzip,
        [ [ id:'no_fai' ], [] ]
    )
    // индексируем картирование ридов
    SAMTOOLS_INDEX (
        BOWTIE2_CONTIGS_ALIGN.out.bam
    )
    // формируем канал формата [[meta], bam, bai]
    mixed_bam = BOWTIE2_CONTIGS_ALIGN.out.bam.combine (
        SAMTOOLS_INDEX.out.bai, by: 0
    )

    ch_mixed_coverage = mixed_bam.combine (
        GUNZIP_CONTIGS.out.gunzip, by: 0
    ).combine (
        SAMTOOLS_FAIDX.out.fai, by: 0
    )

    ch_bams = ch_mixed_coverage.map {
        meta, bam, bai, contigs, fai -> [meta, bam, bai]
    }
    ch_contigs = ch_mixed_coverage.map {
        meta, bam, bai, contigs, fai -> [meta, contigs]
    }
    ch_fasta = ch_mixed_coverage.map {
        meta, bam, bai, contigs, fai -> [meta, fai]
    }

    // запускаем вычисление покрытия контигов
    SAMTOOLS_COVERAGE (
        ch_bams,
        ch_contigs,
        ch_fasta
    )

    ch_pivot_channel_tables = PARSE_MMSEQS_NUCLEOTIDE_RESULTS.out.table.combine(
            PARSE_MMSEQS_PROTEIN_RESULTS.out.table, by: 0
        ).combine(
            SAMTOOLS_COVERAGE.out.coverage, by: 0
        ).filter {
            meta, nucl, prot, coverage ->
                prot.countLines() > 1
        }

    CREATE_EXTENDED_TABLE (
        ch_pivot_channel_tables
    )

    ch_final_grep = ch_potential_viruses.combine(
        CREATE_EXTENDED_TABLE.out.as_list, by: 0
    ).filter {
        it -> it[2].countLines() > 0
    }.multiMap {
        it -> 
            fastx: [it[0], it[1]]
            list: [it[2]]
    }

    SEQKIT_GREP_FOR_EXTENSION (
        ch_final_grep.fastx,
        ch_final_grep.list
    )
    ch_versions = ch_versions.mix(SEQKIT_GREP_FOR_EXTENSION.out.versions.first())

    // тут добавляем увеличение контигов

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_mqc_versions.yaml')
    )

    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        [],
        [],
        [],
        []
    )

    multiqc_report = MULTIQC.out.report.toList()

}
