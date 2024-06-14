/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_mag_pipeline'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { BINNING_PREPARATION             } from '../subworkflows/local/binning_preparation'
include { BINNING                         } from '../subworkflows/local/binning'
include { BINNING_REFINEMENT              } from '../subworkflows/local/binning_refinement'
include { BUSCO_QC                        } from '../subworkflows/local/busco_qc'
include { VIRUS_IDENTIFICATION            } from '../subworkflows/local/virus_identification'
include { CHECKM_QC                       } from '../subworkflows/local/checkm_qc'
include { GUNC_QC                         } from '../subworkflows/local/gunc_qc'
include { GTDBTK                          } from '../subworkflows/local/gtdbtk'
include { ANCIENT_DNA_ASSEMBLY_VALIDATION } from '../subworkflows/local/ancient_dna'
include { DOMAIN_CLASSIFICATION           } from '../subworkflows/local/domain_classification'
include { DEPTHS                          } from '../subworkflows/local/depths'

//
// MODULE: Installed directly from nf-core/modules
//
include { ARIA2 as ARIA2_UNTAR                   } from '../modules/nf-core/aria2/main'
include { FASTQC as FASTQC_RAW                   } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIMMED               } from '../modules/nf-core/fastqc/main'
include { SEQTK_MERGEPE                          } from '../modules/nf-core/seqtk/mergepe/main'
include { BBMAP_BBNORM                           } from '../modules/nf-core/bbmap/bbnorm/main'
include { FASTP                                  } from '../modules/nf-core/fastp/main'
include { ADAPTERREMOVAL as ADAPTERREMOVAL_PE    } from '../modules/nf-core/adapterremoval/main'
include { ADAPTERREMOVAL as ADAPTERREMOVAL_SE    } from '../modules/nf-core/adapterremoval/main'
include { CAT_FASTQ                              } from '../modules/nf-core/cat/fastq/main'
include { PRODIGAL                               } from '../modules/nf-core/prodigal/main'
include { PROKKA                                 } from '../modules/nf-core/prokka/main'
include { MMSEQS_DATABASES                       } from '../modules/nf-core/mmseqs/databases/main'
include { METAEUK_EASYPREDICT                    } from '../modules/nf-core/metaeuk/easypredict/main'

//
// MODULE: Local to the pipeline
//
include { BOWTIE2_REMOVAL_BUILD as BOWTIE2_HOST_REMOVAL_BUILD } from '../modules/local/bowtie2_removal_build'
include { BOWTIE2_REMOVAL_ALIGN as BOWTIE2_HOST_REMOVAL_ALIGN } from '../modules/local/bowtie2_removal_align'
include { BOWTIE2_REMOVAL_BUILD as BOWTIE2_PHIX_REMOVAL_BUILD } from '../modules/local/bowtie2_removal_build'
include { BOWTIE2_REMOVAL_ALIGN as BOWTIE2_PHIX_REMOVAL_ALIGN } from '../modules/local/bowtie2_removal_align'
include { PORECHOP                                            } from '../modules/local/porechop'
include { NANOLYSE                                            } from '../modules/local/nanolyse'
include { FILTLONG                                            } from '../modules/local/filtlong'
include { NANOPLOT as NANOPLOT_RAW                            } from '../modules/local/nanoplot'
include { NANOPLOT as NANOPLOT_FILTERED                       } from '../modules/local/nanoplot'
include { CENTRIFUGE_DB_PREPARATION                           } from '../modules/local/centrifuge_db_preparation'
include { CENTRIFUGE                                          } from '../modules/local/centrifuge'
include { KRAKEN2_DB_PREPARATION                              } from '../modules/local/kraken2_db_preparation'
include { KRAKEN2                                             } from '../modules/local/kraken2'
include { KRONA_DB                                            } from '../modules/local/krona_db'
include { KRONA                                               } from '../modules/local/krona'
include { POOL_SINGLE_READS as POOL_SHORT_SINGLE_READS        } from '../modules/local/pool_single_reads'
include { POOL_PAIRED_READS                                   } from '../modules/local/pool_paired_reads'
include { POOL_SINGLE_READS as POOL_LONG_READS                } from '../modules/local/pool_single_reads'
include { MEGAHIT                                             } from '../modules/local/megahit'
include { SPADES                                              } from '../modules/local/spades'
include { SPADESHYBRID                                        } from '../modules/local/spadeshybrid'
include { GUNZIP as GUNZIP_ASSEMBLIES                         } from '../modules/nf-core/gunzip'
include { QUAST                                               } from '../modules/local/quast'
include { QUAST_BINS                                          } from '../modules/local/quast_bins'
include { QUAST_BINS_SUMMARY                                  } from '../modules/local/quast_bins_summary'
include { CAT_DB                                              } from '../modules/local/cat_db'
include { CAT_DB_GENERATE                                     } from '../modules/local/cat_db_generate'
include { CAT                                                 } from '../modules/local/cat'
include { CAT_SUMMARY                                         } from "../modules/local/cat_summary"
include { BIN_SUMMARY                                         } from '../modules/local/bin_summary'
include { COMBINE_TSV as COMBINE_SUMMARY_TSV                  } from '../modules/local/combine_tsv'

////////////////////////////////////////////////////
/* --  Create channel for reference databases  -- */
////////////////////////////////////////////////////

if ( params.host_genome ) {
    host_fasta = params.genomes[params.host_genome].fasta ?: false
    ch_host_fasta = Channel
        .value(file( "${host_fasta}" ))
    host_bowtie2index = params.genomes[params.host_genome].bowtie2 ?: false
    ch_host_bowtie2index = Channel
        .value(file( "${host_bowtie2index}/*" ))
} else if ( params.host_fasta ) {
    ch_host_fasta = Channel
        .value(file( "${params.host_fasta}" ))
} else {
    ch_host_fasta = Channel.empty()
}

if (params.busco_db) {
    ch_busco_db = file(params.busco_db, checkIfExists: true)
} else {
    ch_busco_db = []
}

if(params.checkm_db) {
    ch_checkm_db = file(params.checkm_db, checkIfExists: true)
}

if (params.gunc_db) {
    ch_gunc_db = file(params.gunc_db, checkIfExists: true)
} else {
    ch_gunc_db = Channel.empty()
}

if(params.centrifuge_db){
    ch_centrifuge_db_file = file(params.centrifuge_db, checkIfExists: true)
} else {
    ch_centrifuge_db_file = []
}

if(params.kraken2_db){
    ch_kraken2_db_file = file(params.kraken2_db, checkIfExists: true)
} else {
    ch_kraken2_db_file = []
}

if(params.cat_db){
    ch_cat_db_file = Channel
        .value(file( "${params.cat_db}" ))
} else {
    ch_cat_db_file = Channel.empty()
}

if(params.krona_db){
    ch_krona_db_file = Channel
        .value(file( "${params.krona_db}" ))
} else {
    ch_krona_db_file = Channel.empty()
}

if(!params.keep_phix) {
    ch_phix_db_file = Channel
        .value(file( "${params.phix_reference}" ))
}

if (!params.keep_lambda) {
    ch_nanolyse_db = Channel
        .value(file( "${params.lambda_reference}" ))
}

if (params.genomad_db){
    ch_genomad_db = file(params.genomad_db, checkIfExists: true)
} else {
    ch_genomad_db = Channel.empty()
}

gtdb = ( params.skip_binqc || params.skip_gtdbtk ) ? false : params.gtdb_db

if (gtdb) {
    gtdb = file( "${gtdb}", checkIfExists: true)
    gtdb_mash = params.gtdb_mash ? file("${params.gtdb_mash}", checkIfExists: true) : []
} else {
    gtdb = []
}

if(params.metaeuk_db && !params.skip_metaeuk) {
    ch_metaeuk_db = Channel.
        value(file("${params.metaeuk_db}", checkIfExists: true))
} else {
    ch_metaeuk_db = Channel.empty()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Additional info for completion email and summary
def busco_failed_bins = [:]

workflow MAG {
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Get mmseqs db for MetaEuk if requested
    if (!params.skip_metaeuk && params.metaeuk_mmseqs_db) {
        MMSEQS_DATABASES(params.metaeuk_mmseqs_db)
        ch_metaeuk_db = MMSEQS_DATABASES.out.database
        ch_versions = ch_versions.mix(MMSEQS_DATABASES.out.versions)
    }

    /*
    * Bin QC subworkflows: for checking bin completeness with either BUSCO, CHECKM, and/or GUNC
    */

    ch_input_bins_for_qc = Channel.fromPath("${params.binput_dir}/*.fna")
        | map { bin ->
            def id = bin.getSimpleName()
            def group = "binput"
            def assembler = "Unknown"
            def domain = "Unclassified"
            def meta = [id: id, group: group, assembler: assembler, domain: domain]
            [ meta, bin ]
        }
        | groupTuple(by: 0)

    ch_input_bins_for_qc = ch_input_for_postbinning_bins_unbins.transpose()

    if (!params.skip_binqc && params.binqc_tool == 'busco'){
        /*
        * BUSCO subworkflow: Quantitative measures for the assessment of genome assembly
        */

        BUSCO_QC (
            ch_busco_db,
            ch_input_bins_for_qc
        )
        ch_busco_summary = BUSCO_QC.out.summary
        ch_versions = ch_versions.mix(BUSCO_QC.out.versions.first())
        // process information if BUSCO analysis failed for individual bins due to no matching genes
        BUSCO_QC.out
            .failed_bin
            .splitCsv(sep: '\t')
            .map { bin, error -> if (!bin.contains(".unbinned.")) busco_failed_bins[bin] = error }
    }

    if (!params.skip_binqc && params.binqc_tool == 'checkm'){
        /*
        * CheckM subworkflow: Quantitative measures for the assessment of genome assembly
        */

        ch_input_bins_for_checkm = ch_input_bins_for_qc
            .filter { meta, bins ->
                meta.domain != "eukarya"
            }

        CHECKM_QC (
            ch_input_bins_for_checkm.groupTuple(),
            ch_checkm_db
        )
        ch_checkm_summary = CHECKM_QC.out.summary

        ch_versions = ch_versions.mix(CHECKM_QC.out.versions)

    }

    if ( params.run_gunc && params.binqc_tool == 'checkm' ) {
        GUNC_QC ( ch_input_bins_for_checkm, ch_gunc_db, CHECKM_QC.out.checkm_tsv )
        ch_versions = ch_versions.mix( GUNC_QC.out.versions )
    } else if ( params.run_gunc ) {
        ch_input_bins_for_gunc = ch_input_for_postbinning_bins_unbins
            .filter { meta, bins ->
                meta.domain != "eukarya"
            }
        GUNC_QC ( ch_input_bins_for_qc, ch_gunc_db, [] )
        ch_versions = ch_versions.mix( GUNC_QC.out.versions )
    }

    ch_quast_bins_summary = Channel.empty()
    if (!params.skip_quast){
        ch_input_for_quast_bins = ch_input_for_postbinning_bins_unbins
                                    .groupTuple()
                                    .map {
                                        meta, bins ->
                                            def new_bins = bins.flatten()
                                            [meta, new_bins]
                                        }

        QUAST_BINS ( ch_input_for_quast_bins )
        ch_versions = ch_versions.mix(QUAST_BINS.out.versions.first())
        ch_quast_bin_summary = QUAST_BINS.out.quast_bin_summaries
            .collectFile(keepHeader: true) {
                meta, summary ->
                ["${meta.id}.tsv", summary]
        }
        QUAST_BINS_SUMMARY ( ch_quast_bin_summary.collect() )
        ch_quast_bins_summary = QUAST_BINS_SUMMARY.out.summary
    }

    /*
        * CAT: Bin Annotation Tool (BAT) are pipelines for the taxonomic classification of long DNA sequences and metagenome assembled genomes (MAGs/bins)
        */
    ch_cat_db = Channel.empty()
    if (params.cat_db){
        CAT_DB ( ch_cat_db_file )
        ch_cat_db = CAT_DB.out.db
    } else if (params.cat_db_generate){
        CAT_DB_GENERATE ()
        ch_cat_db = CAT_DB_GENERATE.out.db
    }
    CAT (
        ch_input_for_postbinning_bins_unbins,
        ch_cat_db
    )
    // Group all classification results for each sample in a single file
    ch_cat_summary = CAT.out.tax_classification_names
        .collectFile(keepHeader: true) {
                meta, classification ->
                ["${meta.id}.txt", classification]
        }   
    // Group all classification results for the whole run in a single file
    CAT_SUMMARY(
        ch_cat_summary.collect()
    )
    ch_versions = ch_versions.mix(CAT.out.versions.first())
    ch_versions = ch_versions.mix(CAT_SUMMARY.out.versions)

    // If CAT is not run, then the CAT global summary should be an empty channel
    if ( params.cat_db_generate || params.cat_db) {
        ch_cat_global_summary = CAT_SUMMARY.out.combined
    } else {
        ch_cat_global_summary = Channel.empty()
    }

    /*
        * GTDB-tk: taxonomic classifications using GTDB reference
        */

    if ( !params.skip_gtdbtk ) {

        ch_gtdbtk_summary = Channel.empty()
        if ( gtdb ){

            ch_gtdb_bins = ch_input_for_postbinning_bins_unbins
                .filter { meta, bins ->
                    meta.domain != "eukarya"
                }

            GTDBTK (
                ch_gtdb_bins,
                ch_busco_summary,
                ch_checkm_summary,
                gtdb,
                gtdb_mash
            )
            ch_versions = ch_versions.mix(GTDBTK.out.versions.first())
            ch_gtdbtk_summary = GTDBTK.out.summary
        }
    } else {
        ch_gtdbtk_summary = Channel.empty()
    }

    if ( ( !params.skip_binqc ) || !params.skip_quast || !params.skip_gtdbtk){
        BIN_SUMMARY (
            ch_input_for_binsummary,
            ch_busco_summary.ifEmpty([]),
            ch_checkm_summary.ifEmpty([]),
            ch_quast_bins_summary.ifEmpty([]),
            ch_gtdbtk_summary.ifEmpty([]),
            ch_cat_global_summary.ifEmpty([])
        )
    }

    /*
        * Prokka: Genome annotation
        */

    if (!params.skip_prokka){
        ch_bins_for_prokka = ch_input_for_postbinning_bins_unbins.transpose()
        .map { meta, bin ->
            def meta_new = meta + [id: bin.getBaseName()]
            [ meta_new, bin ]
        }
        .filter { meta, bin ->
            meta.domain != "eukarya"
        }

        PROKKA (
            ch_bins_for_prokka,
            [],
            []
        )
        ch_versions = ch_versions.mix(PROKKA.out.versions.first())
    }

    if (!params.skip_metaeuk && (params.metaeuk_db || params.metaeuk_mmseqs_db)) {
        ch_bins_for_metaeuk = ch_input_for_postbinning_bins_unbins.transpose()
            .filter { meta, bin ->
                meta.domain in ["eukarya", "unclassified"]
            }
            .map { meta, bin ->
                def meta_new = meta + [id: bin.getBaseName()]
                [ meta_new, bin ]
            }

        METAEUK_EASYPREDICT (ch_bins_for_metaeuk, ch_metaeuk_db)
        ch_versions = ch_versions.mix(METAEUK_EASYPREDICT.out.versions)
    }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
