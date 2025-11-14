#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { HOMOLOGY_SEARCH } from './modules/homology_search.nf'
include { ORTHOLOG_ASSIGN } from './modules/ortholog_assign.nf'
include { DOMAIN_SCAN } from './modules/domain_scan.nf'
include { GENERATE_SUMMARY_REPORT } from './modules/sum_repo.nf'

workflow {
    params.csv_file  = file(params.csv_file)
    params.threads   = params.threads ?: 10
    params.domain_db = params.domain_db ?: ""
    params.eggnog_db = params.eggnog_db ?: ""

    // Read CSV and create a channel with query simpleName and species as keys
    query_db_ch = Channel.fromPath(params.csv_file)
        .splitCsv(header:true)
        .map { row ->
            tuple(
                row.query.tokenize('/')[-1].tokenize('.')[0], // query simpleName
                file(row.database).simpleName,               // species
                file(row.query),
                file(row.database),
                row.kog_id ?: "",
                row.target_domain ?: ""
            )
        }

    // HOMOLOGY_SEARCH
    homology_results = HOMOLOGY_SEARCH(
        query_db_ch.map { q_name, species, q, db, k, td -> tuple(q, db) }
    )

    // ORTHOLOG_ASSIGN
    ortho_input_ch = query_db_ch
        .map { q_name, species, q, db, k, td -> tuple(q_name, species, q, db, k, td) }
        .join(
            homology_results.hits_fasta.map { q, db, hits -> tuple(q.simpleName, db.simpleName, q, db, hits) },
            by: [0, 1]
        )
        .map { q_name, species, q, db, k, td, q2, db2, hits ->
            tuple(q, hits, params.threads, k, params.eggnog_db, species)
        }

    ortholog_results = ORTHOLOG_ASSIGN(ortho_input_ch)

    // DOMAIN_SCAN - FIXED VERSION
    domain_input_ch = ortholog_results.orthologs_fa
        .map { q, fa, species ->
            def fa_file = file(fa)
            tuple(q.simpleName, species, q, fa_file, species, fa_file.exists() && fa_file.size() > 0)
        }
        .join(
            query_db_ch.map { q_name, species, q, db, k, td -> tuple(q_name, species, q, db, k, td) },
            by: [0, 1]
        )
        .map { q_name, species, q1, fa, species2, fa_valid, q2, db, k, td ->
            (q1.name == q2.name && fa_valid && species == species2) ? 
                tuple(q1, fa, params.threads, td, params.domain_db, species) : null
        }
        .filter { it != null }

    domain_results = DOMAIN_SCAN(domain_input_ch)

    // Collect outputs for summary report
    homology_hits_ch = homology_results.hits_list.collect()
    ortholog_fastas_ch = ortholog_results.orthologs_fa.map { q, fa, species -> fa }.collect()
    filtered_orthologs_ch = domain_results.filtered_orthologs.collect()
    domain_files_ch = domain_results.ortholog_domains.collect()

    // Generate summary report
    GENERATE_SUMMARY_REPORT(
        homology_hits_ch,
        ortholog_fastas_ch,
        filtered_orthologs_ch,
        domain_files_ch
    )
}
