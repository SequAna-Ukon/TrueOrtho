#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { HOMOLOGY_SEARCH } from './modules/homology_search.nf'
include { DATABASE_SETUP } from './modules/database_setup.nf'
include { ORTHOLOG_ASSIGN } from './modules/ortholog_assign.nf'
include { DOMAIN_SCAN } from './modules/domain_scan.nf'
include { GENERATE_SUMMARY_REPORT } from './modules/sum_repo.nf'

workflow {
    def input_path = params.input ?: params.csv_file
    
    if (!input_path) {
        error "ERROR: No input CSV specified. Please use --input <file>"
    }

    csv_file_obj = file(input_path)

    params.threads   = params.threads ?: 10
    params.domain_db = params.domain_db ?: ""
    params.eggnog_db = params.eggnog_db ?: ""
    params.outdir = params.outdir ?: "Results"

    // 1. Setup Input Channel
    input_ch = Channel.fromPath(csv_file_obj)
        .splitCsv(header:true)
        .map { row ->
            def q_file = file(row.query)
            def db_file = file(row.database)
            tuple(q_file.simpleName, db_file.simpleName, q_file, db_file, row.kog_id ?: "", row.target_domain ?: "")
        }

    // 2. Setup Databases (Passing workflow.workDir to place them in work/databases/)
    db_results = DATABASE_SETUP(params.eggnog_db, params.domain_db, workflow.workDir)

    // 3. Homology Search
    homology_results = HOMOLOGY_SEARCH(input_ch.map { qid, sp, q, db, k, td -> tuple(q, db) })

    // 4. Ortholog Assignment
    ortho_input_ch = input_ch
        .map { qid, sp, q, db, k, td -> tuple(qid, sp, q, k) }
        .join(homology_results.hits_fasta.map { q, db, fa -> tuple(q.simpleName, db.simpleName, fa) }, by: [0, 1])
        .map { qid, sp, q_file, k_id, hits_fa -> tuple(q_file, hits_fa, params.threads, k_id, sp) }

    // Use the explicit path emitted from DATABASE_SETUP
    ortholog_results = ORTHOLOG_ASSIGN(ortho_input_ch, db_results.egg_dir)

    // 5. Domain Scan
    domain_input_ch = ortholog_results.orthologs_fa
        .filter { q, fa, sp -> fa.size() > 0 }
        .map { q, fa, sp -> tuple(q.simpleName, sp, q, fa) }
        .join(input_ch.map { qid, sp, q, db, k, td -> tuple(qid, sp, td) }, by: [0, 1])
        .map { qid, sp, q_file, ortho_fa, td -> tuple(q_file, ortho_fa, params.threads, td, sp) }

    domain_results = DOMAIN_SCAN(domain_input_ch, db_results.hmm_dir)

    // 6. Generate Summary Report
    GENERATE_SUMMARY_REPORT(
        homology_results.hits_list.collect().ifEmpty([]),
        ortholog_results.orthologs_fa.map { it[1] }.collect().ifEmpty([]),
        domain_results.filtered_orthologs.collect().ifEmpty([]),
        domain_results.ortholog_domains.collect().ifEmpty([])
    )
}
