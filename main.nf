#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { DATABASE_SETUP }            from './modules/database_setup.nf'
include { HOMOLOGY_SEARCH }           from './modules/homology_search.nf'
include { ORTHOLOG_ASSIGN }           from './modules/ortholog_assign.nf'
include { DOMAIN_SCAN }               from './modules/domain_scan.nf'
include { STRUCTURAL_ORTHOLOGY_EVAL } from './modules/structure.nf'
include { GENERATE_SUMMARY_REPORT }   from './modules/sum_repo.nf'


workflow {

    // 1. Parameter Validation & Channel Construction
    def input_path = params.input ?: params.csv_file

    if (!input_path) {
        error "ERROR: No input CSV specified. Please provide --input <file.csv>"
    }

    csv_file_obj = file(input_path)

    no_file_placeholder = file("${projectDir}/assets/NO_FILE")
    if (!no_file_placeholder.exists()) {
        no_file_placeholder.parent.mkdirs()
        no_file_placeholder.text = ''
    }

    input_ch = Channel.fromPath(csv_file_obj)
        .splitCsv(header: true, quote: '"')
        .map { row ->
            def q_file  = file(row.query)
            def db_file = file(row.database)
            tuple(
                q_file.simpleName,
                db_file.simpleName,
                q_file,
                db_file,
                row.kog_id?.trim() ?: "",
                row.target_domain?.trim() ?: ""
            )
        }

    // 2. Setup Databases (eggNOG, HMM, Foldseek ProstT5)
    db_results = DATABASE_SETUP(
        params.eggnog_db ?: "",
        params.domain_db ?: "",
        params.prostt5_db ?: ""
    )

    // 3. Homology Search (jackhmmer)
    homology_results = HOMOLOGY_SEARCH(
        input_ch.map { qid, sp, q, db, k, td -> tuple(q, db) }
    )

    // 4. EggNOG Ortholog Assignment
    ortho_input_ch = input_ch
        .map { qid, sp, q, db, k, td -> tuple(qid, sp, q, k) }
        .join(
            homology_results.hits_fasta.map { q, db, fa -> tuple(q.simpleName, db.simpleName, fa) },
            by: [0, 1]
        )
        .map { qid, sp, q_file, k_id, hits_fa -> tuple(q_file, hits_fa, params.threads, k_id, sp) }

    ortholog_results = ORTHOLOG_ASSIGN(
        ortho_input_ch,
        db_results.egg_dir
    )

    // 5. Domain Scan & Domain Architecture Validation (HMMER)
    domain_input_ch = ortholog_results.orthologs_fa
        .filter { q, fa, sp -> fa.size() > 0 }
        .map { q, fa, sp -> tuple(q.simpleName, sp, q, fa) }
        .join(
            input_ch.map { qid, sp, q, db, k, td -> tuple(qid, sp, td) },
            by: [0, 1]
        )
        .map { qid, sp, q_file, ortho_fa, td -> tuple(q_file, ortho_fa, params.threads, td, sp) }

    domain_results = DOMAIN_SCAN(
        domain_input_ch,
        db_results.hmm_dir
    )

    // 6. Channel Pairing & Grouping
    struct_input_ch = input_ch
        .map { row -> tuple(row[0], row[1], row[2]) }
        .join(
            domain_results.filtered_orthologs.map { q, sp, fa -> tuple(q, sp, fa) },
            by: [0, 1],
            remainder: true
        )
        .map { qid, sp, q_file, fa -> tuple(qid, q_file, sp, fa ?: no_file_placeholder) }
        .groupTuple(by: [0, 1])

    // 7. Structural, Embedding, and Phylogenetic Analysis

    struct_results = STRUCTURAL_ORTHOLOGY_EVAL(
        struct_input_ch,
        db_results.prostt5_dir
    )

    // 8. Integrated Classification & HTML Report Generation
    GENERATE_SUMMARY_REPORT(
        struct_results.forward_m8.collect().ifEmpty([]),
        struct_results.recip_m8.collect().ifEmpty([]),
        struct_results.esm_csv.collect().ifEmpty([]),
        struct_results.treefile.collect().ifEmpty([]),
        homology_results.hits_list.collect().ifEmpty([]),
        ortholog_results.orthologs_fa.map { it[1] }.collect().ifEmpty([]),
        domain_results.filtered_orthologs.map { it[2] }.collect().ifEmpty([]),
        domain_results.ortholog_domains.collect().ifEmpty([]),
        file("${projectDir}/scripts/classify_and_report.py")
    )
}
