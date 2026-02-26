process DOMAIN_SCAN {
    tag "${query.simpleName}_${species}"
    
    publishDir "${params.outdir}/domain_scan/${query.simpleName}_${species}",
        mode: 'copy',
        saveAs: { f -> f.endsWith('.log') ? null : f }

    input:
    tuple path(query), path(orthologs_fa), val(threads), val(target_domain), val(species)
    path hmm_db_dir
    
    output:
    path "${query.simpleName}_${species}_filtered_orthologs.fa", emit: filtered_orthologs
    path "${query.simpleName}_${species}_ortholog_domains.txt", emit: ortholog_domains
    path "${query.simpleName}_${species}_domains.tblout", optional: true, emit: domains_tblout
    path "${query.simpleName}_${species}_query_domains.txt", optional: true, emit: query_domains
    path "${query.simpleName}_${species}_target_domains.txt", optional: true, emit: target_domains

    conda "bioconda::hmmer=3.4 bioconda::seqkit=2.8.0"

    script:
    """
    #!/bin/bash
    set -euo pipefail

    # --- Initialize output files ---
    OUT_FA="${query.simpleName}_${species}_filtered_orthologs.fa"
    OUT_TXT="${query.simpleName}_${species}_ortholog_domains.txt"
    touch "\$OUT_FA" "\$OUT_TXT"

    if [ ! -s "${orthologs_fa}" ]; then
        exit 0
    fi

    # --- FIX: Database path relative to work/databases/ ---
    DB_PATH="${hmm_db_dir}/Pf_Sm"
    
    required_domains_file="required.list"
    > "\$required_domains_file"

    if [ -n "${target_domain}" ] && [ "${target_domain}" != "null" ]; then
        echo "${target_domain}" | tr ',' '\\n' | grep -v '^[[:space:]]*\$' > "\$required_domains_file"
        cp "\$required_domains_file" "${query.simpleName}_${species}_target_domains.txt"
    else
        hmmscan --domtblout query_domains.tblout --noali -E 1e-5 --cpu ${threads} "\$DB_PATH" "${query}" > query_scan.log 2>&1
        if [ -s "query_domains.tblout" ]; then
            awk '\$1 !~ /^#/ {print \$1}' query_domains.tblout | sort -u > "\$required_domains_file"
            cp "\$required_domains_file" "${query.simpleName}_${species}_query_domains.txt"
        fi
    fi

    hmmscan --domtblout ortholog_domains.tblout --noali -E 1e-5 --cpu ${threads} "\$DB_PATH" "${orthologs_fa}" > ortholog_scan.log 2>&1

    required_count=\$(wc -l < "\$required_domains_file")
    matched_ids_file="matched_ids.tmp"
    > "\$matched_ids_file"

    seqkit seq --name --only-id "${orthologs_fa}" > all_ids.txt

    while read -r seq_id; do
        seq_domains=\$(awk -v seq="\$seq_id" '\$4 == seq && \$1 !~ /^#/ {print \$1}' ortholog_domains.tblout | sort -u)
        missing_count=0
        if [ "\$required_count" -gt 0 ]; then
            while read -r req; do
                if ! echo "\$seq_domains" | grep -qx "\$req"; then
                    missing_count=\$((missing_count + 1))
                fi
            done < "\$required_domains_file"
        fi

        if [ "\$missing_count" -eq 0 ]; then
            echo "\$seq_id" >> "\$matched_ids_file"
        fi
    done < all_ids.txt

    if [ -s "\$matched_ids_file" ]; then
        seqkit grep -f "\$matched_ids_file" "${orthologs_fa}" > "\$OUT_FA"
        awk 'FNR==NR {ids[\$1]; next} \$4 in ids && \$1 !~ /^#/ {print \$4 " " \$1}' "\$matched_ids_file" ortholog_domains.tblout | sort -u > "\$OUT_TXT"
        mv ortholog_domains.tblout "${query.simpleName}_${species}_domains.tblout"
    fi
    """
}
