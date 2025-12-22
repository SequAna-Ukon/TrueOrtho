process DOMAIN_SCAN {
    tag "${query.simpleName}_${species}"
    
    publishDir "results/domain_scan/${query.simpleName}_${species}",
        mode: 'copy',
        saveAs: { f -> f.endsWith('.log') ? null : f }

    input:
    tuple path(query), path(orthologs_fa), val(threads), val(target_domain), val(species)
    path hmm_db_dir
    
    output:
    path "${query.simpleName}_${species}_filtered_orthologs.fa", optional: true, emit: filtered_orthologs
    path "${query.simpleName}_${species}_domains.tblout", optional: true, emit: domains_tblout
    path "${query.simpleName}_${species}_query_domains.txt", optional: true, emit: query_domains
    path "${query.simpleName}_${species}_ortholog_domains.txt", optional: true, emit: ortholog_domains
    path "${query.simpleName}_${species}_target_domains.txt", optional: true, emit: target_domains

    conda "bioconda::hmmer=3.4 bioconda::seqkit=2.8.0"

    script:
    """
    #!/bin/bash
    set -euo pipefail

    echo "[INFO] Starting DOMAIN_SCAN for ${query.simpleName}_${species}"
    echo "[INFO] Target domain: '${target_domain}'"
    echo "[INFO] Orthologs: ${orthologs_fa}"
    echo "[INFO] HMM database: ${hmm_db_dir}/Pf_Sm"

    if [ ! -s "${orthologs_fa}" ]; then
        echo "[INFO] Orthologs file is empty or missing. Skipping domain scan."
        exit 0
    fi

    # Set database path
    DB_PATH="${hmm_db_dir}/Pf_Sm"
    
    # === DOMAIN SCAN & FILTER ===
    required_domains_list=""
    if [ -n "${target_domain}" ]; then
        echo "[INFO] Using provided target domain(s): ${target_domain}"
        IFS=',' read -ra arr <<< "${target_domain}"
        for d in "\${arr[@]}"; do
            required_domains_list="\${required_domains_list}\${d}"\$'\\n'
        done
        required_domains_list=\$(echo "\$required_domains_list" | grep -v '^[[:space:]]*\$')
        echo "# Target domains for ${query.simpleName}_${species}" > "${query.simpleName}_${species}_target_domains.txt"
        echo "\$required_domains_list" >> "${query.simpleName}_${species}_target_domains.txt"
    else
        echo "[INFO] Auto-detecting domains in query..."
        hmmscan --domtblout query_domains.tblout --noali -E 1e-5 --cpu ${threads} "\$DB_PATH" "${query}" 2> query_scan.log
        if [ -s "query_domains.tblout" ]; then
            required_domains_list=\$(awk '\$1 !~ /^#/ {print \$1}' query_domains.tblout | sort -u)
            query_domain_count=\$(echo "\$required_domains_list" | wc -l)
            if [ "\$query_domain_count" -gt 0 ]; then
                echo "# Domains found in query ${query.simpleName}" > "${query.simpleName}_${species}_query_domains.txt"
                echo "\$required_domains_list" >> "${query.simpleName}_${species}_query_domains.txt"
            fi
        fi
    fi

    echo "[INFO] Scanning orthologs..."
    hmmscan --domtblout ortholog_domains.tblout --noali -E 1e-5 --cpu ${threads} "\$DB_PATH" "${orthologs_fa}" 2> ortholog_scan.log

    required_count=\$(echo "\$required_domains_list" | wc -l)
    ortholog_ids=\$(seqkit seq --name --only-id "${orthologs_fa}" | sort -u)
    matched_ids_file="matched_ids.tmp"
    > "\$matched_ids_file"

    echo "\$ortholog_ids" | while read seq_id; do
        seq_domains=\$(awk -v seq="\$seq_id" '\$4 == seq && \$1 !~ /^#/ {print \$1}' ortholog_domains.tblout | sort -u)
        missing_count=0
        if [ "\$required_count" -gt 0 ]; then
            while read req; do
                [ -z "\$req" ] && continue
                if ! echo "\$seq_domains" | grep -q "^\$req\$"; then
                    missing_count=\$((missing_count + 1))
                fi
            done <<< "\$required_domains_list"
        fi
        if [ "\$required_count" -eq 0 ] || [ "\$missing_count" -eq 0 ]; then
            echo "\$seq_id" >> "\$matched_ids_file"
        fi
    done

    match_count=\$(wc -l < "\$matched_ids_file" 2>/dev/null || echo 0)
    echo "[INFO] Selected \$match_count orthologs"

    if [ "\$match_count" -gt 0 ]; then
        seqkit grep -f "\$matched_ids_file" "${orthologs_fa}" > "${query.simpleName}_${species}_filtered_orthologs.fa"
        awk '\$1 !~ /^#/ {print \$4 " " \$1}' ortholog_domains.tblout | sort -u > "${query.simpleName}_${species}_ortholog_domains.txt"
        mv ortholog_domains.tblout "${query.simpleName}_${species}_domains.tblout"
    fi

    rm -f "\$matched_ids_file" 2>/dev/null || true
    echo "[INFO] DOMAIN_SCAN completed for ${query.simpleName}_${species}"
    """
}
