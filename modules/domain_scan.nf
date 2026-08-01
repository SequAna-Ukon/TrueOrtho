process DOMAIN_SCAN {
    tag "${query.simpleName}_${species}"
    
    publishDir "${params.outdir}/domain_scan/${query.simpleName}_${species}",
        mode: 'copy',
        saveAs: { f -> f.endsWith('.log') ? null : f }

    input:
    tuple path(query), path(orthologs_fa), val(threads), val(target_domain), val(species)
    path hmm_db_dir
    
    output:
    tuple val("${query.simpleName}"), val(species), path("${query.simpleName}_${species}_filtered_orthologs.fa"), emit: filtered_orthologs
    path "${query.simpleName}_${species}_ortholog_domains.txt", emit: ortholog_domains
    path "${query.simpleName}_${species}_domains.tblout", optional: true, emit: domains_tblout
    path "${query.simpleName}_${species}_query_domains.txt", optional: true, emit: query_domains
    path "${query.simpleName}_${species}_target_domains.txt", optional: true, emit: target_domains


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

    DB_PATH="${hmm_db_dir}/Pf_Sm"
    
    # 1. Prepare required domains 
    required_domains_file="required.list"
    if [ -n "${target_domain}" ] && [ "${target_domain}" != "null" ]; then
        # This part handles "Domain1, Domain2" by splitting at the comma and trimming spaces
        echo "${target_domain}" | tr ',' '\\n' | sed 's/^[[:space:]]*//;s/[[:space:]]*\$//' | grep -v '^[[:space:]]*\$' | sort -u > "\$required_domains_file"
        cp "\$required_domains_file" "${query.simpleName}_${species}_target_domains.txt"
    else
        hmmscan --domtblout query_domains.tblout --noali -E 1e-5 --cpu ${threads} "\$DB_PATH" "${query}" > query_scan.log 2>&1
        if [ -s "query_domains.tblout" ]; then
            awk '\$1 !~ /^#/ {print \$1}' query_domains.tblout | sort -u > "\$required_domains_file"
            cp "\$required_domains_file" "${query.simpleName}_${species}_query_domains.txt"
        fi
    fi

    # 2. Scan the candidate orthologs
    hmmscan --domtblout ortholog_domains.tblout --noali -E 1e-5 --cpu ${threads} "\$DB_PATH" "${orthologs_fa}" > ortholog_scan.log 2>&1

    # 3. Validation Logic
    num_req=\$(wc -l < "\$required_domains_file")
    
    if [ "\$num_req" -eq 0 ]; then
        cp "${orthologs_fa}" "\$OUT_FA"
    else
        awk -v req_count="\$num_req" '
            FNR==NR { req[\$1]; next } 
            (\$1 in req) { matches[\$4][\$1] } 
            END { 
                for (seq in matches) {
                    count = 0; for (d in matches[seq]) count++;
                    if (count == req_count) print seq
                }
            }' "\$required_domains_file" ortholog_domains.tblout > matched_ids.tmp

        if [ -s matched_ids.tmp ]; then
            seqkit grep -f matched_ids.tmp "${orthologs_fa}" > "\$OUT_FA"
            awk 'FNR==NR {ids[\$1]; next} (\$4 in ids) && (\$1 !~ /^#/) {print \$4 "\t" \$1}' matched_ids.tmp ortholog_domains.tblout | sort -u > "\$OUT_TXT"
        fi
    fi

    mv ortholog_domains.tblout "${query.simpleName}_${species}_domains.tblout"
    """
}
