process ORTHOLOG_ASSIGN {
    tag "${query.simpleName}_${species}"
    
    // Publish only the three key files
    publishDir "results/ortholog_assign/${query.simpleName}_${species}", mode: 'copy', pattern: '*_orthologs.fa'
    publishDir "results/ortholog_assign/${query.simpleName}_${species}", mode: 'copy', pattern: '*_hits.txt', optional: true
    publishDir "results/ortholog_assign/${query.simpleName}_${species}", mode: 'copy', pattern: '*_kog_info.txt', optional: true

    input:
    tuple path(query), path(hits_fasta), val(threads), val(kog_id), val(species)
    path eggnog_db_dir
    
    output:
    tuple path(query), path("${query.simpleName}_${species}_orthologs.fa"), val(species), emit: orthologs_fa
    path "${query.simpleName}_${species}_hits.txt", optional: true, emit: hits_to_extract
    path "${query.simpleName}_${species}_kog_info.txt", optional: true, emit: kog_info

    conda "bioconda::eggnog-mapper=2.1.13 bioconda::seqkit=2.8.0"

    script:
    """
    #!/bin/bash
    set -euo pipefail

    orthologs_file="${query.simpleName}_${species}_orthologs.fa"
    hits_file="${query.simpleName}_${species}_hits.txt"
    kog_info_file="${query.simpleName}_${species}_kog_info.txt"
    prefix="${query.simpleName}_${species}"

    exec > "\${prefix}_process.log" 2>&1
    
    echo "=== ORTHOLOG_ASSIGN DEBUG INFO ==="
    echo "[INFO] Query: ${query}"
    echo "[INFO] Hits file: ${hits_fasta}"
    echo "[INFO] Hits file size: \$(wc -l < '${hits_fasta}' 2>/dev/null || echo 0) lines"
    echo "[INFO] Species: ${species}"
    echo "[INFO] KOG ID: '${kog_id}'"
    echo "[INFO] Threads: ${threads}"
    echo "[INFO] EggNOG database directory: ${eggnog_db_dir}"
    echo "=================================="

    # --- FIX: Find database in work/databases/ ---
    EGGNOG_DB_PATH=""
    if [ -f "${eggnog_db_dir}/eggnog.db" ]; then
        EGGNOG_DB_PATH="${eggnog_db_dir}"
    elif [ -f "${eggnog_db_dir}/data/eggnog.db" ]; then
        EGGNOG_DB_PATH="${eggnog_db_dir}/data"
    else
        found_db=\$(find -L "${eggnog_db_dir}" -name "eggnog.db" -type f | head -1)
        if [ -n "\$found_db" ]; then
            EGGNOG_DB_PATH=\$(dirname "\$found_db")
        else
            echo "[ERROR] eggnog.db not found"
            exit 1
        fi
    fi
    export EGGNOG_DATA_DIR="\$EGGNOG_DB_PATH"

    # List contents for debug
    find -L "${eggnog_db_dir}" -type f -name "*.db" -o -name "*.dmnd" -o -name "*.txt" 2>/dev/null | head -20
    
    if [ ! -s "${hits_fasta}" ]; then
        echo "[WARNING] Hits file is empty or missing."
        {
            echo "# Ortholog assignment for \${prefix}"
            echo "Status: No hits found"
        } > "\$kog_info_file"
        touch "\$orthologs_file" "\$hits_file"
        exit 0
    fi

    fasta_seq_count=\$(grep -c '^>' "${hits_fasta}" 2>/dev/null || echo 0)
    if [ "\$fasta_seq_count" -eq 0 ]; then
        touch "\$orthologs_file" "\$hits_file"
        exit 0
    fi

    # Annotate query
    emapper.py -m diamond --cpu ${task.cpus} --data_dir "\$EGGNOG_DATA_DIR" -i "${query}" -o query_emapper
    
    if [ ! -s query_emapper.emapper.annotations ]; then
        touch "\$orthologs_file" "\$hits_file"
        exit 0
    fi

    # Annotate hits
    emapper.py -m diamond --cpu ${task.cpus} --data_dir "\$EGGNOG_DATA_DIR" -i "${hits_fasta}" -o hits_emapper
    
    if [ ! -s hits_emapper.emapper.annotations ]; then
        touch "\$orthologs_file" "\$hits_file"
        exit 0
    fi

    # KOG logic
    target_kog=""
    if [ -n "${kog_id}" ] && [ "${kog_id}" != "null" ] && [ "${kog_id}" != "" ]; then
        target_kog="${kog_id}"
    else
        target_kog=\$(awk -F'\\t' 'NR>1 && \$5 != "-" {split(\$5, a, ","); for(i in a) if(a[i] ~ /KOG[0-9]+@1\\|root\$/) {print a[i]; exit}}' query_emapper.emapper.annotations)
        if [ -z "\$target_kog" ]; then
            target_kog=\$(awk -F'\\t' 'NR>1 && \$5 != "-" {split(\$5, a, ","); for(i in a) if(a[i] ~ /COG[0-9]+@1\\|root\$/) {print a[i]; exit}}' query_emapper.emapper.annotations)
        fi
    fi

    if [ -n "\$target_kog" ]; then
        echo "Used KOG: \$target_kog" >> "\$kog_info_file"
        awk -F'\\t' -v t="\$target_kog" 'NR>1 && \$5 != "-" {split(\$5, a, ","); for(i in a) if(a[i] == t) print \$1}' hits_emapper.emapper.annotations | sort -u > temp_hits.txt
        if [ ! -s temp_hits.txt ]; then
            base_kog=\$(echo "\$target_kog" | cut -d'@' -f1)
            awk -F'\\t' -v b="\$base_kog" 'NR>1 && \$5 != "-" {split(\$5, a, ","); for(i in a) if(a[i] ~ "^" b "@") print \$1}' hits_emapper.emapper.annotations | sort -u > temp_hits.txt
        fi
        mv temp_hits.txt "\$hits_file"
    fi

    if [ -s "\$hits_file" ]; then
        seqkit grep -f "\$hits_file" "${hits_fasta}" > "\$orthologs_file"
    else
        touch "\$orthologs_file"
    fi
    """
}
