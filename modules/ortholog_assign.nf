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

    # List contents of database directory for debugging
    echo "[DEBUG] Contents of eggnog database directory:"
    find "${eggnog_db_dir}" -type f -name "*.db" -o -name "*.dmnd" -o -name "*.txt" 2>/dev/null | head -20
    
    # Check if hits file exists and has content (at least one sequence)
    if [ ! -s "${hits_fasta}" ]; then
        echo "[WARNING] Hits file is empty or missing. No orthologs to process."
        {
            echo "# Ortholog assignment for \${prefix}"
            echo "Query: ${query}"
            echo "Target species: ${species}"
            echo "Provided KOG/COG ID: ${kog_id}"
            echo "Date: \$(date)"
            echo "Status: No hits found in homology search - empty hits file"
            echo "Hits file: ${hits_fasta}"
            echo "Hits file size: \$(wc -c < '${hits_fasta}' 2>/dev/null || echo 0) bytes"
        } > "\$kog_info_file"
        # Create empty output files
        touch "\$orthologs_file" "\$hits_file"
        echo "[INFO] Created empty output files and exiting gracefully"
        exit 0
    fi

    # Check if hits file contains valid FASTA sequences (not just headers)
    fasta_seq_count=\$(grep -c '^>' "${hits_fasta}" 2>/dev/null || echo 0)
    if [ "\$fasta_seq_count" -eq 0 ]; then
        echo "[WARNING] Hits file contains no FASTA sequences (no '>' headers found)"
        {
            echo "# Ortholog assignment for \${prefix}"
            echo "Query: ${query}"
            echo "Target species: ${species}"
            echo "Provided KOG/COG ID: ${kog_id}"
            echo "Date: \$(date)"
            echo "Status: No valid FASTA sequences in hits file"
            echo "FASTA sequence count: \$fasta_seq_count"
            echo "First few lines of hits file:"
            head -n 5 "${hits_fasta}" 2>/dev/null || echo "File not readable"
        } > "\$kog_info_file"
        touch "\$orthologs_file" "\$hits_file"
        echo "[INFO] Created empty output files and exiting gracefully"
        exit 0
    fi

    echo "[INFO] Hits file contains \$fasta_seq_count sequences"

    # Always create KOG info (even if no orthologs)
    {
        echo "# Ortholog assignment for \${prefix}"
        echo "Query: ${query}"
        echo "Target species: ${species}"
        echo "Provided KOG/COG ID: ${kog_id}"
        echo "Date: \$(date)"
        echo "Input hits: \$fasta_seq_count sequences"
    } > "\$kog_info_file"

    # Find the eggnog.db file - it might be in the directory or in a subdirectory
    EGGNOG_DB_PATH=""
    if [ -f "${eggnog_db_dir}/eggnog.db" ]; then
        EGGNOG_DB_PATH="${eggnog_db_dir}"
    elif [ -f "${eggnog_db_dir}/data/eggnog.db" ]; then
        EGGNOG_DB_PATH="${eggnog_db_dir}/data"
    else
        # Search for eggnog.db in any subdirectory
        found_db=\$(find "${eggnog_db_dir}" -name "eggnog.db" -type f | head -1)
        if [ -n "\$found_db" ]; then
            EGGNOG_DB_PATH=\$(dirname "\$found_db")
        else
            echo "[ERROR] eggnog.db not found in ${eggnog_db_dir} or subdirectories"
            echo "[INFO] Searching for any .db files:"
            find "${eggnog_db_dir}" -name "*.db" -type f 2>/dev/null
            exit 1
        fi
    fi
    
    echo "[INFO] Found eggnog.db at: \$EGGNOG_DB_PATH"
    export EGGNOG_DATA_DIR="\$EGGNOG_DB_PATH"
    echo "[INFO] EGGNOG_DATA_DIR set to: \$EGGNOG_DATA_DIR"
    
    # Check for other essential files in the same directory
    if [ ! -f "\$EGGNOG_DATA_DIR/eggnog_proteins.dmnd" ]; then
        # Try to find it in the original directory
        if [ -f "${eggnog_db_dir}/eggnog_proteins.dmnd" ]; then
            echo "[INFO] Found eggnog_proteins.dmnd in ${eggnog_db_dir}"
            # Create symlink or copy
            ln -sf "${eggnog_db_dir}/eggnog_proteins.dmnd" "\$EGGNOG_DATA_DIR/eggnog_proteins.dmnd" 2>/dev/null || \
            cp "${eggnog_db_dir}/eggnog_proteins.dmnd" "\$EGGNOG_DATA_DIR/eggnog_proteins.dmnd"
        else
            echo "[WARNING] eggnog_proteins.dmnd not found in \$EGGNOG_DATA_DIR"
            echo "[INFO] emapper.py might fail or try to download it"
        fi
    fi

    # Annotate query
    echo "[INFO] Annotating query sequence..."
    echo "[DEBUG] Running: emapper.py -m diamond --cpu ${task.cpus} --data_dir \"\$EGGNOG_DATA_DIR\" -i \"${query}\" -o query_emapper"
    emapper.py -m diamond --cpu ${task.cpus} --data_dir "\$EGGNOG_DATA_DIR" -i "${query}" -o query_emapper 2>&1 | tee query_annotation.log
    
    # Check if annotation succeeded
    if [ ! -s query_emapper.emapper.annotations ]; then
        echo "[ERROR] Query annotation failed or produced empty results"
        echo "[INFO] Checking if annotation files were created:"
        ls -la query_emapper.* 2>/dev/null || echo "No query_emapper files found"
        echo "[INFO] Last 20 lines of query_annotation.log:"
        tail -20 query_annotation.log
        echo "Query annotation failed" >> "\$kog_info_file"
        touch "\$orthologs_file" "\$hits_file"
        exit 0
    fi
    
    echo "[INFO] Query annotation completed successfully"

    # Annotate hits
    echo "[INFO] Annotating hit sequences..."
    echo "[DEBUG] Running: emapper.py -m diamond --cpu ${task.cpus} --data_dir \"\$EGGNOG_DATA_DIR\" -i \"${hits_fasta}\" -o hits_emapper"
    emapper.py -m diamond --cpu ${task.cpus} --data_dir "\$EGGNOG_DATA_DIR" -i "${hits_fasta}" -o hits_emapper 2>&1 | tee hits_annotation.log
    
    if [ ! -s hits_emapper.emapper.annotations ]; then
        echo "[WARNING] Hits annotation failed or produced empty results"
        echo "[INFO] Checking if annotation files were created:"
        ls -la hits_emapper.* 2>/dev/null || echo "No hits_emapper files found"
        echo "Hits annotation failed - possibly no significant matches" >> "\$kog_info_file"
        touch "\$orthologs_file" "\$hits_file"
        exit 0
    fi

    echo "[INFO] Hits annotation completed successfully"

    # KOG logic
    target_kog=""
    if [ -n "${kog_id}" ] && [ "${kog_id}" != "null" ] && [ "${kog_id}" != "" ]; then
        target_kog="${kog_id}"
        echo "[INFO] Using provided KOG ID: \$target_kog"
    else
        echo "[INFO] Auto-detecting KOG from query annotation..."
        target_kog=\$(awk -F'\\t' 'NR>1 && \$5 != "-" {split(\$5, a, ","); for(i in a) if(a[i] ~ /KOG[0-9]+@1\\|root\$/) {print a[i]; exit}}' query_emapper.emapper.annotations)
        if [ -z "\$target_kog" ]; then
            # Try COG pattern if KOG not found
            target_kog=\$(awk -F'\\t' 'NR>1 && \$5 != "-" {split(\$5, a, ","); for(i in a) if(a[i] ~ /COG[0-9]+@1\\|root\$/) {print a[i]; exit}}' query_emapper.emapper.annotations)
        fi
        if [ -z "\$target_kog" ]; then
            echo "[INFO] No KOG/COG found in query annotation"
            echo "No KOG/COG in query" >> "\$kog_info_file"
            touch "\$orthologs_file" "\$hits_file"
            exit 0
        fi
        echo "[INFO] Auto-detected KOG: \$target_kog"
    fi
    
    echo "Used KOG: \$target_kog" >> "\$kog_info_file"

    # Find matches
    echo "[INFO] Finding orthologs with KOG: \$target_kog"
    awk -F'\\t' -v t="\$target_kog" 'NR>1 && \$5 != "-" {split(\$5, a, ","); for(i in a) if(a[i] == t) print \$1}' hits_emapper.emapper.annotations | sort -u > temp_hits.txt

    if [ ! -s temp_hits.txt ]; then
        echo "[INFO] No exact KOG matches found, trying base KOG..."
        base_kog=\$(echo "\$target_kog" | cut -d'@' -f1)
        awk -F'\\t' -v b="\$base_kog" 'NR>1 && \$5 != "-" {split(\$5, a, ","); for(i in a) if(a[i] ~ "^" b "@") print \$1}' hits_emapper.emapper.annotations | sort -u > temp_hits.txt
        if [ ! -s temp_hits.txt ]; then
            echo "[INFO] No orthologs found with KOG pattern: \$base_kog"
            echo "No orthologs found" >> "\$kog_info_file"
            touch "\$orthologs_file" "\$hits_file"
            exit 0
        fi
    fi

    mv temp_hits.txt "\$hits_file"
    
    # Extract ortholog sequences
    echo "[INFO] Extracting ortholog sequences..."
    seqkit grep -f "\$hits_file" "${hits_fasta}" > "\$orthologs_file"
    
    if [ ! -s "\$orthologs_file" ]; then
        echo "[WARNING] No sequences extracted despite having hits"
        echo "No sequences extracted from hits file" >> "\$kog_info_file"
    else
        count=\$(seqkit stats -T "\$orthologs_file" | tail -n +2 | awk '{print \$4}')
        echo "Extracted orthologs: \$count" >> "\$kog_info_file"
        echo "[INFO] Successfully extracted \$count orthologs"
    fi

    echo "[INFO] ORTHOLOG_ASSIGN completed successfully"
    """
}
