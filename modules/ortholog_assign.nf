process ORTHOLOG_ASSIGN {
    tag "${query.simpleName}_${species}"

    // Publish only the three key files
    publishDir "results/ortholog_assign/${query.simpleName}_${species}", mode: 'copy', pattern: '*_orthologs.fa'
    publishDir "results/ortholog_assign/${query.simpleName}_${species}", mode: 'copy', pattern: '*_hits.txt', optional: true
    publishDir "results/ortholog_assign/${query.simpleName}_${species}", mode: 'copy', pattern: '*_kog_info.txt', optional: true

    input:
    tuple path(query), path(hits_fasta), val(threads), val(kog_id), val(eggnog_db), val(species)

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
    echo "=================================="

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

    # EggNOG DB setup
    eggnog_db_dir=""
    if [ -n "${eggnog_db}" ] && [ "${eggnog_db}" != "null" ] && [ -d "${eggnog_db}" ]; then
        # Remove trailing slash if present
        eggnog_db_dir=\$(echo "${eggnog_db}" | sed 's:/*\$::')
        echo "[INFO] Using provided eggNOG database: \$eggnog_db_dir"
        
        # Check if database files exist
        if [ ! -f "\$eggnog_db_dir/eggnog_proteins.dmnd" ]; then
            echo "[ERROR] eggNOG database not found at: \$eggnog_db_dir/eggnog_proteins.dmnd"
            exit 1
        fi
    else
        echo "[INFO] No valid eggNOG database provided, using default location"
        eggnog_db_dir="\${workflow.workDir}/eggnog_database"
        mkdir -p "\$eggnog_db_dir"
        echo "[INFO] Downloading eggNOG data..."
        download_eggnog_data.py -y --data_dir "\$eggnog_db_dir" || {
            echo "[ERROR] Failed to download eggNOG data"
            exit 1
        }
    fi
    
    export EGGNOG_DATA_DIR="\$eggnog_db_dir"
    echo "[INFO] EGGNOG_DATA_DIR set to: \$EGGNOG_DATA_DIR"

    # Annotate query
    echo "[INFO] Annotating query sequence..."
    emapper.py -m diamond --cpu ${task.cpus} --data_dir "\$eggnog_db_dir" -i "${query}" -o query_emapper 2>&1 | tee query_annotation.log
    for i in {1..30}; do [ -s query_emapper.emapper.annotations ] && break; sleep 2; done
    if [ ! -s query_emapper.emapper.annotations ]; then
        echo "[ERROR] Query annotation failed or produced empty results"
        echo "Query annotation failed" >> "\$kog_info_file"
        touch "\$orthologs_file" "\$hits_file"
        exit 0
    fi

    # Annotate hits
    echo "[INFO] Annotating hit sequences..."
    emapper.py -m diamond --cpu ${task.cpus} --data_dir "\$eggnog_db_dir" -i "${hits_fasta}" -o hits_emapper 2>&1 | tee hits_annotation.log
    for i in {1..30}; do [ -s hits_emapper.emapper.annotations ] && break; sleep 2; done
    if [ ! -s hits_emapper.emapper.annotations ]; then
        echo "[WARNING] Hits annotation failed or produced empty results"
        echo "Hits annotation failed - possibly no significant matches" >> "\$kog_info_file"
        touch "\$orthologs_file" "\$hits_file"
        exit 0
    fi

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
