process DATABASE_SETUP {
    tag "database_setup"
    
    publishDir "databases", mode: 'copy'
    
    input:
    val eggnog_db
    val domain_db
    
    output:
    tuple path("eggnog_database"), path("hmm_database"), emit: db_dirs
    
    conda "bioconda::hmmer=3.4"
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    echo "=== DATABASE SETUP (eggNOG + HMM) ==="
    
    # =================== eggNOG DATABASE SETUP ===================
    echo "[INFO] Setting up eggNOG database..."
    
    # If a database path is provided, use it
    if [ -n "${eggnog_db}" ] && [ "${eggnog_db}" != "null" ] && [ -d "${eggnog_db}" ]; then
        echo "[INFO] Using provided eggNOG database directory"
        # Create a symlink to the provided database
        ln -sf "${eggnog_db}" eggnog_database
    else
        echo "[INFO] No valid eggNOG database provided, downloading essential eggNOG 5.0.2 database..."
        
        # Create local directory
        mkdir -p eggnog_database
        cd eggnog_database
        
        # Base URL for eggNOG 5.0.2 database
        BASE_URL="http://eggnog5.embl.de/download/emapperdb-5.0.2"
        
        echo "[INFO] Downloading essential eggNOG files from: \$BASE_URL"
        
        # Download the 3 essential files
        ESSENTIAL_FILES=(
            "eggnog.db.gz"
            "eggnog_proteins.dmnd.gz"
            "eggnog.taxa.tar.gz"
        )
        
        for file in "\${ESSENTIAL_FILES[@]}"; do
            echo "[INFO] Downloading \$file..."
            wget -q "\$BASE_URL/\$file" || {
                echo "[ERROR] Failed to download \$file"
                exit 1
            }
        done
        
        # Extract files
        echo "[INFO] Extracting eggnog.db.gz..."
        gunzip eggnog.db.gz || {
            echo "[ERROR] Failed to extract eggnog.db.gz"
            exit 1
        }
        
        echo "[INFO] Extracting eggnog_proteins.dmnd.gz..."
        gunzip eggnog_proteins.dmnd.gz || {
            echo "[ERROR] Failed to extract eggnog_proteins.dmnd.gz"
            exit 1
        }
        
        echo "[INFO] Extracting eggnog.taxa.tar.gz..."
        tar -xzf eggnog.taxa.tar.gz || {
            echo "[ERROR] Failed to extract eggnog.taxa.tar.gz"
            exit 1
        }
        
        # Clean up compressed files
        rm -f *.gz *.tar.gz 2>/dev/null || true
        
        echo "[INFO] eggNOG database downloaded successfully"
        cd ..
    fi
    
    # Verify eggNOG database exists
    echo "[INFO] Verifying eggNOG database files..."
    if [ ! -f "eggnog_database/eggnog.db" ]; then
        echo "[ERROR] Missing essential file: eggnog.db"
        ls -la eggnog_database/ 2>/dev/null || echo "No files found"
        exit 1
    fi
    
    if [ ! -f "eggnog_database/eggnog_proteins.dmnd" ]; then
        echo "[ERROR] Missing essential file: eggnog_proteins.dmnd"
        exit 1
    fi
    
    echo "[INFO] eggNOG database setup completed"
    
    # =================== HMM DATABASE SETUP ===================
    echo "[INFO] Setting up HMM database..."
    
    # If a database path is provided, use it
    if [ -n "${domain_db}" ] && [ "${domain_db}" != "null" ] && [ -f "${domain_db}" ]; then
        echo "[INFO] Using provided HMM database file"
        # Create directory and copy the database
        mkdir -p hmm_database
        
        if [ -f "${domain_db}" ]; then
            cp "${domain_db}" hmm_database/Pf_Sm
            echo "[INFO] Copied ${domain_db} to hmm_database/Pf_Sm"
        elif [ -d "${domain_db}" ]; then
            if [ -f "${domain_db}/Pf_Sm" ]; then
                cp "${domain_db}/Pf_Sm" hmm_database/Pf_Sm
                echo "[INFO] Copied ${domain_db}/Pf_Sm to hmm_database/Pf_Sm"
            else
                found_hmm=\$(find "${domain_db}" -name "*.hmm" -type f | head -1)
                if [ -n "\$found_hmm" ]; then
                    cp "\$found_hmm" hmm_database/Pf_Sm
                    echo "[INFO] Copied \$found_hmm to hmm_database/Pf_Sm"
                else
                    echo "[ERROR] No HMM database file found in ${domain_db}"
                    exit 1
                fi
            fi
        fi
        
        # Press the database if needed
        if [ ! -f "hmm_database/Pf_Sm.h3m" ]; then
            echo "[INFO] Pressing HMM database..."
            cd hmm_database
            hmmpress Pf_Sm || {
                echo "[ERROR] Failed to press HMM database"
                exit 1
            }
            cd ..
        fi
    else
        echo "[INFO] No valid HMM database provided, downloading Pfam/SMART database..."
        
        # Create local directory
        mkdir -p hmm_database
        cd hmm_database
        
        DB_FILE="Pf_Sm"
        if [ ! -f "\$DB_FILE" ]; then
            echo "[INFO] Downloading Pfam/SMART database..."
            wget -q "https://cloud.uni-konstanz.de/index.php/s/MMzHFEBjZLrpmN7/download" -O "\$DB_FILE" || {
                echo "[ERROR] Failed to download Pfam/SMART database"
                exit 1
            }
            echo "[INFO] Database downloaded successfully"
        else
            echo "[INFO] Database already exists, skipping download"
        fi
        
        # Press the database if needed
        if [ ! -f "\$DB_FILE.h3m" ]; then
            echo "[INFO] Pressing HMM database..."
            hmmpress "\$DB_FILE" || {
                echo "[ERROR] Failed to press HMM database"
                exit 1
            }
            echo "[INFO] HMM database pressed successfully"
        else
            echo "[INFO] HMM database already pressed, skipping"
        fi
        
        cd ..
    fi
    
    # Verify HMM database exists
    if [ ! -f "hmm_database/Pf_Sm.h3m" ]; then
        echo "[ERROR] HMM database not found or not pressed properly"
        echo "[INFO] Available files in hmm_database/:"
        ls -la hmm_database/ 2>/dev/null || echo "No files found"
        exit 1
    fi
    
    echo "[INFO] HMM database setup completed"
    
    # =================== FINAL OUTPUT ===================
    echo "[INFO] === DATABASE SETUP COMPLETED ==="
    echo "[INFO] eggNOG database: \$(readlink -f eggnog_database)"
    echo "[INFO] HMM database: \$(readlink -f hmm_database)"
    """
}
