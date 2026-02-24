process DATABASE_SETUP {
    tag "database_setup"
    
    input:
    val eggnog_db
    val domain_db    // Back to val so it doesn't crash on 'Pf_Sm' string
    val work_dir
    
    output:
    path "eggnog_database", emit: egg_dir
    path "hmm_database", emit: hmm_dir
    
    conda "bioconda::hmmer=3.4"
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    FINAL_DB_DIR="${work_dir}/databases"
    mkdir -p "\$FINAL_DB_DIR"
    
    # 1. eggNOG Setup
    if [ -n "${eggnog_db}" ] && [ "${eggnog_db}" != "null" ] && [ -d "${eggnog_db}" ]; then
        echo "[INFO] Linking existing eggNOG DB"
        ln -snf \$(readlink -f "${eggnog_db}") "\$FINAL_DB_DIR/eggnog_database"
    else
        echo "[INFO] Downloading eggNOG DB"
        mkdir -p eggnog_local
        cd eggnog_local
        BASE_URL="http://eggnog5.embl.de/download/emapperdb-5.0.2"
        wget -q "\$BASE_URL/eggnog.db.gz" "\$BASE_URL/eggnog_proteins.dmnd.gz" "\$BASE_URL/eggnog.taxa.tar.gz"
        gunzip *.gz && tar -xzf *.tar.gz && rm -f *.tar.gz
        cd ..
        rm -rf "\$FINAL_DB_DIR/eggnog_database"
        mv eggnog_local "\$FINAL_DB_DIR/eggnog_database"
    fi
    
    # 2. HMM Setup
    mkdir -p hmm_local
    # Check if domain_db is a real file on the system
    if [ -f "${domain_db}" ]; then
        echo "[INFO] Using existing HMM file"
        cp \$(readlink -f "${domain_db}") hmm_local/Pf_Sm
    elif [ -d "${domain_db}" ]; then
        echo "[INFO] Using HMM from directory"
        found=\$(find "${domain_db}" -name "*.hmm" -type f | head -1)
        cp \$(readlink -f "\$found") hmm_local/Pf_Sm
    else
        echo "[INFO] HMM file not found locally. Downloading..."
        wget -q "https://cloud.uni-konstanz.de/index.php/s/MMzHFEBjZLrpmN7/download" -O hmm_local/Pf_Sm
    fi

    # Press the HMM
    hmmpress -f hmm_local/Pf_Sm
    
    # Move to central work/databases/
    rm -rf "\$FINAL_DB_DIR/hmm_database"
    mv hmm_local "\$FINAL_DB_DIR/hmm_database"

    # 3. Satisfy Nextflow outputs
    ln -snf "\$FINAL_DB_DIR/eggnog_database" eggnog_database
    ln -snf "\$FINAL_DB_DIR/hmm_database" hmm_database
    """
}
