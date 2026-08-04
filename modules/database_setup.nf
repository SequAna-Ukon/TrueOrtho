process DATABASE_SETUP {
    tag "database_setup"
    storeDir "${params.outdir}/databases"

    input:
    val eggnog_db
    val domain_db
    val prostt5_db

    output:
    path "eggnog_database",   emit: egg_dir
    path "hmm_database",      emit: hmm_dir
    path "prostt5_model",     emit: prostt5_dir

    script:
    """
    #!/bin/bash
    set -euo pipefail

    # 1. eggNOG Setup
    if [ -n "${eggnog_db}" ] && [ "${eggnog_db}" != "null" ] && [ -d "${eggnog_db}" ]; then
        echo "[INFO] Linking existing eggNOG DB"
        ln -snf \$(readlink -f "${eggnog_db}") eggnog_database
    else
        echo "[INFO] Downloading eggNOG DB"
        mkdir -p eggnog_database
        cd eggnog_database
        BASE_URL="http://eggnog5.embl.de/download/emapperdb-5.0.2"
        wget -q "\$BASE_URL/eggnog.db.gz" "\$BASE_URL/eggnog_proteins.dmnd.gz" "\$BASE_URL/eggnog.taxa.tar.gz"
        gunzip *.gz && tar -xzf *.tar.gz && rm -f *.tar.gz
        cd ..
    fi

    # 2. HMM Setup
    mkdir -p hmm_database

    copy_or_decompress() {
        local src="\$1"
        if [[ "\$src" == *.gz ]]; then
            echo "[INFO] Decompressing \$src -> hmm_database/Pf_Sm"
            gzip -dc "\$src" > hmm_database/Pf_Sm
        else
            echo "[INFO] Copying \$src -> hmm_database/Pf_Sm"
            cp "\$src" hmm_database/Pf_Sm
        fi
    }

    if [ -f "${domain_db}" ]; then
        echo "[INFO] Using user-specified local HMM file: ${domain_db}"
        copy_or_decompress \$(readlink -f "${domain_db}")

    elif [ -d "${domain_db}" ]; then
        echo "[INFO] Searching for HMM file in directory: ${domain_db}"
        found=\$(find "${domain_db}" -maxdepth 2 -type f -name "*.hmm*" -o -name "Pf_Sm*" | head -1)
        if [ -n "\$found" ]; then
            copy_or_decompress \$(readlink -f "\$found")
        else
            echo "[ERROR] No matching HMM file found in directory ${domain_db}" >&2
            exit 1
        fi

    elif [ -f "${projectDir}/databases/Pf_Sm.gz" ]; then
        echo "[INFO] Found repo database at ${projectDir}/databases/Pf_Sm.gz"
        copy_or_decompress "${projectDir}/databases/Pf_Sm.gz"

    elif [ -f "${projectDir}/databases/Pf_Sm" ]; then
        echo "[INFO] Found repo database at ${projectDir}/databases/Pf_Sm"
        copy_or_decompress "${projectDir}/databases/Pf_Sm"

    else
        echo "[ERROR] HMM database not found! Checked --domain_db and ${projectDir}/databases/Pf_Sm.gz" >&2
        exit 1
    fi

    echo "[INFO] Formatting HMM database with hmmpress..."
    hmmpress -f hmm_database/Pf_Sm

    # 3. Foldseek ProstT5 Setup
    if [ -n "${prostt5_db}" ] && [ "${prostt5_db}" != "null" ] && [ -d "${prostt5_db}" ]; then
        echo "[INFO] Linking existing Foldseek ProstT5 model"
        ln -snf \$(readlink -f "${prostt5_db}") prostt5_model
    elif [ ! -d prostt5_model ]; then
        echo "[INFO] Pre-existing ProstT5 model not found. Downloading..."
        mkdir -p tmp_fs
        foldseek databases ProstT5 prostt5_model tmp_fs --threads ${task.cpus}
        rm -rf tmp_fs
    fi
    """
}
