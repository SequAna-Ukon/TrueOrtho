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
    // Update USER, REPO, and BRANCH to match your GitHub details
    def github_user   = "SequAna-Ukon"
    def github_repo   = "TrueOrtho"
    def github_branch = "main"
    def raw_github_url = "https://raw.githubusercontent.com/${github_user}/${github_repo}/${github_branch}/databases/Pf_Sm"

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

    # 2. HMM Setup (GitHub Fallback)
    mkdir -p hmm_database

    if [ -f "${domain_db}" ]; then
        echo "[INFO] Using existing local HMM file: ${domain_db}"
        cp \$(readlink -f "${domain_db}") hmm_database/Pf_Sm

    elif [ -d "${domain_db}" ]; then
        echo "[INFO] Using HMM file from directory: ${domain_db}"
        found=\$(find "${domain_db}" -name "*.hmm" -type f -o -name "Pf_Sm*" -type f | head -1)
        cp \$(readlink -f "\$found") hmm_database/Pf_Sm

    elif [ -f "${projectDir}/databases/Pf_Sm" ]; then
        echo "[INFO] Found repo database at ${projectDir}/databases/Pf_Sm"
        cp "${projectDir}/databases/Pf_Sm" hmm_database/Pf_Sm

    else
        echo "[INFO] HMM file not found locally. Downloading from GitHub..."
        wget -q "${raw_github_url}" -O hmm_database/Pf_Sm
    fi

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
