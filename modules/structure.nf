process STRUCTURAL_ORTHOLOGY_EVAL {
    tag "Foldseek, ESM-2 & Tree Analysis"
    publishDir "${params.outdir}/structural_and_phylo_analysis", mode: 'copy'

    input:
    path query_fastas
    path filtered_ortholog_fastas
    path prostt5_model

    output:
    path "aln.m8",               emit: forward_m8
    path "recip_aln.m8",         emit: recip_m8
    path "esm_sim.csv",          emit: esm_csv
    path "aligned.fa",           emit: alignment
    path "iqtree_out.treefile",  emit: treefile

    script:
    """
    set -euo pipefail

    # ==========================================
    # 0. Combine and Deduplicate Input FASTAs
    # ==========================================
    cat ${query_fastas} | seqkit rmdup -s > query_fa.fa

    if [ -n "${filtered_ortholog_fastas}" ]; then
        cat ${filtered_ortholog_fastas} | seqkit rmdup -s > candidates_fa.fa
    else
        touch candidates_fa.fa
    fi

    # ==========================================
    # 1. Foldseek RBH Search (Structural Alignment)
    # ==========================================
    mkdir -p tmp_fs foldseek_db foldseek_out

    # Create Foldseek DBs using the global ProstT5 model
    foldseek createdb query_fa.fa foldseek_db/queryDB --prostt5-model ${prostt5_model} --threads ${task.cpus}
    foldseek createdb candidates_fa.fa foldseek_db/targetDB --prostt5-model ${prostt5_model} --threads ${task.cpus}

    # Forward Search (Query -> Targets)
    foldseek search foldseek_db/queryDB foldseek_db/targetDB foldseek_out/aln_db tmp_fs -a --threads ${task.cpus}
    foldseek convertalis foldseek_db/queryDB foldseek_db/targetDB foldseek_out/aln_db aln.m8 --format-output "query,target,pident,alnlen,bits,evalue"

    # Reciprocal Search (Targets -> Query)
    foldseek search foldseek_db/targetDB foldseek_db/queryDB foldseek_out/recip_aln_db tmp_fs -a --threads ${task.cpus}
    foldseek convertalis foldseek_db/targetDB foldseek_db/queryDB foldseek_out/recip_aln_db recip_aln.m8 --format-output "query,target,pident,alnlen,bits,evalue"

    # ==========================================
    # 2. ESM-2 Deep Learning Embeddings
    # ==========================================
    
    export HF_HOME=\$PWD/.cache
    mkdir -p \$HF_HOME

    python ${projectDir}/scripts/calculate_esm_embeddings.py \
        --query query_fa.fa \
        --candidates candidates_fa.fa \
        --output esm_sim.csv

    # ==========================================
    # 3. Alignment & Phylogeny (MAFFT & IQ-TREE)
    # ==========================================
    cat query_fa.fa candidates_fa.fa | seqkit rmdup -s > combined_clean.fa
    mafft --auto --thread ${task.cpus} combined_clean.fa > aligned.fa

    iqtree -s aligned.fa \
           -m TEST \
           -nt AUTO \
           -pre iqtree_out \
           -fast \
           -redo
    """
}
