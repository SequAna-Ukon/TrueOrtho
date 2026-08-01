process STRUCTURAL_ORTHOLOGY_EVAL {
    tag "${qid}"
    publishDir "${params.outdir}/structural_and_phylo_analysis/${qid}", mode: 'copy'

    input:
    tuple val(qid), path(query_fasta), val(sp_list), path(fa_list)
    path prostt5_model

    output:
    path "${qid}_*_fwd.m8",   emit: forward_m8
    path "${qid}_*_recip.m8", emit: recip_m8
    path "${qid}_*_esm_sim.csv",   emit: esm_csv
    path "${qid}_iqtree.treefile", optional: true, emit: treefile
    path "${qid}_aligned.fa",      emit: alignment

    script:
    """
    set -euo pipefail

    cat ${query_fasta} | seqkit rmdup -s > query_fa.fa

    sp_arr=(${sp_list.join(' ')})
    fa_arr=(${fa_list.join(' ')})

    touch all_candidates.fa

    for i in "\${!sp_arr[@]}"; do
        sp="\${sp_arr[\$i]}"
        fa_file="\${fa_arr[\$i]}"

        if [ -f "\$fa_file" ] && [ "\$fa_file" != "NO_FILE" ]; then
            cat "\$fa_file" | seqkit rmdup -s > "candidates_\${sp}.fa"
        else
            touch "candidates_\${sp}.fa"
        fi

        cand_count=\$(grep -c "^>" "candidates_\${sp}.fa" 2>/dev/null || true)
        cand_count=\${cand_count:-0}

        if [ "\$cand_count" -eq 0 ]; then
            echo "[WARNING] No candidate orthologs for ${qid} vs \${sp} — skipping Foldseek/ESM for this species."
            touch "${qid}_\${sp}_fwd.m8" "${qid}_\${sp}_recip.m8"
            echo "target,esm2_cosine_sim" > "${qid}_\${sp}_esm_sim.csv"
        else
            mkdir -p "tmp_fs_\${sp}" "foldseek_db_\${sp}" "foldseek_out_\${sp}"

            foldseek createdb query_fa.fa "foldseek_db_\${sp}/queryDB" --prostt5-model ${prostt5_model} --threads ${task.cpus}
            foldseek createdb "candidates_\${sp}.fa" "foldseek_db_\${sp}/targetDB" --prostt5-model ${prostt5_model} --threads ${task.cpus}

            foldseek search "foldseek_db_\${sp}/queryDB" "foldseek_db_\${sp}/targetDB" "foldseek_out_\${sp}/aln_db" "tmp_fs_\${sp}" -a --threads ${task.cpus}
            foldseek convertalis "foldseek_db_\${sp}/queryDB" "foldseek_db_\${sp}/targetDB" "foldseek_out_\${sp}/aln_db" "${qid}_\${sp}_fwd.m8" --format-output "query,target,pident,alnlen,bits,evalue"

            foldseek search "foldseek_db_\${sp}/targetDB" "foldseek_db_\${sp}/queryDB" "foldseek_out_\${sp}/recip_aln_db" "tmp_fs_\${sp}" -a --threads ${task.cpus}
            foldseek convertalis "foldseek_db_\${sp}/targetDB" "foldseek_db_\${sp}/queryDB" "foldseek_out_\${sp}/recip_aln_db" "${qid}_\${sp}_recip.m8" --format-output "query,target,pident,alnlen,bits,evalue"

            export HF_HOME="\$PWD/.cache_\${sp}"
            mkdir -p "\$HF_HOME"
            python ${projectDir}/scripts/calculate_esm_embeddings.py \
                --query query_fa.fa \
                --candidates "candidates_\${sp}.fa" \
                --output "${qid}_\${sp}_esm_sim.csv"
        fi

            cat "candidates_\${sp}.fa" >> all_candidates.fa
        done
   
    # Pool Alignment & Phylogeny 
  
    cat query_fa.fa all_candidates.fa | seqkit rmdup -s > "${qid}_pooled.fa"
    mafft --auto --thread ${task.cpus} "${qid}_pooled.fa" > "${qid}_aligned.fa"

    seq_count=\$(grep -c "^>" "${qid}_aligned.fa" 2>/dev/null || true)
    seq_count=\${seq_count:-0}

    if [ "\$seq_count" -ge 3 ]; then
        echo "[INFO] Running IQ-TREE for query ${qid} with \$seq_count pooled sequences..."
        iqtree -s "${qid}_aligned.fa" \
               -m TEST \
               -nt AUTO \
               -pre "${qid}_iqtree" \
               -fast \
               -redo
    else
        echo "[WARNING] Skipping IQ-TREE for ${qid}: only \$seq_count sequence(s) (minimum 3 required)."
        touch "${qid}_iqtree.treefile"
    fi
    """
}
