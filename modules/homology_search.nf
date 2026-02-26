process HOMOLOGY_SEARCH {
    
    tag "${query.simpleName}_${database.simpleName}"
    publishDir "${params.outdir}/homology_search/${query.simpleName}_${database.simpleName}", 
        mode: 'copy',
        saveAs: { filename ->
            // Only keep hits FASTA and hits list, skip query/database files
            if (filename.endsWith('_hits.fa') || filename.endsWith('_hits.list')) {
                return filename
            }
            return null  // Don't copy query.fsa, database.fasta, etc.
        }

    input:
    tuple path(query), path(database)

    output:
    tuple path(query), path(database), path("${query.simpleName}_vs_${database.simpleName}_hits.fa"), emit: hits_fasta
    path("${query.simpleName}_vs_${database.simpleName}_hits.list"), emit: hits_list

    conda "bioconda::hmmer=3.4"

    script:
    """
    echo "[INFO] Query: $query"
    echo "[INFO] Database: $database"
    echo "[INFO] CPUs: ${task.cpus}"

    # 1. Rename duplicate database FASTA headers
    awk '/^>/ {
        name=\$0; count[name]++;
        if(count[name]>1) name=sprintf("%s_%d", name, count[name]);
        print name; next
    } {print}' $database > db_renamed.fasta

    # 2. Run jackhmmer
    jackhmmer --tblout results.jack --cpu ${task.cpus} -N 10 --noali $query db_renamed.fasta

    # 3. Extract hit IDs (create hits list file)
    grep -v '^#' results.jack | awk '{print \$1}' | sort -u > ${query.simpleName}_vs_${database.simpleName}_hits.list

    # 4. Filter FASTA sequences (hits only)
    awk '
      BEGIN { while ((getline < "${query.simpleName}_vs_${database.simpleName}_hits.list") > 0) ids[\$1]=1 }
      /^>/ { keep=0; header=substr(\$0,2); split(header,a," "); header=a[1]; if(ids[header]){ keep=1; print \$0; next } }
      keep { print \$0 }
    ' db_renamed.fasta > ${query.simpleName}_vs_${database.simpleName}_hits.fa

    echo "[INFO] Hits found: \$(wc -l < ${query.simpleName}_vs_${database.simpleName}_hits.list)"
    echo "[INFO] Done: ${query.simpleName}_vs_${database.simpleName}_hits.fa"
    """
}
