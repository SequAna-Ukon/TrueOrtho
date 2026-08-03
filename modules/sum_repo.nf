process GENERATE_SUMMARY_REPORT {
    tag "Final Analysis & Report Generation"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path forward_m8
    path recip_m8
    path esm_csv
    path treefile
    path homology_hits_lists
    path ortholog_fastas
    path filtered_ortholog_fastas
    path domain_files
    path script_file

    output:
    path "final_orthology_evidence_report.csv", emit: report_csv
    path "orthology_evidence_plot_*.png", optional: true, emit: scatter_plots
    path "annotated_tree_*.png", optional: true, emit: tree_plots
    path "summary_counts.tsv",                   emit: summary_counts
    path "summary_report.html",                  emit: html_report

    script:
    """
    #!/bin/bash
    set -euo pipefail
    shopt -s nullglob

    # Fix Matplotlib and Fontconfig cache directory
    export MPLCONFIGDIR="/tmp/matplotlib_${params.outdir}"
    export FC_CACHEDIR="/tmp/fontconfig_${params.outdir}"
    mkdir -p "\$MPLCONFIGDIR" "\$FC_CACHEDIR"

    echo "[INFO] Running Python orthology classification script..."

    # 1. Execute classification 

    python ${script_file} \
        --forward-m8 *_fwd.m8 \
        --recip-m8 *_recip.m8 \
        --esm-csv *_esm_sim.csv \
        --treefile *.treefile \
        --out-report "final_orthology_evidence_report.csv" \
        --out-plot-prefix "orthology_evidence_plot" \
        --out-tree-prefix "annotated_tree"

    echo "[INFO] Calculating sample statistics across all pipeline stages..."

    # 2. Build summary counts TSV 
    echo -e "Sample\tHomology_Hits\tEggNOG_Orthologs\tDomain_Screened\tFinal_Structural_Orthologs" > summary_counts.tsv

    hits_files=( *.list )

    if [ \${#hits_files[@]} -gt 0 ]; then
        for hits_file in "\${hits_files[@]}"; do
            [ -f "\$hits_file" ] || continue

            base_name=\$(basename "\$hits_file" _hits.list)
            base_name=\$(basename "\$base_name" .list)

            query=\$(echo "\$base_name" | sed 's/_vs_/ /' | cut -d' ' -f1)
            species=\$(echo "\$base_name" | sed 's/_vs_/ /' | cut -d' ' -f2)
            sample="\${query}_\${species}"

            # 1. Initial Homology Hits
            hits_count=\$(wc -l < "\$hits_file" 2>/dev/null || echo "0")

            # 2. Reciprocal / EggNOG Orthologs
            ortho_file="\${query}_\${species}_orthologs.fa"
            ortho_count=0
            if [ -f "\$ortho_file" ] && [ -s "\$ortho_file" ]; then
                ortho_count=\$(seqkit seq --name --only-id "\$ortho_file" | wc -l 2>/dev/null || echo "0")
            fi

            # 3. Domain Screened
            final_file="\${query}_\${species}_filtered_orthologs.fa"
            domain_count=0
            if [ -f "\$final_file" ] && [ -s "\$final_file" ]; then
                domain_count=\$(seqkit seq --name --only-id "\$final_file" | wc -l 2>/dev/null || echo "0")
            fi

            # 4. Final Structural True Orthologs
            #
            # NOTE (fix): final_orthology_evidence_report.csv columns are:
            #   1:target 2:esm2_cosine_sim 3:struct_score 4:orthology_confidence_score
            #   5:is_rbh 6:is_self_hit 7:classification
            # classification is column 7, not 6 (it shifted when is_self_hit was added
            # to classify_and_report.py). Also, self-hits (a query found against its
            # own source species) are trivially the highest-confidence ortholog case,
            # not a failure case, so SELF_HIT rows now count toward the final total
            # alongside PRIMARY/SECONDARY/CO_ORTHOLOG classes.
            final_struct_count=0
            if [ -f "final_orthology_evidence_report.csv" ] && [ -f "\$final_file" ]; then
                seqkit seq --name --only-id "\$final_file" > "sample_ids.tmp"
                final_struct_count=\$(awk -F',' 'NR==FNR {ids[\$1]=1; next} (\$1 in ids) && (\$7 ~ /PRIMARY|SECONDARY|CO_ORTHOLOG|SELF_HIT/) {count++} END {print count+0}' "sample_ids.tmp" final_orthology_evidence_report.csv)
                rm -f "sample_ids.tmp"
            fi

            echo -e "\$sample\t\$hits_count\t\$ortho_count\t\$domain_count\t\$final_struct_count" >> summary_counts.tsv
        done
    fi

    echo "[INFO] Generating integrated HTML report..."

    # 3. Build HTML report template
    cat > summary_report.html << 'EOF'
<html>
<head>
    <title>TrueOrtho Pipeline Summary & Classification Report</title>
    <style>
        body { font-family: 'Segoe UI', Helvetica, Arial, sans-serif; margin: 30px; line-height: 1.6; color: #2c3e50; background: #f8f9fa; }
        .container { max-width: 1150px; margin: 0 auto; background: #ffffff; padding: 35px; border-radius: 8px; box-shadow: 0 4px 10px rgba(0,0,0,0.05); }
        h1 { color: #1a365d; border-bottom: 2px solid #e2e8f0; padding-bottom: 12px; margin-top: 0; }
        h2 { color: #2b6cb0; margin-top: 35px; border-bottom: 1px solid #edf2f7; padding-bottom: 8px; }
        h3 { color: #2d3748; margin-bottom: 8px; }
        table { border-collapse: collapse; width: 100%; margin: 15px 0; font-size: 0.95em; }
        th, td { border: 1px solid #e2e8f0; padding: 10px 12px; text-align: left; }
        th { background-color: #edf2f7; font-weight: 600; color: #2d3748; }
        tr:nth-child(even) { background-color: #f7fafc; }
        .section { margin: 30px 0; }
        .count { font-weight: bold; color: #2b6cb0; }
        .domains { font-size: 0.85em; color: #c53030; background: #fff5f5; padding: 2px 6px; border-radius: 4px; font-family: monospace; }
        .img-card { text-align: center; margin: 25px 0; background: #fafafa; border: 1px solid #e2e8f0; padding: 15px; border-radius: 8px; }
        .img-card img { max-width: 100%; height: auto; border-radius: 4px; }
        .sample-card { background: #ffffff; padding: 20px; border-radius: 6px; margin-bottom: 25px; border: 1px solid #e2e8f0; border-left: 5px solid #3182ce; box-shadow: 0 2px 4px rgba(0,0,0,0.02); }
        .badge { font-weight: bold; padding: 3px 8px; border-radius: 4px; font-size: 0.82em; text-transform: uppercase; display: inline-block; }
        .badge-primary { background: #c6f6d5; color: #22543d; }
        .badge-secondary { background: #ebf8ff; color: #2c5282; }
        .badge-self { background: #fefcbf; color: #744210; }
        .badge-other { background: #edf2f7; color: #4a5568; }
    </style>
</head>
<body>
<div class="container">
    <h1>TrueOrtho Pipeline Summary & Classification Report</h1>
    <p><strong>Execution Date:</strong> SYSTEM_DATE_PLACEHOLDER</p>

    <div class="section">
        <h2>Orthology Filtering Funnel Counts</h2>
        <table>
            <thead>
                <tr>
                    <th>Sample</th>
                    <th>Homology Hits</th>
                    <th>EggNOG Orthologs</th>
                    <th>Domain Screened</th>
                    <th>Final Structural Orthologs</th>
                </tr>
            </thead>
            <tbody>
EOF

    curr_date=\$(date)
    sed -i "s/SYSTEM_DATE_PLACEHOLDER/\$curr_date/g" summary_report.html

    awk -F'\t' 'NR>1 { print "<tr><td><b>" \$1 "</b></td><td>" \$2 "</td><td>" \$3 "</td><td>" \$4 "</td><td><b style=\\"color:#2b6cb0;\\">" \$5 "</b></td></tr>" }' summary_counts.tsv >> summary_report.html

    cat >> summary_report.html << 'EOF'
            </tbody>
        </table>
    </div>
EOF

    # Embed one Scatter Plot per query
    scatter_files=( orthology_evidence_plot_*.png )
    if [ \${#scatter_files[@]} -gt 0 ]; then
        cat >> summary_report.html << 'EOF'
    <div class="section">
        <h2>1. Structural & Sequence Embedding Evidence Mapping</h2>
EOF
        for scatter_file in "\${scatter_files[@]}"; do
            [ -f "\$scatter_file" ] || continue
            scatter_query=\$(basename "\$scatter_file" .png)
            scatter_query=\${scatter_query#orthology_evidence_plot_}

            scatter_b64=\$(base64 -w 0 "\$scatter_file" 2>/dev/null || base64 "\$scatter_file")
            cat >> summary_report.html << EOF
        <div class="img-card">
            <h3>Query: \$scatter_query</h3>
            <img src="data:image/png;base64,\${scatter_b64}" alt="Orthology Evidence Scatter Plot for \$scatter_query" />
        </div>
EOF
        done
        cat >> summary_report.html << 'EOF'
    </div>
EOF
    fi

    # Embed one Phylogenetic Tree per query (not per species-sample)
    tree_files=( annotated_tree_*.png )
    if [ \${#tree_files[@]} -gt 0 ]; then
        cat >> summary_report.html << 'EOF'
    <div class="section">
        <h2>2. Per-Query Phylogenetic Trees (IQ-TREE)</h2>
EOF
        for tree_file in "\${tree_files[@]}"; do
            [ -f "\$tree_file" ] || continue
            tree_query=\$(basename "\$tree_file" .png)
            tree_query=\${tree_query#annotated_tree_}

            tree_b64=\$(base64 -w 0 "\$tree_file" 2>/dev/null || base64 "\$tree_file")
            cat >> summary_report.html << EOF
        <div class="img-card">
            <h3>Query: \$tree_query</h3>
            <img src="data:image/png;base64,\${tree_b64}" alt="Annotated Tree for \$tree_query" />
        </div>
EOF
        done
        cat >> summary_report.html << 'EOF'
    </div>
EOF
    fi

    # Append per-sample detailed sections with full structural table integration
    cat >> summary_report.html << 'EOF'
    <div class="section">
        <h2>3. Detailed Sample Results & Integrated Evidence Classification</h2>
EOF

    if [ \${#hits_files[@]} -gt 0 ]; then
        for hits_file in "\${hits_files[@]}"; do
            [ -f "\$hits_file" ] || continue

            base_name=\$(basename "\$hits_file" _hits.list)
            base_name=\$(basename "\$base_name" .list)

            query=\$(echo "\$base_name" | sed 's/_vs_/ /' | cut -d' ' -f1)
            species=\$(echo "\$base_name" | sed 's/_vs_/ /' | cut -d' ' -f2)
            sample="\${query}_\${species}"

            ortho_file="\${query}_\${species}_orthologs.fa"
            final_file="\${query}_\${species}_filtered_orthologs.fa"
            domain_file="\${query}_\${species}_ortholog_domains.txt"

            hits_count=\$(wc -l < "\$hits_file" 2>/dev/null || echo "0")
            ortho_count=\$(if [ -f "\$ortho_file" ] && [ -s "\$ortho_file" ]; then seqkit seq --name --only-id "\$ortho_file" | wc -l 2>/dev/null || echo "0"; else echo "0"; fi)
            domain_count=\$(if [ -f "\$final_file" ] && [ -s "\$final_file" ]; then seqkit seq --name --only-id "\$final_file" | wc -l 2>/dev/null || echo "0"; else echo "0"; fi)

            cat >> summary_report.html << EOF
        <div class="sample-card">
            <h3>Sample: \$sample</h3>
            <p><strong>Homology Hits:</strong> <span class="count">\$hits_count</span> | 
               <strong>EggNOG Orthologs:</strong> <span class="count">\$ortho_count</span> | 
               <strong>Domain Screened:</strong> <span class="count">\$domain_count</span></p>

            <h4>Verified Orthologs with Structural Evidence & Domain Annotations:</h4>
EOF

            if [ -f "final_orthology_evidence_report.csv" ] && [ -s "final_orthology_evidence_report.csv" ]; then
                cat >> summary_report.html << EOF
            <table>
                <thead>
                    <tr>
                        <th>Target ID</th>
                        <th>ESM-2 Cosine Sim</th>
                        <th>Foldseek Score</th>
                        <th>Confidence Score</th>
                        <th>RBH</th>
                        <th>Annotated Domains</th>
                        <th>Classification</th>
                    </tr>
                </thead>
                <tbody>
EOF
                if [ -f "\$final_file" ] && [ -s "\$final_file" ]; then
                    seqkit seq --name --only-id "\$final_file" > "sample_targets.tmp"
                else
                    touch "sample_targets.tmp"
                fi

                awk -F',' -v domain_file="\$domain_file" '
                BEGIN {
                    if ((getline < domain_file) > 0) {
                        close(domain_file);
                        while ((getline line < domain_file) > 0) {
                            split(line, parts, /[ \t]+/);
                            if (parts[1] != "" && parts[2] != "") {
                                if (dom_map[parts[1]] != "") {
                                    dom_map[parts[1]] = dom_map[parts[1]] "," parts[2];
                                } else {
                                    dom_map[parts[1]] = parts[2];
                                }
                            }
                        }
                    }
                }
                NR==FNR {
                    valid_targets[\$1] = 1;
                    next;
                }
                FNR>1 {
                    target = \$1;
                    if (!(target in valid_targets)) next;

                    esm = \$2;
                    struct = \$3;
                    conf = \$4;
                    rbh = \$5;
                    class = \$7;

                    gsub(/\r/, "", class);

                    badge_class = "badge-other";
                    if (class ~ /SELF_HIT/) { badge_class = "badge-self"; }
                    else if (class ~ /PRIMARY/) { badge_class = "badge-primary"; }
                    else if (class ~ /SECONDARY/) { badge_class = "badge-secondary"; }

                    target_doms = (dom_map[target] != "") ? dom_map[target] : "No domains";

                    print "<tr>" \
                          "<td><b>" target "</b></td>" \
                          "<td>" esm "</td>" \
                          "<td>" struct "</td>" \
                          "<td><b>" conf "</b></td>" \
                          "<td>" rbh "</td>" \
                          "<td><span class=\\"domains\\">" target_doms "</span></td>" \
                          "<td><span class=\\"badge " badge_class "\\">" class "</span></td>" \
                          "</tr>";
                }' "sample_targets.tmp" final_orthology_evidence_report.csv >> summary_report.html

                rm -f "sample_targets.tmp"

                cat >> summary_report.html << EOF
                </tbody>
            </table>
EOF
            else
                cat >> summary_report.html << EOF
            <p><em>No structural evidence data recorded for this sample.</em></p>
EOF
            fi

            cat >> summary_report.html << EOF
        </div>
EOF
        done
    fi

    cat >> summary_report.html << 'EOF'
    </div>
</div>
</body>
</html>
EOF

    echo "[INFO] GENERATE_SUMMARY_REPORT process completed successfully."
    """
}
