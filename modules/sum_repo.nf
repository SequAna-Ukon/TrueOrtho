process GENERATE_SUMMARY_REPORT {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path homology_hits_lists
    path ortholog_fastas
    path filtered_ortholog_fastas
    path domain_files

    output:
    path "summary_report.html"
    path "summary_counts.tsv"

    script:
    """
    #!/bin/bash
    set -euo pipefail

    echo "[INFO] Generating TrueOrtho report..."

    # Debug: list what files we received
    echo "[DEBUG] Domain files received:"
    for f in ${domain_files}; do
        echo "  - \$f"
    done

    # Create simple summary counts
    echo -e "Sample\tHomology_Hits\tOrthologs\tFinal_Orthologs" > summary_counts.tsv

    # Process each query-species combination for the TSV
    for hits_file in *_hits.list; do
        [ -f "\$hits_file" ] || continue
        
        base_name=\$(basename "\$hits_file" _hits.list)
        query=\$(echo "\$base_name" | sed 's/_vs_/ /' | cut -d' ' -f1)
        species=\$(echo "\$base_name" | sed 's/_vs_/ /' | cut -d' ' -f2)
        sample="\${query}_\${species}"
        
        # Count hits
        hits_count=\$(wc -l < "\$hits_file" 2>/dev/null || echo "0")
        
        # Count orthologs
        ortho_file="\${query}_\${species}_orthologs.fa"
        ortho_count=0
        if [ -f "\$ortho_file" ] && [ -s "\$ortho_file" ]; then
            ortho_count=\$(seqkit seq --name --only-id "\$ortho_file" | wc -l 2>/dev/null || echo "0")
        fi
        
        # Count final orthologs
        final_file="\${query}_\${species}_filtered_orthologs.fa"
        final_count=0
        if [ -f "\$final_file" ] && [ -s "\$final_file" ]; then
            final_count=\$(seqkit seq --name --only-id "\$final_file" | wc -l 2>/dev/null || echo "0")
        fi
        
        echo -e "\$sample\t\$hits_count\t\$ortho_count\t\$final_count" >> summary_counts.tsv
    done

    # Generate HTML report header
    cat > summary_report.html << EOF
<html>
<head>
    <title>TrueOrtho Pipeline Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; line-height: 1.6; color: #333; }
        table { border-collapse: collapse; width: 100%; max-width: 800px; margin: 20px 0; }
        th, td { border: 1px solid #ddd; padding: 12px; text-align: left; }
        th { background-color: #f8f9fa; font-weight: bold; }
        tr:nth-child(even) { background-color: #fcfcfc; }
        .section { margin: 40px 0; border-top: 2px solid #eee; padding-top: 20px; }
        .sample-section { background: #f9f9f9; padding: 15px; margin-bottom: 20px; border-radius: 5px; }
        .count { font-weight: bold; color: #2E86AB; }
        .domains { font-size: 0.85em; color: #d63384; background: #fff0f5; padding: 2px 5px; border-radius: 3px; }
        h1 { color: #2c3e50; }
        h2 { color: #34495e; border-bottom: 1px solid #eee; padding-bottom: 10px; }
    </style>
</head>
<body>
    <h1>TrueOrtho Pipeline Report</h1>
    <p><strong>Generated on:</strong> \$(date)</p>
    
    <div class="section">
        <h2>Summary Table</h2>
        <table>
            <thead>
                <tr><th>Sample</th><th>Homology Hits</th><th>Orthologs</th><th>Final Orthologs</th></tr>
            </thead>
            <tbody>
EOF

    # Add summary table rows using awk
    awk -F'\t' '
    NR>1 {
        print "                <tr><td>" \$1 "</td><td>" \$2 "</td><td>" \$3 "</td><td>" \$4 "</td></tr>"
    }
    ' summary_counts.tsv >> summary_report.html

    cat >> summary_report.html << EOF
            </tbody>
        </table>
    </div>
    
    <div class="section">
        <h2>Detailed Sample Results</h2>
EOF

    # Add detailed information for each sample
    for hits_file in *_hits.list; do
        [ -f "\$hits_file" ] || continue
        
        base_name=\$(basename "\$hits_file" _hits.list)
        query=\$(echo "\$base_name" | sed 's/_vs_/ /' | cut -d' ' -f1)
        species=\$(echo "\$base_name" | sed 's/_vs_/ /' | cut -d' ' -f2)
        sample="\${query}_\${species}"
        
        ortho_file="\${query}_\${species}_orthologs.fa"
        final_file="\${query}_\${species}_filtered_orthologs.fa"
        domain_file="\${query}_\${species}_ortholog_domains.txt"
        
        hits_count=\$(wc -l < "\$hits_file" 2>/dev/null || echo "0")
        ortho_count=\$(if [ -f "\$ortho_file" ] && [ -s "\$ortho_file" ]; then seqkit seq --name --only-id "\$ortho_file" | wc -l 2>/dev/null || echo "0"; else echo "0"; fi)
        final_count=\$(if [ -f "\$final_file" ] && [ -s "\$final_file" ]; then seqkit seq --name --only-id "\$final_file" | wc -l 2>/dev/null || echo "0"; else echo "0"; fi)
        
        cat >> summary_report.html << EOF
        <div class="sample-section">
            <h3>Sample: \$sample</h3>
            
            <p><strong>Homology Hits:</strong> <span class="count">\$hits_count</span></p>
            <p><strong>Orthologs Identified:</strong> <span class="count">\$ortho_count</span></p>
            <p><strong>Final Orthologs (Domain Filtered):</strong> <span class="count">\$final_count</span></p>

            <h4>Top 10 Final Orthologs & Domains:</h4>
            <ul>
EOF

        # List final orthologs with domain information
        if [ -f "\$final_file" ] && [ -s "\$final_file" ]; then
            seqkit seq --name --only-id "\$final_file" | head -10 | while read final_id; do
                domains=""
                if [ -f "\$domain_file" ]; then
                    # Extract domains from domain file for this ID
                    domains=\$(grep "^\$final_id" "\$domain_file" | awk '{print \$2}' | sort -u | tr '\n' ',' | sed 's/,\$//' 2>/dev/null || echo "")
                fi
                
                if [ -n "\$domains" ]; then
                    echo "                <li>\$final_id <span class=\"domains\">\$domains</span></li>" >> summary_report.html
                else
                    echo "                <li>\$final_id <span class=\"domains\">No domain data</span></li>" >> summary_report.html
                fi
            done
            if [ \$final_count -gt 10 ]; then
                remaining=\$(( final_count - 10 ))
                echo "                <li><em>... and \$remaining more</em></li>" >> summary_report.html
            fi
        else
            echo "                <li>No filtered orthologs passed domain validation.</li>" >> summary_report.html
        fi
        
        cat >> summary_report.html << EOF
            </ul>
        </div>
EOF
    done

    cat >> summary_report.html << EOF
    </div>
</body>
</html>
EOF

    echo "[INFO] TrueOrtho report generated: summary_report.html"
    echo "[INFO] Summary counts: summary_counts.tsv"
    """
}
