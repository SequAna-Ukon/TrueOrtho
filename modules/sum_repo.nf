process GENERATE_SUMMARY_REPORT {
    publishDir "${workflow.launchDir}/results", mode: 'copy'

    input:
    path homology_hits_lists
    path ortholog_fastas
    path filtered_ortholog_fastas
    path domain_files

    output:
    path "orthofind_report.html"
    path "summary_counts.tsv"

    script:
    """
    #!/bin/bash
    set -euo pipefail

    echo "[INFO] Generating OrthoFind report..."

    # Debug: list what files we received
    echo "[DEBUG] Domain files received:"
    for f in ${domain_files}; do
        echo "  - \$f"
    done

    echo "[DEBUG] Current directory contents:"
    ls -la

    # Create simple summary counts
    echo -e "Sample\tHomology_Hits\tOrthologs\tFinal_Orthologs" > summary_counts.tsv

    # Process each query-species combination
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

    # Generate HTML report
    cat > orthofind_report.html << EOF
<html>
<head>
    <title>OrthoFind Pipeline Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        table { border-collapse: collapse; margin: 20px 0; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
        .section { margin: 30px 0; }
        .count { font-weight: bold; color: #2E86AB; }
        .domains { font-size: 0.9em; color: #666; font-style: italic; }
    </style>
</head>
<body>
    <h1>OrthoFind Pipeline Report</h1>
    <p>Generated on: \$(date)</p>
    
    <div class="section">
        <h2>Summary Counts</h2>
        <table>
            <tr><th>Sample</th><th>Homology Hits</th><th>Orthologs</th><th>Final Orthologs</th></tr>
EOF

    # Add summary table rows using awk
    awk -F'\t' '
    NR>1 {
        print "            <tr><td>" \$1 "</td><td>" \$2 "</td><td>" \$3 "</td><td>" \$4 "</td></tr>"
    }
    ' summary_counts.tsv >> orthofind_report.html

    cat >> orthofind_report.html << EOF
        </table>
    </div>
    
    <div class="section">
        <h2>Detailed Results</h2>
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
        
        # Get counts
        hits_count=\$(wc -l < "\$hits_file" 2>/dev/null || echo "0")
        ortho_count=\$(if [ -f "\$ortho_file" ] && [ -s "\$ortho_file" ]; then seqkit seq --name --only-id "\$ortho_file" | wc -l 2>/dev/null || echo "0"; else echo "0"; fi)
        final_count=\$(if [ -f "\$final_file" ] && [ -s "\$final_file" ]; then seqkit seq --name --only-id "\$final_file" | wc -l 2>/dev/null || echo "0"; else echo "0"; fi)
        
        cat >> orthofind_report.html << EOF
        <div class="sample-section">
            <h3>\$sample</h3>
            
            <h4>Homology Hits (<span class="count">\$hits_count</span>):</h4>
            <ul>
EOF
        # List homology hits (first 10)
        if [ -s "\$hits_file" ]; then
            head -10 "\$hits_file" | while read hit; do
                echo "                <li>\$hit</li>" >> orthofind_report.html
            done
            if [ \$hits_count -gt 10 ]; then
                remaining=\$(( hits_count - 10 ))
                echo "                <li>... and \$remaining more</li>" >> orthofind_report.html
            fi
        else
            echo "                <li>No hits found</li>" >> orthofind_report.html
        fi
        
        cat >> orthofind_report.html << EOF
            </ul>
            
            <h4>Orthologs (<span class="count">\$ortho_count</span>):</h4>
            <ul>
EOF
        # List orthologs
        if [ -f "\$ortho_file" ] && [ -s "\$ortho_file" ]; then
            seqkit seq --name --only-id "\$ortho_file" | head -10 | while read ortho; do
                echo "                <li>\$ortho</li>" >> orthofind_report.html
            done
            if [ \$ortho_count -gt 10 ]; then
                remaining=\$(( ortho_count - 10 ))
                echo "                <li>... and \$remaining more</li>" >> orthofind_report.html
            fi
        else
            echo "                <li>No orthologs found</li>" >> orthofind_report.html
        fi
        
        cat >> orthofind_report.html << EOF
            </ul>
            
            <h4>Final Orthologs with Domains (<span class="count">\$final_count</span>):</h4>
            <ul>
EOF
        # List final orthologs with domain information
        if [ -f "\$final_file" ] && [ -s "\$final_file" ]; then
            seqkit seq --name --only-id "\$final_file" | head -10 | while read final_id; do
                # Look for domain information
                domains=""
                if [ -f "\$domain_file" ]; then
                    domains=\$(grep "^\$final_id" "\$domain_file" | awk '{print \$2}' | sort -u | tr '\n' ',' | sed 's/,\$//' 2>/dev/null || echo "")
                    echo "[DEBUG] For \$final_id, found domains: \$domains" >&2
                else
                    echo "[DEBUG] Domain file not found: \$domain_file" >&2
                fi
                if [ -n "\$domains" ]; then
                    echo "                <li>\$final_id <span class=\"domains\">(Domains: \$domains)</span></li>" >> orthofind_report.html
                else
                    echo "                <li>\$final_id <span class=\"domains\">(No domain information)</span></li>" >> orthofind_report.html
                fi
            done
            if [ \$final_count -gt 10 ]; then
                remaining=\$(( final_count - 10 ))
                echo "                <li>... and \$remaining more</li>" >> orthofind_report.html
            fi
        else
            echo "                <li>No final orthologs found</li>" >> orthofind_report.html
        fi
        
        cat >> orthofind_report.html << EOF
            </ul>
        </div>
        <hr>
EOF
    done

    cat >> orthofind_report.html << EOF
    </div>
</body>
</html>
EOF

    echo "[INFO] OrthoFind report generated: orthofind_report.html"
    echo "[INFO] Summary counts: summary_counts.tsv"
    """
}
