# TrueOrtho 2.0
TrueOrtho 2.0 is an expanded Nextflow pipeline for automated, structural- and embedding-aware orthology identification across multiple species. In addition to traditional sequence homology and domain architecture filtering, version 2.0 integrates deep-learning structural embeddings (ESM-2), 3D structural alignment scores (Foldseek), and automated phylogenetic trees (IQ-TREE) into a unified evidence matrix.

## 🌟 What's New in Version 2.0?
- **Structural & Embedding Evidence:** Calculates ESM-2 cosine similarity and Foldseek structural scores for candidate orthologs.
- **Combined Confidence Scoring:** Assigns an integrated confidence metric and reciprocal best hit (RBH) verification.
- **Phylogenetic Context:** Constructs and annotates phylogenetic trees (IQ-TREE) automatically within the workflow.
- **Enhanced 4-Stage Filtering Funnel:** Tracks target progression through Homology $\rightarrow$ Reciprocal/EggNOG $\rightarrow$ Domain Screen $\rightarrow$ Final Structural Orthologs.
- **Upgraded HTML Reporting:** Interactive HTML reports with embedded scatter plots, trees, and detailed domain annotations.

## Workflow Overview

[Query + Target DB]
        │
        ▼
 1. HOMOLOGY_SEARCH (jackhmmer / DIAMOND)
        │
        ▼
 2. ORTHOLOG_ASSIGN (eggNOG-mapper / KOG / RBH)
        │
        ▼
 3. DOMAIN_SCAN (HMMER hmmscan / Pfam / SMART)
        │
        ▼
 4. STRUCTURAL & EMBEDDING EVALUATION (ESM-2 + Foldseek)
        │
        ▼
 5. PHYLOGENETIC TREE (IQ-TREE)
        │
        ▼
 6. GENERATE_SUMMARY_REPORT (Interactive HTML + TSV Funnel)

Quick StartPrerequisites

Nextflow ($\ge$ 20.07.1)Conda / Mamba, Docker, or Singularity for dependency management.Basic UsageBash

nextflow run main.nf \
  --input input.csv \
  --eggnog_db /path/to/eggnog_db \
  --domain_db /path/to/pfam_smart.hmm \
  --threads 12

📋 Input CSV FormatCreate an input CSV defining your query-database pairs:Code snippetquery,database,kog_id,target_domain
````bash
../GS_q.fsa,../Smic.fasta,,
../GOGAT_q.fsa,../Smic.fasta,,
../GDH_q.fsa,../Smic.fasta,,
../NIR_q.fsa,../Smic.fasta,,
../NR_q.fsa,../Smic.fasta,,"PF00069,PF00070"
````
Column Specificationsquery: Path to query protein sequence FASTA file (Required).database: Path to target species proteome FASTA file (Required).kog_id: Optional KOG/COG ID for explicit target assignment.target_domain: Optional comma-separated domain IDs (e.g., Pfam/SMART accessions) for domain conservation filtering. Enclose multiple domains in quotes ("Dom1,Dom2").⚙️ Pipeline ParametersRequired ParametersParameterDescription--inputPath to input CSV file containing query-database pairs.Database & Resource ParametersParameterDefaultDescription--eggnog_dbAuto-download (v5.0.2)Directory path to pre-downloaded eggNOG database.--domain_dbAuto-download (Pfam/SMART)Path to HMM database file for domain scanning.--threads10Total CPU threads allocated for parallel processing steps.--outdir./resultsDirectory where output files and summary reports are saved.

📁 Output StructurePlaintextresults/
├── homology_search/
│   └── {query}_{species}/
│       ├── {query}_vs_{species}_hits.fa
│       └── {query}_vs_{species}_hits.list
├── ortholog_assign/
│   └── {query}_{species}/
│       ├── {query}_{species}_orthologs.fa
│       ├── {query}_{species}_hits.txt
│       └── {query}_{species}_kog_info.txt
├── domain_scan/
│   └── {query}_{species}/
│       ├── {query}_{species}_filtered_orthologs.fa
│       ├── {query}_{species}_domains.tblout
│       ├── {query}_{species}_ortholog_domains.txt
│       └── {query}_{species}_target_domains.txt
├── structural_analysis/
│   ├── esm_embeddings.csv
│   └── foldseek_alignments.tsv
├── summary_counts.tsv
├── final_orthology_evidence_report.csv
├── orthology_evidence_plot.png
├── annotated_orthology_tree.png
└── summary_report.html
📊 Output File Descriptions1. summary_counts.tsvTracks the exact progression of candidate targets across all 4 pipeline stages:Code snippetSample	Homology_Hits	EggNOG_Orthologs	Domain_Screened	Final_Structural_Orthologs
NR_new_q_Smic	154	12	4	3

2. final_orthology_evidence_report.csvDetailed candidate scoring matrix including deep learning similarity, structural score, and final decision:Target IDESM2_Cosine_SimFoldseek_ScoreConfidence_ScoreRBHClassificationOLP86534.10.97940.60160.8324YESPRIMARY_TRUE_ORTHOLOGOLP86515.10.97780.55650.8137YESPRIMARY_TRUE_ORTHOLOGOLP89482.10.97740.50970.7948YESPRIMARY_TRUE_ORTHOLOG3. Interactive summary_report.htmlA single, self-contained HTML file containing:Filtering Funnel Summary Table: Direct comparison of candidate reduction counts.Evidence Scatter Plot: Base64-embedded visualization comparing ESM-2 Cosine Similarity vs. Foldseek Structural Scores.Annotated Phylogenetic Tree: High-resolution IQ-TREE output tree diagram.Per-Sample Evidence Cards: Complete tables detailing sequence scores, domain assignments (e.g., PF00069), RBH status, and classification badges.

🛠️ Pipeline StagesHOMOLOGY_SEARCH (jackhmmer / DIAMOND)Identifies candidate homologous sequences across target species databases with automatic header sanitization.ORTHOLOG_ASSIGN (eggNOG-mapper)Maps hits against eggNOG orthology groups and extracts reciprocal best hits (RBH).DOMAIN_SCAN (HMMER hmmscan)Scans candidate sequences against Pfam/SMART HMM profiles to enforce target domain architecture preservation.STRUCTURAL_EVALUATION (ESM-2 + Foldseek)Computes sequence embedding similarities via ESM-2 and 3D structural alignment metrics via Foldseek.PHYLOGENETIC_INFERENCE (IQ-TREE)Aligns domain-verified candidates and builds annotated phylogenetic trees to confirm orthology clades.GENERATE_SUMMARY_REPORTSynthesizes all data streams into summary_counts.tsv and summary_report.html.

💡 Key Features & TipsNextflow DSL2 Support: Fully modular process design built for scalability.Resume Capability: Run with -resume to restart from interrupted steps without recalculating complete stages.Performance Optimization: Provide local paths for --eggnog_db and --domain_db to bypass download overhead on large clusters.

📜 Citation & LicenseLicenseThis pipeline is open-source and released under the MIT License. Third-party dependencies maintain their respective open-source licenses (located under /third_party_licenses/).CitationIf you use TrueOrtho in your research, please cite:Sharaf, A., & Voolstra, C. R. (2025). SequAna-Ukon/TrueOrtho (Version 2.0.0). Zenodo. https://doi.org/10.5281/zenodo.17867442Acknowledgments & SupportSupported by the Sequencing Analysis (SequAna) Core Facility at the University of Konstanz (biologie.uni-konstanz.de/sequana).

For inquiries, issues, or feature requests, contact abdoallah.sharaf@uni-konstanz.de.
