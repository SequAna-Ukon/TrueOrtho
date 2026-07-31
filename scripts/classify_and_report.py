#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Phylo

def main():
    parser = argparse.ArgumentParser(description="Classify orthologs based on Foldseek, ESM-2, and IQ-TREE inputs.")
    parser.add_argument("--forward-m8", required=True, help="Path to forward Foldseek m8 file")
    parser.add_argument("--recip-m8", required=True, help="Path to reciprocal Foldseek m8 file")
    parser.add_argument("--esm-csv", required=True, help="Path to ESM-2 similarity CSV")
    parser.add_argument("--treefile", required=True, help="Path to IQ-TREE newick file")
    parser.add_argument("--out-report", default="final_orthology_evidence_report.csv", help="Output CSV path")
    parser.add_argument("--out-plot", default="orthology_evidence_plot.png", help="Output scatter plot path")
    parser.add_argument("--out-tree-plot", default="annotated_orthology_tree.png", help="Output annotated tree plot path")
    args = parser.parse_args()

    # 1. Load Forward Foldseek Hits
    if os.path.exists(args.forward_m8) and os.path.getsize(args.forward_m8) > 0:
        df_fs = pd.read_csv(args.forward_m8, sep="\t", names=["query", "target", "pident", "alnlen", "bits", "evalue"])
        df_fs["target"] = df_fs["target"].str.replace(".pdb", "", regex=False)
        max_bits = df_fs["bits"].max() if df_fs["bits"].max() > 0 else 1.0
        df_fs["struct_score"] = (df_fs["bits"] / max_bits).round(4)
    else:
        df_fs = pd.DataFrame(columns=["target", "struct_score"])

    # 2. Load Reciprocal Foldseek Best Hits
    rbh_targets = set()
    if os.path.exists(args.recip_m8) and os.path.getsize(args.recip_m8) > 0:
        df_recip = pd.read_csv(args.recip_m8, sep="\t", names=["query", "target", "pident", "alnlen", "bits", "evalue"])
        top_recip = df_recip.sort_values(by=["query", "bits"], ascending=[True, False]).groupby("query").first().reset_index()
        rbh_targets = set(top_recip["query"])

    # 3. Load ESM-2 Cosine Similarities
    df_esm = pd.read_csv(args.esm_csv) if os.path.exists(args.esm_csv) else pd.DataFrame(columns=["target", "esm2_cosine_sim"])

    # 4. Merge Metrics & Calculate Score
    df_final = pd.merge(df_esm, df_fs[["target", "struct_score"]], on="target", how="left").fillna(0)
    df_final["orthology_confidence_score"] = ((df_final["struct_score"] * 0.40) + (df_final["esm2_cosine_sim"] * 0.40) + 0.20).round(4)
    df_final["is_rbh"] = df_final["target"].apply(lambda x: "YES" if x in rbh_targets else "NO")

    # Exclude Self-Hits
    candidates = df_final[df_final["orthology_confidence_score"] < 0.999].copy()

    def classify_ortholog(row):
        if row["is_rbh"] == "YES" and row["esm2_cosine_sim"] >= 0.95 and row["struct_score"] >= 0.50:
            return "PRIMARY_TRUE_ORTHOLOG"
        elif row["esm2_cosine_sim"] >= 0.95 and row["struct_score"] >= 0.45:
            return "CO_ORTHOLOG_PARALOG"
        else:
            return "DISTANT_HOMOLOG"

    candidates["classification"] = candidates.apply(classify_ortholog, axis=1)
    candidates = candidates.sort_values(by="orthology_confidence_score", ascending=False)
    candidates.to_csv(args.out_report, index=False)

    # 5. Generate Evidence Scatter Plot
    plt.figure(figsize=(9, 6))
    scatter = plt.scatter(candidates["esm2_cosine_sim"], candidates["struct_score"], c=candidates["orthology_confidence_score"], cmap="viridis", s=140, edgecolors="k")
    plt.colorbar(scatter, label="Orthology Confidence Score")

    for _, r in candidates.iterrows():
        label = f"{r['target']} ({'RBH' if r['is_rbh'] == 'YES' else 'No-RBH'})"
        plt.annotate(label, (r["esm2_cosine_sim"], r["struct_score"]), textcoords="offset points", xytext=(6, 6), fontsize=8)

    plt.xlabel("ESM-2 Cosine Similarity")
    plt.ylabel("Foldseek ProstT5 Structure Score")
    plt.title("Orthology Evidence Mapping (RBH + Sequence + Structure)")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(args.out_plot, dpi=300)
    plt.close()

    # 6. Generate Annotated Tree Plot
    if os.path.exists(args.treefile):
        tree = Phylo.read(args.treefile, "newick")
        df_idx = candidates.set_index("target")

        for leaf in tree.get_terminals():
            name = leaf.name.split()[0]
            if name in df_idx.index:
                class_type = df_idx.loc[name, "classification"]
                score = df_idx.loc[name, "orthology_confidence_score"]
                leaf.name = f"{name} [{class_type}|Score:{score}]"

        fig = plt.figure(figsize=(12, 7), dpi=300)
        axes = fig.add_subplot(1, 1, 1)
        Phylo.draw(tree, axes=axes, do_show=False)
        plt.title("IQ-TREE Annotated with Automated Orthology Classification")
        plt.tight_layout()
        plt.savefig(args.out_tree_plot, dpi=300)
        plt.close()

if __name__ == "__main__":
    main()
