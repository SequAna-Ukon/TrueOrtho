#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Phylo

def extract_qid(filepath):
    """Extract the query id from a filename of the form '<qid>_vs_<sp>_...' """
    base = os.path.basename(filepath)
    if "_vs_" in base:
        return base.split("_vs_")[0]
    return base

def load_and_concat_csvs(file_list, sep=",", names=None):
    dfs = []
    for f in file_list:
        if os.path.exists(f) and os.path.getsize(f) > 0:
            try:
                if names:
                    df = pd.read_csv(f, sep=sep, names=names)
                else:
                    df = pd.read_csv(f, sep=sep)
                if not df.empty:
                    df["qid"] = extract_qid(f)
                    dfs.append(df)
            except Exception as e:
                print(f"[WARNING] Could not read {f}: {e}")
    if dfs:
        return pd.concat(dfs, ignore_index=True)
    return pd.DataFrame()

def main():
    parser = argparse.ArgumentParser(description="Classify orthologs across all samples.")
    parser.add_argument("--forward-m8", nargs="+", required=True, help="Path(s) to forward Foldseek m8 file(s)")
    parser.add_argument("--recip-m8", nargs="+", required=True, help="Path(s) to reciprocal Foldseek m8 file(s)")
    parser.add_argument("--esm-csv", nargs="+", required=True, help="Path(s) to ESM-2 similarity CSV(s)")
    parser.add_argument("--treefile", nargs="+", required=True, help="Path(s) to IQ-TREE newick file(s), one per query")
    parser.add_argument("--out-report", default="final_orthology_evidence_report.csv", help="Output CSV path")
    parser.add_argument("--out-plot-prefix", default="orthology_evidence_plot", help="Filename prefix for per-query scatter plots")
    parser.add_argument("--out-tree-prefix", default="annotated_tree", help="Filename prefix for per-query annotated tree plots")
    args = parser.parse_args()

    # 1. Load Forward Foldseek Hits (tagged with qid from filename)
    df_fs = load_and_concat_csvs(args.forward_m8, sep="\t", names=["query", "target", "pident", "alnlen", "bits", "evalue"])
    if not df_fs.empty:
        df_fs["target"] = df_fs["target"].astype(str).str.split().str[0].str.replace(".pdb", "", regex=False)
        max_bits = df_fs["bits"].max() if df_fs["bits"].max() > 0 else 1.0
        df_fs["struct_score"] = (df_fs["bits"] / max_bits).round(4)
    else:
        df_fs = pd.DataFrame(columns=["query", "target", "struct_score", "qid"])

    # 2. Load Reciprocal Foldseek Hits
    rbh_targets = set()
    df_recip = load_and_concat_csvs(args.recip_m8, sep="\t", names=["query", "target", "pident", "alnlen", "bits", "evalue"])
    if not df_recip.empty:
        df_recip["query"] = df_recip["query"].astype(str).str.split().str[0]
        top_recip = df_recip.sort_values(by=["query", "bits"], ascending=[True, False]).groupby("query").first().reset_index()
        rbh_targets = set(top_recip["query"])

    # 3. Load ESM-2 Cosine Similarities (tagged with qid from filename)
    df_esm = load_and_concat_csvs(args.esm_csv, sep=",")
    if not df_esm.empty and "target" in df_esm.columns:
        df_esm["target"] = df_esm["target"].astype(str).str.split().str[0]
    else:
        df_esm = pd.DataFrame(columns=["target", "esm2_cosine_sim", "qid"])

    # 4. Merge Metrics (join on target AND qid, so per-query separation is preserved)
    self_hit_targets = set(df_fs.loc[df_fs["query"] == df_fs["target"], "target"].astype(str)) if not df_fs.empty else set()
    df_final = pd.merge(df_esm, df_fs[["target", "qid", "struct_score"]], on=["target", "qid"], how="outer").fillna(0)

    if not df_final.empty:
        df_final["target"] = df_final["target"].astype(str).str.split().str[0]
        df_final = df_final.drop_duplicates(subset=["qid", "target"])

        df_final["struct_score"] = pd.to_numeric(df_final["struct_score"], errors="coerce").fillna(0.0)
        df_final["esm2_cosine_sim"] = pd.to_numeric(df_final["esm2_cosine_sim"], errors="coerce").fillna(0.0)

        df_final["orthology_confidence_score"] = ((df_final["struct_score"] * 0.40) + (df_final["esm2_cosine_sim"] * 0.40) + 0.20).round(4)
        df_final["is_rbh"] = df_final["target"].apply(lambda x: "YES" if str(x) in rbh_targets else "NO")

        if self_hit_targets:
            candidates = df_final[~df_final["target"].isin(self_hit_targets)].copy()
        else:
            candidates = df_final.copy()

        def classify_ortholog(row):
            if row["is_rbh"] == "YES" and row["esm2_cosine_sim"] >= 0.90 and row["struct_score"] >= 0.40:
                return "PRIMARY_TRUE_ORTHOLOG"
            elif row["esm2_cosine_sim"] >= 0.85 and row["struct_score"] >= 0.35:
                return "SECONDARY_ORTHOLOG"
            elif row["esm2_cosine_sim"] >= 0.70 or row["struct_score"] >= 0.30:
                return "CO_ORTHOLOG_PARALOG"
            else:
                return "DISTANT_HOMOLOG"

        if not candidates.empty:
            candidates["classification"] = candidates.apply(classify_ortholog, axis=1)
            candidates = candidates.sort_values(by="orthology_confidence_score", ascending=False)
        else:
            candidates["classification"] = []
    else:
        candidates = pd.DataFrame(columns=["target", "qid", "esm2_cosine_sim", "struct_score", "orthology_confidence_score", "is_rbh", "classification"])

    candidates.drop(columns=["qid"]).to_csv(args.out_report, index=False)

    # 5. Per-Query Scatter Plots
    if not candidates.empty:
        unique_qids = sorted(candidates["qid"].unique())
    else:
        unique_qids = []

    for qid in unique_qids:
        sub = candidates[candidates["qid"] == qid]
        out_path = f"{args.out_plot_prefix}_{qid}.png"

        plt.figure(figsize=(9, 6))
        if not sub.empty:
            scatter = plt.scatter(sub["esm2_cosine_sim"], sub["struct_score"], c=sub["orthology_confidence_score"], cmap="viridis", s=140, edgecolors="k")
            plt.colorbar(scatter, label="Orthology Confidence Score")

            for _, r in sub.iterrows():
                label = f"{r['target']} ({'RBH' if r['is_rbh'] == 'YES' else 'No-RBH'})"
                plt.annotate(label, (r["esm2_cosine_sim"], r["struct_score"]), textcoords="offset points", xytext=(6, 6), fontsize=8)

        plt.xlabel("ESM-2 Cosine Similarity")
        plt.ylabel("Foldseek ProstT5 Structure Score")
        plt.title(f"Orthology Evidence Mapping — {qid}")
        plt.grid(True, linestyle="--", alpha=0.5)
        plt.tight_layout()
        plt.savefig(out_path, dpi=300)
        plt.close()

    # 6. Per-Query Tree Plots
    df_idx_full = candidates.set_index("target") if not candidates.empty else pd.DataFrame()

    for tf in args.treefile:
        query_name = os.path.basename(tf)
        if query_name.endswith("_iqtree.treefile"):
            query_name = query_name[: -len("_iqtree.treefile")]

        out_path = f"{args.out_tree_prefix}_{query_name}.png"

        tree_obj = None
        tree_ok = os.path.exists(tf) and os.path.getsize(tf) > 0
        if tree_ok:
            try:
                tree_obj = Phylo.read(tf, "newick")
            except Exception as e:
                print(f"[WARNING] Could not parse treefile {tf}: {e}")
                tree_ok = False

        if tree_ok and tree_obj is not None:
            for leaf in tree_obj.get_terminals():
                name = leaf.name.split()[0]
                if not df_idx_full.empty and name in df_idx_full.index:
                    class_type = df_idx_full.loc[name, "classification"]
                    score = df_idx_full.loc[name, "orthology_confidence_score"]
                    leaf.name = f"{name} [{class_type}|Score:{score}]"

            fig = plt.figure(figsize=(12, 7), dpi=300)
            axes = fig.add_subplot(1, 1, 1)
            Phylo.draw(tree_obj, axes=axes, do_show=False)
            plt.title(f"IQ-TREE Annotated Orthology Classification — {query_name}")
            plt.tight_layout()
            plt.savefig(out_path, dpi=300)
            plt.close()
        else:
            fig, ax = plt.subplots(figsize=(6, 2))
            ax.text(0.5, 0.5, f"Phylogenetic tree omitted for {query_name} (<3 sequences)", ha="center", va="center")
            ax.axis("off")
            plt.savefig(out_path, dpi=300)
            plt.close()

if __name__ == "__main__":
    main()
