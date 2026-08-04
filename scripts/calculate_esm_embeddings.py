#!/usr/bin/env python3
import argparse
import torch
import pandas as pd
from Bio import SeqIO
from transformers import AutoTokenizer, EsmModel
from sklearn.metrics.pairwise import cosine_similarity

def main():
    parser = argparse.ArgumentParser(description="Calculate ESM-2 cosine similarity between query and candidate sequences.")
    parser.add_argument("--query", required=True, help="Path to query FASTA")
    parser.add_argument("--candidates", required=True, help="Path to candidates FASTA")
    parser.add_argument("--output", default="esm_sim.csv", help="Output CSV path")
    parser.add_argument("--model", default="facebook/esm2_t6_8M_UR50D", help="ESM model name")
    args = parser.parse_args()

    device = "cuda" if torch.cuda.is_available() else "cpu"
    tokenizer = AutoTokenizer.from_pretrained(args.model)
    model = EsmModel.from_pretrained(args.model).to(device)
    model.eval()

    def get_embedding(seq):
        inputs = tokenizer(seq, return_tensors="pt", padding=True, truncation=True, max_length=1024).to(device)
        with torch.no_grad():
            outputs = model(**inputs)
        return outputs.last_hidden_state.mean(dim=1).squeeze().cpu().numpy()

    query_rec = next(SeqIO.parse(args.query, "fasta"))
    query_emb = get_embedding(str(query_rec.seq)).reshape(1, -1)

    results = []
    for record in SeqIO.parse(args.candidates, "fasta"):
        target_emb = get_embedding(str(record.seq)).reshape(1, -1)
        sim = cosine_similarity(query_emb, target_emb)[0][0]
        results.append({
            "target": record.id.split()[0],
            "esm2_cosine_sim": round(float(sim), 4)
        })

    pd.DataFrame(results).to_csv(args.output, index=False)

if __name__ == "__main__":
    main()
