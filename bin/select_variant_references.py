#!/usr/bin/env python
"""
Pick the reference genomes a sample's reads will be mapped to for variant calling.

Input is the MMseqs2 nucleotide top-hit table of the sample's viral contigs
(MMSEQS_NUCLEOTIDES `*_tophit_aln`, --format-mode 4, i.e. with a header line)
and, optionally, the parsed taxonomy table (`*_filtered_taxa.tsv`).

For every contig the best hit (max bits) is taken; hits that are too short or too
divergent are dropped. Hits are then aggregated per target: a reference supported by
more aligned contig sequence ranks higher. Output is one row per candidate reference,
best supported first, capped at --max-refs.

Outputs:
  <prefix>_ref_candidates.tsv   candidate references with support and taxonomy
  <prefix>_ref_candidates.txt   their ids, one per line (for extraction)
"""

import argparse
import sys

import pandas as pd

CANDIDATE_COLUMNS = [
    "reference", "ref_length", "n_contigs", "contig_aln_bp", "contig_fident",
    "best_contig", "taxid", "taxon_name", "lineage", "contig_believe",
]


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tophit", required=True, help="MMseqs2 *_tophit_aln table (with header)")
    parser.add_argument("--taxa", default=None, help="*_filtered_taxa.tsv from PARSE_MMSEQS_NUCLEOTIDE_RESULTS")
    parser.add_argument("--min-aln-len", type=int, default=200)
    parser.add_argument("--min-fident", type=float, default=0.7)
    parser.add_argument("--max-refs", type=int, default=50)
    parser.add_argument("--prefix", required=True)
    return parser.parse_args()


def write_outputs(df: pd.DataFrame, prefix: str) -> None:
    df = df.reindex(columns=CANDIDATE_COLUMNS)
    df.to_csv(f"{prefix}_ref_candidates.tsv", sep="\t", index=False)
    df["reference"].to_csv(f"{prefix}_ref_candidates.txt", index=False, header=False)


def load_taxa(path):
    if not path:
        return None
    try:
        taxa = pd.read_csv(path, sep="\t")
    except (pd.errors.EmptyDataError, FileNotFoundError):
        return None
    if "query" not in taxa.columns:
        return None
    keep = [c for c in ("query", "taxid", "name", "lineage", "believe") if c in taxa.columns]
    return taxa[keep].drop_duplicates("query")


def main():
    args = parse_arguments()

    try:
        aln = pd.read_csv(args.tophit, sep="\t")
    except pd.errors.EmptyDataError:
        aln = pd.DataFrame()

    required = {"query", "target", "tlen", "fident", "alnlen", "bits"}
    if aln.empty or not required.issubset(aln.columns):
        if not aln.empty:
            sys.stderr.write(f"WARNING: {args.tophit} lacks columns {sorted(required - set(aln.columns))}\n")
        write_outputs(pd.DataFrame(columns=CANDIDATE_COLUMNS), args.prefix)
        return

    # fident from MMseqs2 is a fraction (0..1); guard against a percent scale
    if aln["fident"].max() > 1.0:
        aln["fident"] = aln["fident"] / 100.0

    taxa = load_taxa(args.taxa)
    if taxa is not None:
        # only contigs that survived the nucleotide taxonomy parsing
        aln = aln[aln["query"].isin(taxa["query"])]

    # best hit per contig
    best = (aln.sort_values("bits", ascending=False)
               .drop_duplicates("query", keep="first"))
    best = best[(best["alnlen"] >= args.min_aln_len) & (best["fident"] >= args.min_fident)]

    if best.empty:
        write_outputs(pd.DataFrame(columns=CANDIDATE_COLUMNS), args.prefix)
        return

    if taxa is not None:
        best = best.merge(taxa, on="query", how="left")
    for col in ("taxid", "name", "lineage", "believe"):
        if col not in best.columns:
            best[col] = pd.NA

    best["aln_x_fident"] = best["alnlen"] * best["fident"]

    rows = []
    for ref, grp in best.groupby("target", sort=False):
        grp = grp.sort_values("alnlen", ascending=False)
        top = grp.iloc[0]
        aln_bp = int(grp["alnlen"].sum())
        rows.append({
            "reference": ref,
            "ref_length": int(top["tlen"]),
            "n_contigs": int(len(grp)),
            "contig_aln_bp": aln_bp,
            "contig_fident": round(float(grp["aln_x_fident"].sum() / aln_bp), 4),
            "best_contig": top["query"],
            "taxid": top["taxid"],
            "taxon_name": top["name"],
            "lineage": top["lineage"],
            "contig_believe": top["believe"],
        })

    cand = (pd.DataFrame(rows)
              .sort_values(["contig_aln_bp", "contig_fident"], ascending=False)
              .head(args.max_refs))
    write_outputs(cand, args.prefix)


if __name__ == "__main__":
    main()
