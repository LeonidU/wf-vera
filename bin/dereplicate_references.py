#!/usr/bin/env python
"""
Collapse near-identical candidate references to one per cluster.

Reads are mapped competitively to all references of a sample at once. If two strains
of the same virus are both in the reference set, reads split between them, get MAPQ 0
and are lost to variant calling. So candidates are clustered first (MMseqs2
easy-cluster, `*_cluster.tsv`: representative<TAB>member) and from every cluster we
keep the member with the strongest support from the sample's own contigs (not the
MMseqs2 representative, which is chosen without that knowledge).

Without a cluster table every candidate is kept.

Outputs:
  <prefix>_references.tsv   final references (+ `collapsed_refs`: the near-identical
                            candidates this reference stands for)
  <prefix>_references.txt   their ids
"""

import argparse

import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--clusters", default=None)
    parser.add_argument("--present", default=None,
                        help="ids actually extracted from the reference source; others are dropped")
    parser.add_argument("--prefix", required=True)
    return parser.parse_args()


def main():
    args = parse_arguments()
    cand = pd.read_csv(args.candidates, sep="\t", dtype={"reference": str})
    cand["reference"] = cand["reference"].astype(str)
    if args.present:
        with open(args.present) as fh:
            present = {line.strip() for line in fh if line.strip()}
        cand = cand[cand["reference"].isin(present)]

    if args.clusters:
        clu = pd.read_csv(args.clusters, sep="\t", header=None,
                          names=["representative", "reference"], dtype=str)
        cand = cand.merge(clu, on="reference", how="left")
    else:
        cand["representative"] = pd.NA
    # candidates missing from the cluster table form their own cluster
    cand["representative"] = cand["representative"].fillna(cand["reference"])

    def numeric(col):
        if col not in cand.columns:          # fasta_all mode has no contig support
            return pd.Series(0, index=cand.index)
        return pd.to_numeric(cand[col], errors="coerce").fillna(0)

    cand["_support"] = numeric("contig_aln_bp")
    cand["_length"] = numeric("ref_length")

    keep = []
    for _, grp in cand.groupby("representative", sort=False):
        grp = grp.sort_values(["_support", "_length"], ascending=False)
        top = grp.iloc[0].copy()
        others = grp["reference"].iloc[1:].tolist()
        top["collapsed_refs"] = ",".join(others)
        if len(grp) > 1 and "n_contigs" in grp.columns:
            # the kept reference now also stands for the contigs of the collapsed ones
            top["n_contigs"] = int(pd.to_numeric(grp["n_contigs"], errors="coerce").fillna(0).sum())
            top["contig_aln_bp"] = int(grp["_support"].sum())
        keep.append(top)

    final = (pd.DataFrame(keep)
               .sort_values(["_support", "_length"], ascending=False)
               .drop(columns=["representative", "_support", "_length"]))
    final.to_csv(f"{args.prefix}_references.tsv", sep="\t", index=False)
    final["reference"].to_csv(f"{args.prefix}_references.txt", index=False, header=False)


if __name__ == "__main__":
    main()
