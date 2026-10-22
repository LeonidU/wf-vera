#!/usr/bin/env python

import argparse
import re

import numpy as np
import pandas as pd

from mmseqs_lca import read_lca


def what_taxa(x):
    # x is the 4-rank lineage, its first field the top-level taxon (see mmseqs_lca.py)
    if not isinstance(x, str):
        return "Undefinied"
    top = x.split(";")[0]
    if top in ("Eukaryota", "Bacteria", "Viruses"):
        return top
    return "Other"


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('-in_lca', '--input_lca',
                        type=str, required=True,
                        help="Input kraken report table")
    
    parser.add_argument('-out_p', '--out_pattern',
                        type=str, required=True, 
                        help="Patter for output files")
    
    parser.add_argument('-v', '--version',
                        action='version', version="0.0.1")

    return parser.parse_args()


if __name__ == '__main__':

    args = parse_arguments()

    df_lca = read_lca(args.input_lca)
    df_lca["contig_length"] = df_lca["query"].apply(lambda x: int(re.search(r"(?<=_length_)(\d{1,9}?)(?=_cov_)", x).group()))
    df_lca["taxa"] = df_lca["lineage"].apply(lambda x: what_taxa(x))
    df_lca["contig_name"] = df_lca["query"].apply(lambda x: x.split("_frame")[0])

    viruses = df_lca[df_lca["taxa"] == "Viruses"]
    unique_viral_contigs = viruses["contig_name"].unique()

    keeper = list()

    for contig in unique_viral_contigs:
        sub = df_lca[df_lca["contig_name"] == contig]
        unique_alignments = sub["lineage"].dropna().unique()
        if len(unique_alignments) == 1:
            if len(set(unique_alignments[0].split(";"))) > 2:
                conclusion = "clean lineage"
            else:
                conclusion = "unstable lineage"
            rdf = sub[sub["lineage"] == unique_alignments[0]]
            lineage = rdf.iloc[0]["lineage"]
            rank = rdf.iloc[0]["rank"]
            name = rdf.iloc[0]["name"]
            taxid = rdf.iloc[0]["taxid"]
            contig_length = rdf.iloc[0]["contig_length"]
        else:
            personal_df = pd.DataFrame([x.split(";") for x in unique_alignments],
                                       columns=["superkingdom","order","family","genus"])
            conclusion = ""
            if len(personal_df["superkingdom"].unique()) > 1:
                conclusion = "mixed lineage"
                top_index = personal_df[personal_df['superkingdom'] == 'Viruses'].apply(lambda x: len(set(x)), axis=1).idxmax()
            else:
                conclusion = "unstable lineage"
                top_index = personal_df.apply(lambda x: len(set(x)), axis=1).idxmax()
            lineage = ";".join(personal_df.loc[top_index].values)
            rdf = sub[sub["lineage"] == lineage]
            rank = rdf.iloc[0]["rank"]
            name = rdf.iloc[0]["name"]
            taxid = rdf.iloc[0]["taxid"]
            contig_length = rdf.iloc[0]["contig_length"]
        keeper.append([contig, contig_length, taxid, rank, name, lineage, conclusion])

    result = pd.DataFrame(keeper,
                          columns=["contig_name", "contig_length", "taxid", "rank", "name", "lineage", "conclusion"])
    result.to_csv(f"{args.out_pattern}_selected.tsv", sep="\t", index=False)
    result["contig_name"].to_csv(f"{args.out_pattern}_selected.txt", index=False, header=False)
