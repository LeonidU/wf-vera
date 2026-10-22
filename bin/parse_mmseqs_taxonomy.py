#!/usr/bin/env python

import argparse

import numpy as np
import pandas as pd

from mmseqs_lca import read_lca


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('-in_lca', '--input_lca',
                        type=str, required=True,
                        help="Input kraken report table")
    
    parser.add_argument('-in_aln', '--input_aln',
                        type=str, required=True,
                        help="Input kraken report table")
    
    parser.add_argument('-hot_hit_size_match_cutoff', '--hot_hit_size_match_cutoff',
                        type=float, default=0.8,
                        help="Minimal cutoff for close hit size fraction")

    parser.add_argument('-hot_hit_length_cutoff', '--hot_hit_length_cutoff',
                        type=int, default=100,
                        help="Minimal length for close hit")
    
    parser.add_argument('-query_coverage_cutoff', '--query_coverage_cutoff',
                        type=float, default=0.5,
                        help="Minimal query coverage")
    
    parser.add_argument('-out_p', '--out_pattern',
                        type=str, required=True, 
                        help="Patter for output files")
    
    parser.add_argument('-v', '--version',
                        action='version', version="0.0.1")
        
    return parser.parse_args()


if __name__ == '__main__':

    args = parse_arguments()

    df_lca = read_lca(args.input_lca)
    df_lca.dropna(subset=['lineage'], axis=0, inplace=True)

    df_lca["believe"] = "NP"
    df_lca["record_len"] = 0
    df_lca["len_coverage"] = 0.0

    df_aln = pd.read_csv(args.input_aln, 
                         sep="\t"
                         )
    
    for f_idx, record in df_lca.iterrows():
        sub_df = df_aln[df_aln["query"] == record['query']]

        df_lca.loc[f_idx, 'record_len'] = sub_df.iloc[0]["qlen"]
        coverage = np.zeros(sub_df.iloc[0]["qlen"])

        at_least_one_close = False
        
        for idx, row in sub_df.iterrows():
            # если хотя бы один из хитов составляет как минимум 50% длины, то мы уже считаем его хорошей классификацией
            q2t = min(row['tlen'], row['qlen']) / max(row['tlen'], row['qlen'])
            if q2t >= args.hot_hit_size_match_cutoff and row['alnlen'] >= args.hot_hit_length_cutoff:
                at_least_one_close = True
            start = min(row['qstart'], row['qend'])
            end = max(row['qstart'], row['qend'])
            # в любом случае смотрим на хит
            coverage[start-1:end] = 1

        len_coverage = np.count_nonzero(coverage) / sub_df.iloc[0]["qlen"]
        if len_coverage >= args.query_coverage_cutoff and at_least_one_close:
            df_lca.loc[f_idx, 'believe'] = "complete close"
        elif len_coverage >= args.query_coverage_cutoff and not at_least_one_close:
            df_lca.loc[f_idx, 'believe'] = "fragment"
        elif len_coverage < args.query_coverage_cutoff and at_least_one_close:
            df_lca.loc[f_idx, 'believe'] = "complete far"
        else:
            df_lca.loc[f_idx, 'believe'] = "unknown"
        df_lca.loc[f_idx, "len_coverage"] = round(len_coverage, 3)

    df_lca.sort_values(by=['believe', 'record_len'], ascending=[True, False]).to_csv(args.out_pattern + "_filtered_taxa.tsv", index=False, sep="\t")
