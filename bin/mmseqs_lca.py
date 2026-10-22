"""
Reading MMseqs2 `easy-taxonomy` LCA tables independently of the NCBI taxonomy version.

The pipeline runs easy-taxonomy with `--lca-ranks superkingdom,order,family,genus` and
`--tax-lineage 1`, so `*_lca.tsv` has the columns

    query  taxid  rank  name  lineage  full_lineage

`lineage` is the fixed 4-rank string the parsers work with. Since the March 2025 NCBI
taxonomy update the rank "superkingdom" no longer exists (Viruses became an
"acellular root", Bacteria / Archaea / Eukaryota a "domain"), and MMseqs2 prints
"unknown" for it, e.g. `unknown;Picornavirales;Picornaviridae;Enterovirus`.
`full_lineage` (`-_Viruses;-_Riboviria;...;g_Enterovirus;s_Enterovirus A`) still
names the top-level taxon, so the first field of `lineage` is restored from it.
Tables written without `--tax-lineage` (5 columns) are read as before.
"""

import pandas as pd

LCA_COLUMNS = ["query", "taxid", "rank", "name", "lineage", "full_lineage"]
TOP_LEVEL = ("Viruses", "Bacteria", "Archaea", "Eukaryota")


def top_level_taxon(full_lineage):
    """'Viruses' / 'Bacteria' / 'Archaea' / 'Eukaryota' from a full lineage, else None."""
    if not isinstance(full_lineage, str):
        return None
    for item in full_lineage.split(";"):
        # items are "<rank code>_<name>"; the rank code never contains "_"
        name = item.split("_", 1)[1] if "_" in item else item
        if name in TOP_LEVEL:
            return name
    return None


def fix_lineage(lineage, full_lineage):
    """Put the top-level taxon into the first field of a 4-rank lineage."""
    if not isinstance(lineage, str):
        return lineage
    top = top_level_taxon(full_lineage)
    if top is None:
        return lineage
    fields = lineage.split(";")
    fields[0] = top
    return ";".join(fields)


def read_lca(path):
    """LCA table with `lineage` normalised; `full_lineage` is NaN for 5-column tables."""
    try:
        df = pd.read_csv(path, sep="\t", header=None, names=LCA_COLUMNS,
                         dtype={"query": str}, keep_default_na=True)
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=LCA_COLUMNS)
    df["lineage"] = [fix_lineage(l, f) for l, f in zip(df["lineage"], df["full_lineage"])]
    return df
