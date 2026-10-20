#!/usr/bin/env python
"""
Turn raw variant calls (LoFreq VCF or iVar TSV) of one sample into the tables that
feed the LLM report.

`caller_qual` is LoFreq's QUAL (phred) or iVar's p-value, depending on the caller.

Every call is normalised to VCF-style alleles, classified (SNV / MNV / INS / DEL,
consensus vs minor) and filtered with the same rules for both callers, so switching
the caller does not change the report format.

Outputs (<prefix> = sample id):
  <prefix>_variants.tsv          every call with its filter status (PASS or reasons)
  <prefix>_variant_summary.tsv   one row per reference: coverage, callable bases,
                                 counts of consensus / minor variants, divergence
  <prefix>_variants.json         compact LLM payload: per-reference summary plus the
                                 strongest PASS variants (capped per reference)
"""

import argparse
import json
import math
import sys
from collections import defaultdict

import pandas as pd

VARIANT_COLUMNS = [
    "sample", "reference", "taxon_name", "pos", "ref", "alt", "type", "class",
    "af", "depth", "alt_depth", "ref_fwd", "ref_rev", "alt_fwd", "alt_rev",
    "homopolymer", "caller_qual", "caller_filter", "filter",
]


# ─────────────────────────────────────────────────────────────────────────────
# inputs
# ─────────────────────────────────────────────────────────────────────────────
def parse_arguments():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample", required=True)
    p.add_argument("--caller", required=True, choices=["lofreq", "ivar"])
    p.add_argument("--platform", default="illumina", choices=["illumina", "iontorrent"])
    p.add_argument("--calls", required=True, help="LoFreq VCF or iVar variants TSV")
    p.add_argument("--references", required=True, help="*_references.tsv")
    p.add_argument("--ref-fasta", required=True)
    p.add_argument("--coverage", required=True, help="samtools coverage output")
    p.add_argument("--callable", required=True, help="reference<TAB>bases with depth >= min-depth")
    p.add_argument("--min-depth", type=int, default=10)
    p.add_argument("--min-af", type=float, default=0.03)
    p.add_argument("--min-alt-reads", type=int, default=3)
    p.add_argument("--consensus-af", type=float, default=0.5)
    p.add_argument("--min-callable-frac", type=float, default=0.5,
                   help="below this fraction of callable bases a reference is reported as low_coverage")
    p.add_argument("--strand-min-frac", type=float, default=0.1)
    p.add_argument("--iontorrent-max-hp", type=int, default=4,
                   help="Ion Torrent indels in homopolymers of at least this length are filtered")
    p.add_argument("--max-report-per-ref", type=int, default=100)
    p.add_argument("--prefix", required=True)
    return p.parse_args()


def read_fasta(path):
    seqs, desc, name, chunks = {}, {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                header = line[1:].strip()
                name = header.split()[0] if header else ""
                desc[name] = header[len(name):].strip()
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks).upper()
    return seqs, desc


def read_lofreq_vcf(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            chrom, pos, _id, ref, alt, qual, filt, info = f[:8]
            tags = {}
            for item in info.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    tags[k] = v
                else:
                    tags[item] = True
            dp4 = [int(x) for x in tags.get("DP4", "0,0,0,0").split(",")]
            for a in alt.split(","):
                rows.append({
                    "reference": chrom, "pos": int(pos), "ref": ref.upper(), "alt": a.upper(),
                    "af": float(tags.get("AF", "nan")), "depth": int(tags.get("DP", sum(dp4))),
                    "ref_fwd": dp4[0], "ref_rev": dp4[1], "alt_fwd": dp4[2], "alt_rev": dp4[3],
                    "alt_depth": dp4[2] + dp4[3],
                    "caller_qual": None if qual == "." else float(qual),
                    "caller_filter": filt,
                })
    return rows


def read_ivar_tsv(path, seqs):
    try:
        df = pd.read_csv(path, sep="\t", dtype={"REGION": str})
    except pd.errors.EmptyDataError:
        return []
    if df.empty:
        return []
    # iVar prints one row per overlapping GFF feature; one call per allele is enough
    df = df.drop_duplicates(["REGION", "POS", "REF", "ALT"])
    rows = []
    for r in df.itertuples(index=False):
        ref_base = str(r.REF).upper()
        alt = str(r.ALT).upper()
        is_indel = alt[:1] in ("+", "-")
        # iVar indels: ALT is "+SEQ" / "-SEQ" anchored on the base at POS (VCF convention)
        if alt.startswith("+"):
            ref, alt = ref_base, ref_base + alt[1:]
        elif alt.startswith("-"):
            ref, alt = ref_base + alt[1:], ref_base
        else:
            ref = ref_base
        alt_dp = int(r.ALT_DP)
        alt_rv = int(r.ALT_RV)
        ref_dp = int(r.REF_DP)
        ref_rv = int(r.REF_RV)
        # iVar does not track strands for indels (ALT_RV is always 0): leave them unknown
        strands = (None, None, None, None) if is_indel else \
                  (ref_dp - ref_rv, ref_rv, alt_dp - alt_rv, alt_rv)
        rows.append({
            "reference": r.REGION, "pos": int(r.POS), "ref": ref, "alt": alt,
            "af": float(r.ALT_FREQ), "depth": int(r.TOTAL_DP),
            "ref_fwd": strands[0], "ref_rev": strands[1],
            "alt_fwd": strands[2], "alt_rev": strands[3], "alt_depth": alt_dp,
            "caller_qual": float(r.PVAL),
            "caller_filter": "PASS" if str(r.PASS).upper() == "TRUE" else "ivar_pval",
        })
    return rows


# ─────────────────────────────────────────────────────────────────────────────
# classification and filtering
# ─────────────────────────────────────────────────────────────────────────────
def variant_type(ref, alt):
    if len(ref) == len(alt):
        return "SNV" if len(ref) == 1 else "MNV"
    return "INS" if len(alt) > len(ref) else "DEL"


def run_length(seq, i):
    """Length of the run of identical bases covering 0-based position i."""
    if not seq or i < 0 or i >= len(seq):
        return 0
    b = seq[i]
    lo = i
    while lo > 0 and seq[lo - 1] == b:
        lo -= 1
    hi = i
    while hi < len(seq) - 1 and seq[hi + 1] == b:
        hi += 1
    return hi - lo + 1


def homopolymer_len(seq, pos, vtype):
    """Longest reference homopolymer touching the indel (anchor base is at 1-based pos)."""
    if vtype not in ("INS", "DEL") or seq is None:
        return 0
    # 0-based: pos-1 is the anchor, pos is the first base after it
    return max(run_length(seq, pos - 1), run_length(seq, pos))


def assign_filters(v, a):
    reasons = []
    if v["caller_filter"] not in ("PASS", "."):
        reasons.append(v["caller_filter"])
    if v["depth"] < a.min_depth:
        reasons.append("low_depth")
    if not (v["af"] >= a.min_af):          # also catches NaN
        reasons.append("low_af")
    if v["alt_depth"] < a.min_alt_reads:
        reasons.append("few_alt_reads")
    # strand bias only when strands are known and the locus is covered from both strands
    # (amplicon ends are not)
    strands_known = None not in (v["ref_fwd"], v["ref_rev"], v["alt_fwd"], v["alt_rev"])
    if strands_known and v["alt_depth"] >= 10 and v["ref_fwd"] > 0 and v["ref_rev"] > 0:
        minor_strand = min(v["alt_fwd"], v["alt_rev"]) / v["alt_depth"]
        if minor_strand < a.strand_min_frac:
            reasons.append("strand_bias")
    if a.platform == "iontorrent" and v["type"] in ("INS", "DEL") \
            and v["homopolymer"] >= a.iontorrent_max_hp:
        reasons.append("homopolymer_indel")
    return ";".join(dict.fromkeys(reasons)) or "PASS"


# ─────────────────────────────────────────────────────────────────────────────
# per-reference summary
# ─────────────────────────────────────────────────────────────────────────────
def read_coverage(path):
    try:
        cov = pd.read_csv(path, sep="\t", dtype={"#rname": str})
    except pd.errors.EmptyDataError:
        return {}
    cov = cov.rename(columns={"#rname": "reference"})
    return {r["reference"]: r for r in cov.to_dict("records")}


def read_callable(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0]:
                out[f[0]] = int(f[1])
    return out


def nan_to_none(x):
    if x is None:
        return None
    if isinstance(x, float) and math.isnan(x):
        return None
    try:
        if pd.isna(x):
            return None
    except (TypeError, ValueError):
        pass
    return x


def main():
    a = parse_arguments()

    seqs, desc = read_fasta(a.ref_fasta)
    refs = pd.read_csv(a.references, sep="\t", dtype={"reference": str})
    refs_meta = {r["reference"]: r for r in refs.to_dict("records")}
    cov = read_coverage(a.coverage)
    callable_bp = read_callable(a.callable)

    raw = read_lofreq_vcf(a.calls) if a.caller == "lofreq" else read_ivar_tsv(a.calls, seqs)

    variants = []
    for v in raw:
        v["sample"] = a.sample
        v["type"] = variant_type(v["ref"], v["alt"])
        v["class"] = "consensus" if v["af"] >= a.consensus_af else "minor"
        v["homopolymer"] = homopolymer_len(seqs.get(v["reference"]), v["pos"], v["type"])
        v["taxon_name"] = nan_to_none(refs_meta.get(v["reference"], {}).get("taxon_name"))
        v["filter"] = assign_filters(v, a)
        variants.append(v)

    vdf = pd.DataFrame(variants, columns=VARIANT_COLUMNS)
    if not vdf.empty:
        vdf = vdf.sort_values(["reference", "pos", "alt"])
        vdf["af"] = vdf["af"].round(4)
    vdf.to_csv(f"{a.prefix}_variants.tsv", sep="\t", index=False)

    passed = vdf[vdf["filter"] == "PASS"] if not vdf.empty else vdf
    by_ref = defaultdict(list)
    for v in passed.to_dict("records"):
        by_ref[v["reference"]].append(v)

    summary_rows, json_refs = [], []
    for ref in refs["reference"]:
        meta = refs_meta[ref]
        c = cov.get(ref, {})
        ref_len = len(seqs.get(ref, "")) or int(c.get("endpos", 0) or 0)
        n_callable = callable_bp.get(ref, 0)
        vs = by_ref.get(ref, [])

        def count(cls, types):
            return sum(1 for v in vs if v["class"] == cls and v["type"] in types)

        n_cons_snv = count("consensus", ("SNV", "MNV"))
        n_cons_indel = count("consensus", ("INS", "DEL"))
        n_minor_snv = count("minor", ("SNV", "MNV"))
        n_minor_indel = count("minor", ("INS", "DEL"))
        minor_afs = [v["af"] for v in vs if v["class"] == "minor"]
        mapped = int(c.get("numreads", 0) or 0)
        callable_frac = n_callable / ref_len if ref_len else 0.0

        if mapped == 0:
            status = "no_reads"
        elif callable_frac < a.min_callable_frac:
            status = "low_coverage"
        else:
            status = "ok"

        row = {
            "sample": a.sample,
            "reference": ref,
            "reference_description": desc.get(ref, ""),
            "taxid": nan_to_none(meta.get("taxid")),
            "taxon_name": nan_to_none(meta.get("taxon_name")),
            "lineage": nan_to_none(meta.get("lineage")),
            "ref_length": ref_len,
            "n_contigs": nan_to_none(meta.get("n_contigs")),
            "contig_aln_bp": nan_to_none(meta.get("contig_aln_bp")),
            "contig_fident": nan_to_none(meta.get("contig_fident")),
            "collapsed_refs": nan_to_none(meta.get("collapsed_refs")) or "",
            "mapped_reads": mapped,
            "breadth_pct": round(float(c.get("coverage", 0.0) or 0.0), 2),
            "mean_depth": round(float(c.get("meandepth", 0.0) or 0.0), 2),
            "mean_mapq": round(float(c.get("meanmapq", 0.0) or 0.0), 1),
            "callable_bases": n_callable,
            "callable_pct": round(100.0 * callable_frac, 2),
            "n_consensus_snv": n_cons_snv,
            "n_consensus_indel": n_cons_indel,
            "n_minor_snv": n_minor_snv,
            "n_minor_indel": n_minor_indel,
            "consensus_snv_per_kb": round(1000.0 * n_cons_snv / n_callable, 3) if n_callable else None,
            "minor_snv_per_kb": round(1000.0 * n_minor_snv / n_callable, 3) if n_callable else None,
            "consensus_identity_est": round(1.0 - (n_cons_snv + n_cons_indel) / n_callable, 5)
                                       if n_callable else None,
            "max_minor_af": round(max(minor_afs), 4) if minor_afs else None,
            "mean_minor_af": round(sum(minor_afs) / len(minor_afs), 4) if minor_afs else None,
            "status": status,
        }
        summary_rows.append(row)

        ranked = sorted(vs, key=lambda v: (v["class"] != "consensus", -v["af"], v["pos"]))
        shown = ranked[:a.max_report_per_ref]
        json_refs.append({
            **{k: nan_to_none(v) for k, v in row.items() if k != "sample"},
            "variants": [
                {k: nan_to_none(v[k]) for k in ("pos", "ref", "alt", "type", "class", "af",
                                                "depth", "alt_depth")}
                for v in shown
            ],
            "variants_not_shown": max(0, len(ranked) - len(shown)),
        })

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(f"{a.prefix}_variant_summary.tsv", sep="\t", index=False)

    payload = {
        "sample": a.sample,
        "caller": a.caller,
        "platform": a.platform,
        "thresholds": {
            "min_depth": a.min_depth, "min_af": a.min_af, "min_alt_reads": a.min_alt_reads,
            "consensus_af": a.consensus_af,
        },
        "definitions": {
            "consensus": f"allele frequency >= {a.consensus_af}: the sample's majority base differs from the reference",
            "minor": f"allele frequency between {a.min_af} and {a.consensus_af}: intra-sample (sub-population) variant",
            "callable_bases": f"reference positions with read depth >= {a.min_depth}; only these can yield calls",
            "consensus_identity_est": "1 - (consensus SNV+indel events) / callable bases",
            "status": "ok | low_coverage (callable < "
                      f"{int(100 * a.min_callable_frac)}% of reference) | no_reads",
        },
        "n_variants_total": int(len(vdf)),
        "n_variants_pass": int(len(passed)),
        "references": json_refs,
    }
    with open(f"{a.prefix}_variants.json", "w") as fh:
        json.dump(payload, fh, indent=2, default=str)

    sys.stderr.write(f"[summarize_variants] {a.sample}: {len(refs)} reference(s), "
                     f"{len(vdf)} call(s), {len(passed)} PASS\n")


if __name__ == "__main__":
    main()
