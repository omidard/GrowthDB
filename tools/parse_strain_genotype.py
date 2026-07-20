#!/usr/bin/env python3
"""Structure every GrowthDB strain string into {parent, genotype, is_wt, simulable} so records can
be organised as a hierarchy (parent strain -> wild type + derivatives) and the autocurator can
run genotype-aware simulations: knock the deleted gene(s) out of the model, load the record's
medium, simulate, compare to the measured value.

parent      : the base strain designation (MG1655, BW25113, ATCC 14266, N16961, H26, ...)
genotype    : list of ops, each {op, genes|desc}. op in:
                'del'  gene deletion(s)          -> simulable as an FBA knockout
                'other' point mutant / expression / relocation / complementation -> NOT a simple KO
is_wt       : the record is the wild-type / parental strain (no genotype ops)
simulable   : True if wild type OR every op is a gene deletion (the autocurator can reproduce it)

Deletion gene tokens use operon shorthand: 'lutABC' -> lutA,lutB,lutC ; 'ldhA' -> ldhA ;
'pcaIJ' -> pcaI,pcaJ. Writes record['strain_parsed'] in place (dry-run unless --write).
"""
import json, os, re, sys, collections

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GR = os.path.join(ROOT, "data", "growth_records.json")
WRITE = "--write" in sys.argv

CC = r"ATCC|DSMZ|DSM|NCTC|CCUG|JCM|NBRC|IFO|CECT|LMG|CIP|NCIMB|BCRC|KCTC|NRRL|CGMCC|MCCC|PCC|UTEX|CBS|KACC|VPI|NCDO|NCFB|BCCM|KCCM"
UNID = re.compile(r"not (identified|stated|specified|given|reported)|unknown|unspecified|^\s*n/?a\s*$|^\s*-\s*$", re.I)
WT = re.compile(r"\bwild[- ]?type\b|\bwt\b|\bparent(al)?\b|\breference\b", re.I)
# a background strain: "BW25113 background", "H26 background", "in MG1655" — the genetic parent to group under
BG = re.compile(r"([A-Za-z0-9][A-Za-z0-9\-]{2,})\s+background\b", re.I)
GENE_STOP = {"the", "and", "both", "gene", "genes", "region", "operon", "its", "all", "native",
             "ectopic", "each", "for", "with", "from", "type", "wild", "empty", "strain", "mutant",
             "double", "triple", "single", "delta", "evolved", "derived", "based", "grown", "tested",
             "carrying", "expressing", "encoding", "using", "containing"}
# deletion phrasings (gene token kept case-sensitive so operon shorthand splits correctly)
DEL = re.compile(r"(?:[Dd]elta[-\s]?|Δ\s?|deletion of\s+(?:the\s+|a\s+)?|[Kk]nock(?:ed)?[- ]?out of\s+)([a-z]{2,4}[A-Za-z]*\d*)")
DEL2 = re.compile(r"\b([a-z]{2,4}[A-Za-z]*\d*)\s+(?:deletion|knockout|knocked out|deleted)\b")
# non-KO modifications -> 'other'
OTHER = re.compile(r"point mutant|[A-Z]\d{1,4}[A-Z]\b|relocat|ectopic|expressing|carrying|overexpress|complement|empty vector|::|\bp[A-Z]{2,}\d|transposon|insertion|engineered|evolved|\bALE\b|reporter|catalytically dead|amber[- ]mutant", re.I)


def split_genes(tok):
    """operon shorthand -> individual genes. lutABC -> [lutA,lutB,lutC]; ldhA -> [ldhA]; pcaIJ -> [pcaI,pcaJ]."""
    m = re.match(r"^([a-z]{2,4})([A-Z][A-Za-z]*)$", tok)
    if not m:
        return [tok]
    prefix, suffix = m.group(1), m.group(2)
    # a run of single capitals (ABC) is an operon; a capital followed by lowercase is one gene name
    if re.fullmatch(r"[A-Z]{2,}", suffix):
        return [prefix + c for c in suffix]
    return [prefix + suffix]


def base_strain(text):
    """the parent strain to GROUP under. An explicit genetic background wins (so a derivative like
    'WM1659 (ftsA*, BW25113 background)' groups under BW25113, not its own name)."""
    t = text
    m = BG.search(t)                                     # explicit "X background" -> that IS the parent
    if m and not re.match(r"^(same|this|identical|its)$", m.group(1), re.I):
        return m.group(1)
    m = re.search(r"\b(" + CC + r")\s*[-: ]?\s*(\d+[A-Za-z]?)\b", t, re.I)
    if m:
        return (m.group(1).upper() + " " + m.group(2))
    # common lab designations: letters+digits token (MG1655, BW25113, N16961, MIT9301, WH8102, TS559, H26)
    m = re.search(r"\b([A-Z]{1,4}\d{2,6}[A-Za-z]?)\b", t)
    if m:
        return m.group(1)
    m = re.search(r"\bK-?12\b", t, re.I)
    if m:
        return "K-12"
    # first clean word
    w = re.split(r"[ (,;]", t.strip())[0]
    return w or None


def parse(text):
    if not text or UNID.search(text):
        return {"parent": None, "genotype": [], "is_wt": False, "simulable": False, "unidentified": True}
    genotype = []
    dels = []
    for tok in DEL.findall(text) + DEL2.findall(text):
        if tok.lower() in GENE_STOP:
            continue
        for g in split_genes(tok):
            if g.lower() not in [d.lower() for d in dels]:
                dels.append(g)
    if dels:
        genotype.append({"op": "del", "genes": dels})
    # non-KO modifications (only count if not already fully explained by deletions + wt)
    stripped = DEL.sub(" ", text)
    if OTHER.search(stripped) and not (WT.search(text) and not dels):
        genotype.append({"op": "other", "desc": text.strip()[:120]})
    is_wt = bool(WT.search(text)) and not dels and not any(o["op"] == "other" for o in genotype)
    simulable = is_wt or (bool(genotype) and all(o["op"] == "del" for o in genotype))
    return {"parent": base_strain(text), "genotype": genotype, "is_wt": is_wt,
            "simulable": simulable, "unidentified": False}


def main():
    gr = json.load(open(GR))
    st = collections.Counter(); parents = collections.Counter(); allgenes = collections.Counter()
    rel = 0
    for r in gr:
        relevant = r.get("growth_rate_per_h") is not None or r.get("uptake_rates") or r.get("secretion_rates")
        p = parse(r.get("strain") if isinstance(r.get("strain"), str) else None)
        r["strain_parsed"] = p
        if not relevant:
            continue
        rel += 1
        if p["unidentified"]:
            st["unidentified"] += 1
        elif p["is_wt"]:
            st["wild_type"] += 1
        elif any(o["op"] == "del" for o in p["genotype"]):
            st["deletion_derivative"] += 1
            for o in p["genotype"]:
                if o["op"] == "del":
                    for g in o["genes"]:
                        allgenes[g] += 1
        elif p["genotype"]:
            st["other_modification"] += 1
        else:
            st["plain_named"] += 1
        if p["parent"]:
            parents[p["parent"]] += 1
    print("validation-relevant records: %d" % rel)
    for k, v in st.most_common():
        print(f"  {v:5d}  {k}")
    print("\ndistinct parent strains:", len(parents))
    print("top parents:", parents.most_common(12))
    print("\ntop deletion genes (simulable KO targets):", allgenes.most_common(25))
    if WRITE:
        json.dump(gr, open(GR, "w"), separators=(",", ":"))
        print("\nWROTE strain_parsed into", GR)
    else:
        print("\n(dry run — pass --write to persist)")


if __name__ == "__main__":
    main()
