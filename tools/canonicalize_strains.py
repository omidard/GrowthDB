#!/usr/bin/env python3
"""Canonicalise every GrowthDB parent strain so a strain and its synonyms/deposits/backgrounds
collapse to ONE grouping key — no more 'MG1655' vs 'K-12 MG1655' vs 'ATCC 47076' split, and an
engineered derivative 'ASC662 (MG1655 lacZ)' groups under its reference background MG1655, not
its construct name.

Two mechanisms:
  1. ORGANISM-SCOPED reference registry: per organism, a canonical strain -> all its surface forms
     (name variants, culture-collection deposits, known-derivative names). Scoping stops numeric
     synonyms ('168', '824', 'S2') cross-matching between organisms.
  2. Structural rules: a known reference found ANYWHERE in the string wins as the parent (so the
     background in parentheses beats the construct name); culture-collection ids are normalised
     (ATCC13032 -> ATCC 13032; DSMZ -> DSM); junk tokens (str., wild, the, single letters) are
     never a parent.

Writes strain_parsed.parent (canonicalised) + strain_parsed.parent_raw (the pre-canonical guess)
in place. Dry-run unless --write.
"""
import json, os, re, sys, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import parse_strain_genotype as PSG

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GR = os.path.join(ROOT, "data", "growth_records.json")
WRITE = "--write" in sys.argv

# organism (lowercase species substring) -> { canonical parent : [surface-form regexes] }.
# Deposits (ATCC/DSM/CGSC numbers) are strain-unique, so they safely map a deposit id to the strain.
REGISTRY = {
    "escherichia coli": {
        "MG1655": [r"mg[\s\-]?1655", r"atcc[\s\-]?47076", r"atcc[\s\-]?700926", r"atcc[\s\-]?700925",
                   r"cgsc[\s\-]?7740", r"dsm[\s\-]?18039", r"\bsij488\b", r"\bmds42\b", r"\bju15\b", r"\becl45\b", r"\bmgcsc"],
        "W3110": [r"\bw3110\b", r"atcc[\s\-]?27325", r"dsmz?[\s\-]?5911", r"cgsc[\s\-]?4474", r"\blj110\b", r"\bwl3110\b"],
        "BW25113": [r"\bbw25113\b", r"cgsc[\s\-]?7636", r"\bjw\d{3,4}\b", r"\bkeio\b"],
        "NCM3722": [r"\bncm3722\b", r"cgsc[\s\-]?12355", r"\bnq\d{3,4}\b"],
        "BW2952": [r"\bbw2952\b"],
        "MC4100": [r"\bmc4100\b"],
        "MG1063/REL606": [r"\brel606\b"],
        "Nissle 1917": [r"nissle", r"dsm[\s\-]?6601"],
        "BL21": [r"\bbl21\b"],
        "MC1061": [r"\bmc1061\b"],
        "DH5alpha": [r"dh5[\s]?[αa]"],
        "DH10B": [r"\bdh10b\b"],
        "JM109": [r"\bjm109\b"],
        "JM108": [r"\bjm108\b"],
        "TOP10": [r"\btop10\b"],
        "CFT073": [r"\bcft073\b"], "UTI89": [r"\buti89\b"], "O157:H7 Sakai": [r"\bsakai\b"],
        "Crooks ATCC 8739": [r"atcc[\s\-]?8739", r"\bcrooks\b"],
        "E. coli W": [r"atcc[\s\-]?9637", r"\bkj134\b"],
    },
    "corynebacterium glutamicum": {"ATCC 13032": [r"atcc[\s\-]?13032", r"\b13032\b"],
                                   "ATCC 13869": [r"atcc[\s\-]?13869", r"\b13869\b"], "R": [r"\bstrain r\b"]},
    "bacillus subtilis": {"168": [r"\b168\b", r"\bw168\b", r"\bbsb168\b", r"\b1a1\b"], "PY79": [r"\bpy79\b"],
                          "JH642": [r"\bjh642\b"]},
    "pseudomonas putida": {"KT2440": [r"\bkt2440\b", r"dsm[\s\-]?6125", r"atcc[\s\-]?47054"],
                           "KT2442": [r"\bkt2442\b"], "S12": [r"\bs12\b"]},
    "zymomonas mobilis": {"ZM4": [r"\bzm4\b", r"atcc[\s\-]?31821"], "CP4": [r"\bcp4\b"], "8b": [r"\b8b\b"]},
    "synechocystis": {"PCC 6803": [r"pcc[\s\-]?6803"]},
    "synechococcus": {"PCC 7002": [r"pcc[\s\-]?7002"], "PCC 7942": [r"pcc[\s\-]?7942"],
                      "WH8102": [r"\bwh8102\b"], "WH7803": [r"\bwh7803\b"]},
    "vibrio natriegens": {"ATCC 14048": [r"atcc[\s\-]?14048", r"dsm[\s\-]?759", r"\bvmax\b"]},
    "staphylococcus aureus": {"NCTC 8325": [r"nctc[\s\-]?8325", r"\b8325[\s\-]?4?\b", r"\bhg00[13]\b", r"\bsh1000\b"],
                              "USA300": [r"usa\s?300", r"\blac\b", r"\bfpr3757\b", r"\bje2\b"],
                              "Newman": [r"\bnewman\b"], "COL": [r"\bcol\b"], "N315": [r"\bn315\b"], "Mu50": [r"\bmu50\b"]},
    "pseudomonas aeruginosa": {"PAO1": [r"\bpao1\b"], "PA14": [r"\bpa14\b", r"\bucbpp[\s\-]?pa14\b"], "PAK": [r"\bpak\b"]},
    "clostridium acetobutylicum": {"ATCC 824": [r"atcc[\s\-]?824", r"\b824\b", r"dsm[\s\-]?792"]},
    "clostridium ljungdahlii": {"DSM 13528": [r"dsm[\s\-]?13528", r"atcc[\s\-]?55383", r"\bpetc\b"]},
    "methanococcus maripaludis": {"S2": [r"\bs2\b", r"dsm[\s\-]?14266", r"\bs0001\b", r"\bmm901\b", r"\bll\b"]},
    "geobacter sulfurreducens": {"PCA": [r"\bpca\b", r"dsmz?[\s\-]?12127", r"atcc[\s\-]?51573", r"\bdl1\b"]},
    "geobacter metallireducens": {"GS-15": [r"gs[\s\-]?15", r"atcc[\s\-]?53774", r"dsm[\s\-]?7210"]},
    "mycobacterium smegmatis": {"mc2 155": [r"mc.?2.?155", r"mc²155"]},
    "mycolicibacterium smegmatis": {"mc2 155": [r"mc.?2.?155", r"mc²155"]},
    "cupriavidus necator": {"H16": [r"\bh16\b", r"atcc[\s\-]?17699", r"dsm[\s\-]?428", r"eutropha h16"]},
    "ralstonia eutropha": {"H16": [r"\bh16\b", r"atcc[\s\-]?17699", r"dsm[\s\-]?428"]},
    "saccharomyces cerevisiae": {"S288C": [r"\bs288c\b", r"\bby4741\b", r"\bby4742\b", r"\bcen\.?pk"], "CEN.PK": [r"cen\.?pk"]},
    "streptomyces coelicolor": {"A3(2)": [r"a3\(2\)", r"\bm145\b", r"\bm600\b", r"\bj1501\b"]},
    "lactococcus lactis": {"MG1363": [r"\bmg1363\b", r"\bnz9000\b"], "IL1403": [r"\bil1403\b"]},
    "streptococcus thermophilus": {"LMD-9": [r"lmd[\s\-]?9", r"\bst18311\b"], "LMG 18311": [r"lmg[\s\-]?18311"]},
}

CC = r"ATCC|DSMZ|DSM|NCTC|CCUG|JCM|NBRC|IFO|CECT|LMG|CIP|NCIMB|BCRC|KCTC|NRRL|CGMCC|MCCC|PCC|UTEX|CBS|KACC|VPI|NCDO|NCFB|BCCM|KCCM|CGSC"
JUNK = re.compile(r"^(wild|type|wild-?type|derivative|derivatives|mutant|strain|parent|parental|deletion|the|and|of|a|an|E|sp|K-?12|K12|type|control|reference|delta|background|clone|isolate|evolved|engineered|str|substr|in|not|uncultured|enrichment|community|recombinant|variant|producer|host)$", re.I)


def reg_for(org):
    o = (org or "").lower()
    for key, reg in REGISTRY.items():
        if key in o:
            return reg
    return None


def normalize_cc(s):
    """Extract + normalise a culture-collection id (ATCC13032 -> ATCC 13032; DSMZ -> DSM)."""
    m = re.search(r"\b(" + CC + r")\s*[\-:# ]?\s*(\d+[A-Za-z]?)\b", s, re.I)
    if not m:
        return None
    coll = m.group(1).upper().replace("DSMZ", "DSM").replace("NCFB", "NCDO")
    return coll + " " + m.group(2)


def canonical_parent(strain_str, organism):
    if not strain_str:
        return None
    s = str(strain_str)
    low = s.lower()
    reg = reg_for(organism)
    if reg:
        for canon, pats in reg.items():
            for p in pats:
                if re.search(p, low):
                    return canon
    # no registry hit -> a normalised culture-collection id is a good stable parent
    cc = normalize_cc(s)
    if cc:
        return cc
    # else the base designation from the parser, but reject junk tokens
    base = PSG.base_strain(s)
    if base and not JUNK.match(base) and len(base) > 1 and not base[0].islower():
        # normalise MG-1655 -> MG1655 style (drop an internal hyphen between letters and digits)
        base = re.sub(r"^([A-Za-z]+)-(\d)", r"\1\2", base)
        return base
    # a distinctive leading construct/strain token (the head before the first paren/comma), the strain's own identity
    head = re.split(r"\s*[(,]", s.strip())[0].strip()
    htoks = head.split()
    if 1 <= len(htoks) <= 2 and len(head) <= 24:
        cand = htoks[-1] if (len(htoks) == 2 and re.match(r"^K[\s\-]?12$", htoks[0], re.I)) else head
        if not JUNK.match(cand.split()[0]) and (re.search(r"\d", cand) or re.search(r"[a-z][A-Z]", cand)) and not re.match(r"^(19|20)\d\d$", cand):
            return re.sub(r"^([A-Za-z]+)-(\d)", r"\1\2", cand)
    # last resort: an explicit lineage (K-12 / B / C) so it's grouped by lineage rather than left blank
    if re.search(r"\bK[\s\-]?12\b", s, re.I):
        return "K-12 (unspecified derivative)"
    if re.search(r"\bB\s?strain\b|\bBL\b", s, re.I):
        return "B lineage"
    return None


def parent_key(name):
    """A GENERAL, registry-independent grouping key: two surface forms of the SAME strain collapse
    regardless of spacing / hyphens / dots / case / 'str.'/'substr.'/'K-12' prefixes. This is what
    stops 'NCM3722' vs 'NCM 3722', 'str. SolV' vs 'SolV', 'EA2' vs 'Ea2', 'DH5-alpha' vs 'DH5alpha'
    from ever being two strains, even for strains not in the registry."""
    s = (name or "").lower()
    s = re.sub(r"^(str\.?|substr\.?|strain|substrain)\s+", "", s)
    s = re.sub(r"\bk[\s\-]?12\b", "", s)
    s = re.sub(r"\(unspecified derivative\)", "", s)
    s = re.sub(r"\balpha\b", "a", s)                     # DH5alpha == DH5a
    return re.sub(r"[^a-z0-9]", "", s)                   # drop spaces/hyphens/dots/^ etc.


def main():
    gr = json.load(open(GR))
    rel = [r for r in gr if r.get("growth_rate_per_h") is not None or r.get("uptake_rates") or r.get("secretion_rates")]
    # pass 1: raw canonical parent per record
    for r in gr:
        spar = r.get("strain_parsed") or {}
        org = r.get("gtdb_species") or r.get("species") or r.get("organism")
        spar["parent_raw"] = spar.get("parent")
        spar["parent"] = canonical_parent(r.get("strain") if isinstance(r.get("strain"), str) else None, org)
        r["strain_parsed"] = spar
    # pass 2: unify every parent that shares a key to ONE canonical display (registry name wins,
    # else the most common surface form) — so equivalent forms are byte-identical everywhere
    reg_names = set()
    for reg in REGISTRY.values():
        reg_names.update(reg.keys())
    key_forms = collections.defaultdict(collections.Counter)
    for r in gr:
        p = (r.get("strain_parsed") or {}).get("parent")
        if p:
            key_forms[parent_key(p)][p] += 1
    canon_by_key = {}
    for k, forms in key_forms.items():
        reg_hit = [f for f in forms if f in reg_names]
        canon_by_key[k] = reg_hit[0] if reg_hit else forms.most_common(1)[0][0]
    before = collections.Counter(); after = collections.Counter()
    ecoli_before = set(); ecoli_after = set()
    for r in gr:
        spar = r.get("strain_parsed") or {}
        p = spar.get("parent")
        if p:
            k = parent_key(p); spar["pkey"] = k; spar["parent"] = canon_by_key[k]
        r["strain_parsed"] = spar
        relevant = r.get("growth_rate_per_h") is not None or r.get("uptake_rates") or r.get("secretion_rates")
        if relevant:
            if spar.get("parent_raw"):
                before[spar["parent_raw"]] += 1
            if spar.get("parent"):
                after[spar["parent"]] += 1
            org = (r.get("gtdb_species") or r.get("species") or r.get("organism") or "")
            if "coli" in org.lower():
                if spar.get("parent_raw"):
                    ecoli_before.add(spar["parent_raw"])
                if spar.get("parent"):
                    ecoli_after.add(spar["parent"])
    # verify: no two distinct parents share a key
    coll = collections.defaultdict(set)
    for r in rel:
        p = (r.get("strain_parsed") or {}).get("parent")
        if p:
            coll[parent_key(p)].add(p)
    residual = {k: v for k, v in coll.items() if len(v) > 1}
    print("distinct parents (all organisms):  before %d -> after %d" % (len(before), len(after)))
    print("distinct E. coli parents:          before %d -> after %d" % (len(ecoli_before), len(ecoli_after)))
    print("residual key COLLISIONS (should be 0):", len(residual), residual if residual else "")
    print("\ntop canonical parents after:")
    for p, n in after.most_common(15):
        print(f"  {n:4d}  {p}")
    if WRITE:
        json.dump(gr, open(GR, "w"), separators=(",", ":"))
        print("\nWROTE", GR)
    else:
        print("\n(dry run — pass --write to persist)")


if __name__ == "__main__":
    main()
