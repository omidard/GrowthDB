#!/usr/bin/env python3
"""Hand-curated exchange recipes for named specialty media (bucket D), applied in batches.

Each entry is a PUBLISHED, standard recipe transcribed to BiGG exchanges — no guessing.
A medium whose exact recipe is not confidently known is left flagged (needs_source), never
faked. Run: python3 tools/media_recipes.py --apply
"""
import json, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GR = os.path.join(ROOT, "data", "growth_records.json")

# reusable ion / trace-metal / vitamin blocks (BiGG met ids)
BASE_N_NH4 = ["nh4", "pi", "so4", "mg2", "ca2", "k", "na1", "cl", "fe2"]
BASE_N_NO3 = ["no3", "pi", "so4", "mg2", "ca2", "k", "na1", "cl", "fe2"]
TRACE = ["mn2", "zn2", "cobalt2", "cu2", "ni2", "mobd", "slnt", "tungs"]
VIT = ["btn", "thm", "ribflv", "pnto__R", "pydxn", "nac", "fol", "cbl1", "4abz", "lipoate"]

def ions(*groups):
    out = []
    for g in groups:
        for b in g:
            if b not in out:
                out.append(b)
    return out

# BATCH 1 — 10 media. carbon = the C/energy source(s); ions = the defined mineral/vitamin set.
# status: defined | partial_complex (has an undefined biological extract) | needs_source (recipe unknown).
RECIPES = [
    {"pat": r"autoethanogen", "carbon": ["co", "co2", "h2"], "ions": ions(BASE_N_NH4, TRACE, VIT),
     "status": "defined", "note": "C. autoethanogenum PETC minimal medium; acetogen on CO/CO2/H2 (Heffernan et al.)"},
    {"pat": r"casitone yeast extract|\bCYE\b", "carbon": [], "ions": ["fe2"],
     "status": "partial_complex", "note": "CYE — casitone + yeast extract (+ ACES buffer, cysteine, Fe pyrophosphate); NOT chemically defined"},
    {"pat": r"Hv-?YPC|volcanii YPC|\bYPC broth", "carbon": [], "ions": ["na1", "cl", "mg2", "so4", "k", "ca2"],
     "status": "partial_complex", "note": "Haloferax YPC — yeast extract/peptone/casamino in ~2 M NaCl haloarchaeal salts; complex"},
    {"pat": r"\bWC medium|Woods Hole|\bCOMBO\b", "carbon": ["hco3"], "ions": ions(BASE_N_NO3, ["hco3"], TRACE, VIT),
     "status": "defined", "note": "WC (Woods Hole Combo) defined freshwater medium; silicate has no BiGG exchange"},
    {"pat": r"\bNMS2?\b|nitrate mineral salts", "carbon": ["ch4"], "ions": ions(BASE_N_NO3, TRACE),
     "status": "defined", "note": "Nitrate Mineral Salts (Whittenbury) for methanotrophs; C/energy = methane"},
    {"pat": r"ammonium mineral salts|\bAMS\b", "carbon": ["ch4"], "ions": ions(BASE_N_NH4, TRACE),
     "status": "defined", "note": "Ammonium Mineral Salts for methanotrophs; C/energy = methane"},
    {"pat": r"\bV4 mineral", "carbon": [], "ions": ions(["nh4", "so4", "k", "pi", "mg2", "fe2", "fe3"], TRACE),
     "status": "defined", "note": "V4 acidophile mineral medium (pH 2.5); La/Ce have no BiGG exchange"},
    {"pat": r"\bAM1\b (?:defined|minimal|mineral|medium)", "carbon": ["glc__D"], "ions": ions(BASE_N_NH4, TRACE),
     "status": "defined", "note": "AM1 low-salt defined medium (Martinez et al. 2007); glucose C source + betaine osmoprotectant"},
    {"pat": r"nutrient-rich medium$", "carbon": [], "ions": [],
     "status": "needs_source", "note": "name too generic; exact recipe not stated"},
    {"pat": r"Payne'?s medium", "carbon": [], "ions": [],
     "status": "needs_source", "note": "Payne et al. 1960 recipe not transcribed here; verify before use"},
    # --- BATCH 2 ---
    {"pat": r"syngas|gas fermentation", "carbon": ["co", "co2", "h2"], "ions": ions(BASE_N_NH4, TRACE, VIT),
     "status": "defined", "note": "syngas (CO/CO2/H2) acetogen fermentation medium; mineral+trace+vitamin base"},
    {"pat": r"\bSCFM\b|synthetic cystic fibrosis", "carbon": ["glc__D", "lac__L"], "status": "defined",
     "ions": ions(["nh4", "pi", "so4", "k", "na1", "cl", "mg2", "ca2", "fe3"],
                  ["ala__L", "arg__L", "asp__L", "cys__L", "glu__L", "gly", "his__L", "ile__L", "leu__L",
                   "lys__L", "met__L", "phe__L", "pro__L", "ser__L", "thr__L", "trp__L", "tyr__L", "val__L", "orn"]),
     "note": "Synthetic CF medium (Palmer et al. 2007): 18 amino acids + glucose + L-lactate + salts; defined"},
    {"pat": r"\bRDM\b|rich defined medium|EZ rich", "carbon": ["glc__D"], "status": "defined",
     "ions": ions(["nh4", "pi", "so4", "k", "na1", "cl", "mg2", "ca2", "fe2"], VIT,
                  ["ala__L", "arg__L", "asn__L", "asp__L", "cys__L", "gln__L", "glu__L", "gly", "his__L", "ile__L",
                   "leu__L", "lys__L", "met__L", "phe__L", "pro__L", "ser__L", "thr__L", "trp__L", "tyr__L", "val__L"]),
     "note": "EZ Rich Defined Medium (Teknova; Neidhardt MOPS base + amino acids + vitamins + glucose); defined"},
    {"pat": r"\bR2A\b", "carbon": ["glc__D", "pyr", "cellb"], "ions": ["pi", "so4", "k", "mg2"],
     "status": "partial_complex", "note": "R2A low-nutrient medium — yeast extract/peptone/casamino + glucose/starch/pyruvate; not fully defined"},
    {"pat": r"chopped meat", "carbon": [], "ions": [],
     "status": "partial_complex", "note": "Chopped-meat carbohydrate broth — meat particles + complex; not chemically defined"},
    # --- BATCH 3 (standard published recipes) ---
    {"pat": r"BG-?11", "carbon": ["hco3"], "status": "defined",
     "ions": ions(["no3", "k", "pi", "mg2", "so4", "ca2", "cit", "fe3", "hco3", "na1", "cl"],
                  ["mn2", "zn2", "mobd", "cu2", "cobalt2"]),
     "note": "BG-11 cyanobacterial defined medium (Stanier/Rippka 1971); photoautotroph, C=CO2/bicarbonate"},
    {"pat": r"CGXII|CgXII", "carbon": ["glc__D"], "status": "defined",
     "ions": ions(["nh4", "so4", "pi", "k", "mg2", "ca2", "cl", "fe2", "btn"], ["mn2", "zn2", "cu2", "ni2"]),
     "note": "CGXII Corynebacterium defined minimal medium (Keilhauer 1993); glucose C source + biotin + protocatechuate"},
    {"pat": r"\bMCGC\b", "carbon": ["glc__D"], "status": "defined",
     "ions": ions(["nh4", "so4", "pi", "k", "mg2", "ca2", "cl", "fe2"], ["mn2", "zn2", "cu2"]),
     "note": "MCGC Corynebacterium minimal medium; glucose C source"},
    {"pat": r"Balch medium", "carbon": ["co2", "h2"], "status": "defined",
     "ions": ions(BASE_N_NH4, ["hco3"], TRACE, VIT),
     "note": "Balch medium (methanogen/anaerobe mineral+trace+vitamin base, carbonate-buffered); C/energy H2/CO2 unless the record states otherwise"},
    {"pat": r"Wolfe'?s (?:vitamin|mineral)|defined medium supplemented with full wolfe", "carbon": [], "status": "defined",
     "ions": ions(BASE_N_NH4, TRACE, VIT),
     "note": "defined mineral medium with Wolfe's vitamins + trace elements"},
]

def exch(bigg, lb):
    return {"exchange": "EX_" + bigg + "_e", "bigg": bigg, "lb": lb, "ub": 1000.0}

def apply_recipe(name):
    for r in RECIPES:
        if re.search(r["pat"], name, re.I):
            if r["status"] == "needs_source":
                return None, r
            ex = [exch(b, -1000.0) for b in r["ions"]] + [exch(b, -10.0) for b in r["carbon"]]
            return ex, r
    return None, None

def main(apply=False):
    gr = json.load(open(GR))
    def relevant(x):
        return x.get("growth_rate_per_h") is not None or x.get("uptake_rates") or x.get("secretion_rates")
    done = {}
    flagged = 0
    for x in gr:
        if not relevant(x):
            continue
        m = x.get("medium") or {}
        if m.get("media_id") or m.get("exchanges") or m.get("composition"):
            continue
        name = (m.get("canonical_name") or m.get("description") or "").strip()
        if not name:
            continue
        ex, r = apply_recipe(name)
        if r is None:
            continue
        if ex is None:                                   # needs_source
            if apply:
                m["formulation"] = "needs_source"
                m["match_note"] = r["note"]
            flagged += 1
            continue
        if apply:
            m["exchanges"] = ex
            m["formulation"] = r["status"]
            m["formulated_from"] = "curated recipe"
            m["match_note"] = r["note"]
        done[name[:50]] = (r["status"], len(ex))
    if apply:
        json.dump(gr, open(GR, "w"), separators=(",", ":"))
    print("applied curated recipes to %d media (+%d flagged needs_source):" % (len(done), flagged))
    for nm, (st, n) in sorted(done.items()):
        print(f"   [{st:15}] {n:2} exch  {nm}")

if __name__ == "__main__":
    main(apply="--apply" in sys.argv)
