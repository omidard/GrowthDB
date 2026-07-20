#!/usr/bin/env python3
"""Precisely formulate the exchange composition of GrowthDB media from the paper's own
recipe (medium.composition), so records can be used to set a GEM's medium for validation.

Every ingredient is classified and mapped:
  - inorganic salt (formula) -> dissociated to its BiGG ion exchanges (curated + parser)
  - organic C/N source       -> its BiGG exchange (via the Media DB Mapper)
  - trace element / vitamin solution -> the standard trace-metal / vitamin exchanges opened
  - irrelevant (antibiotic, inducer, dye, buffer, reductant, agar) -> ignored (not a nutrient)
  - complex (yeast extract, peptone, casamino, tryptone) -> flags the medium as NOT fully defined
  - reference to another base medium not itemised here -> flags 'needs base recipe'

Nothing is guessed: an unrecognised token makes the salt/ingredient UNMAPPED (reported), never
mis-assigned. Writes medium.exchanges + medium.formulation to data/growth_records.json.
"""
import json, os, re, sys, collections

sys.path.insert(0, "/data/media_curate/tools")
from map_metabolite import Mapper
MAP = Mapper()
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GR = os.path.join(ROOT, "data", "growth_records.json")

# --- inorganic-salt FORMULA parser -> BiGG ion ids (None = element BiGG lacks; skip, still a valid salt) ---
# Cations are matched FIRST (multi-char metals before single-letter K/W) so a single-letter anion like
# F (fluoride) can never mis-match the F in Fe. Se/Mo/W/B appear only in their oxyanions here.
CATIONS = [("NH4", "nh4"), ("Na", "na1"), ("Mg", "mg2"), ("Ca", "ca2"), ("Fe", "fe2"), ("Zn", "zn2"),
           ("Mn", "mn2"), ("Co", "cobalt2"), ("Cu", "cu2"), ("Ni", "ni2"), ("Al", None), ("Ba", None),
           ("Sr", None), ("Li", None), ("K", "k")]
ANIONS = [("H2PO4", "pi"), ("HPO4", "pi"), ("PO4", "pi"), ("HCO3", "hco3"), ("CO3", "hco3"),
          ("S2O3", "tsul"), ("SO4", "so4"), ("SO3", "so3"), ("NO3", "no3"), ("NO2", "no2"),
          ("MoO4", "mobd"), ("WO4", "tungs"), ("SeO4", "slnt"), ("SeO3", "slnt"), ("SiO3", None),
          ("B4O7", None), ("BO3", None), ("OH", "oh1"), ("Cl", "cl"), ("Br", None), ("I", "iodine"), ("F", None)]

def dehydrate(f):
    f = re.sub(r"\([^)]*\)$", "", f).strip()                      # trailing note in parens
    f = re.split(r"\s*[·.*xX]\s*\d*\s*H2?O.*$", f)[0]             # ·7H2O / .2H2O / *7H2O / x 12H2O tail
    f = re.sub(r"[·.*]\s*\d*\s*H2?O.*$", "", f)
    return f.strip()

def parse_formula_salt(name):
    """Dissociate a salt written as a chemical FORMULA (KH2PO4, (NH4)2SO4, FeCl3·6H2O) into BiGG
    ion exchanges, [] for a recognised salt whose ions BiGG lacks (e.g. H3BO3), or None if it is
    not a clean recognised inorganic salt. Fe valence inferred from FeCl3 / ferric."""
    f = dehydrate(name)
    if not f or " " in f or re.search(r"[a-z]{3,}", f):           # words -> not a formula
        return None
    ferric = bool(re.search(r"Fe.?Cl3|Fe2\(|Fe2O3|Fe\(III\)|ferric", name, re.I))
    rest = f
    ions, matched, skipped = [], False, False
    for cat, ion in CATIONS:
        if cat in rest:
            matched = True
            i = "fe3" if (cat == "Fe" and ferric) else ion
            if i is None:
                skipped = True
            elif i not in ions:
                ions.append(i)
            rest = rest.replace(cat, " ")
    for grp, ion in ANIONS:
        if grp in rest:
            matched = True
            if ion is None:
                skipped = True
            elif ion not in ions:
                ions.append(ion)
            rest = rest.replace(grp, " ")
    leftover = re.sub(r"[0-9()·.\sH]", "", rest)                  # H (acid/water protons) is benign
    if leftover:                                                  # unrecognised element -> refuse (no guessing)
        return None
    if ions:
        return ions
    return [] if (matched and skipped) else None

CARBON_CAP = 10.0
COMPLEX = re.compile(r"yeast extract|peptone|tryptone|casitone|casamino|casein|beef|meat extract|"
                     r"\bbroth\b|hydrolysate|rumen|serum|blood|infusion|lysate|digest|corn steep|malt extract", re.I)
IRRELEVANT = re.compile(r"\bIPTG\b|kanamycin|ampicillin|tetracyclin|chlorampheni|streptomyc|gentamic|spectinom|"
                        r"carbenicill|resazurin|phenol red|bromothymol|bromocresol|neutral red|indicator|"
                        r"HEPES|MOPS buffer|\bMOPS\b|PIPES|\bTris\b|\bMES\b|\bBis-Tris\b|antifoam|agar\b|"
                        r"reductant|reducing agent|dithion|titanium|\bNa2S\b|sodium sulf?ide|cysteine.?HCl|"
                        r"resorufin|X-gal|inducer|selection|Tween|antibiotic", re.I)
BASE_REF = re.compile(r"salt.* not itemis|not itemized|composition not (?:stated|given)|concentrations not|"
                      r"per (?:the )?(?:cited|original|ref)|full composition not|M199|\bDMEM\b|\bM17\b|"
                      r"(?:AM1|Gutnick|CGXII|CgXII|M9|MOPS|NBS|BG11).{0,30}(?:salts|mineral|base|medium)$", re.I)
TRACE = re.compile(r"trace (?:element|metal)|SL-?\d|Wolfe.?s (?:mineral|metal)", re.I)
VITAMIN = re.compile(r"vitamin", re.I)
GAS = {"co2": "co2", "carbon dioxide": "co2", "h2": "h2", "hydrogen": "h2", "co": "co", "n2": None, "o2": "o2"}
TRACE_METALS = ["fe2", "mn2", "zn2", "cu2", "cobalt2", "mobd", "ni2", "ca2", "mg2", "cl", "so4", "slnt", "tungs"]
VITAMINS_STD = ["btn", "thm", "ribflv", "pnto__R", "pydxn", "nac", "fol", "cbl1", "4abz", "lipoate"]

def carbon_c(bigg):
    props = getattr(MAP, "_bigg", None)
    return None

try:
    from enrich_coverage import parse_salt as word_salt          # word-based salt dissociation ("sodium bicarbonate")
except Exception:
    word_salt = lambda x: None

def ingredient_exchanges(name):
    """Return (list of (bigg, is_carbon), classification) for one ingredient. Nothing is guessed."""
    n = name.strip()
    if not n:
        return [], "empty"
    if IRRELEVANT.search(n):
        return [], "irrelevant"
    if COMPLEX.search(n):
        return [], "complex"
    if BASE_REF.search(n):
        return [], "base_ref"
    if TRACE.search(n):
        return [(b, False) for b in TRACE_METALS], "trace"
    if VITAMIN.search(n):
        return [(b, False) for b in VITAMINS_STD], "vitamin"
    low = n.lower().strip()
    if re.search(r"\b(gas|atmosphere|sparge|headspace|head-space)\b", low) or low in GAS:
        got = []
        if re.search(r"\bco2\b|carbon dioxide", low): got.append(("co2", False))
        if re.search(r"\bh2\b|hydrogen", low): got.append(("h2", True))
        if re.search(r"\bco\b(?!2)", low): got.append(("co", True))
        if got or low in GAS:
            return got, "gas"
    # inorganic salt written as a formula — try the token, and any formula inside parentheses
    for cand in [n] + re.findall(r"\(([^)]*)\)", n):
        ions = parse_formula_salt(cand)
        if ions is not None:
            return [(b, False) for b in ions], "salt"
    # inorganic salt written as words ("sodium bicarbonate", "ferric ammonium citrate")
    ws = word_salt(n)
    if ws:
        return [(b, False) for b in ws], "salt_word"
    # organic / named compound via the Media DB mapper (parentheticals stripped)
    h = MAP.map(name=re.sub(r"\([^)]*\)", "", n).strip())
    if h and h.get("in_biggr") and h.get("bigg_metabolite"):
        b = h["bigg_metabolite"]
        f = (BIGG.get(b, {}).get("xrefs", {}) or {}).get("formula", "")
        mc = re.search(r"C(\d*)", f)
        is_c = bool(mc) and (int(mc.group(1) or 1) >= 2) and b not in ("co2", "hco3")
        return [(b, is_c)], "organic"
    return [], "unmapped"

BIGG = json.load(open("/data/media_curate/tools/bigg_metabolite_dict.json"))
INORG_ION = {"pi", "so4", "cl", "hco3", "no3", "no2", "mobd", "tungs", "slnt", "so3", "tsul", "oh1", "iodine",
             "na1", "k", "nh4", "mg2", "ca2", "fe2", "fe3", "zn2", "mn2", "cobalt2", "cu2", "ni2", "h2o", "h"}

def formulate(comp):
    ex = {}                       # bigg -> lb (most-open wins for ions; carbon capped)
    classes = collections.Counter()
    unmapped = []
    for c in comp:
        nm = (c.get("name") or "").strip()
        got, cls = ingredient_exchanges(nm)
        classes[cls] += 1
        if cls in ("unmapped", "vitamin_unlisted") and nm:
            unmapped.append(nm)
        for bigg, is_c in got:
            lb = -CARBON_CAP if (is_c and bigg not in INORG_ION) else -1000.0
            # keep the most permissive uptake bound for a metabolite
            if bigg not in ex or lb < ex[bigg]:
                ex[bigg] = min(ex.get(bigg, 0), lb)
    exch = [{"exchange": "EX_" + b + "_e", "bigg": b, "lb": ex[b], "ub": 1000.0} for b in ex]
    return exch, classes, unmapped

def status(classes, exch):
    if classes.get("complex"):
        return "partial_complex"                # contains undefined biological extract
    if classes.get("base_ref") and not exch:
        return "needs_base_recipe"
    if classes.get("unmapped") or classes.get("vitamin_unlisted"):
        return "partial_unmapped"
    if any(x["lb"] == -CARBON_CAP for x in exch):
        return "defined"                        # has a defined carbon source
    return "defined_nocarbon"

def main(apply=False):
    gr = json.load(open(GR))
    done = collections.Counter()
    unmapped_all = collections.Counter()
    n = 0
    for r in gr:
        m = r.get("medium") or {}
        if m.get("media_id") or not m.get("composition"):
            continue
        exch, classes, unmapped = formulate(m["composition"])
        st = status(classes, exch)
        done[st] += 1
        for u in unmapped:
            unmapped_all[u] += 1
        if apply:
            m["exchanges"] = exch
            m["formulation"] = st
            m["formulated_from"] = "composition"
        n += 1
    if apply:
        json.dump(gr, open(GR, "w"), separators=(",", ":"))
    print("formulated %d media-with-composition:" % n, dict(done))
    print("top unmapped ingredients (need a table entry or are genuinely complex):")
    for u, c in unmapped_all.most_common(25):
        print(f"   {c:4}  {u}")

if __name__ == "__main__":
    main(apply="--apply" in sys.argv)
