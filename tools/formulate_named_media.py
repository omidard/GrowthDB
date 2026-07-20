#!/usr/bin/env python3
"""Resolve the media that have a NAME but no composition list, precisely and conservatively:

  A. UNUSABLE (in-vivo, environmental, 'not stated') -> flagged not_defined; never faked.
  B. generic minimal + a named carbon ('glucose minimal medium', 'malate minimal-salt medium')
     -> the carbon's BiGG exchange (the standard mineral base is supplied by the autocurator);
        formulation = defined_from_name.
  C. standard medium whose recipe is fixed (MRS, Nutrient broth, DSM<n>, Marine broth ...)
     -> matched to a Media DB entry by name/number (conservative, carbon-consistent).
  D. everything else (specialty defined media: autoethanogenum, AM1, NMS2, WC, CSBK ...)
     -> left for hand-formulation in batches of 10 (written to media_todo_batches.json).

Nothing is guessed; an unclear medium goes to bucket D for manual precision work.
"""
import json, os, re, sys, collections

sys.path.insert(0, "/data/media_curate/tools")
from map_metabolite import Mapper
MAP = Mapper()
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GR = os.path.join(ROOT, "data", "growth_records.json")
BIGG = json.load(open("/data/media_curate/tools/bigg_metabolite_dict.json"))

# genuinely NOT formulatable (no defined recipe anywhere): host/in-vivo, environmental samples, or unstated.
# Media whose recipe merely lives in a citation/SI go to bucket D (findable), not here.
UNUSABLE = re.compile(r"^(?:not stated|not specified|not explicitly|unspecified|n/s)\b|"
                      r"in vivo|in-vivo|intracellular|within (?:the )?cyto|bladder|urinary|murine|"
                      r"\bmouse\b|human (?:bladder|gut|host)|\bhost\b|"
                      r"reservoir|seawater|sea water|water column|epilimnion|sediment|\bsoil\b|"
                      r"rumen fluid|feces|fecal|gut content|sputum|natural seawater|environmental sample", re.I)
MINIMAL = re.compile(r"minimal|defined|mineral|chemically defined|synthetic", re.I)
NOT_CARBON = {"pi", "so4", "cl", "na1", "k", "nh4", "mg2", "ca2", "fe2", "fe3", "mn2", "zn2", "cu2",
              "cobalt2", "mobd", "ni2", "h2o", "h", "co2", "hco3", "o2", "no3", "no2"}

def carbon_c(b):
    f = (BIGG.get(b, {}).get("xrefs", {}) or {}).get("formula", "")
    m = re.search(r"C(\d*)", f)
    return int(m.group(1) or 1) if m else 0

def named_carbon(name):
    """The carbon source stated in a medium name, mapped to (bigg, exchange), or None."""
    cands = re.findall(r"\(([^)]+)\)", name)
    cands += re.findall(r"([A-Za-z][A-Za-z\- ]{2,20})\s+(?:as (?:the )?(?:sole )?(?:carbon|c|energy) source|minimal|as carbon)", name, re.I)
    cands += re.findall(r"([A-Za-z][A-Za-z\-]{2,18})[- ](?:limited|fed|grown|restricted|based)\b", name, re.I)   # 'glucose-limited'
    cands += re.findall(r"(?:with|on|,|\+)\s+([A-Za-z][A-Za-z\- ]{2,20})\b", name)
    cands = [name.split()[0]] + cands
    for c in cands:
        c = re.sub(r"[- ]?(limited|fed|grown|restricted|based|only|medium|minimal)\b", "", c.strip(), flags=re.I).strip()
        if len(c) < 3:
            continue
        h = MAP.map(name=c)
        if h and h.get("in_biggr") and h.get("bigg_metabolite"):
            b = h["bigg_metabolite"]
            if b not in NOT_CARBON and carbon_c(b) >= 2:
                return b, h["exchange"], c
    return None

# --- Media DB matcher for standard fixed-recipe media ---
_MEDIA = None
def media_db():
    global _MEDIA
    if _MEDIA is None:
        idx = json.load(open("/data/media_curate/data/index.json"))
        _MEDIA = idx["media"]
    return _MEDIA

STD_MATCH = {  # record-name pattern -> a distinctive Media DB name token to match (recipe is fixed)
    r"^MRS\b|de man rogosa": "mrs",
    r"nutrient broth\b": "nutrient broth",
    r"marine broth|zobell": "marine broth",
    r"\bPYE\b|peptone yeast": "pye",
    r"\bBHI\b|brain heart": "brain heart",
    r"\bTSB\b|tryptic soy|trypticase soy": "tryptic soy",
    r"\bLB\b|luria|lysogeny": "lb ",
    r"\bM17\b": "m17",
}
def dsm_number(name):
    m = re.search(r"\b(?:DSM|DSMZ|mediadive)\s*#?\s*(\d+[a-z]?)\b", name, re.I)
    return m.group(1) if m else None

def media_db_match(name):
    md = media_db()
    dn = dsm_number(name)
    if dn:
        for m in md:
            if m.get("id", "").lower() == "mediadive_" + dn.lower():
                return m["id"], "DSMZ MediaDive " + dn
    for pat, tok in STD_MATCH.items():
        if re.search(pat, name, re.I):
            # prefer a defined/std entry, else a complexlit with a real composition, that names this medium
            best = None
            for m in md:
                nm = (m.get("name") or "").lower()
                if tok in nm and (m.get("n_mapped", 0) or 0) > 0:
                    if best is None or (m.get("defined") and not best.get("defined")):
                        best = m
            if best:
                return best["id"], "name match: " + tok.strip()
    return None

def main(apply=False):
    gr = json.load(open(GR))
    def relevant(r):
        return r.get("growth_rate_per_h") is not None or r.get("uptake_rates") or r.get("secretion_rates")
    buckets = collections.Counter()
    todo = collections.Counter()
    for r in gr:
        if not relevant(r):
            continue
        m = r.get("medium") or {}
        if m.get("media_id") or m.get("exchanges") or m.get("composition"):
            continue
        name = (m.get("canonical_name") or m.get("description") or "").strip()
        if not name:
            continue
        if UNUSABLE.search(name):
            buckets["A_unusable"] += 1
            if apply:
                m["formulation"] = "not_defined"
            continue
        mm = media_db_match(name)
        if mm:
            buckets["C_mediadb"] += 1
            if apply:
                m["media_id"] = mm[0]
                m["media_url"] = "https://omidard.github.io/Media/?medium=" + mm[0]
                m["formulation"] = "matched_media_db"
                m["match_note"] = mm[1]
            continue
        if MINIMAL.search(name):
            nc = named_carbon(name)
            if nc:
                buckets["B_minimal_carbon"] += 1
                if apply:
                    m["exchanges"] = [{"exchange": nc[1], "bigg": nc[0], "lb": -10.0, "ub": 1000.0}]
                    m["formulation"] = "defined_from_name"
                    m["formulated_from"] = "name(minimal+carbon)"
                continue
        buckets["D_needs_recipe"] += 1
        todo[name] += 1
    if apply:
        json.dump(gr, open(GR, "w"), separators=(",", ":"))
    print("buckets:", dict(buckets))
    # write the batch-of-10 to-do for bucket D
    items = [{"medium": nm, "n_records": c} for nm, c in todo.most_common()]
    batches = [items[i:i + 10] for i in range(0, len(items), 10)]
    json.dump({"n_media": len(items), "n_batches": len(batches), "batches": batches},
              open(os.path.join(ROOT, "data", "media_todo_batches.json"), "w"), indent=1)
    print("bucket D (needs recipe): %d media -> %d batches of 10 (data/media_todo_batches.json)" % (len(items), len(batches)))

if __name__ == "__main__":
    main(apply="--apply" in sys.argv)
