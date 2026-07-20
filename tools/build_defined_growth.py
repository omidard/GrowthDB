#!/usr/bin/env python3
"""Extract MEDIA-GROUNDED substrate-growth evidence from the growth records.

A distinct data type from Biolog PM plates and from qS/qP rate measurements: when a paper
grows an organism on a DEFINED medium whose carbon/energy source is a specific compound
(e.g. "M9 minimal (maltose)" — M9 salts, no glucose, maltose as the sole C source), that is
direct evidence the organism GROWS ON that substrate on a fully-formulatable medium — exactly
what you need to validate a GEM (set the medium, check biomass feasibility).

For each record on a defined medium we emit one growth-on-substrate call:
  {organism, gtdb_species, strain, strain_std, substrate, bigg, exchange, base_family,
   base_medium, media_id, growth_rate_per_h, doubling_time_h, oxygen, temperature_C, pH,
   citation, doi, record_id}
Output: data/defined_medium_growth.json (+ per-species shard section, + index counters).
"""
import json, os, re, sys, collections

sys.path.insert(0, "/data/media_curate/tools")
from map_metabolite import Mapper

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GR = "/data/modelseed_cache/gdb_growth_records.json"
BIGG = json.load(open("/data/media_curate/tools/bigg_metabolite_dict.json"))
MAP = Mapper()

# complex-medium markers: presence => NOT a defined medium (skip)
COMPLEX = re.compile(r"yeast extract|peptone|tryptone|casitone|casamino|casein|beef|meat extract|"
                     r"\bbroth\b|\bLB\b|\bBHI\b|\bTSB\b|\bMRS\b|rumen|serum|blood|hydrolysate|digest|"
                     r"lysate|infusion|molasses|corn steep|malt extract|proteose", re.I)
# ingredients that map to BiGG but are NOT the carbon/energy source (ions, buffers, vitamins,
# and common reducing agents / supplements that are present at trace levels, not as the C source)
NOT_CARBON = {"pi", "so4", "cl", "na1", "k", "nh4", "mg2", "ca2", "fe2", "fe3", "mn2", "zn2", "cu2",
              "cobalt2", "mobd", "ni2", "sel", "h2o", "h", "co2", "hco3", "o2", "no3", "no2", "n2",
              "cbl1", "btn", "thm", "ribflv", "pnto__R", "fol", "pydxn", "nac", "4abz",
              "lipoate", "edta", "hepes", "mops", "pipes", "tris", "cynt", "tungs",
              "cys__L", "cys__D", "Lcystin", "gthrd", "gthox", "ascb__L", "so3", "tsul", "dtt"}
REDUCING = {"cys__L", "cys__D", "gthrd", "ascb__L", "so3", "tsul"}
DEFINED_HINT = re.compile(r"\bM9\b|\bM63\b|\bMOPS\b|minimal|defined|synthetic|mineral|Gutnick|"
                          r"\bNBS\b|CGXII|Davis|Neidhardt|chemically|BG11|Wolfe", re.I)

def base_family(name):
    for pat, fam in [(r"\bM9\b", "M9"), (r"\bM63\b", "M63"), (r"\bMOPS\b", "MOPS"), (r"Gutnick", "Gutnick"),
                     (r"CGXII", "CGXII"), (r"\bNBS\b", "NBS"), (r"BG11", "BG11"), (r"Davis", "Davis"),
                     (r"Wolfe|mineral", "mineral"), (r"minimal|defined|synthetic", "defined-minimal")]:
        if re.search(pat, name or "", re.I):
            return fam
    return "defined"

def carbon_c(bigg):
    f = (BIGG.get(bigg, {}).get("xrefs", {}) or {}).get("formula") or BIGG.get(bigg, {}).get("formula") or ""
    m = re.search(r"C(\d*)", f)
    if not m:
        return 0
    return int(m.group(1) or 1) if f.startswith("C") or "C" in f[:2] else (int(m.group(1) or 1))

def name_carbon(medium_name):
    """Carbon named in the medium string, e.g. 'M9 minimal (maltose)' or 'MOPS + 4 g/L xylose'."""
    out = []
    for m in re.findall(r"\(([^)]+)\)", medium_name or ""):
        out.append(m)
    for m in re.findall(r"(?:\+|with|on)\s+(?:[\d.]+\s*(?:g/l|g/L|%|mm|mM)?\s*)?([A-Za-z][A-Za-z\- ]{2,30})", medium_name or ""):
        out.append(m)
    return out

def is_carbon_hit(h):
    if not h or not h.get("in_biggr") or not h.get("bigg_metabolite"):
        return None
    b = h["bigg_metabolite"]
    if b in NOT_CARBON:
        return None
    if carbon_c(b) < 2:      # need an actual carbon skeleton (>=2 C) to be a C/energy source
        return None
    return b

def to_gl(amount, unit):
    """Rough normalisation of a concentration to g/L for picking the dominant C source."""
    try:
        v = float(re.findall(r"[\d.]+", str(amount))[0])
    except Exception:
        return None
    u = (unit or "").lower()
    if "%" in u:
        return v * 10.0            # % w/v -> g/L
    if u.startswith("g/l") or u == "g/l" or "g l" in u or u == "g/l":
        return v
    if "mg/l" in u or "mg/ml" in u:
        return v / 1000.0 if "mg/l" in u else v
    if "mm" in u:
        return v * 0.15            # ~molar mass proxy for a small sugar/acid (approx, for ranking only)
    if "g/kg" in u:
        return v
    return None

def extract(rec):
    med = rec.get("medium") or {}
    name = (med.get("canonical_name") or med.get("description") or "")
    comp = med.get("composition") or []
    comp_text = " ".join((c.get("name") or "") for c in comp)
    if COMPLEX.search(name) or COMPLEX.search(comp_text):
        return []
    if not (DEFINED_HINT.search(name) or comp):
        return []
    carbons = []           # (bigg, display_name, exchange, confidence)
    seen = set()
    def add(bigg, disp, ex, conf):
        if bigg and bigg not in seen:
            seen.add(bigg); carbons.append((bigg, disp, ex, conf))
    # 1) carbon named explicitly in the medium string -> highest confidence (e.g. "M9 minimal (maltose)")
    for cand in name_carbon(name):
        h = MAP.map(name=cand.strip()); b = is_carbon_hit(h)
        if b:
            add(b, cand.strip(), h["exchange"], "named")
    # 2) composition organics, ranked by concentration; drop trace reducing agents/supplements
    conc = []
    for c in comp:
        nm = (c.get("name") or "").strip()
        if not nm:
            continue
        h = MAP.map(name=nm); b = is_carbon_hit(h)
        if not b or b in seen:
            continue
        gl = to_gl(c.get("amount"), c.get("unit"))
        conc.append((b, nm, h["exchange"], gl if gl is not None else -1))
    if conc:
        mx = max((g for *_ , g in conc if g >= 0), default=-1)
        for b, nm, ex, g in conc:
            if b in REDUCING and (mx > 0 and (g < 0 or g < 0.5 * mx)):
                continue                       # reducing agent at trace level -> not the C source
            if mx > 0 and g >= 0 and g < 0.15 * mx and not any(x[0] == b for x in carbons):
                continue                       # trace organic vs the bulk substrate -> skip
            add(b, nm, ex, "dominant" if (mx <= 0 or g >= 0.5 * mx) else "listed")
    # 3) fall back to the Madin/species carbon_substrates only if nothing else and it's clearly defined
    if not carbons and DEFINED_HINT.search(name):
        for n in (rec.get("carbon_substrates") or []):
            h = MAP.map(name=str(n).replace("_", " ")); b = is_carbon_hit(h)
            if b:
                add(b, str(n), h["exchange"], "trait")
    if not carbons:
        return []
    cond = rec.get("conditions") or {}
    prov = rec.get("provenance") or {}
    out = []
    for b, cand, ex, conf in carbons[:6]:
        out.append({
            "organism": rec.get("organism"), "gtdb_species": rec.get("gtdb_species"),
            "strain": rec.get("strain"), "substrate": cand, "bigg": b, "exchange": ex, "confidence": conf,
            "base_family": base_family(name), "base_medium": name[:80], "media_id": med.get("media_id"),
            "media_url": med.get("media_url"),
            "growth_rate_per_h": rec.get("growth_rate_per_h"), "doubling_time_h": rec.get("doubling_time_h"),
            "oxygen": cond.get("oxygen"), "temperature_C": cond.get("temperature_C"), "pH": cond.get("pH"),
            "citation": (prov.get("citation") or "")[:180], "doi": prov.get("doi"), "record_id": rec.get("id"),
        })
    return out

def main():
    gr = json.load(open(GR))
    records = []
    for rec in gr:
        records.extend(extract(rec))
    # dedupe identical (organism, substrate, base_family, record_id)
    seen, uniq = set(), []
    for r in records:
        k = (r["organism"], r["bigg"], r["base_family"], r["record_id"])
        if k in seen:
            continue
        seen.add(k)
        uniq.append(r)
    OUT = os.path.join(ROOT, "data", "defined_medium_growth.json")
    n_mapped_media = sum(1 for r in uniq if r["media_id"])
    by_sub = collections.Counter(r["bigg"] for r in uniq)
    by_base = collections.Counter(r["base_family"] for r in uniq)
    summary = {"count": len(uniq), "n_species": len({r["gtdb_species"] or r["organism"] for r in uniq}),
               "n_substrates": len(by_sub), "n_with_media_id": n_mapped_media,
               "top_substrates": by_sub.most_common(25), "by_base_family": dict(by_base)}
    json.dump({"summary": summary, "records": uniq}, open(OUT, "w"), separators=(",", ":"))
    print("defined_medium_growth: %d growth-on-substrate calls | %d species | %d substrates | %d with a linked media_id"
          % (len(uniq), summary["n_species"], summary["n_substrates"], n_mapped_media))
    print("top substrates:", by_sub.most_common(12))
    print("by base family:", dict(by_base))
    print("media formulation gap: %d of %d need an exact medium formulated/linked" % (len(uniq) - n_mapped_media, len(uniq)))

    # ---- inject a per-species `carbon_growth` section into the species shards (presented separately) ----
    def norm_sp(r):
        s = r.get("gtdb_species") or r.get("organism") or ""
        p = s.replace("[", "").replace("]", "").split()
        return (p[0] + " " + p[1]) if len(p) >= 2 else s
    bysp = collections.defaultdict(list)
    for r in uniq:
        bysp[norm_sp(r)].append(r)
    idx_path = os.path.join(ROOT, "data", "species_index.json")
    sidx = json.load(open(idx_path)) if os.path.exists(idx_path) else None
    rows = (sidx.get("species") if isinstance(sidx, dict) else sidx) or []
    slug_of = {r.get("s"): r.get("slug") for r in rows if isinstance(r, dict)}
    touched = 0
    for sp, recs in bysp.items():
        slug = slug_of.get(sp)
        f = os.path.join(ROOT, "data", "species", (slug or "") + ".json")
        if not slug or not os.path.exists(f):
            continue
        d = json.load(open(f))
        d["carbon_growth"] = [{"substrate": r["substrate"], "bigg": r["bigg"], "exchange": r["exchange"],
                                "base": r["base_family"], "medium": r["base_medium"], "media_id": r["media_id"],
                                "media_url": r["media_url"], "mu": r["growth_rate_per_h"], "oxygen": r["oxygen"],
                                "confidence": r["confidence"], "citation": r["citation"], "doi": r["doi"]} for r in recs]
        json.dump(d, open(f, "w"), separators=(",", ":"))
        touched += 1
    print("injected carbon_growth into %d species shards" % touched)
    # index counters
    ri_path = os.path.join(ROOT, "data", "records_index.json")
    if os.path.exists(ri_path):
        ri = json.load(open(ri_path))
        ri["n_defined_growth"] = len(uniq)
        ri["n_defined_growth_species"] = summary["n_species"]
        json.dump(ri, open(ri_path, "w"))
        print("records_index: n_defined_growth=%d" % len(uniq))

if __name__ == "__main__":
    main()
