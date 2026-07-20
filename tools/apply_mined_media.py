#!/usr/bin/env python3
"""Formulate the mined medium compositions -> BiGG exchanges, link them to the GrowthDB
records (by PMC + medium name), and emit a media_to_add.json for the Media DB.

The LLM extracted the ingredient list from each paper; here the deterministic (validated)
formula-salt parser + Media DB Mapper map it to exchanges — so precision stays with code.
"""
import json, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import formulate_media as FM

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GR = os.path.join(ROOT, "data", "growth_records.json")
MINED = json.load(open("/tmp/mined_media.json"))

def norm(s):
    return re.sub(r"[^a-z0-9]", "", (s or "").lower())

def toks(s):
    return set(re.findall(r"[a-z0-9]{3,}", (s or "").lower()))

def formulate(comp):
    ex = {}
    unmapped = 0
    for c in comp:
        got, cls = FM.ingredient_exchanges((c.get("name") or "").strip())
        if cls in ("unmapped", "vitamin_unlisted", "base_ref"):
            unmapped += 1
        for bigg, is_c in got:
            lb = -FM.CARBON_CAP if (is_c and bigg not in FM.INORG_ION) else -1000.0
            if bigg not in ex or lb < ex[bigg]:
                ex[bigg] = lb
    return [{"exchange": "EX_" + b + "_e", "bigg": b, "lb": ex[b], "ub": 1000.0} for b in ex], unmapped

gr = json.load(open(GR))
# index records by pmc
bypmc = {}
for r in gr:
    m = re.search(r"(PMC\d+)", r.get("id", ""))
    if m:
        bypmc.setdefault(m.group(1), []).append(r)

linked = 0
to_add = []
for mm in MINED:
    pmc = mm["pmc"]
    exch, unmapped = formulate(mm["composition"])
    if not exch:
        continue
    is_def = mm.get("is_defined")
    status = "partial_complex" if not is_def else ("partial_unmapped" if unmapped else "defined")
    mt = toks(mm["name"])
    # link to the pmc's records whose medium name best-matches and is not yet resolved
    cand = []
    for r in bypmc.get(pmc, []):
        med = r.get("medium") or {}
        if med.get("media_id") or med.get("exchanges") or med.get("composition"):
            continue
        rn = (med.get("canonical_name") or med.get("description") or "")
        ov = len(mt & toks(rn))
        if ov >= 2 or norm(mm["name"])[:12] and norm(mm["name"])[:12] in norm(rn):
            cand.append((ov, r))
    for ov, r in cand:
        med = r["medium"]
        med["exchanges"] = exch
        med["formulation"] = status
        med["formulated_from"] = "mined recipe (paper Methods)"
        med["composition"] = mm["composition"]
        linked += 1
    # add to Media DB regardless (reusable), tagged with the paper
    to_add.append({"medium_name": mm["name"], "pmcid": pmc, "composition": mm["composition"],
                   "is_defined": bool(is_def), "status": status})

json.dump(gr, open(GR, "w"), separators=(",", ":"))
json.dump({"pending_media_for_media_repo": to_add}, open(os.path.join(ROOT, "data", "mined_media_to_add.json"), "w"), indent=1)
print("mined media formulated: %d | GrowthDB records linked: %d | media_to_add for Media DB: %d"
      % (len(MINED), linked, len(to_add)))
