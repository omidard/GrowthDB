#!/usr/bin/env python3
"""Add a standard strain identifier (`strain_std`) to every GrowthDB record.

GrowthDB stores the strain as free text ("BW25113 (Keio collection wild-type
background)", "K-12 MG1655", "ATCC 14266 (wild type)"). Validation and cross-
referencing need a *comparable* token, so this derives a normalised standard strain
id from that prose (and from the organism name when the strain field is empty):

  priority: culture-collection accession (ATCC/DSM/NCTC/…) > str./substr. designation
            > lab designation (MG1655, BW25113, N16961, …) > K-12

It writes `strain_std` into:
  data/growth_records.json      (canonical records)
  data/species/*.json           (per-species record copies)
  data/strains.json             (NEW: standard-strain-id catalog: token -> species, counts)
and refreshes the strain counters in data/records_index.json.

Idempotent — safe to re-run. Run from the repo root: python3 tools/standardize_strains.py
"""
import json, os, re, glob

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
_CC = r"ATCC|DSMZ|DSM|NCTC|CCUG|JCM|NBRC|IFO|CECT|LMG|CIP|NCIMB|BCRC|KCTC|NRRL|CGMCC|MCCC|PCC|UTEX|CBS|KACC|VPI|NCDO|NCFB|BCCM|KCCM"

def strain_std(text):
    if not text:
        return None
    t = str(text)
    m = re.search(r"\b(" + _CC + r")\s*[-: ]?\s*(\d+[A-Za-z]?)\b", t, re.I)
    if m:
        return (m.group(1) + m.group(2)).upper()
    m = re.search(r"substr\.?\s+([A-Za-z0-9][A-Za-z0-9\-]{1,})", t, re.I) or re.search(r"(?:str\.?|strain)\s+([A-Za-z0-9][A-Za-z0-9\-]{1,})", t, re.I)
    if m and re.search(r"\d", m.group(1)):
        return re.sub(r"[^A-Z0-9]", "", m.group(1).upper())
    m = re.search(r"\b([A-Z]{1,4}\d{2,6}[A-Za-z]?)\b", t)
    if m:
        return m.group(1).upper()
    if re.search(r"\bK-?12\b", t, re.I):
        return "K12"
    return None

def load(p):
    with open(p) as fh:
        return json.load(fh)

def dump(obj, p):
    with open(p, "w") as fh:
        json.dump(obj, fh, separators=(",", ":"))

# 1) canonical growth_records.json -> strain_std + id map + catalog
gr_path = os.path.join(DATA, "growth_records.json")
records = load(gr_path)
id2std, catalog = {}, {}
n_std = 0
for r in records:
    std = strain_std(r.get("strain")) or strain_std(r.get("organism"))
    r["strain_std"] = std
    if std:
        n_std += 1
        rid = r.get("id")
        if rid:
            id2std[rid] = std
        sp = r.get("gtdb_species") or r.get("species") or r.get("organism") or "unknown"
        c = catalog.setdefault(std, {"n_records": 0, "species": set()})
        c["n_records"] += 1
        c["species"].add(sp)
dump(records, gr_path)
print("growth_records.json: %d / %d records now carry strain_std" % (n_std, len(records)))

# 2) per-species files
sp_files = glob.glob(os.path.join(DATA, "species", "*.json"))
touched = 0
for f in sp_files:
    d = load(f)
    grs = d.get("growth_records")
    if not isinstance(grs, list):
        continue
    changed = False
    for r in grs:
        std = id2std.get(r.get("id")) or strain_std(r.get("strain"))
        if r.get("strain_std") != std:
            r["strain_std"] = std
            changed = True
    if changed:
        dump(d, f)
        touched += 1
print("species files updated: %d / %d" % (touched, len(sp_files)))

# 3) strains.json catalog
cat_out = {k: {"n_records": v["n_records"], "n_species": len(v["species"]),
               "species": sorted(v["species"])[:12]} for k, v in catalog.items()}
dump({"count": len(cat_out), "note": "Standard strain identifiers derived from GrowthDB strain prose; token -> record/species counts.",
      "strains": cat_out}, os.path.join(DATA, "strains.json"))
print("strains.json: %d distinct standard strain ids" % len(cat_out))

# 4) refresh records_index.json counters
idx_path = os.path.join(DATA, "records_index.json")
if os.path.exists(idx_path):
    idx = load(idx_path)
    idx["n_with_strain_std"] = n_std
    idx["n_standard_strains"] = len(cat_out)
    dump(idx, idx_path)
    print("records_index.json: n_with_strain_std=%d, n_standard_strains=%d" % (n_std, len(cat_out)))
