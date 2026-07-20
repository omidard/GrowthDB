#!/usr/bin/env python3
"""Separate substrate phenotypes by EVIDENCE TYPE so a modeller uses the right ones.

GrowthDB conflated three very different assays into one consensus per (organism, substrate):
  * GROWTH / assimilation / Biolog PM   -> the organism GROWS on the substrate as a sole
    C/N/energy source  => validates model biomass feasibility with that substrate.
  * "builds acid from" / fermentation    -> the organism produces ACID from the substrate
    in a rich base  => this is a SECRETION/fermentation phenotype, NOT growth.
  * hydrolysis / reduction / degradation -> an enzyme/activity assay, not growth at all.

We tag every phenotype with `ptype` (growth / acid_production / enzyme / growth+acid / other)
across phenotypes_full.json.gz, the compact phenotypes.json.gz, phenotypes_summary.json and
the per-species shards, so the autocurator can filter to ptype=='growth' for growth validation
and read ptype=='acid_production' as fermentation/secretion evidence. Idempotent.
"""
import json, os, gzip, glob, collections

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
GROWTH = ("assimil", "utiliz", "carbon source", "energy source", "nitrogen source", "phosphorus source",
          "sulfur source", "sole", "respiration", "oxid", "biolog", "growth", "grows", "substrate util")
ACID = ("builds acid", "acid from", "gas from", "ferment")
ENZ = ("hydrol", "degrad", "reduc", "electron", "activity", "assay")

def ptype(kinds):
    ks = [(k or "").lower() for k in (kinds or [])]
    g = any(any(s in k for s in GROWTH) for k in ks)
    a = any(any(s in k for s in ACID) for k in ks)
    e = any(any(s in k for s in ENZ) for k in ks)
    if g and a:
        return "growth+acid"
    if a:
        return "acid_production"
    if g:
        return "growth"
    if e:
        return "enzyme"
    return "other"

# 1) full (has kinds) -> add ptype
full = json.load(gzip.open(os.path.join(DATA, "phenotypes_full.json.gz"), "rt"))
pt_by_key = {}
cnt = collections.Counter()
for r in full:
    p = ptype(r.get("kinds"))
    r["ptype"] = p
    cnt[p] += 1
    pt_by_key[(r["organism"], (r["substrate"] or "").lower(), r.get("exchange") or "")] = p
with gzip.open(os.path.join(DATA, "phenotypes_full.json.gz"), "wt") as f:
    json.dump(full, f, separators=(",", ":"))
print("phenotypes_full: ptype", dict(cnt))

# 2) compact -> add 't'
comp = json.load(gzip.open(os.path.join(DATA, "phenotypes.json.gz"), "rt"))
for r in comp:
    r["t"] = pt_by_key.get((r["o"], (r["s"] or "").lower(), r.get("e") or ""), "other")
with gzip.open(os.path.join(DATA, "phenotypes.json.gz"), "wt") as f:
    json.dump(comp, f, separators=(",", ":"))

# 3) summary -> by_ptype and a growth-only mapped count
summ = json.load(open(os.path.join(DATA, "phenotypes_summary.json")))
summ["by_ptype"] = dict(cnt)
summ["n_growth_phenotypes"] = cnt["growth"]
summ["n_acid_production"] = cnt["acid_production"]
json.dump(summ, open(os.path.join(DATA, "phenotypes_summary.json"), "w"), separators=(",", ":"))

# 4) species shards -> add ptype to each phenotype (browser reads these)
touched = 0
for fp in glob.glob(os.path.join(DATA, "species", "*.json")):
    d = json.load(open(fp)); changed = False
    for r in (d.get("phenotypes") or []):
        p = ptype(r.get("kinds"))
        if r.get("ptype") != p:
            r["ptype"] = p; changed = True
    # regroup phenotype_groups by ptype too, if present
    for g in (d.get("phenotype_groups") or []):
        for r in (g.get("calls") or []):
            r["ptype"] = ptype(r.get("kinds"))
    if changed:
        json.dump(d, open(fp, "w"), separators=(",", ":")); touched += 1
print("added ptype to %d species shards" % touched)
