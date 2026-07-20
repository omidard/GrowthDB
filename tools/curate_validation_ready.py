#!/usr/bin/env python3
"""Curate GrowthDB so a genome-scale modeller can trust it for GEM validation.

Two record-level problems break validation and are fixed here:

1. RATE UNITS. Uptake/secretion rates are only usable as FBA flux bounds when they are
   biomass-SPECIFIC (mmol gDW-1 h-1). GrowthDB mixes those with volumetric (g/L/h, mM/h),
   per-cell, per-area and qualitative rates in the same list. We tag every rate with
   `unit_class` and a boolean `flux_usable` (True only for gDW-specific), and normalise the
   specific-family unit string to 'mmol/gDW/h'. A modeller / the autocurator can then filter
   to flux_usable rates and use them directly as bounds.

2. GROWTH-RATE SANITY. Some µ are impossible (>3.5 h-1 ~ <12 min doubling), non-positive, or
   inconsistent with the reported doubling time. We tag each record with `mu_qc`
   (ok / high_check / implausible_high / nonpositive / doubling_mismatch) and `mu_usable`.

Writes the flags into data/growth_records.json (canonical) AND the compact records in
data/species/*.json, and refreshes counters in data/records_index.json. Idempotent.
"""
import json, os, re, math, glob, collections

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GR = os.path.join(ROOT, "data", "growth_records.json")

def unit_class(u):
    if not u:
        return "unknown"
    s = u.lower(); c = s.replace(" ", "")
    if any(g in c for g in ("mmol/gdw/h", "mmol/gdcw/h", "mmol/gcdw/h", "mmol/gprotein/h", "mmol/gproteinh")):
        return "specific"
    if any(g in c for g in ("g/l/h", "mmol/l/h", "mm/h", "mmol/l/h")) or c in ("mm", "g/l", "mmol/l"):
        return "volumetric"
    if any(g in s for g in ("/cell", "per cell", "cell-1", "biovolume", "mg protein", "/mg protein")):
        return "per_cell"
    if any(g in s for g in ("sediment", "cm-3", "cm3", "mm3")):
        return "per_area"
    if "day" in s or "d-1" in s or "/d " in s or s.endswith("/d"):
        return "per_day"
    if any(g in s for g in ("qualitative", "not quantified", "fermented", "remineral", "sole ", "consumed over", "monod", "km_")):
        return "qualitative"
    return "other"

def curate_rate(x):
    uc = unit_class(x.get("units"))
    x["unit_class"] = uc
    x["flux_usable"] = (uc == "specific")
    if uc == "specific":
        x["units_norm"] = "mmol/gDW/h"
    return x

def mu_qc(mu, td):
    if mu is None:
        return None, None
    if mu <= 0:
        return "nonpositive", False
    if mu > 3.5:                       # fastest known prokaryotes ~4 h-1; >3.5 is almost always a unit error
        return "implausible_high", False
    if td and td > 0:
        exp = math.log(2) / td
        if abs(exp - mu) / max(exp, mu) > 0.35:
            return "doubling_mismatch", True
    if mu > 2.5:
        return "high_check", True
    return "ok", True

def main():
    gr = json.load(open(GR))
    idflags = {}
    n_up = n_up_use = n_sec = n_sec_use = 0
    muc = collections.Counter()
    for r in gr:
        for x in (r.get("uptake_rates") or []):
            curate_rate(x); n_up += 1; n_up_use += x["flux_usable"]
        for x in (r.get("secretion_rates") or []):
            curate_rate(x); n_sec += 1; n_sec_use += x["flux_usable"]
        qc, usable = mu_qc(r.get("growth_rate_per_h"), r.get("doubling_time_h"))
        r["mu_qc"] = qc
        r["mu_usable"] = usable
        if qc:
            muc[qc] += 1
        # record-level convenience flag: can this record constrain/validate a GEM?
        has_flux = any(x.get("flux_usable") for x in (r.get("uptake_rates") or []) + (r.get("secretion_rates") or []))
        r["validation_usable"] = bool(usable) or has_flux
        if r.get("id"):
            idflags[r["id"]] = {
                "mu_qc": qc, "mu_usable": usable, "validation_usable": r["validation_usable"],
                "up": {x.get("exchange"): (x.get("unit_class"), x.get("flux_usable")) for x in (r.get("uptake_rates") or [])},
                "sec": {x.get("exchange"): (x.get("unit_class"), x.get("flux_usable")) for x in (r.get("secretion_rates") or [])},
            }
    json.dump(gr, open(GR, "w"), separators=(",", ":"))
    print("uptake rates: %d, flux-usable %d (%.0f%%)" % (n_up, n_up_use, 100 * n_up_use / max(1, n_up)))
    print("secretion rates: %d, flux-usable %d (%.0f%%)" % (n_sec, n_sec_use, 100 * n_sec_use / max(1, n_sec)))
    print("mu_qc:", dict(muc))

    # propagate flags into the compact species-shard records (by id)
    touched = 0
    for f in glob.glob(os.path.join(ROOT, "data", "species", "*.json")):
        d = json.load(open(f)); changed = False
        for r in (d.get("growth_records") or []):
            fl = idflags.get(r.get("id"))
            if not fl:
                continue
            r["mu_qc"] = fl["mu_qc"]; r["mu_usable"] = fl["mu_usable"]; r["validation_usable"] = fl["validation_usable"]
            for x in (r.get("up") or []):
                uc = unit_class(x.get("units")); x["unit_class"] = uc; x["flux_usable"] = (uc == "specific")
            for x in (r.get("sec") or []):
                uc = unit_class(x.get("units")); x["unit_class"] = uc; x["flux_usable"] = (uc == "specific")
            changed = True
        if changed:
            json.dump(d, open(f, "w"), separators=(",", ":")); touched += 1
    print("propagated flags into %d species shards" % touched)

    ri_path = os.path.join(ROOT, "data", "records_index.json")
    if os.path.exists(ri_path):
        ri = json.load(open(ri_path))
        ri["n_flux_usable_uptake"] = n_up_use
        ri["n_flux_usable_secretion"] = n_sec_use
        ri["n_mu_implausible"] = muc.get("implausible_high", 0) + muc.get("nonpositive", 0)
        json.dump(ri, open(ri_path, "w"))
        print("records_index updated")

if __name__ == "__main__":
    main()
