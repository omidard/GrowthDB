#!/usr/bin/env python3
"""Formulate a small set of WELL-KNOWN named standard media whose published composition is
unambiguous, to BiGG exchanges. These are textbook/standard recipes (not paper-specific), so
they can be formulated precisely without mining. Each recipe below is transcribed from the
canonical source; only recipes I can state with confidence are included — anything ambiguous
is deliberately omitted and left for source mining.

For FBA, the exchange set is what matters: the defined salt background (unlimited) + the
carbon/energy source (organic capped at -10; inorganic/gas unlimited). Exact trace-metal
grammage does not change feasibility.

Writes medium.exchanges + medium.formulation='curated_named_recipe' in place (dry-run unless --write).
"""
import json, os, re, sys, collections

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GR = os.path.join(ROOT, "data", "growth_records.json")
WRITE = "--write" in sys.argv
CAP = 10.0

# standard defined mineral background (unlimited uptake), reused by several recipes
MB = ["pi", "so4", "cl", "na1", "k", "mg2", "ca2", "fe2", "fe3", "mn2", "zn2", "cu2", "cobalt2", "mobd", "h2o", "h"]


def ex(bigg_lbs):
    """dict bigg->lb  ->  exchange list"""
    return [{"exchange": "EX_" + b + "_e", "bigg": b, "lb": float(lb), "ub": 1000.0} for b, lb in bigg_lbs.items()]


def recipe_bg11():
    # Rippka et al. 1979 BG-11: nitrate as N, CO2/bicarbonate as carbon (photoautotroph), no organic C.
    d = {b: -1000.0 for b in MB}
    d.update({"no3": -1000.0, "co2": -1000.0, "hco3": -1000.0})   # N = nitrate; C = CO2/HCO3
    return d, "photoautotroph"


def recipe_pro99():
    # Moore et al. Pro99 (marine cyanobacterium): artificial-seawater base, NH4Cl as N, phosphate,
    # trace metals; carbon = CO2/HCO3 (photoautotroph).
    d = {b: -1000.0 for b in MB}
    d.update({"nh4": -1000.0, "co2": -1000.0, "hco3": -1000.0})
    return d, "photoautotroph"


def recipe_sistrom():
    # Sistrom 1962 minimal (Rhodobacter): succinate as principal carbon; (NH4)2SO4 as N; small
    # glutamate + aspartate; phosphate; Mg/Ca; NTA chelator (omitted); trace metals + vitamins.
    d = {b: -1000.0 for b in MB}
    d.update({"nh4": -1000.0, "succ": -CAP, "glu__L": -1.0, "asp__L": -1.0})
    return d, "defined"


def recipe_atgn():
    # AT minimal salts + Glucose + NH4Cl (Agrobacterium ATGN): glucose carbon, ammonium N.
    d = {b: -1000.0 for b in MB}
    d.update({"nh4": -1000.0, "glc__D": -CAP})
    return d, "defined"


# name pattern -> recipe builder. Order matters (first match wins).
RECIPES = [
    (re.compile(r"\bBG[- ]?11\b", re.I), recipe_bg11, "BG-11 (Rippka 1979)"),
    (re.compile(r"pro[- ]?99", re.I), recipe_pro99, "Pro99 marine cyanobacterial medium (Moore et al.)"),
    (re.compile(r"sistrom", re.I), recipe_sistrom, "Sistrom's minimal medium (Sistrom 1962)"),
    (re.compile(r"\bATGN\b", re.I), recipe_atgn, "ATGN — AT salts + glucose + NH4Cl"),
]


def main():
    gr = json.load(open(GR))
    done = collections.Counter(); by = collections.Counter()
    for r in gr:
        m = r.get("medium") or {}
        if m.get("media_id") or m.get("exchanges") or m.get("medium_type") != "defined":
            continue
        if not (r.get("growth_rate_per_h") is not None or r.get("uptake_rates") or r.get("secretion_rates")):
            continue
        name = m.get("canonical_name") or m.get("description") or ""
        for pat, build, label in RECIPES:
            if pat.search(name):
                d, status = build()
                m["exchanges"] = ex(d)
                m["formulation"] = "curated_named_recipe"
                m["formulated_from"] = label
                if status == "photoautotroph":
                    m["formulation"] = "curated_named_recipe_photoautotroph"
                done["formulated"] += 1; by[label] += 1
                break
    print("curated named-recipe formulation:")
    for k, v in done.most_common():
        print(f"  {v:3d}  {k}")
    print("\nby recipe:")
    for k, v in by.most_common():
        print(f"  {v:3d}  {k}")
    if WRITE:
        json.dump(gr, open(GR, "w"), separators=(",", ":"))
        print("\nWROTE exchanges into", GR)
    else:
        print("\n(dry run — pass --write to persist)")


if __name__ == "__main__":
    main()
