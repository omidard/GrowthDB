#!/usr/bin/env python3
"""Formulate the UNRESOLVED, chemically-DEFINED minimal media to BiGG exchanges.

Insight for FBA validation: a defined minimal medium's exchange set is deterministic — the
standard mineral base (unlimited) + the named carbon/nitrogen/energy source(s) (capped). The
exact trace-mineral grammage does not change FBA feasibility, so we do NOT need the paper's
weigh-out; we need the correct exchange set, which the medium's own description gives us.

PRECISION RULE (validation-critical): only formulate when a carbon/energy source is confidently
identified from the medium's name or its enumerated composition. If none is found (e.g. an
amino-acid-supplemented CDM whose AA list isn't given, or a bare "defined medium"), the record
is LEFT UNRESOLVED and flagged for source mining — never faked.

Writes medium.exchanges + medium.formulation='defined_minimal_reconstructed' in place
(dry-run unless --write).
"""
import json, os, re, sys, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import formulate_media as FM

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GR = os.path.join(ROOT, "data", "growth_records.json")
WRITE = "--write" in sys.argv

# standard defined-minimal mineral base (unlimited uptake) — the FBA-relevant salt background
MINERAL_BASE = ["pi", "so4", "cl", "na1", "k", "nh4", "mg2", "ca2", "fe2", "fe3", "mn2", "zn2",
                "cu2", "cobalt2", "mobd", "ni2", "h2o", "h"]

# carbon/energy sources we will look for by name (mapped through the Media DB mapper for the BiGG id)
CARBON_WORDS = [
    "glucose", "fructose", "galactose", "mannose", "xylose", "arabinose", "ribose", "rhamnose",
    "sucrose", "maltose", "lactose", "cellobiose", "trehalose", "raffinose", "melibiose",
    "glycerol", "sorbitol", "mannitol", "inositol", "gluconate", "glucuronate",
    "acetate", "lactate", "pyruvate", "malate", "fumarate", "succinate", "citrate", "oxaloacetate",
    "propionate", "butyrate", "formate", "benzoate", "phenol", "phenanthrene", "dibenzofuran",
    "methanol", "ethanol", "propanol", "butanol", "methane", "2,3-butanediol", "butanediol",
    "glutamate", "aspartate", "alanine", "serine", "glycine", "proline", "arginine",
    "n-acetyl-d-glucosamine", "n-acetylglucosamine", "acetylglucosamine", "glucosamine",
    "urea", "methylamine", "trimethylamine", "glycolate", "glyoxylate", "tartrate", "mannose",
]
NITRO_WORDS = {"nitrate": "no3", "nitrite": "no2", "ammonia": "nh4", "ammonium": "nh4", "urea": None}
GAS_WORDS = {"h2": "h2", "hydrogen": "h2", "co2": "co2", "carbon dioxide": "co2", "co": "co",
             "syngas": None, "methane": "ch4", "n2": None}


def carbon_from_text(t):
    """Identified (bigg, is_carbon) exchanges from free text, precise-only."""
    low = " " + t.lower() + " "
    found = {}
    # longest-first, consuming each match so a shorter substring (e.g. 'glucosamine' inside
    # 'n-acetyl-d-glucosamine') can't double-map
    for w in sorted(CARBON_WORDS, key=len, reverse=True):
        pat = r"[^a-z]" + re.escape(w) + r"[^a-z]"
        if re.search(pat, low):
            got, cls = FM.ingredient_exchanges(w)
            for b, isc in got:
                if isc:
                    found[b] = True
            low = re.sub(pat, "  ", low)
    # explicit "X + Y + Z" carbon mixtures
    for chunk in re.findall(r"([a-z]+(?:\s*\+\s*[a-z]+)+)", low):
        for w in re.split(r"\+", chunk):
            w = w.strip()
            got, cls = FM.ingredient_exchanges(w)
            for b, isc in got:
                if isc:
                    found[b] = True
    # gases (autotrophic / methanogenic / gas-fermentation); (?![a-z0-9]) stops 'co' matching 'co2'
    for g, b in GAS_WORDS.items():
        if b and re.search(r"[^a-z]" + re.escape(g) + r"(?![a-z0-9])", low):
            found[b] = (b in ("h2", "co", "ch4"))       # H2/CO/CH4 are energy/C sources; CO2 is inorganic C
    return found


def formulate_record(m):
    name = m.get("canonical_name") or m.get("description") or ""
    carbons = {}
    # 1) enumerated composition (if any)
    for c in (m.get("composition") or []):
        got, cls = FM.ingredient_exchanges((c.get("name") or "").strip())
        for b, isc in got:
            if isc and b not in FM.INORG_ION:
                carbons[b] = True
            elif cls in ("salt", "salt_word") and b not in carbons:
                carbons.setdefault(b, False)   # a defined salt named in composition -> add unlimited
    # 2) the medium NAME
    carbons.update({b: v for b, v in carbon_from_text(name).items()})
    # need at least one genuine carbon/energy source, else do not formulate
    real_c = [b for b, isc in carbons.items() if isc]
    if not real_c:
        return None, "no_carbon_identified"
    ex = {b: -1000.0 for b in MINERAL_BASE}
    # nitrate/nitrite alt-N from name
    for w, b in NITRO_WORDS.items():
        if b and re.search(r"[^a-z]" + re.escape(w) + r"[^a-z]", " " + name.lower() + " "):
            ex[b] = -1000.0
    for b, isc in carbons.items():
        ex[b] = -FM.CARBON_CAP if isc else -1000.0
    exch = [{"exchange": "EX_" + b + "_e", "bigg": b, "lb": ex[b], "ub": 1000.0} for b in ex]
    return exch, "defined_minimal_reconstructed"


def main():
    gr = json.load(open(GR))
    done = collections.Counter(); skipped_names = []
    for r in gr:
        m = r.get("medium") or {}
        if m.get("media_id") or m.get("exchanges") or m.get("medium_type") != "defined":
            continue
        if not (r.get("growth_rate_per_h") is not None or r.get("uptake_rates") or r.get("secretion_rates")):
            continue
        exch, status = formulate_record(m)
        if exch:
            m["exchanges"] = exch
            m["formulation"] = status
            m["formulated_from"] = "defined-minimal reconstruction (standard mineral base + named C/N/energy source)"
            done["formulated"] += 1
        else:
            done["left_for_mining"] += 1
            nm = m.get("canonical_name") or m.get("description") or ""
            if nm:
                skipped_names.append(nm[:70])
    print("defined-minimal formulation:")
    for k, v in done.most_common():
        print(f"  {v:4d}  {k}")
    print("\n--- left for mining (no carbon confidently identified), distinct ---")
    for nm, n in collections.Counter(skipped_names).most_common(25):
        print(f"  {n:2d}  {nm}")
    if WRITE:
        json.dump(gr, open(GR, "w"), separators=(",", ":"))
        print("\nWROTE exchanges into", GR)
    else:
        print("\n(dry run — pass --write to persist)")


if __name__ == "__main__":
    main()
