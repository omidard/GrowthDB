#!/usr/bin/env python3
"""Classify every GrowthDB medium by TYPE so GEM validation is honest about what a record can
support. A defined-medium FBA validation is only meaningful for a chemically defined medium;
an in-vivo/environmental/complex-undefined medium cannot give an exact exchange set, and the
autocurator must SAY so rather than silently fail or fabricate a recipe.

medium_type:
  defined              chemically defined / minimal / mineral-salts (formulable -> exchanges)
  complex_undefined    peptone/yeast-extract/broth/commercial (undefinable pools; partial at best)
  in_vivo              grown in a host / tissue / intracellular (not a lab medium)
  environmental        seawater / soil / freshwater / field sample (not a lab medium)
  unstated             no medium reported

`formulable` is True only for `defined`. Writes medium.medium_type + medium.formulable in place
(dry-run prints the breakdown without writing unless --write).
"""
import json, os, re, sys, collections

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GR = os.path.join(ROOT, "data", "growth_records.json")
WRITE = "--write" in sys.argv

IN_VIVO = re.compile(r"\b(in[ -]?vivo|intracellular|within (the )?cytoplasm|cytosol of|murine|mouse|rat|human|host|xylem sap|phyllosphere|rhizosphere|infection|epithelial cell|macrophage|bladder|urinary tract|spleen|liver|gut (lumen|content)|caecal|cecal|rumen fluid|tissue|in planta|nodule)\b", re.I)
ENVIRON = re.compile(r"\b(seawater|sea water|reservoir|water column|epilimnion|hypolimnion|\bsoil\b|biocrust|freshwater|fresh water|sediment|\blake\b|ocean|sargasso|pond|river|groundwater|aquifer|field sample|natural (water|sample)|rain(water)?|brine|hot spring|hydrothermal|estuar|wastewater|activated sludge)\b", re.I)
COMPLEX = re.compile(r"\b(peptone|yeast extract|tryptone|casitone|casamino|casein|meat extract|beef extract|brain[- ]heart|\bBHI\b|\bLB\b|luria|nutrient broth|nutrient[- ]rich|rich medium|rich broth|chopped meat|marine broth|difco|BD BBL|oxoid|sputum|\bblood\b|\bserum\b|\bmilk\b|\bwort\b|molasses|corn steep|\bTSB\b|tryptic soy|MRS|\bYEP\b|\bYPD\b|complex medium|undefined medium|\bYCFA\b|gifu|GAM broth|reinforced clostridial|todd[- ]hewitt|mueller[- ]hinton|\bBBL\b)\b", re.I)
DEFINED = re.compile(r"\b(defined|minimal|mineral salts|mineral medium|mineral base|\bM9\b|\bM63\b|\bMOPS\b|chemically defined|synthetic medium|CDM\b|basal salts|salts medium|\bBG-?11\b|\bAMS\b|nitrate mineral|\bDSMZ? \d|artificial seawater medium|\bASW\b)\b", re.I)
UNSTATED = re.compile(r"^\s*(not stated|not specified|unspecified|unknown|n/s|not explicitly|not reported|n/a|routine|standard (liquid )?medium)\b", re.I)
# variable-composition feedstocks/mixtures — chemically NOT defined even if called a "medium"
FEEDSTOCK = re.compile(r"\b(crude glycerol|plant oil|lignocellulos|hydrolysate|pyrolysis oil|waste|sludge|molasses|biomass|straw|bagasse|manure|digestate|hydrolyzate|distillers|whey|steep liquor)\b", re.I)


def classify(name, comp):
    t = (name or "").strip()
    if not t or UNSTATED.match(t):
        return "unstated", False
    # "no/not defined" and "undefined" must NOT count as a defined medium
    has_defined = bool(DEFINED.search(t)) and not re.search(r"undefined|\b(?:no|not)\s+defined", t, re.I)
    # in-vivo take precedence — not a lab medium at all
    if IN_VIVO.search(t):
        return "in_vivo", False
    # variable-composition feedstock -> not chemically defined
    if FEEDSTOCK.search(t):
        return "complex_undefined", False
    # environmental: natural seawater/soil/field samples aren't formulable; ARTIFICIAL/synthetic seawater is
    if ENVIRON.search(t):
        if re.search(r"\bartificial\b|\bsynthetic\b", t, re.I):
            return "defined", True
        return "environmental", False
    # a medium that names peptone/yeast-extract/commercial base is complex, UNLESS it is explicitly
    # "defined medium without yeast extract" etc.
    without = re.search(r"\b(without|no|free of|minus|lacking)\b[^,;.]{0,25}(yeast extract|peptone|casein|complex)", t, re.I)
    if COMPLEX.search(t) and not without:
        return "complex_undefined", False
    if has_defined or without:
        return "defined", True
    # composition-based fallback: if the enumerated composition is all simple chemicals, treat as defined
    if isinstance(comp, list) and comp:
        undef = sum(1 for c in comp if COMPLEX.search(c.get("name", "") or ""))
        if undef == 0 and all(c.get("name") for c in comp):
            return "defined", True
        if undef:
            return "complex_undefined", False
    return "named_other", False           # a named medium we can't confidently type -> review


def main():
    gr = json.load(open(GR))
    st = collections.Counter(); rel = collections.Counter(); examples = collections.defaultdict(list)
    for r in gr:
        m = r.get("medium") or {}
        if m.get("media_id") or m.get("exchanges"):
            st["already_resolved"] += 1
            m.setdefault("medium_type", "defined"); m.setdefault("formulable", True)
            continue
        name = m.get("canonical_name") or m.get("description") or ""
        mt, formul = classify(name, m.get("composition"))
        m["medium_type"] = mt; m["formulable"] = formul
        st[mt] += 1
        relevant = r.get("growth_rate_per_h") is not None or r.get("uptake_rates") or r.get("secretion_rates")
        if relevant:
            rel[mt] += 1
        if len(examples[mt]) < 6 and name:
            examples[mt].append(name[:65])
    print("ALL records by medium_type:")
    for k, v in st.most_common():
        print(f"  {v:6d}  {k}")
    print("\nVALIDATION-RELEVANT (has mu/rates) unresolved by medium_type:")
    for k, v in rel.most_common():
        print(f"  {v:5d}  {k}")
    print("\nExamples:")
    for k in ("defined", "complex_undefined", "in_vivo", "environmental", "named_other"):
        if examples[k]:
            print(f"  [{k}]")
            for e in examples[k]:
                print(f"      {e}")
    if WRITE:
        json.dump(gr, open(GR, "w"), separators=(",", ":"))
        print("\nWROTE medium_type + formulable into", GR)
    else:
        print("\n(dry run — pass --write to persist)")


if __name__ == "__main__":
    main()
