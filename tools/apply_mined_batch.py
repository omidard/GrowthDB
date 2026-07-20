#!/usr/bin/env python3
"""Formulate the mined medium compositions (/tmp/mined/*.json, one file per paper) to BiGG
exchanges with the VALIDATED formula-salt parser, link them to the GrowthDB growth records by
PMC + medium-name, and stage a media-DB payload. Precision stays in code: the agents only
transcribed ingredient lists from the papers; the deterministic parser does the chemistry.

Only media the agent actually found (non-empty composition) are formulated; a 'not in paper'
result is skipped and left flagged.
"""
import json, os, re, sys, glob, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import formulate_media as FM

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GR = os.path.join(ROOT, "data", "growth_records.json")
MINED_DIR = "/tmp/mined"
WRITE = "--write" in sys.argv

# a defined MINIMAL medium provides an unlimited mineral background; add any of these the recipe omits
MINERAL_BASE = ["pi", "so4", "cl", "na1", "k", "nh4", "mg2", "ca2", "fe2", "fe3", "mn2", "zn2", "cu2", "cobalt2", "mobd", "h2o", "h"]


def toks(s):
    return set(re.findall(r"[a-z0-9]{3,}", (s or "").lower()))


# reducing agents / redox poisers — NEVER a carbon source
REDUCING = re.compile(r"cystein|\bna2s\b|sodium sulfide|sulfide|dithion|dithiothreitol|\bdtt\b|titanium.?(iii|citrate)|ascorb|thioglycol|resazurin|mercaptoethanol", re.I)
# variable-composition feedstock named as the substrate -> the medium is not chemically defined
FEEDSTOCK = re.compile(r"crude oil|\boil\b|lignocellulos|hydrolysate|pyrolysis|molasses|\bwaste\b|sludge|distillers|whey|steep liquor|biomass\b", re.I)
# growth-factor vitamins — organic (so _is_carbon is True) but NEVER the growth carbon source
VITAMINISH = re.compile(r"aminobenzo|\bpaba\b|nicotin|thiamin|riboflavin|biotin|folate|folic|pantothen|pyridox|cobalamin|\bb12\b|lipoic|menaquinone|\bvitamin", re.I)


def _is_carbon(b):
    """Carbon source iff the BiGG formula contains carbon and it isn't CO2/bicarbonate or an inorganic ion.
    This correctly treats C1 substrates (methanol/formate) and salt-derived organic anions (acetate from
    sodium acetate) as carbon, which the salt parser's hardcoded is_c=False misses."""
    if b in ("co2", "hco3") or b in FM.INORG_ION:
        return False
    f = (FM.BIGG.get(b, {}).get("xrefs", {}) or {}).get("formula", "") or ""
    return bool(re.search(r"C\d*", f))


def carbon_ex(names):
    """Authoritative carbon/energy sources -> {bigg: is_carbon}. Parentheticals stripped first (so a
    description like 'CO2 (reduced to methane)' can't leak methane); reducing agents excluded; H2 kept
    as an energy source."""
    out = {}
    for nm in names or []:
        base = re.sub(r"\([^)]*\)", "", str(nm)).strip()
        if not base or REDUCING.search(base):
            continue
        low = base.lower()
        if re.search(r"\bh2\b|hydrogen", low):
            out["h2"] = True                             # energy source (no carbon, but drives autotrophy)
        if re.search(r"\bco2\b|carbon dioxide|bicarbonate", low):
            out.setdefault("co2", False)
        vit = bool(VITAMINISH.search(base))              # a vitamin: map it, but never as the carbon
        got, cls = FM.ingredient_exchanges(base)
        for b, is_c in got:
            out[b] = out.get(b, False) or (_is_carbon(b) and not vit)
    return out


def formulate(comp, is_defined, carbon_sources=None):
    """carbon_sources is AUTHORITATIVE. When a source is named but unmappable (e.g. a polysaccharide),
    we do NOT let a random micronutrient in the composition become the carbon — we skip instead."""
    ex = {}; unmapped = 0; complex_n = 0
    src_named = bool(carbon_sources)
    cfromsrc = carbon_ex(carbon_sources or [])
    for c in comp:
        got, cls = FM.ingredient_exchanges((c.get("name") or "").strip())
        if cls == "complex":
            complex_n += 1
        elif cls in ("unmapped", "vitamin_unlisted", "base_ref"):
            unmapped += 1
        for bigg, is_c in got:
            # a composition compound may be THE carbon only if the paper named no carbon source at all
            treat_c = _is_carbon(bigg) and not src_named
            lb = -FM.CARBON_CAP if treat_c else -1000.0
            if bigg not in ex or lb < ex[bigg]:
                ex[bigg] = lb
    for b, isc in cfromsrc.items():                      # authoritative carbon overrides
        if isc:
            ex[b] = -FM.CARBON_CAP                        # force the designated carbon to the carbon cap
        elif b not in ex or -1000.0 < ex[b]:
            ex[b] = -1000.0
    has_carbon = any(v == -FM.CARBON_CAP for v in ex.values())
    if is_defined:
        for b in MINERAL_BASE:
            ex.setdefault(b, -1000.0)
    exch = [{"exchange": "EX_" + b + "_e", "bigg": b, "lb": ex[b], "ub": 1000.0} for b in ex]
    return exch, unmapped, complex_n, has_carbon


def main():
    gr = json.load(open(GR))
    bypmc = collections.defaultdict(list)
    for r in gr:
        m = re.search(r"(PMC\d+)", r.get("id", ""))
        if m:
            bypmc[m.group(1)].append(r)

    mined = []
    for f in glob.glob(MINED_DIR + "/*.json"):
        try:
            d = json.load(open(f))
        except Exception as e:
            print("skip", f, e); continue
        pmc = d.get("pmc") or os.path.basename(f)[:-5]
        for md in d.get("media", []):
            mined.append((pmc, md))

    linked = 0; form = 0; skipped = 0; to_add = []
    for pmc, md in mined:
        comp = md.get("composition") or []
        carbon_sources = md.get("carbon_sources") or []
        # a variable-composition feedstock (crude oil, hydrolysate…) is not chemically defined — skip
        if any(FEEDSTOCK.search(str(s)) for s in carbon_sources):
            skipped += 1; continue
        # formulate from the enumerated composition; fall back to the confirmed carbon source(s) when
        # the paper only named the carbon source (a defined-minimal medium is {mineral base + carbon})
        exch, unmapped, complex_n, has_c = formulate(comp, md.get("is_defined", True), carbon_sources)
        if not exch or not has_c:
            skipped += 1; continue
        form += 1
        if complex_n:
            status = "mined_partial_complex"          # contains an undefined pool (e.g. yeast extract)
        elif not comp and carbon_sources:
            status = "mined_defined_minimal"           # reconstructed from confirmed carbon source + mineral base
        elif unmapped:
            status = "mined_partial_unmapped"
        else:
            status = "mined_defined"
        mt = toks(md.get("name"))
        # link to this paper's still-unresolved records whose medium name best-matches
        for r in bypmc.get(pmc, []):
            med = r.get("medium") or {}
            if med.get("media_id") or med.get("exchanges"):
                continue
            rn = med.get("canonical_name") or med.get("description") or ""
            if len(mt & toks(rn)) >= 2 or (md.get("name", "")[:15].lower() in rn.lower()):
                med["exchanges"] = exch
                med["formulation"] = status
                med["formulated_from"] = "mined recipe (paper Methods/SI)"
                med["composition"] = comp
                linked += 1
        to_add.append({"medium_name": md.get("name"), "pmcid": pmc, "composition": comp,
                       "is_defined": bool(md.get("is_defined")), "status": status})

    if WRITE:
        json.dump(gr, open(GR, "w"), separators=(",", ":"))
        json.dump({"pending_media_for_media_repo": to_add}, open(os.path.join(ROOT, "data", "mined_batch_to_add.json"), "w"), indent=1)
    print("mined media files: %d | formulated: %d | skipped (no comp / no carbon): %d | records linked: %d"
          % (len(mined), form, skipped, linked))
    if not WRITE:
        print("(dry run — pass --write to persist)")


if __name__ == "__main__":
    main()
