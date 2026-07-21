#!/usr/bin/env python3
"""Formulate the BASE MEDIUM of every substrate-utilisation phenotype to a BiGG exchange set, and
separate the three assay types that were being conflated as 'grows-on-X':

  growth   (base_type 'defined-minimal'): the substrate is the SOLE carbon/N/S/P source on a defined
           minimal base -> a real growth-on-X test. Medium = defined-minimal base + the substrate.
  acid     (base_type 'complex'): fermentative ACID production from the substrate on a COMPLEX peptone
           base (API 50CH style). This is NOT growth on the substrate; it validates the CATABOLIC
           pathway for it. Base is complex -> not exactly formulable; flagged, exported separately.
  enzyme   (base_type 'assay'): a hydrolysis/reduction/activity assay (gelatinase, urease, nitrate
           reductase...). NOT a growth test; validates that a specific REACTION exists in the model.

For growth phenotypes the base is formulated: a defined-minimal mineral base (C/N/S/P-free as
appropriate) + trace metals, plus a vitamin-supplemented variant for fastidious gram-positives
(Biolog IF-10 style). The substrate exchange is added per phenotype at validation time.

Writes phenotype.assay_type + phenotype.base_media_id + phenotype.base_formulable into the species
shards, and the base media into the Media DB. Dry-run unless --write.
"""
import json, os, re, sys, glob, collections

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SPECIES_DIR = os.path.join(ROOT, "data", "species")
MEDIA_REPO = "/data/media_curate"
WRITE = "--write" in sys.argv

MINERALS = ["pi", "so4", "cl", "na1", "k", "nh4", "mg2", "ca2", "fe2", "fe3", "mn2", "zn2", "cu2", "cobalt2", "mobd", "ni2", "h2o", "h"]
VITAMINS = ["thm", "ribflv", "nac", "pnto__R", "pydxn", "fol", "btn", "4abz", "cbl1", "lipoate", "5mthf"]

# defined-minimal base media (exchanges; the tested substrate is added at validation time).
# lb -1000 = unlimited (minerals/vitamins are never limiting for FBA feasibility).
def _ex(biggs, o2=True):
    d = {b: -1000.0 for b in biggs}
    if o2:
        d["o2"] = -1000.0
    return [{"exchange": "EX_" + b + "_e", "bigg": b, "lb": lb, "ub": 1000.0} for b, lb in d.items()]

BASES = {
    "def_min_c": {"name": "Defined-minimal carbon-utilisation base (C-free)",
                  "desc": "Mineral salts + trace metals, no carbon source; the tested substrate is supplied as the sole carbon/energy source (Biolog GN IF-0 / BacDive assimilation).",
                  "ex": _ex(MINERALS)},
    "def_min_c_vit": {"name": "Defined-minimal carbon base + vitamins (gram-positive / fastidious, IF-10 style)",
                      "desc": "Mineral salts + trace metals + B-vitamins, no carbon; for fastidious gram-positives whose Biolog base (IF-10) is supplemented, so vitamin auxotrophy is not scored as a carbon-utilisation failure.",
                      "ex": _ex(MINERALS + VITAMINS)},
    "def_min_n": {"name": "Defined-minimal nitrogen-utilisation base (N-free, +glucose)",
                  "desc": "Mineral salts + trace + D-glucose as carbon, no ammonium; the tested substrate is the sole nitrogen source (Biolog PM3).",
                  "ex": _ex([m for m in MINERALS if m != "nh4"]) + [{"exchange": "EX_glc__D_e", "bigg": "glc__D", "lb": -10.0, "ub": 1000.0}]},
    "def_min_ps": {"name": "Defined-minimal P/S-utilisation base (+glucose, +NH4)",
                   "desc": "Mineral salts + trace + D-glucose + ammonium; the tested substrate is the sole phosphorus or sulfur source (Biolog PM4). The corresponding mineral (Pi or SO4) is withheld at validation time.",
                   "ex": _ex(MINERALS) + [{"exchange": "EX_glc__D_e", "bigg": "glc__D", "lb": -10.0, "ub": 1000.0}]},
}

# gram-positive genera (Bacillota / Actinomycetota) -> use the vitamin-supplemented base
GP_GENERA = set("""bacillus clostridium staphylococcus streptococcus lactobacillus lactococcus enterococcus
listeria corynebacterium mycobacterium mycolicibacterium streptomyces bifidobacterium actinomyces micrococcus
arthrobacter paenarthrobacter rhodococcus nocardia geobacillus paenibacillus lysinibacillus brevibacillus
weissella leuconostoc pediococcus oenococcus carnobacterium propionibacterium cutibacterium eubacterium
faecalibacterium ruminococcus roseburia blautia coprococcus dorea anaerostipes butyrivibrio
gardnerella kocuria brevibacterium microbacterium cellulomonas curtobacterium gordonia tsukamurella
sarcina peptostreptococcus peptococcus finegoldia anaerococcus veillonella megasphaera selenomonas
heliobacterium sporosarcina oceanobacillus halobacillus virgibacillus salinicoccus thermoanaerobacter
moorella caldicellulosiruptor thermoanaerobacterium desulfotomaculum symbiobacterium""".split())


def base_for(cat, gram_pos):
    if cat == "nitrogen":
        return "def_min_n"
    if cat in ("phosphorus", "sulfur"):
        return "def_min_ps"
    return "def_min_c_vit" if gram_pos else "def_min_c"   # carbon (default)


def assay_type_of(bt, kinds):
    ks = set((k or "").lower() for k in (kinds or []))
    if bt == "defined-minimal":
        return "growth"
    if bt == "complex":
        return "acid_production"
    if bt == "assay":
        return "enzyme"
    # unknown: infer from kinds
    if ks & {"builds acid from", "fermentation", "acidification"}:
        return "acid_production"
    if ks & {"hydrolysis", "reduction", "activity"}:
        return "enzyme"
    if ks & {"carbon source", "assimilation", "growth", "respiration", "utilization", "utilisation"}:
        return "growth"
    return "growth_uncertain"


def write_media_db():
    made = 0
    for mid, b in BASES.items():
        comps = [{"name": e["bigg"], "bigg_metabolite": e["bigg"], "exchange": e["exchange"],
                  "lower_bound": e["lb"], "upper_bound": e["ub"], "exchange_source": "biolog_base_formulation",
                  "mapping_method": "curated", "mapping_confidence": "curated"} for e in b["ex"]]
        doc = {"id": "biolog_base_" + mid, "name": b["name"], "category": "growth_medium", "organism_scope": "prokaryote",
               "aerobic": True, "namespace": "bigg", "description": b["desc"],
               "provenance": {"source_type": "standard", "citation": "Biolog Phenotype MicroArray / BacDive substrate-utilisation base medium, formulated to BiGG exchanges for GEM validation.", "doi": "", "url": ""},
               "components": comps, "n_components": len(comps), "n_mapped": len(comps), "n_in_biggr": len(comps),
               "version": 1, "coverage": {"n_compounds": len(comps), "n_covered": len(comps), "n_uncovered": 0, "pct_covered": 100.0, "by_source": {"biggr": len(comps)}},
               "uncovered": [], "oxygen": "aerobic", "curation": "curated"}
        if WRITE:
            with open(os.path.join(MEDIA_REPO, "data", "media", "biolog_base_" + mid + ".json"), "w") as fh:
                json.dump(doc, fh, indent=1)
        made += 1
    return made


def main():
    at = collections.Counter(); baseassign = collections.Counter(); touched = 0
    for f in glob.glob(SPECIES_DIR + "/*.json"):
        try:
            d = json.load(open(f))
        except Exception:
            continue
        genus = (d.get("species") or "").split()[0].lower()
        gram_pos = genus in GP_GENERA
        changed = False
        for p in (d.get("phenotypes") or []):
            a = assay_type_of(p.get("base_medium_type") or p.get("base_type"), p.get("kinds"))
            p["assay_type"] = a
            if a in ("growth", "growth_uncertain"):
                bm = base_for(p.get("category"), gram_pos)
                p["base_media_id"] = "biolog_base_" + bm
                p["base_formulable"] = True
                baseassign[bm] += 1
            elif a == "acid_production":
                p["base_media_id"] = None; p["base_formulable"] = False   # complex peptone base
            else:
                p["base_media_id"] = None; p["base_formulable"] = False   # enzyme assay, not a growth medium
            at[a] += 1
            changed = True
        if changed and WRITE:
            json.dump(d, open(f, "w"), separators=(",", ":"))
            touched += 1
    n_media = write_media_db()
    print("assay_type assigned:")
    for k, v in at.most_common():
        print(f"  {v:7d}  {k}")
    print("\ngrowth-phenotype base assignment:")
    for k, v in baseassign.most_common():
        print(f"  {v:7d}  {k}")
    print("\nbase media formulated:", n_media, "->", list(BASES.keys()))
    if WRITE:
        print("shards patched:", touched)
    else:
        print("(dry run — pass --write to persist)")


if __name__ == "__main__":
    main()
