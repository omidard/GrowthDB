#!/usr/bin/env python3
"""Apply mined strain designations (/tmp/strain_mined/*.json) to the unidentified paper-backed
GrowthDB records, then re-parse the genotype, and classify every record's strain LEVEL so the
autocurator can be honest about resolution:

  level='strain'   a specific strain is known (culture-collection id / lab designation / genotype)
  level='species'  a species-level reference (Madin trait synthesis and similar DB compilations —
                   these never tracked a strain; still usable as a species-typical growth rate)
  level='unknown'  paper-backed but the strain genuinely couldn't be found (should be rare)

Records are matched to a mined strain by PMC + organism/species. Writes in place (dry-run unless --write).
"""
import json, os, re, sys, glob, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import parse_strain_genotype as PSG

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GR = os.path.join(ROOT, "data", "growth_records.json")
MINED_DIR = "/tmp/strain_mined"
WRITE = "--write" in sys.argv


def norm_org(s):
    return re.sub(r"[^a-z]", "", (s or "").lower())[:20]


def main():
    gr = json.load(open(GR))
    # index mined strains by (pmc, normalised organism)
    mined = {}
    for f in glob.glob(MINED_DIR + "/*.json"):
        try:
            d = json.load(open(f))
        except Exception:
            continue
        pmc = d.get("pmc") or os.path.basename(f)[:-5]
        for s in d.get("strains", []):
            strain = (s.get("strain") or "").strip()
            if not strain or re.search(r"not stated|not (given|reported|specified)|unknown", strain, re.I):
                continue
            mined[(pmc, norm_org(s.get("organism")))] = strain

    applied = 0; level = collections.Counter()
    for r in gr:
        relevant = r.get("growth_rate_per_h") is not None or r.get("uptake_rates") or r.get("secretion_rates")
        spar = r.get("strain_parsed") or {}
        # 1) apply a recovered strain to an unidentified paper-backed record
        if spar.get("unidentified"):
            m = re.search(r"(PMC\d+)", r.get("id", ""))
            org = r.get("gtdb_species") or r.get("species") or r.get("organism")
            key = (m.group(1), norm_org(org)) if m else None
            if key and key in mined:
                r["strain"] = mined[key]
                r["strain_parsed"] = PSG.parse(mined[key])
                spar = r["strain_parsed"]
                applied += 1
        # 2) classify the strain level
        st_raw = r.get("strain") if isinstance(r.get("strain"), str) else ""
        blank = (not st_raw.strip()) or bool(re.match(r"^\s*(not stated|not specified|unspecified|unknown|not given|not reported|n/?a|-)\b", st_raw, re.I))
        if spar.get("parent") or (st_raw.strip() and not blank):
            lvl = "strain"; spar.pop("strain_reason", None)   # has a canonical parent OR a real strain description
        elif re.match(r"^madin_", r.get("id", "")) or (r.get("provenance") or {}).get("source_type") == "database":
            lvl = "species"                              # species-level trait compilation, never had a strain
            spar["strain_reason"] = "species-level trait compilation (Madin) — no strain recorded in source"
        else:
            lvl = "unknown"
            # WHY there is no strain — verified by mining the source paper (it genuinely names none)
            st = (r.get("strain") or "").lower()
            if re.search(r"agent-based|generic|theoretical|in silico|parameter", st):
                spar["strain_reason"] = "source is a theoretical/generic model — no specific strain"
            elif re.search(r"multiple isolates|isolates|panel|community", st):
                spar["strain_reason"] = "a panel of multiple isolates — not one strain"
            elif re.search(r"sender|reporter|synthetic|engineered .*(strain|construct)", st):
                spar["strain_reason"] = "an engineered construct with no standard base strain named"
            else:
                spar["strain_reason"] = "the source paper (often a review/modeling reuse) states no strain"
        spar["level"] = lvl
        r["strain_parsed"] = spar
        if relevant:
            level[lvl] += 1

    print("strains recovered from papers and applied:", applied)
    print("validation-relevant records by strain level:")
    for k, v in level.most_common():
        print(f"  {v:5d}  {k}")
    if WRITE:
        json.dump(gr, open(GR, "w"), separators=(",", ":"))
        print("WROTE", GR)
    else:
        print("(dry run — pass --write to persist)")


if __name__ == "__main__":
    main()
