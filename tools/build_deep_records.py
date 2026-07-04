#!/usr/bin/env python3
"""Consolidate the deep Opus re-mine into rich growth records; SUPERSEDE the shallow lit records
from the same (flux) papers. Map every uptake/secretion rate to BiGG, keep composition/yields/
calculations. Emits media_to_add2.json for Media. curate.py is run AFTER this to normalize vocab
+ link media."""
import json, os, re, glob, sys, math
sys.path.insert(0,"/tmp/claude-1000/-data-Brilliant-genomics-department/eb8d91f3-1707-45de-a10d-2de68fef6627/scratchpad/media_work/repo/tools")
from map_metabolite import Mapper
GW="/tmp/claude-1000/-data-Brilliant-genomics-department/eb8d91f3-1707-45de-a10d-2de68fef6627/scratchpad/growthdb_work"
REPO=os.path.join(GW,"repo"); RE=os.path.join(GW,"lit","reextract"); GEMMINE=os.path.join(GW,"lit","gemmine")
m=Mapper()
def num(v):
    try: return float(re.search(r"[-+]?\d*\.?\d+([eE][-+]?\d+)?", str(v)).group(0))
    except: return None
def mapflux(fl, sign):
    out=[]
    for r in (fl or []):
        name=(r.get("compound") or "").strip()
        if not name: continue
        hit=m.map(name=name); ex=hit["exchange"] if (hit and hit["in_biggr"]) else None
        rate=num(r.get("rate"))
        if rate is not None and sign=="uptake" and rate>0: rate=-rate
        if rate is not None and sign=="secretion": rate=abs(rate)
        out.append({"compound":name,"exchange":ex,"bigg_metabolite":(ex[3:-2] if ex else None),
                    "rate":rate,"units":r.get("units") or "mmol/gDW/h","method":(r.get("method") or "reported")})
    return out

deep=[]; deep_pmcids=set(); add_media=[]
for f in glob.glob(os.path.join(RE,"batch_*.json"))+glob.glob(os.path.join(GEMMINE,"batch_*.json")):
    try: data=json.load(open(f))
    except: continue
    for r in data.get("records",[]):
        if not r.get("source_snippet"): continue
        org=r.get("organism")
        if not org: continue
        pmc=r.get("pmcid","PMC"); deep_pmcids.add(pmc)
        mu=num(r.get("growth_rate_per_h")); td=num(r.get("doubling_time_h")); dil=num(r.get("dilution_rate_per_h"))
        if mu is None and td and td>0: mu=round(math.log(2)/td,4)
        if mu is None and dil and (r.get("culture_mode","").lower().startswith(("chemo","turbi","conti"))): mu=dil
        up=mapflux(r.get("uptake_rates"),"uptake"); sec=mapflux(r.get("secretion_rates"),"secretion")
        if mu is None and not up and not sec: continue
        comp=r.get("medium_composition") or []
        if comp and r.get("medium_composition_stated",True):
            add_media.append({"pmcid":pmc,"medium_name":r.get("medium_name"),"composition":comp})
        au=r.get("first_author") or ""; yr=r.get("year") or ""
        rec={"id":f"flux_{pmc}_{len(deep)}","organism":org,"strain":r.get("strain"),
             "growth_rate_per_h":mu,"doubling_time_h":td,
             "conditions":{"culture_mode":r.get("culture_mode") or "unspecified","culture_detail":r.get("culture_detail"),
                 "oxygen":r.get("oxygen"),"aeration_detail":r.get("aeration_detail"),
                 "temperature_C":num(r.get("temperature_C")),"pH":num(r.get("pH")),
                 "optimum_temperature_C":None,"optimum_pH":None,"dilution_rate_per_h":dil},
             "medium":{"media_id":None,"description":r.get("medium_name"),
                 "composition":comp,"composition_stated":r.get("medium_composition_stated",bool(comp)),
                 "note":"medium composition extracted from paper"},
             "uptake_rates":up,"secretion_rates":sec,"mfa_fluxes":r.get("mfa_fluxes") or [],"yields":r.get("yields") or {},"carbon_substrates":[],
             "provenance":{"source_type":"literature",
                 "citation":f"{au} et al. {r.get('paper_title','')}. {r.get('journal','')} ({yr})."+(f" DOI:{r.get('doi')}." if r.get('doi') else "")+f" [{pmc}]",
                 "doi":r.get("doi"),"pmcid":pmc,"url":f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmc}/",
                 "method":r.get("growth_method") or "reported","snippet":re.sub(r'\s+',' ',r.get('source_snippet',''))[:240]},
             "confidence":r.get("confidence"),"curation_notes":r.get("curation_notes")}
        deep.append(rec)
json.dump({"pending_media_for_media_repo":add_media}, open(os.path.join(GW,"lit","media_to_add2.json"),"w"), indent=1)
# supersede: drop old lit records from the re-mined (flux) papers; keep Madin + other lit
cur=json.load(open(os.path.join(REPO,"data","growth_records.json")))
def oldpmc(r): return (r.get("provenance",{}) or {}).get("pmcid")
kept=[r for r in cur if not (r["id"].startswith("flux_") or (r["id"].startswith("lit_") and oldpmc(r) in deep_pmcids))]
dropped=len(cur)-len(kept)
allrec=kept+deep
json.dump(allrec, open(os.path.join(REPO,"data","growth_records.json"),"w"), separators=(",",":"))
print(f"deep records built: {len(deep)} | superseded {dropped} shallow records from {len(deep_pmcids)} flux papers")
print(f"  with uptake: {sum(1 for r in deep if r['uptake_rates'])} | secretion: {sum(1 for r in deep if r['secretion_rates'])} | composition: {sum(1 for r in deep if r['medium']['composition'])}")
print(f"  media compositions queued for Media: {len(add_media)}")
print(f"total growth records now: {len(allrec)}")
