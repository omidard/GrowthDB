#!/usr/bin/env python3
"""Consolidate GrowthDB literature extractions -> growth records. Map uptake/secretion compounds
to BiGG (deterministic), link medium to the Media repo by name, keep method/snippet/citation.
Media not found in Media repo are flagged for addition."""
import json, os, re, glob, sys, math
sys.path.insert(0,"/tmp/claude-1000/-data-Brilliant-genomics-department/eb8d91f3-1707-45de-a10d-2de68fef6627/scratchpad/media_work/repo/tools")
from map_metabolite import Mapper
GW="/tmp/claude-1000/-data-Brilliant-genomics-department/eb8d91f3-1707-45de-a10d-2de68fef6627/scratchpad/growthdb_work"
REPO=os.path.join(GW,"repo"); LIT=os.path.join(GW,"lit","extractions")
m=Mapper()

# Media repo index for linkage (name -> media_id)
MEDIA_IDX=json.load(open("/tmp/claude-1000/-data-Brilliant-genomics-department/eb8d91f3-1707-45de-a10d-2de68fef6627/scratchpad/media_work/repo/data/index.json"))["media"]
def mnorm(s): return re.sub(r"[^a-z0-9]","",(s or "").lower())
MEDIA_BY_NAME={}
for md in MEDIA_IDX:
    MEDIA_BY_NAME.setdefault(mnorm(md["name"]), md["id"])
    # also index common short names
COMMON={"m9":"m9_glucose_aerobic","m9minimalmedium":"m9_glucose_aerobic","lb":"lb_lennox","lbbroth":"lb_lennox",
        "lubriabertani":"lb_lennox","tsb":"tsb","bhi":"bhi","mops":"mops_minimal_glucose","mopsminimal":"mops_minimal_glucose"}
BACKMAP={}
try: BACKMAP=json.load(open(os.path.join(GW,"lit","media_backmap.json")))
except: pass
def link_media(name):
    if not name: return None
    k=mnorm(name)
    if k in MEDIA_BY_NAME: return MEDIA_BY_NAME[k]
    for c,mid in COMMON.items():
        if c in k: return mid
    # substring match against media names
    for md in MEDIA_IDX:
        if k and mnorm(md["name"]).startswith(k) and len(k)>=4: return md["id"]
    return None

def num(v):
    try: return float(re.search(r"[-+]?\d*\.?\d+", str(v)).group(0))
    except: return None
def mapflux(fl, sign):
    out=[]
    for r in (fl or []):
        name=r.get("compound","");
        if not name: continue
        hit=m.map(name=name); ex=hit["exchange"] if (hit and hit["in_biggr"]) else None
        rate=num(r.get("rate"))
        if rate is not None and sign=="uptake" and rate>0: rate=-rate   # BiGG: uptake negative
        if rate is not None and sign=="secretion": rate=abs(rate)
        out.append({"compound":name,"exchange":ex,"bigg_metabolite":(ex[3:-2] if ex else None),
                    "rate":rate,"units":r.get("units") or "mmol/gDW/h","method":r.get("method") or "reported"})
    return out

records=[]; seen=set(); total=0; add_media=[]
for fp in glob.glob(os.path.join(LIT,"batch_*.json")):
    try: data=json.load(open(fp))
    except: continue
    for r in data.get("records",[]):
        total+=1
        if not r.get("source_snippet"): continue
        org=r.get("organism");
        if not org: continue
        mu=num(r.get("growth_rate_per_h")); td=num(r.get("doubling_time_h"))
        if mu is None and td and td>0: mu=round(math.log(2)/td,4)
        dil=num(r.get("dilution_rate_per_h"))
        if mu is None and dil and (r.get("culture_mode","").lower().startswith("chemo")): mu=dil
        up=mapflux(r.get("uptake_rates"),"uptake"); sec=mapflux(r.get("secretion_rates"),"secretion")
        if mu is None and not up and not sec: continue    # need at least one rate
        mid=link_media(r.get("medium_name")) or BACKMAP.get(f"{r.get('pmcid','')}|{r.get('medium_name')}")
        comp=r.get("medium_composition") or []
        if mid is None and comp:
            add_media.append({"pmcid":r.get("pmcid"),"medium_name":r.get("medium_name"),"composition":comp})
        pmcid=r.get("pmcid","PMC"); rid=f"lit_{pmcid}_{len(records)}"
        if rid in seen: continue
        seen.add(rid)
        au=r.get("first_author") or ""; yr=r.get("year") or ""
        rec={"id":rid,"organism":org,"strain":r.get("strain"),
             "growth_rate_per_h":mu,"doubling_time_h":td,
             "conditions":{"culture_mode":r.get("culture_mode") or "unspecified","oxygen":r.get("oxygen"),
                 "temperature_C":num(r.get("temperature_C")),"pH":num(r.get("pH")),"optimum_temperature_C":None,"optimum_pH":None,
                 "dilution_rate_per_h":dil},
             "medium":{"media_id":mid,"description":r.get("medium_name"),
                 "composition":comp if mid is None else None,
                 "note":("linked to Media repo" if mid else "not in Media repo — candidate to add")},
             "uptake_rates":up,"secretion_rates":sec,"carbon_substrates":[],
             "provenance":{"source_type":"literature",
                 "citation":f"{au} et al. {r.get('paper_title','')}. {r.get('journal','')} ({yr})."+(f" DOI:{r.get('doi')}." if r.get('doi') else "")+f" [{pmcid}]",
                 "doi":r.get("doi"),"pmcid":pmcid,"url":f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/",
                 "method":r.get("growth_method") or "reported","snippet":re.sub(r'\s+',' ',r.get('source_snippet',''))[:220]},
             "confidence":r.get("confidence"),"curation_notes":r.get("curation_notes")}
        records.append(rec)
json.dump({"pending_media_for_media_repo":add_media}, open(os.path.join(GW,"lit","media_to_add.json"),"w"), indent=1)
# merge into growth_records.json (Madin seed + literature)
seed=json.load(open(os.path.join(REPO,"data","growth_records.json")))
allrec=seed+records
json.dump(allrec, open(os.path.join(REPO,"data","growth_records.json"),"w"), separators=(",",":"))
print(f"extraction records in: {total} | literature growth records kept: {len(records)}")
print(f"  with growth rate: {sum(1 for r in records if r['growth_rate_per_h'] is not None)} | with uptake: {sum(1 for r in records if r['uptake_rates'])} | with secretion: {sum(1 for r in records if r['secretion_rates'])}")
print(f"  media linked to Media repo: {sum(1 for r in records if r['medium']['media_id'])} | media to ADD: {len(add_media)}")
print(f"total growth records (seed+lit): {len(allrec)}")
