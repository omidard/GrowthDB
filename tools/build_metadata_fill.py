#!/usr/bin/env python3
"""Fill missing temperature/pH/oxygen/medium metadata on growth records from the Methods-section
extractions. Only fills fields that are currently empty; marks what was completed + queues any
newly-found medium compositions for the Media repo."""
import json, glob, re, os
from collections import Counter
GW="/tmp/claude-1000/-data-Brilliant-genomics-department/eb8d91f3-1707-45de-a10d-2de68fef6627/scratchpad/growthdb_work"
def num(v):
    try: return float(re.search(r"[-+]?\d*\.?\d+", str(v)).group(0))
    except: return None
by={}
for f in glob.glob(os.path.join(GW,"lit","metadata","ex","batch_*.json")):
    try: data=json.load(open(f))
    except: continue
    for p in data.get("papers",[]):
        if p.get("pmcid"): by[p["pmcid"]]=p
print("paper metadata records:",len(by))

recs=json.load(open(os.path.join(GW,"repo","data","growth_records.json")))
filled=Counter(); add_media=[]
OXOK={"aerobic","anaerobic","microaerophilic","facultative"}
for r in recs:
    if not r["id"].startswith(("lit_","flux_")): continue
    pmc=(r.get("provenance") or {}).get("pmcid"); p=by.get(pmc)
    if not p: continue
    c=r.setdefault("conditions",{}); med=r.setdefault("medium",{}); done=[]
    if c.get("temperature_C") is None and num(p.get("temperature_C")) is not None:
        c["temperature_C"]=num(p["temperature_C"]); filled["temperature"]+=1; done.append("temperature")
    if c.get("pH") is None and num(p.get("pH")) is not None:
        c["pH"]=num(p["pH"]); filled["pH"]+=1; done.append("pH")
    if not c.get("oxygen") and (p.get("oxygen") in OXOK):
        c["oxygen"]=p["oxygen"]; c["aeration_detail"]=c.get("aeration_detail") or p.get("aeration_detail"); filled["oxygen"]+=1; done.append("oxygen")
    if not c.get("aeration_detail") and p.get("aeration_detail"): c["aeration_detail"]=p["aeration_detail"]
    if not c.get("agitation") and p.get("agitation"): c["agitation"]=p.get("agitation")
    if not (med.get("description") or med.get("media_id")) and p.get("medium_name"):
        med["description"]=p["medium_name"]; filled["medium"]+=1; done.append("medium")
        if p.get("medium_composition"): med["composition"]=p["medium_composition"]
    comp=med.get("composition")
    if comp and med.get("description"):
        add_media.append({"pmcid":pmc,"medium_name":med.get("description"),"composition":comp})
    if done:
        r["curation_notes"]=(r.get("curation_notes") or "")+f" [metadata completed from paper Methods: {', '.join(done)}]"
        r.setdefault("provenance",{})["metadata_from_methods"]=done
json.dump(recs,open(os.path.join(GW,"repo","data","growth_records.json"),"w"),separators=(",",":"))
json.dump({"pending_media_for_media_repo":add_media},open(os.path.join(GW,"lit","media_to_add_meta.json"),"w"))
print("FILLED from Methods:",dict(filled))
lit=[r for r in recs if r["id"].startswith(("lit_","flux_"))]
print("literature records now missing: pH",sum(1 for r in lit if r["conditions"].get("pH") is None),
      "| T",sum(1 for r in lit if r["conditions"].get("temperature_C") is None),
      "| oxygen",sum(1 for r in lit if not r["conditions"].get("oxygen")),
      "| medium",sum(1 for r in lit if not (r["medium"].get("description") or r["medium"].get("media_id"))))
