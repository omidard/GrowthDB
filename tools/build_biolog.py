import json, glob, os, sys, re
sys.path.insert(0,"/tmp/claude-1000/-data-Brilliant-genomics-department/eb8d91f3-1707-45de-a10d-2de68fef6627/scratchpad/media_work/repo/tools")
from map_metabolite import Mapper
m=Mapper()
REPO="/tmp/claude-1000/-data-Brilliant-genomics-department/eb8d91f3-1707-45de-a10d-2de68fef6627/scratchpad/growthdb_work/repo"
rows=[]; seen=set()
for f in glob.glob("lit/biolog_ex/batch_*.json"):
    try: d=json.load(open(f))
    except: continue
    for b in d.get("biolog",[]):
        org=b.get("organism"); sub=b.get("substrate"); ph=(b.get("phenotype") or "").lower()
        if not org or not sub or not b.get("source_snippet"): continue
        key=(org.lower().strip(),(b.get("strain") or "").lower(),sub.lower().strip())
        if key in seen: continue
        seen.add(key)
        hit=m.map(name=sub); ex=hit["exchange"] if (hit and hit["in_biggr"]) else None
        pn={"positive":"positive","negative":"negative","borderline":"borderline","variable":"variable","+":"positive","-":"negative"}.get(ph,ph or None)
        au=b.get("first_author") or ""; yr=b.get("year") or ""; pmc=b.get("pmcid","")
        rows.append({"organism":org,"strain":b.get("strain"),"plate":b.get("plate"),
            "substrate":sub,"exchange":ex,"category":b.get("category"),"phenotype":pn,
            "value":b.get("value"),"units":b.get("units"),
            "citation":f"{au} et al. ({yr}) [{pmc}]","pmcid":pmc,"doi":b.get("doi"),
            "snippet":re.sub(r'\s+',' ',b.get('source_snippet',''))[:180]})
from collections import Counter
json.dump({"count":len(rows),"n_organisms":len(set(r['organism'] for r in rows)),
    "n_mapped":sum(1 for r in rows if r['exchange']),
    "by_phenotype":dict(Counter(r['phenotype'] for r in rows)),
    "by_category":dict(Counter(r['category'] for r in rows)),
    "phenotypes":rows}, open(os.path.join(REPO,"data","biolog_phenotypes.json"),"w"), separators=(",",":"))
print(f"biolog phenotypes: {len(rows)} (deduped) | organisms: {len(set(r['organism'] for r in rows))} | substrate->BiGG mapped: {sum(1 for r in rows if r['exchange'])}")
print("by phenotype:",dict(Counter(r['phenotype'] for r in rows)),"| by category:",dict(Counter(r['category'] for r in rows)))
