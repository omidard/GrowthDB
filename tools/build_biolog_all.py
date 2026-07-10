import json, glob, sys, re
sys.path.insert(0,"/tmp/claude-1000/-data-Brilliant-genomics-department/eb8d91f3-1707-45de-a10d-2de68fef6627/scratchpad/media_work/repo/tools")
from map_metabolite import Mapper
m=Mapper()
PN={'positive':'positive','negative':'negative','borderline':'borderline','variable':'variable','+':'positive','-':'negative','no growth':'negative','growth':'positive'}
out=open("lit/biolog_all_phenotypes.jsonl","w"); n=0; seen=set(); mapped=0
for f in glob.glob("lit/biolog_all_ex/batch_*.json"):
    try: d=json.load(open(f))
    except: continue
    for b in d.get("biolog",[]):
        org=b.get("organism"); sub=b.get("substrate")
        if not org or not sub or not b.get("source_snippet"): continue
        key=(org.lower().strip(),(b.get("strain") or "").lower(),sub.lower().strip())
        if key in seen: continue
        seen.add(key)
        hit=m.map(name=sub); ex=hit["exchange"] if hit else None
        if ex: mapped+=1
        au=b.get("first_author") or ""; yr=b.get("year") or ""; pmc=b.get("pmcid","")
        out.write(json.dumps({"source":"biolog_pm" if b.get("is_biolog") else "sub_util","organism":org,"strain":b.get("strain"),
            "substrate":sub,"exchange":ex,"category":b.get("category") or "carbon","phenotype":PN.get((b.get("phenotype") or "").lower(),b.get("phenotype")),
            "plate":b.get("plate"),"kind":"Biolog PM" if b.get("is_biolog") else "substrate utilization","chebi":None,
            "citation":f"{au} et al. ({yr}) [{pmc}]"})+"\n"); n+=1
out.close()
print(f"biolog_all consolidated: {n} phenotypes (deduped) | substrate->BiGG mapped: {mapped}")
