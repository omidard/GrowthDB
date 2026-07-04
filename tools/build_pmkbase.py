import csv, json, sys, re
sys.path.insert(0,"/tmp/claude-1000/-data-Brilliant-genomics-department/eb8d91f3-1707-45de-a10d-2de68fef6627/scratchpad/media_work/repo/tools")
from map_metabolite import Mapper
m=Mapper()
cite="Pan-Putida: experimentally-validated GEMs + Biolog PM for 24 P. putida strains (mSystems 2025); Zenodo 10.5281/zenodo.17382094"
out=open("lit/pmkbase/pmkbase_phenotypes.jsonl","w"); n=0; mapped=0
for x in csv.DictReader(open("lit/pmkbase/data/growth_summary.csv")):
    sub=(x.get("Compound") or "").strip()
    if not sub or "negative control" in sub.lower(): continue
    try: g=float(x.get("Growth"))
    except: continue
    kegg=(x.get("KEGG ID") or "").strip() or None
    hit=m.map(name=sub, kegg=kegg); ex=hit["exchange"] if (hit and hit["in_biggr"]) else None
    if ex: mapped+=1
    out.write(json.dumps({"source":"pmkbase","organism":"Pseudomonas putida","strain":x.get("Strain"),
        "substrate":sub,"exchange":ex,"kegg":kegg,"category":"carbon",
        "phenotype":("positive" if g>=0.5 else "negative"),"kind":"Biolog PM","citation":cite})+"\n"); n+=1
out.close()
print(f"pmkbase phenotypes: {n} | substrate->BiGG mapped: {mapped}")
