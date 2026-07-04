#!/usr/bin/env python3
"""Pull BacDive (free API, 2026): per strain -> growth conditions + MediaDive-linked media
(-> Media repo) + metabolite-utilization / API phenotypes (ChEBI -> BiGG). Deterministic."""
import json, os, sys, time, urllib.request
sys.path.insert(0,"/tmp/claude-1000/-data-Brilliant-genomics-department/eb8d91f3-1707-45de-a10d-2de68fef6627/scratchpad/media_work/repo/tools")
from map_metabolite import Mapper
m=Mapper()
GW="/tmp/claude-1000/-data-Brilliant-genomics-department/eb8d91f3-1707-45de-a10d-2de68fef6627/scratchpad/growthdb_work"
OUT=os.path.join(GW,"lit","bacdive"); os.makedirs(OUT,exist_ok=True)
MEDIA_IDS={mm["id"] for mm in json.load(open("/tmp/claude-1000/-data-Brilliant-genomics-department/eb8d91f3-1707-45de-a10d-2de68fef6627/scratchpad/media_work/repo/data/index.json"))["media"]}
MAXID=int(sys.argv[1]) if len(sys.argv)>1 else 190000
def get(url):
    for _ in range(3):
        try:
            req=urllib.request.Request(url,headers={"Accept":"application/json","User-Agent":"GrowthDB/1.0"})
            return json.load(urllib.request.urlopen(req,timeout=60))
        except Exception: time.sleep(2)
    return None
def aslist(x): return x if isinstance(x,list) else ([x] if x else [])
def norm_ox(s):
    s=(s or "").lower()
    if "facultative" in s: return "facultative"
    if "microaero" in s: return "microaerophilic"
    if "anaerob" in s or "anoxic" in s: return "anaerobic"
    if "aerob" in s: return "aerobic"
    return s or None
def num(v):
    try: return float(str(v).split("-")[0])
    except: return None

gf=open(os.path.join(OUT,"bacdive_growth.jsonl"),"w")
pf=open(os.path.join(OUT,"bacdive_phenotypes.jsonl"),"w")
ng=0; npg=0; nstrain=0
for start in range(1,MAXID+1,100):
    ids=";".join(str(i) for i in range(start,start+100))
    d=get(f"https://api.bacdive.dsmz.de/fetch/{ids}")
    if not d: continue
    res=d.get("results",{})
    items=res.items() if isinstance(res,dict) else [(None,r) for r in res]
    for sid,rec in items:
        if not isinstance(rec,dict): continue
        nstrain+=1
        tax=rec.get("Name and taxonomic classification",{}) or {}
        sp=tax.get("species") or tax.get("full scientific name") or ((tax.get("genus","")+" "+tax.get("species epithet","")).strip())
        if not sp: continue
        strain=tax.get("strain designation") or (rec.get("General",{}) or {}).get("BacDive-ID")
        cgc=rec.get("Culture and growth conditions",{}) or {}
        pm=rec.get("Physiology and metabolism",{}) or {}
        ox=None
        oxt=pm.get("oxygen tolerance")
        for o in aslist(oxt):
            if isinstance(o,dict): ox=norm_ox(o.get("oxygen tolerance")); break
        temps=[num(t.get("temperature")) for t in aslist(cgc.get("culture temp")) if isinstance(t,dict) and (t.get("growth") in ("positive","yes")) and num(t.get("temperature")) is not None]
        # media
        media=[]
        for md in aslist(cgc.get("culture medium")):
            if not isinstance(md,dict): continue
            link=md.get("link") or ""; mid=None
            if "mediadive.dsmz.de/medium/" in link:
                cand="mediadive_"+link.rstrip("/").split("/")[-1]
                if cand in MEDIA_IDS: mid=cand
            media.append({"name":md.get("name"),"growth":md.get("growth"),"media_id":mid,
                          "media_url":(f"https://omidard.github.io/Media/?medium={mid}" if mid else None)})
        cite=f"BacDive strain {sid} (Reimer LC et al., BacDive 2022, NAR); https://bacdive.dsmz.de/strain/{sid}"
        gf.write(json.dumps({"source":"bacdive","bacdive_id":sid,"organism":sp,"strain":strain,
            "oxygen":ox,"growth_temperatures_C":sorted(set(temps)),
            "optimum_temperature_C":(max(set(temps),key=temps.count) if temps else None),
            "media":media,"citation":cite})+"\n"); ng+=1
        # phenotypes: metabolite utilization + production
        for mu in aslist(pm.get("metabolite utilization")):
            if not isinstance(mu,dict): continue
            met=mu.get("metabolite"); act=(mu.get("utilization activity") or "").strip()
            if not met: continue
            ph="positive" if act=="+" else ("negative" if act=="-" else ("borderline" if act in ("+/-","+ / -") else None))
            chebi=mu.get("Chebi-ID"); hit=m.map(name=met, chebi=str(chebi) if chebi else None)
            ex=hit["exchange"] if (hit and hit["in_biggr"]) else None
            pf.write(json.dumps({"source":"bacdive","organism":sp,"strain":strain,"substrate":met,
                "exchange":ex,"chebi":chebi,"category":"carbon","phenotype":ph,
                "kind":mu.get("kind of utilization tested"),"citation":cite})+"\n"); npg+=1
    if (start//100)%50==0:
        print(f"...{start} ids scanned | strains {nstrain} | growth {ng} | phenotypes {npg}",flush=True)
gf.close(); pf.close()
print(f"DONE: strains {nstrain} | bacdive growth records {ng} | bacdive phenotypes {npg}")
