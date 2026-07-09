#!/usr/bin/env python3
"""Species-first reorganization of GrowthDB: dedup, add a normalized species key, group all
growth records + phenotypes by species into per-species shards + a species index for the browser."""
import json, gzip, os, re
from collections import defaultdict
REPO="repo"
def species_key(o):
    o=(o or "").strip()
    o=re.sub(r'\s+(strain|str\.?|subsp\.?|serovar|biovar|pv\.?|ATCC|DSM|NCTC|NCIMB|CECT|LMG|JCM|KCTC|CIP|NBRC|BCRC|PCC|MG1655|W3110|BL21|K-?12|O157).*','',o,flags=re.I)
    m=re.match(r'^([A-Z][a-z]+ [a-z]+)',o)
    if m: return m.group(1)
    m2=re.match(r'^(Candidatus [A-Z][a-z]+ [a-z]+)',o)
    return m2.group(1) if m2 else o
def slug(s): return re.sub(r'[^a-z0-9]+','_',(s or '').lower()).strip('_')[:80] or 'unknown'

# superkingdom lookup
skmap={}
for o in json.load(open(os.path.join(REPO,"data","organisms.json"))): skmap[o["species"]]=o["superkingdom"]

recs=json.load(open(os.path.join(REPO,"data","growth_records.json")))
# dedup literature records + attach species
seen=set(); out=[]
for r in recs:
    if r["id"].startswith(("lit_","flux_")):
        k=(r["organism"],(r.get("medium") or {}).get("description"),r.get("growth_rate_per_h"),
           (r.get("provenance") or {}).get("pmcid"),
           tuple(sorted(f"{u.get('compound')}|{u.get('rate')}" for u in (r.get("uptake_rates") or []))))
        if k in seen: continue
        seen.add(k)
    r["species"]=species_key(r["organism"])
    out.append(r)
print(f"records {len(recs)} -> {len(out)} after dedup ({len(recs)-len(out)} removed)")
json.dump(out,open(os.path.join(REPO,"data","growth_records.json"),"w"),separators=(",",":"))

bysp_rec=defaultdict(list)
for r in out: bysp_rec[r["species"]].append(r)
# phenotypes by species (full consensus set, with counts + citation)
bysp_ph=defaultdict(list)
with gzip.open(os.path.join(REPO,"data","phenotypes_full.json.gz"),"rt") as f:
    for row in json.load(f):
        bysp_ph[species_key(row["organism"])].append(row)

os.makedirs(os.path.join(REPO,"data","species"),exist_ok=True)
def comp_rec(r):
    c=r.get("conditions",{}); med=r.get("medium",{}); p=r.get("provenance",{})
    return {"id":r["id"],"strain":r.get("strain"),"mu":r.get("growth_rate_per_h"),"td":r.get("doubling_time_h"),
        "mode":c.get("culture_mode"),"cdet":c.get("culture_detail"),"ox":c.get("oxygen"),"T":c.get("temperature_C"),
        "pH":c.get("pH"),"optT":c.get("optimum_temperature_C"),"optpH":c.get("optimum_pH"),"dil":c.get("dilution_rate_per_h"),
        "oxdet":c.get("aeration_detail"),"agit":c.get("agitation"),"meta_methods":(r.get("provenance") or {}).get("metadata_from_methods"),
        "medium":med.get("canonical_name") or med.get("description"),"media_id":med.get("media_id"),"media_url":med.get("media_url"),
        "up":r.get("uptake_rates"),"sec":r.get("secretion_rates"),"mfa":r.get("mfa_fluxes"),"yields":r.get("yields"),
        "cs":r.get("carbon_substrates"),"src":("Literature" if r["id"].startswith(("lit_","flux_")) else "Madin traits"),
        "cite":p.get("citation"),"doi":p.get("doi"),"pmcid":p.get("pmcid"),"method":p.get("method"),
        "snippet":p.get("snippet"),"notes":r.get("curation_notes"),"conf":r.get("confidence")}
def comp_ph(r):
    return {"substrate":r["substrate"],"exchange":r.get("exchange"),"category":r.get("category"),
        "phenotype":r.get("phenotype"),"n_strains":r.get("n_strains"),"n_positive":r.get("n_positive"),
        "n_negative":r.get("n_negative"),"sources":r.get("sources"),"kinds":r.get("kinds"),
        "base_type":r.get("base_medium_type"),"base_medium":r.get("base_medium"),"citation":r.get("citation")}

index=[]; used=set()
for sp in sorted(set(bysp_rec)|set(bysp_ph)):
    grecs=bysp_rec.get(sp,[]); phs=bysp_ph.get(sp,[])
    mus=[r["growth_rate_per_h"] for r in grecs if r.get("growth_rate_per_h") is not None]
    sl=slug(sp); base=sl; k=2
    while sl in used: sl=f"{base}_{k}"; k+=1
    used.add(sl)
    phc=[comp_ph(r) for r in phs]
    # package phenotypes BY assay medium (all substrate calls tested on the same base medium)
    from collections import defaultdict as _dd
    _g=_dd(list)
    for p in phc: _g[(p.get("base_type") or "unknown", p.get("base_medium") or "")].append(p)
    def _lab(bt):
        return {"defined-minimal":"Defined minimal medium (sole-substrate growth)","complex":"Complex/peptone base (acid-production test)",
                "assay":"Enzyme/activity assay medium","unknown":"Unspecified base medium"}.get(bt,bt)
    pheno_groups=[{"base_type":k[0],"label":_lab(k[0]),"base_medium":k[1],"n":len(v),
                   "n_pos":sum(1 for x in v if x["phenotype"]=="positive"),"n_neg":sum(1 for x in v if x["phenotype"]=="negative"),
                   "calls":sorted(v,key=lambda x:(x["phenotype"]!="positive",x["substrate"]))}
                  for k,v in sorted(_g.items(),key=lambda kv:-len(kv[1]))]
    shard={"species":sp,"superkingdom":skmap.get(sp,grecs[0].get("superkingdom") if grecs else ""),
           "n_records":len(grecs),"n_mu":len(mus),"n_phenotypes":len(phs),
           "growth_records":[comp_rec(r) for r in grecs],"phenotypes":phc,"phenotype_groups":pheno_groups}
    json.dump(shard,open(os.path.join(REPO,"data","species",sl+".json"),"w"),separators=(",",":"))
    index.append({"s":sp,"sk":shard["superkingdom"],"ng":len(grecs),"nmu":len(mus),
                  "mumax":(round(max(mus),3) if mus else None),"nu":sum(1 for r in grecs if r.get("uptake_rates")),
                  "np":len(phs),"slug":sl})
index.sort(key=lambda x:(-(x["nmu"] or 0),-(x["ng"] or 0)))
from collections import Counter
json.dump({"count":len(index),"n_with_rate":sum(1 for x in index if x["nmu"]),
    "n_with_pheno":sum(1 for x in index if x["np"]),"by_superkingdom":dict(Counter(x["sk"] for x in index)),
    "species":index},open(os.path.join(REPO,"data","species_index.json"),"w"),separators=(",",":"))
print(f"species shards: {len(index)} | with growth rate: {sum(1 for x in index if x['nmu'])} | with phenotypes: {sum(1 for x in index if x['np'])}")
print("species_index.json:",round(os.path.getsize(os.path.join(REPO,'data','species_index.json'))/1e6,2),"MB")
