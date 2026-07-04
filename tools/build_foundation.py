#!/usr/bin/env python3
"""GrowthDB foundation: prokaryote backbone (all GTDB species) + Madin-derived growth records.
Growth records = species-level growth rate (mu from doubling time) + culture conditions, cited.
Media/uptake/secretion (from literature) come in a later phase."""
import csv, json, os, math, re
W="/tmp/claude-1000/-data-Brilliant-genomics-department/eb8d91f3-1707-45de-a10d-2de68fef6627/scratchpad/growthdb_work"
REPO=os.path.join(W,"repo"); os.makedirs(os.path.join(REPO,"data"),exist_ok=True)

# ---- prokaryote backbone (GTDB species) ----
orgs={}
def add_species(line, sk):
    s=line.strip()[3:] if line.startswith("s__") else line.strip()
    if not s: return
    genus=s.split()[0]
    named = not re.search(r"\bsp\d{5,}$", s) and not re.match(r"^[A-Z0-9-]+ ", s)  # crude 'named' flag
    orgs[s]={"species":s,"genus":genus,"superkingdom":sk,"gtdb":True,"named":named,"has_growth":False}
for line in open(os.path.join(W,"gtdb_species_bac.txt")): add_species(line,"Bacteria")
for line in open(os.path.join(W,"gtdb_species_arc.txt")): add_species(line,"Archaea")
print("GTDB backbone species:", len(orgs), "| bacteria:", sum(1 for o in orgs.values() if o['superkingdom']=='Bacteria'), "| archaea:", sum(1 for o in orgs.values() if o['superkingdom']=='Archaea'))

# ---- Madin growth records ----
def num(v):
    try: return float(v)
    except: return None
records=[]; n_dbl=0; matched=0; extra_orgs=0
csv.field_size_limit(10**7)
with open(os.path.join(W,"madin_traits.csv")) as f:
    for row in csv.DictReader(f):
        sp=row.get("species","").strip()
        if not sp or sp in ("NA","sp."): continue
        dbl=num(row.get("doubling_h")); gt=num(row.get("growth_tmp")); ot=num(row.get("optimum_tmp"))
        oph=num(row.get("optimum_ph")); metab=row.get("metabolism","")
        cs=row.get("carbon_substrates","")
        has=any(x not in (None,"") for x in [dbl,gt,ot,oph]) or metab not in ("","NA") or cs not in ("","NA")
        if not has: continue
        # link to backbone
        if sp in orgs: orgs[sp]["has_growth"]=True; matched+=1
        else:
            orgs[sp]={"species":sp,"genus":row.get("genus",""),"superkingdom":row.get("superkingdom","") or "Bacteria",
                      "gtdb":False,"named":True,"has_growth":True}; extra_orgs+=1
        mu = round(math.log(2)/dbl,4) if (dbl and dbl>0) else None
        if mu is not None: n_dbl+=1
        oxy={"aerobic":"aerobic","obligate aerobic":"aerobic","anaerobic":"anaerobic","obligate anaerobic":"anaerobic",
             "facultative":"facultative","microaerophilic":"microaerophilic"}.get((metab or "").lower(), (metab if metab not in ("","NA") else None))
        rec={"id":f"madin_{row.get('tax_id')}_{len(records)}",
             "organism":sp,"genus":row.get("genus",""),"family":row.get("family",""),"phylum":row.get("phylum",""),
             "superkingdom":row.get("superkingdom",""),"ncbi_tax_id":row.get("tax_id"),"gtdb_species":(sp if sp in orgs and orgs[sp]["gtdb"] else None),
             "growth_rate_per_h":mu,"doubling_time_h":dbl,
             "conditions":{"culture_mode":"unspecified","oxygen":oxy,
                 "temperature_C":gt,"optimum_temperature_C":ot,"optimum_pH":oph,"pH":None,
                 "temperature_range":(row.get("range_tmp") if row.get("range_tmp") not in ("","NA") else None)},
             "carbon_substrates":[c.strip() for c in cs.split(",") if c.strip() and c!="NA"] if cs not in ("","NA") else [],
             "medium":{"media_id":None,"description":None,"note":"Madin trait synthesis does not record the full medium; to be linked/added via literature curation."},
             "uptake_rates":[],"secretion_rates":[],
             "provenance":{"source_type":"database",
                 "citation":f"Madin JS et al. A synthesis of bacterial and archaeal phenotypic trait data. Sci Data 7:170 (2020); source dataset: {row.get('data_source')}.",
                 "doi":"10.1038/s41597-020-0497-4","data_source":row.get("data_source"),
                 "url":"https://github.com/bacteria-archaea-traits/bacteria-archaea-traits","method":("reported" )},
             "curation_notes":("growth_rate_per_h computed as ln(2)/doubling_time_h. " if mu is not None else "")+"Species-level trait record; culture mode/pH often unspecified; medium + uptake/secretion pending literature curation."}
        records.append(rec)
print(f"Madin growth records: {len(records)} | with doubling->mu: {n_dbl} | matched to GTDB: {matched} | NCBI-only orgs added: {extra_orgs}")
print("total organisms (backbone + NCBI-only):", len(orgs), "| with growth data:", sum(1 for o in orgs.values() if o['has_growth']))

json.dump(list(orgs.values()), open(os.path.join(REPO,"data","organisms.json"),"w"), separators=(",",":"))
json.dump(records, open(os.path.join(REPO,"data","growth_records.json"),"w"), separators=(",",":"))
from collections import Counter
sk=Counter(o["superkingdom"] for o in orgs.values())
idx={"n_organisms":len(orgs),"n_with_growth":sum(1 for o in orgs.values() if o['has_growth']),
     "n_growth_records":len(records),"n_with_rate":n_dbl,"by_superkingdom":dict(sk),
     "sources":{"backbone":"GTDB r220 (bac120 + ar53)","growth":"Madin et al. Sci Data 2020"}}
json.dump(idx, open(os.path.join(REPO,"data","index.json"),"w"), indent=1)
print("wrote data/organisms.json, growth_records.json, index.json")
print("sizes:", {f:round(os.path.getsize(os.path.join(REPO,'data',f))/1e6,1) for f in ['organisms.json','growth_records.json']})
