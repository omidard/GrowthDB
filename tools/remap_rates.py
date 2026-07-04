import json, sys
sys.path.insert(0,"/tmp/claude-1000/-data-Brilliant-genomics-department/eb8d91f3-1707-45de-a10d-2de68fef6627/scratchpad/media_work/repo/tools")
from map_metabolite import Mapper
m=Mapper()
recs=json.load(open('repo/data/growth_records.json'))
n=0; before=0; after=0
for r in recs:
    if not r['id'].startswith(('lit_','flux_')): continue
    for key in ('uptake_rates','secretion_rates'):
        for f in (r.get(key) or []):
            if f.get('exchange'): after+=1; continue
            before+=1
            hit=m.map(name=f.get('compound'))
            if hit and hit['in_biggr']:
                f['exchange']=hit['exchange']; f['bigg_metabolite']=hit['bigg_metabolite']; n+=1; after+=1
json.dump(recs,open('repo/data/growth_records.json','w'),separators=(',',':'))
print(f'rates re-mapped: newly mapped {n} | now mapped {after} of {before+after-n+n} total')
