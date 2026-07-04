import json
from collections import Counter
recs=json.load(open('data/growth_records.json'))
orgs=json.load(open('data/organisms.json'))
def src(r): return 'Literature' if r['id'].startswith('lit_') else 'Madin traits'
out=[]
QUANT_ONLY=True
for r in recs:
    lit=r['id'].startswith('lit_'); c=r.get('conditions',{})
    row={'id':r['id'],'org':r['organism'],'strain':r.get('strain') if lit else None,'sk':r.get('superkingdom') or '',
         'mu':r.get('growth_rate_per_h'),'td':r.get('doubling_time_h'),
         'mode':c.get('culture_mode') or 'unspecified','ox':c.get('oxygen'),
         'T':c.get('temperature_C'),'pH':c.get('pH'),'optT':c.get('optimum_temperature_C'),'optpH':c.get('optimum_pH'),
         'dil':c.get('dilution_rate_per_h'),
         'nu':len(r.get('uptake_rates') or []),'ns':len(r.get('secretion_rates') or []),
         'mid':r.get('medium',{}).get('media_id'),'med':r.get('medium',{}).get('description'),
         'src':src(r),'cite':(r.get('provenance',{}).get('citation') or '')[:150],
         'doi':r.get('provenance',{}).get('doi'),'pmcid':r.get('provenance',{}).get('pmcid'),
         'conf':r.get('confidence')}
    if lit:  # include full detail for literature records
        row['up']=r.get('uptake_rates'); row['sec']=r.get('secretion_rates')
        row['snip']=r.get('provenance',{}).get('snippet'); row['notes']=r.get('curation_notes')
        row['method']=r.get('provenance',{}).get('method'); row['cs']=r.get('carbon_substrates')
    if QUANT_ONLY and (r.get("growth_rate_per_h") is None and not r.get("uptake_rates") and not r.get("secretion_rates")): continue
    out.append(row)
idx={'count':len(out),'n_organisms':len(orgs),'n_with_growth':sum(1 for o in orgs if o['has_growth']),
     'by_superkingdom':dict(Counter(o['superkingdom'] for o in orgs)),
     'by_source':dict(Counter(src(r) for r in recs)),
     'n_with_rate':sum(1 for r in recs if r.get('growth_rate_per_h') is not None),
     'n_uptake':sum(1 for r in recs if r.get('uptake_rates')),'n_secretion':sum(1 for r in recs if r.get('secretion_rates')),
     'n_media_linked':sum(1 for r in recs if r.get('medium',{}).get('media_id')),
     'records':out}
json.dump(idx,open('data/records_index.json','w'),separators=(',',':'))
import os;print('records_index.json:',round(os.path.getsize('data/records_index.json')/1e6,1),'MB |',len(out),'records')
print('with rate:',idx['n_with_rate'],'| uptake:',idx['n_uptake'],'| secretion:',idx['n_secretion'],'| by source:',idx['by_source'])
