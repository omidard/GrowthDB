import json, os, gzip, sys
from collections import defaultdict
sys.path.insert(0,"/tmp/claude-1000/-data-Brilliant-genomics-department/eb8d91f3-1707-45de-a10d-2de68fef6627/scratchpad/media_work/repo/tools")
from map_metabolite import Mapper
_m=Mapper(); _remap_cache={}
def _resolve(sub, chebi=None, kegg=None):
    key=(sub,chebi,kegg)
    if key in _remap_cache: return _remap_cache[key]
    h=_m.map(name=sub)
    if h: r=(h["exchange"],"bigg")   # accept ANY valid BiGG exchange (not only local-reactome subset)
    else:
        fb=_m.fallback_exchange(chebi=chebi,kegg=kegg)
        r=(fb["exchange"],fb["namespace"]) if fb else (None,None)
    _remap_cache[key]=r; return r
REPO="repo"
try: _ASSAY=json.load(open("/tmp/claude-1000/-data-Brilliant-genomics-department/eb8d91f3-1707-45de-a10d-2de68fef6627/scratchpad/media_work/repo/tools/assay_base_media.json"))
except Exception: _ASSAY={}
agg=defaultdict(lambda:{'pos':0,'neg':0,'other':0,'n':0,'ex':None,'ns':None,'cat':None,'cite':None,'src':set(),'kinds':set()})
def add(org,sub,ex,cat,ph,cite,src,kind=None,chebi=None,kegg=None):
    ns='bigg' if ex else None
    if not ex and sub: ex,ns=_resolve(sub,chebi,kegg)   # BiGG, else ModelSEED/KEGG fallback id
    k=(org.strip(), (sub or '').strip().lower(), ex or '')
    a=agg[k]; a['n']+=1; a['ex']=ex; a['ns']=a['ns'] or ns; a['cat']=a['cat'] or cat; a['cite']=a['cite'] or cite; a['src'].add(src)
    if kind: a['kinds'].add(kind)
    if ph=='positive': a['pos']+=1
    elif ph=='negative': a['neg']+=1
    else: a['other']+=1
    a['_org']=org; a['_sub']=sub
# paper biolog
for r in json.load(open(os.path.join(REPO,'data','biolog_phenotypes.json')))['phenotypes']:
    add(r['organism'],r['substrate'],r.get('exchange'),r.get('category'),r.get('phenotype'),r.get('citation'),'paper')
# bacdive (has ChEBI -> ModelSEED fallback for non-BiGG substrates)
for l in open('lit/bacdive/bacdive_phenotypes.jsonl'):
    d=json.loads(l); add(d['organism'],d['substrate'],d.get('exchange'),d.get('category'),d.get('phenotype'),d.get('citation'),'bacdive',d.get('kind'),chebi=d.get('chebi'))
# pmkbase (P. putida Biolog PM; has KEGG)
import os as _os
if _os.path.exists('lit/pmkbase/pmkbase_phenotypes.jsonl'):
    for l in open('lit/pmkbase/pmkbase_phenotypes.jsonl'):
        d=json.loads(l); add(d['organism'],d['substrate'],d.get('exchange'),d.get('category'),d.get('phenotype'),d.get('citation'),'pmkbase',d.get('kind'),kegg=d.get('kegg'))
# all-paper Biolog / substrate-utilisation mine (825 papers)
if _os.path.exists('lit/biolog_all_phenotypes.jsonl'):
    for l in open('lit/biolog_all_phenotypes.jsonl'):
        d=json.loads(l); add(d['organism'],d['substrate'],d.get('exchange'),d.get('category'),d.get('phenotype'),d.get('citation'),d.get('source','sub_util'),d.get('kind'))
def assay_context(kinds, sources):
    ks=set((k or '').lower() for k in (kinds or [])); src=set(sources or [])
    def has(*subs): return any(any(s in k for s in subs) for k in ks)
    if 'pmkbase' in src: return ('defined-minimal','Biolog Phenotype MicroArray — defined minimal medium (IF-0) with the substrate as the sole carbon source; positive = growth')
    # growth-relevant sole-substrate utilisation takes priority (the meaningful signal for GEMs)
    if has('assimil','carbon source','energy source','nitrogen source','respiration','oxid','required for growth','sole') or ('growth' in ks) or has('growth'):
        return ('defined-minimal','Sole-substrate growth / assimilation on a defined minimal medium (substrate as C/N/energy source) — positive = the organism grows on it')
    if has('acid','gas from','ferment'):
        return ('complex','Acid/gas-production test: peptone/complex base + substrate + pH indicator — positive = fermentative acid production, NOT growth on the substrate as a sole carbon source (interpret with care for model validation)')
    if has('hydrol','degrad','reduc','electron'):
        return ('assay','Enzyme / activity assay (hydrolysis, degradation, reduction, electron transfer) in a test medium — not a sole-substrate growth test')
    if 'paper' in src: return ('defined-minimal','Phenotype-microarray / substrate-utilisation assay (defined minimal + sole substrate)')
    return ('unknown','Substrate-utilisation assay; base medium per the original source')

rows=[]
for (org,sub,ex),a in agg.items():
    tot=a['pos']+a['neg']+a['other']
    cons='positive' if a['pos']>a['neg'] and a['pos']>0 else ('negative' if a['neg']>a['pos'] else 'variable')
    bt,bd=assay_context(a['kinds'],a['src'])
    bmid=_ASSAY.get(bt); bmurl=(f"https://omidard.github.io/Media/?medium={bmid}" if bmid else None)
    _s=a['src']
    platform=('Biolog PM' if ('pmkbase' in _s or 'biolog_pm' in _s) else ('Biolog/PM (literature)' if 'paper' in _s else ('Substrate utilisation (literature)' if 'sub_util' in _s else 'BacDive metabolite utilization')))
    rows.append({'organism':a['_org'],'substrate':a['_sub'],'exchange':a['ex'] or None,'exchange_ns':a.get('ns'),'category':a['cat'],
        'platform':platform,'base_media_id':bmid,'base_media_url':bmurl,
        'phenotype':cons,'n_strains':a['n'],'n_positive':a['pos'],'n_negative':a['neg'],
        'sources':sorted(a['src']),'kinds':sorted(a['kinds']),'base_medium_type':bt,'base_medium':bd,'citation':a['cite']})
# browser index: mapped-to-BiGG subset (most useful), sorted
from collections import Counter
# compact rows for drill-down (gzipped): {o,s,e,c,p,n}
compact=[{'o':r['organism'],'s':r['substrate'],'e':r['exchange'],'c':r['category'],'p':r['phenotype'],'n':r['n_strains']} for r in rows]
with gzip.open(os.path.join(REPO,'data','phenotypes.json.gz'),'wt') as f: json.dump(compact,f,separators=(',',':'))
# small summary: stats + per-organism counts (datalist + drill-down)
byorg=defaultdict(lambda:{'n':0,'pos':0})
for r in rows:
    byorg[r['organism']]['n']+=1
    if r['phenotype']=='positive': byorg[r['organism']]['pos']+=1
orglist=sorted(([o,v['n'],v['pos']] for o,v in byorg.items()), key=lambda x:-x[1])
summary={'count':len(rows),'n_mapped':sum(1 for r in rows if r['exchange']),'n_organisms':len(byorg),
     'n_raw_datapoints':881958+706,'by_phenotype':dict(Counter(r['phenotype'] for r in rows)),
     'by_category':dict(Counter(r['category'] for r in rows)),'by_source':dict(Counter(s for r in rows for s in r['sources'])),
     'top_organisms':orglist[:2000]}
json.dump(summary,open(os.path.join(REPO,'data','phenotypes_summary.json'),'w'),separators=(',',':'))
# full (with n_pos/n_neg/sources/citation) gzipped for download/programmatic use
with gzip.open(os.path.join(REPO,'data','phenotypes_full.json.gz'),'wt') as f: json.dump(rows,f,separators=(',',':'))
print('consensus phenotypes:',len(rows),'| BiGG-mapped:',summary['n_mapped'],'| organisms:',summary['n_organisms'])
print('sizes: phenotypes.json.gz',round(os.path.getsize(os.path.join(REPO,'data','phenotypes.json.gz'))/1e6,1),
      'MB | summary',round(os.path.getsize(os.path.join(REPO,'data','phenotypes_summary.json'))/1e6,2),
      'MB | full.gz',round(os.path.getsize(os.path.join(REPO,'data','phenotypes_full.json.gz'))/1e6,1),'MB')
print('by phenotype:',summary['by_phenotype'],'| by source:',summary['by_source'])
