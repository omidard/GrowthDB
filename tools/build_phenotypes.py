import json, os, gzip
from collections import defaultdict
REPO="repo"
agg=defaultdict(lambda:{'pos':0,'neg':0,'other':0,'n':0,'ex':None,'cat':None,'cite':None,'src':set(),'kinds':set()})
def add(org,sub,ex,cat,ph,cite,src,kind=None):
    k=(org.strip(), (sub or '').strip().lower(), ex or '')
    a=agg[k]; a['n']+=1; a['ex']=ex; a['cat']=a['cat'] or cat; a['cite']=a['cite'] or cite; a['src'].add(src)
    if kind: a['kinds'].add(kind)
    if ph=='positive': a['pos']+=1
    elif ph=='negative': a['neg']+=1
    else: a['other']+=1
    a['_org']=org; a['_sub']=sub
# paper biolog
for r in json.load(open(os.path.join(REPO,'data','biolog_phenotypes.json')))['phenotypes']:
    add(r['organism'],r['substrate'],r.get('exchange'),r.get('category'),r.get('phenotype'),r.get('citation'),'paper')
# bacdive
for l in open('lit/bacdive/bacdive_phenotypes.jsonl'):
    d=json.loads(l); add(d['organism'],d['substrate'],d.get('exchange'),d.get('category'),d.get('phenotype'),d.get('citation'),'bacdive',d.get('kind'))
rows=[]
for (org,sub,ex),a in agg.items():
    tot=a['pos']+a['neg']+a['other']
    cons='positive' if a['pos']>a['neg'] and a['pos']>0 else ('negative' if a['neg']>a['pos'] else 'variable')
    rows.append({'organism':a['_org'],'substrate':a['_sub'],'exchange':a['ex'] or None,'category':a['cat'],
        'phenotype':cons,'n_strains':a['n'],'n_positive':a['pos'],'n_negative':a['neg'],
        'sources':sorted(a['src']),'kinds':sorted(a['kinds']),'citation':a['cite']})
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
