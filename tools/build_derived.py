#!/usr/bin/env python3
"""Data-science / GEM-modeller layer: derive per-species metrics from the raw records and store
them in each species shard + the index. Computes: biomass yield & maintenance (NGAM) from the
qS-vs-µ chemostat relation, overflow-metabolism byproducts (qP vs µ), cardinal/optimum growth
temperature, oxygen preference, and the substrate-utilization spectrum for GEM validation."""
import json, os, glob
from collections import defaultdict, Counter
REPO="/tmp/claude-1000/-data-Brilliant-genomics-department/eb8d91f3-1707-45de-a10d-2de68fef6627/scratchpad/growthdb_work/repo"
CARBON={"EX_glc__D_e":"glucose","EX_glyc_e":"glycerol","EX_ac_e":"acetate","EX_succ_e":"succinate",
 "EX_lac__L_e":"lactate","EX_xyl__D_e":"xylose","EX_fru_e":"fructose","EX_gal_e":"galactose",
 "EX_malt_e":"maltose","EX_glcn_e":"gluconate","EX_pyr_e":"pyruvate","EX_cit_e":"citrate",
 "EX_etoh_e":"ethanol","EX_meoh_e":"methanol","EX_arab__L_e":"arabinose","EX_man_e":"mannose"}
def lsq(pts):
    n=len(pts)
    if n<3: return None
    sx=sum(p[0] for p in pts); sy=sum(p[1] for p in pts)
    sxx=sum(p[0]**2 for p in pts); sxy=sum(p[0]*p[1] for p in pts)
    d=n*sxx-sx*sx
    if abs(d)<1e-9: return None
    a=(n*sxy-sx*sy)/d; b=(sy-a*sx)/n
    ybar=sy/n; sstot=sum((p[1]-ybar)**2 for p in pts); ssres=sum((p[1]-(a*p[0]+b))**2 for p in pts)
    r2=1-ssres/sstot if sstot>1e-9 else None
    return {"slope":a,"intercept":b,"n":n,"r2":r2}

idx=json.load(open(os.path.join(REPO,"data","species_index.json")))
byslug={x["slug"]:x for x in idx["species"]}
n_maint=n_over=n_card=0
for f in glob.glob(os.path.join(REPO,"data","species","*.json")):
    sh=json.load(open(f)); grecs=sh["growth_records"]; phs=sh["phenotypes"]
    mus=[r for r in grecs if r.get("mu") is not None]
    d={}
    if mus:
        d["mu_max"]=round(max(r["mu"] for r in mus),4)
        d["td_min_h"]=round(0.6931/d["mu_max"],3) if d["mu_max"]>0 else None
    # cardinal / optimum temperature
    tmus=[(r["T"],r["mu"]) for r in mus if r.get("T") is not None]
    if len({t for t,_ in tmus})>=3:
        topt=max(tmus,key=lambda x:x[1])
        d["temp_optimum_C"]=topt[0]; d["mu_at_topt"]=round(topt[1],4)
        d["temp_range_C"]=[min(t for t,_ in tmus),max(t for t,_ in tmus)]; n_card+=1
    optT=[r["optT"] for r in grecs if r.get("optT") is not None]
    if optT: d["reported_optimum_temp_C"]=round(sum(optT)/len(optT),1)
    optpH=[r["optpH"] for r in grecs if r.get("optpH") is not None]
    if optpH: d["reported_optimum_pH"]=round(sum(optpH)/len(optpH),1)
    # oxygen preference
    oxc=Counter(r["ox"] for r in grecs if r.get("ox"))
    if oxc: d["oxygen_preference"]=oxc.most_common(1)[0][0]
    # biomass yield + maintenance from qS vs mu (carbon uptake)
    qs=[]
    for r in mus:
        for u in (r.get("up") or []):
            if u.get("exchange") in CARBON and u.get("rate") is not None:
                qs.append((r["mu"],abs(u["rate"]))); break
    fit=lsq(qs)
    if len(qs)>=3:
        d["qs_mu_points"]=[{"mu":round(m,4),"qS":round(q,3)} for m,q in qs]
        # only report a maintenance/yield MODEL if the qS-vs-µ relation is genuinely linear
        # (a valid derivation needs a coherent chemostat series, not pooled cross-study points)
        if fit and fit["slope"]>1e-6 and fit["r2"] is not None and fit["r2"]>=0.6:
            d["biomass_yield_gDW_per_mmol"]=round(1/fit["slope"],4)
            d["maintenance_NGAM_mmol_gDW_h"]=round(max(0,fit["intercept"]),4)
            d["yield_fit"]={"n":fit["n"],"r2":round(fit["r2"],3),"slope":round(fit["slope"],4),"intercept":round(fit["intercept"],4),
                            "note":"qS = (1/Yxs)·µ + m_s  — linear fit of carbon-source uptake vs growth rate"}
            n_maint+=1
        else:
            d["qs_mu_note"]="pooled across studies/media — not a single chemostat series, so no maintenance model fitted"
    # overflow metabolism: secreted products increasing with mu
    prod=defaultdict(list)
    for r in mus:
        for s in (r.get("sec") or []):
            if s.get("exchange") and s.get("rate") is not None: prod[s["exchange"]].append((r["mu"],abs(s["rate"])))
    overflow=[]
    for ex,pts in prod.items():
        fp=lsq(pts)
        if fp and fp["slope"]>0.05 and len(pts)>=3:
            overflow.append({"exchange":ex,"slope":round(fp["slope"],3),"n":len(pts),
                "points":[{"mu":round(m,4),"qP":round(q,3)} for m,q in pts]})
    if overflow: d["overflow_products"]=sorted(overflow,key=lambda x:-x["slope"]); n_over+=1
    # substrate utilization spectrum (for GEM validation)
    pos=[p for p in phs if p.get("phenotype")=="positive"]; neg=[p for p in phs if p.get("phenotype")=="negative"]
    d["n_substrates_positive"]=len(pos); d["n_substrates_negative"]=len(neg)
    d["utilizable_exchanges"]=sorted({p["exchange"] for p in pos if p.get("exchange")})
    if pos: d["substrate_categories"]=dict(Counter(p.get("category") for p in pos))
    sh["derived"]=d
    json.dump(sh,open(f,"w"),separators=(",",":"))
    # index summary
    e=byslug.get(sh.get("species") and next((x["slug"] for x in idx["species"] if x["s"]==sh["species"]),None))
    sl=os.path.basename(f)[:-5]
    if sl in byslug:
        byslug[sl]["mumax"]=d.get("mu_max"); byslug[sl]["topt"]=d.get("temp_optimum_C") or d.get("reported_optimum_temp_C")
        byslug[sl]["yxs"]=d.get("biomass_yield_gDW_per_mmol"); byslug[sl]["ngam"]=d.get("maintenance_NGAM_mmol_gDW_h")
        byslug[sl]["over"]=bool(d.get("overflow_products")); byslug[sl]["nsub"]=d.get("n_substrates_positive",0)
idx["n_with_maintenance"]=n_maint; idx["n_with_overflow"]=n_over; idx["n_with_cardinalT"]=n_card
json.dump(idx,open(os.path.join(REPO,"data","species_index.json"),"w"),separators=(",",":"))
print(f"derived metrics: maintenance/yield {n_maint} species | overflow {n_over} | cardinal-T {n_card}")
