#!/usr/bin/env python3
"""GrowthDB curation engine: controlled vocabularies (culture_mode, oxygen) + medium resolver
(canonical name, synonym collapse, dedup, deep link to a specific Media entry). Applied to every
record; free-text detail preserved. Reports before/after and any standard media missing from Media."""
import json, os, re
GW="/tmp/claude-1000/-data-Brilliant-genomics-department/eb8d91f3-1707-45de-a10d-2de68fef6627/scratchpad/growthdb_work"
REPO=os.path.join(GW,"repo")
MEDIA_REPO="/tmp/claude-1000/-data-Brilliant-genomics-department/eb8d91f3-1707-45de-a10d-2de68fef6627/scratchpad/media_work/repo"
MEDIA_IDS={m["id"] for m in json.load(open(os.path.join(MEDIA_REPO,"data","index.json")))["media"]}
BACKMAP=json.load(open(os.path.join(GW,"lit","media_backmap.json"))) if os.path.exists(os.path.join(GW,"lit","media_backmap.json")) else {}

# ---------- controlled culture_mode ----------
def culture_mode(raw):
    s=(raw or "").lower().strip()
    if not s or s in ("na","none","unspecified","other","unknown"): return "unspecified", (raw or None)
    order=[("turbidostat","turbidostat"),("retentostat","retentostat"),
           ("mother machine|mother-machine|microfluidic|single-cell|single cell|time-lapse microscop","single-cell"),
           ("chemostat|continuous cultur|continuous exponential|constant dilution|steady.?state.*dilution|\\bd\\s*=","chemostat"),
           ("fed-batch|fed batch|fedbatch","fed-batch"),
           ("co-cultur|cocultur|community|consortium","co-culture"),
           ("in vivo|in-host|in host|infection|coloniz|gnotobiotic|host cell|intracellular|cecal|gut lumen","in vivo"),
           ("in situ|environment|metagenom|epilimnion|sediment|field|enrichment","environmental"),
           ("batch|flask|shake|tube|bottle|microcosm|deep-well|deep well|microtiter|microplate|well plate|serum|plate reader|agar|colony|solid medium","batch")]
    for pat,mode in order:
        if re.search(pat,s): return mode, (raw if raw.strip().lower()!=mode else None)
    return "other", raw

# ---------- controlled oxygen ----------
def oxygen(raw):
    s=(raw or "").lower().strip()
    if not s or s in ("na","none"): return None, None
    if "not stated" in s or "not specified" in s or "unknown" in s: return "not_stated", (raw or None)
    if "facultative" in s: return "facultative", (raw if s!="facultative" else None)
    if "microaero" in s: return "microaerophilic", (raw if s not in("microaerophilic","microaerobic") else None)
    if "anoxic" in s or "anaerob" in s or "n2 headspace" in s or "oxygen-free" in s: return "anaerobic", (raw if s!="anaerobic" else None)
    if "aerobic" in s or "oxic" in s or "shaking" in s or "rpm" in s or "ambient" in s or "aerated" in s or "aeration" in s: return "aerobic", (raw if s!="aerobic" else None)
    return "other", raw

# ---------- medium resolver ----------
STD=[  # (regex, canonical_base, media_id or None, is_minimal)
 (r"lysogeny|luria[\s\-]*bertani|luria broth|lennox|miller.?s? lb|\blb\b|\blbk\b|\blb0?\b","LB","lb_lennox",False),
 (r"brain[\s\-]*heart|\bbhis?\b","BHI","bhi",False),
 (r"tryptic[\s\-]*soy|trypticase[\s\-]*soy|\btsb\b|\btsa\b","TSB / TSA","tsb",False),
 (r"\bmops\b","MOPS minimal","mops_minimal_glucose",True),
 (r"\bm9\b","M9 minimal",None,True),
 (r"\bm63\b","M63 minimal",None,True),
 (r"davis","Davis minimal","davis_minimal_glucose",True),
 (r"gutnick","Gutnick minimal",None,True),
 (r"nutrient[\s\-]*(broth|agar|medium|bouillon)","Nutrient broth",None,False),
 (r"terrific","Terrific broth",None,False),
 (r"2\s*[x×]\s*yt|2xyt|double.?strength yt","2xYT",None,False),
 (r"\br2a\b","R2A",None,False),
 (r"columbia","Columbia blood agar","blood_agar_columbia",False),
 (r"marine[\s\-]*broth|zobell","Marine broth (ZoBell)",None,False),
 (r"reinforced clostridial|\brcm\b","Reinforced clostridial medium",None,False),
 (r"\bycfa\b|\byca\b","YCFA",None,False),
 (r"\bmrs\b|de man.?rogosa","MRS",None,False),
]
CARBON=[("glucose|dextrose|d-glc","glucose"),("glycerol","glycerol"),("acetate|acetic","acetate"),
 ("succinate|succinic","succinate"),("lactate|lactic","lactate"),("pyruvate","pyruvate"),("fructose","fructose"),
 ("xylose","xylose"),("galactose","galactose"),("mannitol","mannitol"),("sucrose","sucrose"),("citrate","citrate"),
 ("gluconate","gluconate"),("ethanol","ethanol"),("methanol","methanol"),("arabinose","arabinose"),("maltose","maltose"),
 ("mannose","mannose"),("sorbitol","sorbitol"),("glutamate","glutamate"),("casamino|casein hydrolys","casamino acids")]
def carbon_of(s):
    for pat,c in CARBON:
        if re.search(pat,s): return c
    return None
def norm_key(name):
    n=(name or "").lower()
    n=re.sub(r"\(.*?\)","",n)                        # drop parentheticals like (LB)
    n=re.sub(r"\b(medium|media|broth|agar|liquid|defined|chemically|conventional|standard|modified|supplemented|with|the)\b"," ",n)
    n=re.sub(r"[^a-z0-9]+"," ",n).strip()
    return re.sub(r"\s+"," ",n)
def resolve_media(name, pmcid=""):
    if not name: return None, None, None
    s=name.lower()
    for pat,base,mid,minimal in STD:
        if re.search(pat,s):
            if minimal:
                c=carbon_of(s)
                canon=f"{base} ({c})" if c else base
                # M9/M63 carbon-specific link
                if base.startswith("M9") and c:
                    cand=f"m9_{c}_aerobic";
                    if cand in MEDIA_IDS: mid=cand
                    elif f"m9_{c}" in MEDIA_IDS: mid=f"m9_{c}"
                return canon, norm_key(canon), (mid if mid in MEDIA_IDS else None)
            return base, norm_key(base), (mid if mid and mid in MEDIA_IDS else None)
    # non-standard: canonical = cleaned name; link via backmap (composition-added media)
    bm=BACKMAP.get(f"{pmcid}|{name}")
    return name.strip(), norm_key(name), (bm if bm in MEDIA_IDS else None)

# ---------- apply ----------
recs=json.load(open(os.path.join(REPO,"data","growth_records.json")))
from collections import Counter
before_cm=Counter(); after_cm=Counter(); canon_media=set(); linked=0; namednolink=0; missing_std=Counter()
STD_IDS={mid for _,_,mid,_ in STD if mid}
for r in recs:
    c=r.setdefault("conditions",{})
    lit=r["id"].startswith("lit_")
    if lit: before_cm[c.get("culture_mode") or "NULL"]+=1
    cm,cmd=culture_mode(c.get("culture_mode"))
    c["culture_mode"]=cm; c["culture_detail"]=cmd
    ox,oxd=oxygen(c.get("oxygen"))
    c["oxygen"]=ox; c["aeration_detail"]=oxd
    if lit: after_cm[cm]+=1
    med=r.setdefault("medium",{})
    nm=med.get("description")
    canon,key,mid=resolve_media(nm, r.get("provenance",{}).get("pmcid",""))
    med["canonical_name"]=canon; med["canonical_key"]=key
    if mid: med["media_id"]=mid
    med["media_url"]=(f"https://omidard.github.io/Media/?medium={med['media_id']}" if med.get("media_id") else None)
    if key: canon_media.add(key)
    if lit and nm:
        if med.get("media_id"): linked+=1
        else: namednolink+=1
        # track standard media families that SHOULD have a Media entry but don't link
        cbase=canon.split(" (")[0]
        if med.get("media_id") is None and any(re.search(p,nm.lower()) for p,_,mid2,_ in STD if mid2 is None):
            missing_std[cbase]+=1
json.dump(recs, open(os.path.join(REPO,"data","growth_records.json"),"w"), separators=(",",":"))
print("=== culture_mode: before",len(before_cm),"distinct -> after",len(after_cm),"controlled ===")
for k,v in after_cm.most_common(): print(f"   {v:5d}  {k}")
print("\n=== medium: 695 raw names -> ",len(canon_media),"canonical keys ===")
print(f"literature media LINKED to a Media entry: {linked} | named but still unlinked: {namednolink}")
print("\n=== standard-media families reported but MISSING a Media entry (should add):")
for k,v in missing_std.most_common(12): print(f"   {v:4d}  {k}")
