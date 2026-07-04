# GrowthDB — design & data model

Experimental **growth data for prokaryotes**, curated for constraint-based metabolic
modelling: measured **growth rates**, substrate **uptake rates** and product **secretion
rates**, together with the exact **medium** (linked to the [`Media`](https://github.com/omidard/Media)
repo) and the **culture conditions** (batch vs chemostat, aeration, temperature, pH) — every
value referenced to its primary source.

> **Why this exists.** GEM validation and flux-constraint work needs measured rates under
> known media and conditions. That data is scattered across thousands of papers in
> inconsistent units and rarely machine-readable. GrowthDB collects it once, with the medium
> resolved to BiGG exchanges (via the Media repo) so a rate can be dropped straight onto an
> exchange reaction as a constraint.

## Two record types

### 1) Organism (the backbone) — `data/organisms.json`
Every prokaryote, whether or not growth data exists yet:
```jsonc
{ "species":"Escherichia coli", "genus":"Escherichia", "superkingdom":"Bacteria",
  "gtdb":true, "named":true, "has_growth":true }
```
Backbone = **GTDB r220** species clusters (bacteria `bac120` + archaea `ar53`), plus named
species that appear in growth data but aren't GTDB clusters.

### 2) Growth record — `data/growth_records.json`
One measured/curated growth observation:
```jsonc
{
  "id":"...", "organism":"Vibrio natriegens", "genus":..., "phylum":..., "superkingdom":"Bacteria",
  "ncbi_tax_id":"...", "gtdb_species":"...",
  "growth_rate_per_h": 4.33,          // µ (h⁻¹); computed as ln(2)/doubling_time when not reported
  "doubling_time_h": 0.16,
  "conditions": {
    "culture_mode": "batch | chemostat | fed-batch | unspecified",
    "oxygen": "aerobic | anaerobic | facultative | microaerophilic",
    "temperature_C": 37, "optimum_temperature_C": null, "pH": 7.0, "optimum_pH": null
  },
  "carbon_substrates": ["D-glucose", ...],
  "medium": { "media_id":"m9_glucose_aerobic", "description":"...", "note":"..." },   // -> Media repo id
  "uptake_rates":   [ { "compound":"D-glucose", "exchange":"EX_glc__D_e", "rate":-10.5, "units":"mmol/gDW/h", "method":"reported|calculated" } ],
  "secretion_rates":[ { "compound":"acetate",   "exchange":"EX_ac_e",     "rate":  4.2, "units":"mmol/gDW/h", "method":"reported|calculated" } ],
  "provenance": { "source_type":"literature|database", "citation":"...", "doi":"...", "pmcid":"...", "method":"reported|calculated", "snippet":"..." },
  "curation_notes": "how rates were derived (e.g. µ from OD slope; qS from yield × µ)"
}
```

## Rules

- **Media are pointers into the `Media` repo** (`media_id`). If a paper's medium is not yet in
  `Media`, it is **added there** first (so `Media` becomes exhaustive) and then linked here.
- **Rates carry a `method`**: `reported` (stated in the paper) or `calculated` (derived by
  curation from the data the paper provides — e.g. µ from an OD/CFU time-course, specific uptake
  qS = (µ · Yxs⁻¹) or from residual-substrate + biomass, secretion from product-yield). The
  derivation is recorded in `curation_notes` — never invented.
- **Units are explicit** (`mmol·gDW⁻¹·h⁻¹` for fluxes, `h⁻¹` for µ). Rates are signed to match
  BiGG exchange convention (uptake negative, secretion positive) where a mapping exists.
- **Every value is referenced** to its primary source (DOI/PMCID) — database-derived seeds cite
  the database *and* its underlying source.

## Sourcing roadmap
1. **Backbone** — all GTDB prokaryote species. *(done)*
2. **Trait seed** — Madin *et al.* (Sci Data 2020) synthesis: growth rate (from doubling time),
   temperature, pH, oxygen relationship, carbon substrates. *(done)*
3. **BacDive** (DSMZ) — culture temperature/pH/oxygen and media for ~90k strains (needs API key).
4. **Literature mining** — measured & derivable growth / uptake / secretion rates + exact medium
   + conditions from the primary literature, curated (rates computed from raw data where needed),
   each snippet-verified and cited; media resolved to the `Media` repo.
