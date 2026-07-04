# GrowthDB

**Experimental prokaryote growth data for constraint-based metabolic modelling** — growth
rates, substrate uptake rates and product secretion rates, with the exact **medium** (linked to
the [`Media`](https://github.com/omidard/Media) repo) and **culture conditions** (batch vs
chemostat, aeration, temperature, pH), every value **referenced** to its primary source.

Genome-scale metabolic models are only as good as the data that constrains and validates them.
Measured rates under a defined medium and condition are exactly what you need to bound an
exchange reaction or check a predicted growth rate — but that data is scattered across thousands
of papers in inconsistent units. GrowthDB collects it once, machine-readable, with the medium
resolved to BiGG exchanges so a measured uptake rate drops straight onto `EX_glc__D_e`.

## Explore it online

Interactive browser → **https://omidard.github.io/GrowthDB/** — filter the quantitative growth
records by domain, culture mode, oxygen and flux; open any record for its conditions, uptake /
secretion rate tables, the medium (linked to the Media library), and the citation + snippet.

## What's here

```
GrowthDB/
├── data/
│   ├── organisms.json        # every prokaryote (backbone), with a has_growth flag
│   ├── growth_records.json   # growth rate / uptake / secretion records + conditions + medium link
│   ├── records_index.json    # compact index of the quantitative records (browser)
│   └── index.json            # summary counts
├── index.html                # interactive browser (GitHub Pages)
├── DESIGN.md                 # schema, curation & referencing rules
├── README.md
└── LICENSE
```

## Current contents

- **Backbone: 124,666 prokaryotes** — every **GTDB r220** species (107,235 bacteria + 10,122
  archaea) plus 7,309 named species seen in growth data. Listed whether or not growth data exists.
- **51,881 growth records** in total:
  - **50,524** from the **Madin *et al.* (Sci Data 2020)** prokaryote-traits synthesis — culture
    temperature, optimum temperature/pH, oxygen relationship, carbon substrates.
  - **1,357** curated from **1,201 primary papers** (multi-agent mining) — growth rate, specific
    **uptake** and **secretion** rates, and the exact medium + culture mode (batch/chemostat),
    aeration, T and pH; each snippet-verified and cited.
- **2,367 measured growth rates**, **294** records with uptake rates, **194** with secretion rates.
- **664 records linked to a medium** in the [`Media`](https://github.com/omidard/Media) library
  (105 new media were added to `Media` from these papers, so it stays exhaustive).

Still to come: DSMZ **BacDive** conditions (needs an API key) and continued literature curation.

## Every rate is referenced and method-tagged

Each rate records whether it was **`reported`** (stated in the paper) or **`calculated`** (derived
by curation from the paper's data — e.g. µ from an OD/CFU time-course; specific uptake from
yield × µ; secretion from product yield), with the derivation in `curation_notes`. Media are
pointers into the `Media` repo; a medium not yet there is **added to `Media`** first, then linked.

## Related

Part of a set of GEM data resources: [`Media`](https://github.com/omidard/Media) (growth &
simulation media, BiGG-mapped), [`panGEMs`](https://github.com/omidard/panGEMs),
[`EcopanGEM`](https://github.com/omidard/EcopanGEM).

## License

Data: **CC-BY-4.0** (cite each record's original source, in its `provenance`). Code: **MIT**.
Backbone from [GTDB](https://gtdb.ecogenomic.org); trait seed from Madin *et al.* Sci Data 2020.
