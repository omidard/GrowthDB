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

## What's here

```
GrowthDB/
├── data/
│   ├── organisms.json        # every prokaryote (backbone), with a has_growth flag
│   ├── growth_records.json   # growth rate / uptake / secretion records + conditions + medium link
│   └── index.json            # summary counts
├── DESIGN.md                 # schema, curation & referencing rules
├── README.md
└── LICENSE
```

## Current contents

- **Backbone: 124,666 prokaryotes** — every **GTDB r220** species (107,235 bacteria + 10,122
  archaea) plus 7,309 named species seen in growth data. Listed whether or not growth data
  exists yet.
- **50,524 growth records** across **14,613 organisms** from the **Madin *et al.* (Sci Data 2020)**
  prokaryote-traits synthesis — culture temperature, optimum temperature/pH, oxygen relationship,
  carbon substrates, and **1,134 measured growth rates** (µ = ln2 / doubling-time), each cited.

This is the seed. The bulk of quantitative **growth / uptake / secretion rates** under defined
media come from ongoing curation of the primary literature (see [DESIGN.md](DESIGN.md) roadmap),
plus DSMZ **BacDive** conditions.

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
