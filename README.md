# RDBES ecoregion landings shares

This TAF repository estimates, for each ICES stock and catch year, the share of uniquely assigned RDBES landings reported from statistical rectangles associated with each ICES ecoregion.

The workflow uses RDBES commercial landings for **2021–2025** and interprets them against a fixed **2026 ICES stock catalogue**. It first matches landings to stocks, then calculates ecoregion percentages from the matched records.

> **In one sentence:** first decide which ICES stock each RDBES landings record can belong to; then calculate what share of the clearly assigned landings was reported from rectangles linked to each ecoregion.

## What the result means

For each stock, year and ecoregion:

```text
ShareLandingsInEcoregion (%) = 100 × Inside_kg / Assigned_kg
```

- `Assigned_kg` is all RDBES weight uniquely assigned to that stock and year.
- `Inside_kg` is the part of that assigned weight whose statistical rectangle belongs to the selected ecoregion rectangle list.
- A missing share means **unavailable**, not zero real-world landings.
- The estimate refers to the reported fishing location, not the port of landing.
- Shares are based only on the RDBES subset that received a unique stock assignment. They are not verified percentages of all landings from the biological stock.

A statistical rectangle can occur in more than one ecoregion list. Shares calculated independently across ecoregions are therefore not expected to add to exactly 100%.

## Repository structure

The main folders and files are:

```text
boot/
  initial/data/                         Input data supplied with the repository
    RDBES Landings CL 2021 2025/       RDBES CommercialLanding input
    ICES_StatRec_mapto_ICES_Areas/     Official 2017 rectangle-to-area map
    ICES_nephrops_FU/                   Nephrops Functional Unit shapefile

data/                                  Prepared data and cached ICES references
  extendedRDBES.RData                   Prepared RDBES landings object
  rdbes_landings_compact.rds            Compact landings table
  ices_cache_simple/                    Cached SID, Vocabulary and GIS results

model/                                 Saved stock matching and diagnostics
  rdbes_stock_matches.RData             Global stock-matching checkpoint
  stock_matching_rectangle_first/       Matching summaries and audit CSV files

output/                                Ecoregion-specific results

data.R                                Prepares the RDBES input and species lookup
model_1_global_stock_matching.R        Global stock matching
model_2_*.R                            Ecoregion-share calculation
model.R / output.R / report.R          Optional TAF orchestration/reporting scripts
```

The repository may contain older or comparison scripts retained for validation. The current method is the hybrid, rectangle-first workflow described below.

## Workflow

The method has two separate stages.

### Stage 1: match RDBES landings to ICES stocks

1. **Prepare the landings**  
   Read RDBES CommercialLanding data, keep landings records, standardise year, area, rectangle and species fields, and aggregate identical records while preserving total scientific weight.

2. **Identify the species**  
   Use the WoRMS/Aphia or FAO-ASFIS species code supplied in RDBES. A missing code is filled only when the ASFIS–WoRMS relationship is unambiguous. Conflicting taxonomic identifiers remain unresolved.

3. **List possible stocks**  
   Use the fixed 2026 SID catalogue and ICES Vocabulary relationships to find stocks for the same species and their ICES-area definitions. Same-species alternatives remain candidates until geography rules them out.

4. **Resolve location for ordinary stocks**  
   Translate each statistical rectangle to ICES area(s) using the official 2017 rectangle-to-area map. Use the reported RDBES ICES area (`CLarea`) to narrow boundary rectangles or as a fallback when the rectangle is unavailable.

5. **Handle rectangle/CLarea disagreements conservatively**  
   The disagreement remains visible for auditing. Each source is checked separately against each candidate stock:

   - if both sources lead to the same stock-level conclusion, retain that conclusion;
   - if the two sources disagree about whether the record belongs to the stock, keep the evidence `unknown`.

   This avoids choosing one source as universally correct. In the audit, 232 of 233 affected ordinary stock-location relationships were outside under both sources. One 1 kg `lin.27.5b` record at rectangle `49E4` remains genuinely source-dependent.

6. **Handle Nephrops separately**  
   Nephrops Functional Units cannot be represented reliably by ordinary ICES-area codes. For Nephrops, the method intersects the reported ICES area and statistical rectangle to create the best available geographic support, then compares that support with Functional Unit polygons.

7. **Make a unique assignment**  
   A record is assigned only when exactly one candidate stock remains and its geographic evidence is `full`. If several stocks fit, or the location evidence is `partial`, `unknown` or held for `review`, the record remains unresolved.

8. **Apply review and special rules**  
   The workflow preserves `raw_evidence` before applying the review policy. This makes it possible to distinguish a real geographic limitation from a deliberate decision to hold a stock for review. Documented sardine, anchovy and sandeel rules are applied separately.

9. **Save the matching checkpoint**  
   The stock assignments, unresolved candidates, fixed stock catalogue and matching metadata are saved in `model/rdbes_stock_matches.RData`.

   A **checkpoint** is simply a saved intermediate result. It allows the expensive stock-matching stage to be reused for all ecoregions instead of being repeated eleven times.

### Stage 2: calculate ecoregion shares

For each ecoregion, the second model script:

1. loads the completed stock-matching checkpoint;
2. reads or reuses the ecoregion statistical-rectangle list;
3. marks uniquely assigned landings as inside or outside the ecoregion;
4. calculates `Inside_kg / Assigned_kg` for each stock and year;
5. applies the clearly labelled Division 9.a anchovy proxy where required;
6. writes ecoregion-specific summaries and audit files.

## Supported ecoregions

The current global matching scope covers:

- Baltic Sea
- Bay of Biscay and the Iberian Coast
- Celtic Seas
- Greater North Sea
- Norwegian Sea
- Icelandic Waters
- Barents Sea
- Greenland Sea
- Faroes
- Oceanic Northeast Atlantic
- Azores

## Data sources

| Source | Role in the workflow |
|---|---|
| RDBES CommercialLanding, 2021–2025 | Year, scientific weight, species identifiers, reported ICES area and statistical rectangle. |
| ASFIS species table | Current ASFIS species-code reference. |
| ICES Vocabulary: `SpecASFIS` ↔ `SpecWoRMS` | Connects FAO-ASFIS species codes to WoRMS/Aphia identifiers. |
| ICES Stock Information Database (SID), 2026 | Fixed stock codes, descriptions and ecoregion labels used for all five catch years. |
| ICES Vocabulary: stock ↔ species | Links ICES stock codes to species. |
| ICES Vocabulary: stock ↔ ICES area | Defines the ordinary ICES-area coverage of each stock. |
| Official 2017 ICES statistical rectangle → area map | Converts rectangles to ICES areas for ordinary-stock matching. |
| ICES GIS layers | Area and rectangle geometry used by the Nephrops branch. |
| Nephrops Functional Unit shapefile | Nephrops stock boundaries. |
| `icesVMS` ecoregion rectangle lists | Defines which rectangles count as inside each ecoregion for the final percentage. |

### Important coverage limitation

The method can use only records present in the RDBES extract. The extract used here contains no submissions from **Norway, the Faroe Islands or Iceland**. Stocks strongly associated with those countries may therefore have missing or incomplete landings even when the matching method is working correctly. Completeness for other countries, fleets and reporting strata has not been assessed here.

## Evidence labels

| Label | Meaning |
|---|---|
| `full` | The available location is entirely compatible with the stock. This can support a unique assignment. |
| `outside` | The location does not fit the stock. The stock is removed as a candidate. |
| `partial` | The available location partly overlaps the stock or is broader than the stock boundary. This is not enough for a unique assignment. |
| `unknown` | There is not enough usable information, or the geographic sources disagree about the stock. |
| `review` | The geographic result was not outside, but the stock is deliberately held unresolved pending a stock-specific rule or decision. |

## Special cases

### Nephrops

Nephrops retains the geometry-based Functional Unit method. The accepted FU25 five-rectangle correction uses `15E0`, `15E1`, `16E0`, `16E1` and `17E1`.

`nep.fu.30` remains unresolved because the relevant RDBES locations only partly overlap FU30. The remaining partial exposure is about 549,280 kg across six locations. About 90% is in rectangle `02E2`, where only 37.7% of the available spatial support lies inside FU30. The method does not assign the whole record to FU30 and does not split kilograms in proportion to area, because both choices would add unsupported assumptions about where catches occurred within the reported location.

### Sardine

PIL records in rectangles `25E4` and `25E5` that would otherwise be assigned to `pil.27.7` are reassigned to `pil.27.8abd` under the documented production rule. The compact matching table does not retain country, so the country condition behind the original operational rule cannot be verified in this workflow.

### Anchovy in Division 9.a

RDBES geography cannot uniquely split `ane.27.9aS` from `ane.27.9aW`. The workflow therefore does not duplicate or force kilograms into either stock. For the Bay of Biscay and Iberian Coast, both stocks receive the same pooled annual Division 9.a share. This is a regional percentage proxy, not a unique assignment of kilograms.

### Sandeel

The seven `san.sa.*` assessment stocks are not ordinary ICES-area stocks in this workflow. They remain under a separate conservative review rule until a maintainable assessment-area mapping is available.

### Salmon

Salmon is outside the current fisheriesXplorer advice scope. The workflow does not attempt to repair the following salmon assignments: `sal.nac.all`, `sal.neac.all`, `sal.wgc.all`, `sal.27.22-31`, and `sal.27.32`.

## Review decisions

The original review list was intentionally conservative, but it also blocked some full, uncontested records. The following stocks were released after stock-specific diagnostics showed that the normal matching logic could recover a clear subset while leaving contested records unresolved:

<details>
<summary>Released stocks</summary>

- `cod.27.46a7d20`
- `cod.27.22-24`
- `cod.27.24-32`
- `pra.27.3a4a`
- `her.27.3a47d`
- `her.27.20-24`
- `cyo.27.nea`
- `her.27.irls`
- `her.27.6aS7bc`

</details>

Stocks are not released merely to remove a zero. A release requires full, uncontested evidence that can be recovered without forcing the contested part of the fishery into the stock.

## Why some stocks still have no assigned landings

A stock with no assigned landings is not automatically a code failure. The audit distinguishes several causes.

### No RDBES records for the stock taxon

The current matching scope contains no RDBES records for:

- `bwp.27.2729-32`
- `bzq.27.2425`
- `bzq.27.2628`
- `seh.27.1`
- `seh.27.125a14`

This may reflect RDBES coverage, including the absence of Norway, Faroe Islands and Iceland submissions.

### The species is present, but observed geography is outside the stock

<details>
<summary>Stocks in this category</summary>

- `aru.27.5a14`
- `cod.21.1.isc`
- `cod.21.1.osc`
- `cod.21.27.1.14`
- `ghl.27.1-2`
- `jaa.27.10a2`
- `raj.27.1012`
- `reb.2127.dp`
- `reb.2127.sp`
- `reb.27.1-2`
- `reb.27.14b`
- `reb.27.5a14`
- `reg.27.1-2`
- `rng.27.5a10b12ac14b`
- `sbr.27.10`
- `thr.27.nea`

A second group was previously labelled unknown because rectangle and `CLarea` conflicted. Stock-level comparison showed that both sources nevertheless excluded the target stock in 232 of 233 tested relationships. After the refined rule, these relationships are treated as outside. The only genuine source-dependent ordinary-stock case is the 1 kg `lin.27.5b` record at `49E4`.

</details>

### Shared or incomplete geographic evidence

Some records fit more than one same-species stock, or only non-full evidence exists. Examples include:

- `cod.27.21`
- `cod.27.1-2`
- `cod.27.1-2.coastN`
- `cod.27.2.coastS`
- `her.27.nirs`
- `her.27.6aN`
- `her.27.1-24a514a`
- `her.27.5a`
- `pra.27.4a`
- `guq.27.nea`

The workflow keeps these unresolved rather than copying the same weight into several stocks.

### Capelin

RDBES contains a small amount of capelin data, but all observed locations are outside the definitions of `cap.27.1-2` and `cap.27.2a514`. Removing the review flag would therefore not create an assignment.

### Other potential causes

- the fixed 2026 stock catalogue may not reflect historical stock codes or boundaries for all catch years;
- reported rectangle or `CLarea` may be missing, inconsistent or too broad;
- stock origin or management-area rules may not exist in the current ICES Vocabulary relationships;
- country and fleet coverage may be incomplete.

## Running the workflow

The repository uses a TAF-style directory structure, but the two model stages are kept separate so the global matching can be reused.

### 1. Install/load required R packages

The scripts use packages including:

```r
icesTAF
dplyr
tidyr
stringr
data.table
sf
httr
jsonlite
icesSD
icesVocab
icesVMS
```

Use the repository lockfile or project-specific dependency management if available.

### 2. Prepare the data

From the repository root:

```r
source("data.R")
```

This creates or refreshes `data/extendedRDBES.RData` and related species-reference diagnostics.

### 3. Run global stock matching

```r
source("model_1_global_stock_matching.R")
```

This creates:

```text
model/rdbes_stock_matches.RData
```

and diagnostic CSV files under:

```text
model/stock_matching_rectangle_first/
```

Rebuild this checkpoint when the RDBES input, species references, SID reference year, stock-area rules, review decisions or matching logic changes.

### 4. Run the ecoregion-share workflow

Run the repository's `model_2` ecoregion-share script after the checkpoint exists. The script reads the checkpoint, processes the supported ecoregions, and writes regional outputs under `output/`.

Do not rerun global matching separately for each ecoregion unless the checkpoint scope or reference data changed.

## Main outputs

| Output/field | Meaning |
|---|---|
| `model/rdbes_stock_matches.RData` | Saved global stock assignments and candidate information. |
| `stock_assignment_summary.csv` | Assigned kilograms by stock and year. |
| `matching_summary.csv` | Weight and row counts by matching status. |
| `unresolved_landings.csv` | Records that did not receive a unique stock assignment. |
| `Assigned_kg` | Weight uniquely assigned to a stock-year. |
| `Inside_kg` | Assigned weight whose rectangle is in the ecoregion list. |
| `Outside_kg` | Assigned weight outside the ecoregion list. |
| `Missing_rectangle_kg` | Assigned weight with a missing rectangle. |
| `ShareLandingsInEcoregion` | `100 × Inside_kg / Assigned_kg`. |
| `Unresolved_candidate_kg` | Overlapping diagnostic exposure for unresolved records. **Do not sum this across stocks.** |
| `raw_evidence` | Geographic evidence before the review override. |
| `evidence` | Final evidence after the review policy. |
| `share_method`, `share_source`, `special_rule` | Identifies proxy and special-rule results. |

## Validation

Before adoption, the hybrid method was compared with the older polygon-heavy workflow:

- Nephrops matched the baseline exactly after the validated FU treatment was restored: **137/137 stock-years**.
- For non-Nephrops, **1,649** stock-year-ecoregion shares were comparable.
- **1,556** were exactly equal.
- **1,640 of 1,649 (99.45%)** were within 0.01 percentage points.
- Only three comparisons were above 0.1 percentage points.
- None was above 1 percentage point.
- The maximum difference was approximately 0.566 percentage points.

The later rectangle/`CLarea` conflict refinement is deliberately narrow. A refreshed full comparison should be recorded after the final production rerun.

## Reproducibility and maintenance

- Keep the exact RDBES input and static 2017 rectangle map under TAF/bootstrap control.
- Keep cached ICES reference files under `data/ices_cache_simple/` unless deliberately refreshing references.
- Record the fixed SID reference year in outputs and run notes.
- Rebuild the matching checkpoint after any change to inputs, taxonomy, stock definitions, review decisions or special rules.
- Preserve unresolved records and diagnostics; do not convert unavailable estimates to zero.
- Do not sum candidate exposure across stocks because one unresolved record can appear against several candidates.

## Documentation

The fuller plain-language method and decision record is maintained in:

```text
RDBES_ecoregion_landings_shares_method_annex_v3.docx
```

That annex contains the detailed unresolved-stock categories, FU30 investigation, review decisions and interpretation notes.

