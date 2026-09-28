# regmorbidity

Condition indicators, with dates, from Danish register data - defined in
**editable CSV code lists** rather than in code. It computes no score and is
not an implementation of the Danish Multimorbidity Index or of any other
published index. The bundled lists take their starting point in Prior et al.
(2016) and are a starting point to revise, not an instrument to cite.

Authors: Jie Zhang, Sara Schwartz (saras@clin.au.dk)

## Install

```r
remotes::install_github("sara-schwartz/regmorbidity")
library(regmorbidity)
```

Local clone (run from the parent of `regmorbidity/`):
`install.packages("regmorbidity", repos = NULL, type = "source")`.

## What you supply

1. **LPR** - **one** combined diagnosis table with person id, ICD code,
   contact date (defaults: `pnr`, `C_DIAG`, `D_INDDTO`). The package does
   not merge LPR2 / LPR3 / psychiatric tables for you.
2. **LMDB** - person id, **full** ATC, dispensing date (defaults: `pnr`,
   `atc`, `eksd`). For LMDB, pass `from = as.Date("1997-01-01")` unless you
   have handled the pre-1997 mother-CPR issue another way.
3. Code lists - `mb_codelist()` or your own CSVs.

## Quick start

Extractors take the **full code list** and write one onset file per
condition in one pass. Prefer that. Use a sequential, one-condition-at-a-
time run only when the cohort or year range is so large that a full-list
pass is awkward (see vignette).

```r
library(regmorbidity)

codes <- mb_codelist()

# Diagnoses (one combined LPR table)
mb_extract_diagnosis(
  lpr,
  codes  = codes,
  outdir = "data/dx",
  from   = as.Date("1995-01-01")   # example window; set to your study
)

# Medications (batch: all ATC conditions in two queries)
mb_extract_medication_batch(
  lmdb,
  codes  = codes,
  outdir = "data/rx",
  from   = as.Date("1997-01-01")
)

# Combine dx + rx per condition, then reshape / count
mb_merge_all("data/dx", "data/rx", outdir = "data/conditions")
long <- mb_load_conditions("data/conditions")
wide <- mb_to_wide(long)
mb_count_conditions(wide, as_of = "2015-01-01")
```

One condition at a time (sequential medication extract): pass a code list
filtered to that condition into `mb_extract_medication()`, or the same
idea with `mb_extract_diagnosis()`. Details and when to prefer it:
`vignette("regmorbidity")`.

## More detail

`vignette("regmorbidity")` - full walkthrough (sequential extract,
`keep_events`, prevalence, checks).

`?regmorbidity` - all 13 exported functions.

`ASSUMPTIONS_AND_LIMITATIONS.txt` and `DECISIONS.md` in the source repo -
read before reporting a number. Output is onset-only by default.
