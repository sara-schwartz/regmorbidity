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

## What you supply

1. **LPR** - **one** combined diagnosis table with person id, ICD code,
   contact date (defaults: `pnr`, `C_DIAG`, `D_INDDTO`). The package does
   not merge LPR2 / LPR3 / psychiatric tables for you.
2. **LMDB** - person id, **full** ATC, dispensing date (defaults: `pnr`,
   `atc`, `eksd`). Medication extractors **require** `from` (typically
   `from = as.Date("1997-01-01")` unless you have handled the pre-1997
   mother-CPR issue another way).
3. Code lists - `mb_codelist()` or your own CSVs.

## Real register data: parquet + DuckDB

For real register work, use **parquet** and open it through raw DuckDB + DBI +
`dbplyr`; do not load a SAS register into R's memory. In particular,
`mb_extract_medication_batch()` and the medication DuckDB path need a `dbplyr`
`tbl_lazy` made from a plain `DBI::dbConnect(duckdb::duckdb())` connection and
`dplyr::tbl(con, ...)`:

```r
library(DBI)
library(duckdb)
library(dplyr)

con <- dbConnect(duckdb())
# path = folder of parquet files for LMDB (or LPR)
dbExecute(con, "CREATE VIEW lmdb AS SELECT * FROM read_parquet('path/to/lmdb/**/*.parquet')")
lmdb <- tbl(con, "lmdb")

# Then use the lazy table, for example:
# mb_extract_medication_batch(lmdb, codes = rx_codes, ...)
# dbDisconnect(con, shutdown = TRUE) when done
```

For LPR, create an `lpr` view from the LPR parquet folder and use
`tbl(con, "lpr")` in the same way. Diagnosis has no `window_order()` rule, but
this opening pattern is still recommended for consistency and RAM use. For a
real run, install the packages used here (`duckdb`, `DBI`, `dbplyr`, and
`arrow`). Small in-memory data frames still work for tiny/toy runs, as in the
vignette; use DuckDB + parquet for a real LMDB.

Do **not** feed medication extractors with
`duckplyr::read_parquet_duckdb()` or `fastreg::read_register()`; those are not
the supported path for medication extracts. `fastreg` is useful for converting
SAS to parquet and setting up a project layout: [fastreg](https://dp-next.github.io/fastreg/),
specifically [Getting started](https://dp-next.github.io/fastreg/articles/fastreg.html).

## Quick start

Extractors take the **full code list** and write one onset file per
condition in one pass. Prefer that. Use a sequential, one-condition-at-a-
time run only when the cohort or year range is so large that a full-list
pass is awkward (see vignette).

```r
library(regmorbidity)

# For real parquet-backed registers, define lpr/lmdb as shown above.
# The vignette uses tiny in-memory toy tables instead.

# Bundled lists are ATC + first-pass ICD (2026-09-28). distress archived out of default set.
rx_codes <- mb_codelist()
dx_codes <- mb_codelist(conditions = "hypertension", vocab = "ICD10")

# Diagnoses (one combined LPR table; ICD code list, not the ATC bundle)
mb_extract_diagnosis(
  lpr,
  codes  = dx_codes,
  outdir = "data/dx",
  from   = as.Date("1995-01-01")   # example window; set to your study
)

# Medications (batch: all ATC conditions in two queries)
mb_extract_medication_batch(
  lmdb,
  codes  = rx_codes,
  outdir = "data/rx",
  from   = as.Date("1997-01-01")
)

# Merge returns a long data frame (no outdir). Load is for extract outdirs.
long <- mb_merge_all("data/dx", "data/rx")
# long <- mb_load_conditions("data/rx")   # one extract dir, no merge
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

`?regmorbidity` - all 16 exported functions.

`ASSUMPTIONS_AND_LIMITATIONS.txt` and `DECISIONS.md` in the source repo -
read before reporting a number. Output is onset-only by default.
