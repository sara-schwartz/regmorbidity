# regmorbidity

Condition indicators, with dates, from Danish register data — defined in
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

1. **LPR** — **one** combined diagnosis table with person id, ICD code,
   contact date (defaults: `pnr`, `C_DIAG`, `D_INDDTO`). The package does
   not merge LPR2 / LPR3 / psychiatric tables for you.
2. **LMDB** — person id, **full** ATC, dispensing date (defaults: `pnr`,
   `atc`, `eksd`). Medication extractors **require** `from` (recommend
   `from = as.Date("1997-01-01")` unless you have handled the pre-1997
   mother-CPR issue another way).
3. Code lists — `mb_codelist()` (bundled: **ATC + first-pass ICD**; distress
   archived out of the default set) or your own CSVs.

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

Fifteen exports. Prefer this path; see **Functions** below for the rest.

```r
library(regmorbidity)

# For real parquet-backed registers, define lpr/lmdb as shown above.
# The vignette uses tiny in-memory toy tables instead.

rx_codes <- mb_codelist(vocab = "ATC")
dx_codes <- mb_codelist(conditions = "hypertension", vocab = "ICD10")

mb_extract_diagnosis(
  lpr,
  codes  = dx_codes,
  outdir = "data/dx",
  from   = as.Date("1995-01-01")   # example window; set to your study
)

mb_extract_medication_batch(
  lmdb,
  codes  = rx_codes,
  outdir = "data/rx",
  from   = as.Date("1997-01-01")
)

long <- mb_merge_all("data/dx", "data/rx")
wide <- mb_to_wide(long)
mb_count_conditions(wide, as_of = "2015-01-01")
```

Sequential extract, `keep_events`, and prevalence:
`vignette("regmorbidity")`.

## Functions

**Happy path (loud)**

1. `mb_codelist` — load/filter lists (`conditions=`, `vocab=`)
2. `mb_extract_diagnosis` — LPR onset
3. `mb_extract_medication_batch` — preferred medication onset (SQL, low RAM)
4. `mb_merge_all` — combine dx+rx extract directories → long
5. `mb_to_wide` → `mb_count_conditions` — ever-after counts as of a date

**Backup / special**

- `mb_extract_medication` — sequential; use when you need `keep_events=TRUE` or one-condition debug
- `mb_merge_conditions` — merge one condition’s two data frames (in memory)
- `mb_load_conditions` — load one extract outdir to long (skips `*_all_events.rds`); medication-only studies
- `mb_prevalence` — lookback prevalence on **raw events**, not onset `.rds`

**QA (optional, before long DST runs)**

- `mb_check_codes` — does each list code match anything in the register?
- `mb_lookup` — which conditions claim this code?
- `mb_overlap` — codes shared between conditions

**Advanced / demoted**

- `mb_inspect_codes` — reports **code column string lengths** in the register (e.g. is `atc2` 3 chars?). Only needed before sequential extract / prefilter debugging. Not list review.

**Provisional**

- `mb_apply_exclusions` — optional stage-2; incomplete vs Prior; not happy path. HTN still has C03 / HF / CKD gaps.

## More detail

`vignette("regmorbidity")` — sequential extract, `keep_events`, prevalence, QA.

`?regmorbidity` — all 15 exported functions.

Design notes (`ASSUMPTIONS_AND_LIMITATIONS.txt`, `DECISIONS.md`) live in the
source working tree if present; otherwise use the vignette and `?regmorbidity`.
Output is onset-only by default.
