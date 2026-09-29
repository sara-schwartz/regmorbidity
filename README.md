# regmorbidity

Condition indicators, with dates, from Danish register data - defined in
**editable CSV code lists** rather than in code.

**What it is:** load your lists, extract onset dates from LPR (diagnoses) and
LMDB (medications), merge the halves, and summarise ever-after counts (or
lookback prevalence on raw events).

**What it is not:** it computes no score and is not an implementation of the
Danish Multimorbidity Index or of any other published index. The bundled lists
take their starting point in Prior et al. (2016) and are a starting point to
revise, not an instrument to cite.

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
   `atc`, `eksd`). Medication extractors **require** `from` (recommend
   `from = as.Date("1997-01-01")` unless you have handled the pre-1997
   mother-CPR issue another way).
3. **Code lists** - `mb_codelist()` (bundled: **ATC + first-pass ICD**;
   distress archived out of the default set) or your own CSVs.

## On DST: parquet + DuckDB

On Statistics Denmark (DST), registers should be **parquet**, opened through
raw DuckDB + DBI + `dbplyr` - not loaded as SAS into R's memory. Medication
batch needs a `dbplyr` `tbl_lazy` from a plain
`DBI::dbConnect(duckdb::duckdb())` connection.

Replace the folder path below with **your** LMDB parquet folder on the project
drive (ask a colleague if unsure). The `**` is a DuckDB glob, not something you
fill in by hand: it means “this folder and every subfolder”. `*.parquet` means
“every file ending in `.parquet`”. Together, `'…/lmdb/**/*.parquet'` reads all
parquet files under that LMDB folder (typical when years are split into
subfolders). If everything sits in one flat folder with no subfolders, you can
use `'…/lmdb/*.parquet'` instead.

```r
library(DBI)
library(duckdb)
library(dplyr)

# 1. Open an empty DuckDB database in memory
con <- dbConnect(duckdb())

# 2. Point DuckDB at your LMDB parquet files (edit THIS path only)
dbExecute(
  con,
  "CREATE VIEW lmdb AS
   SELECT * FROM read_parquet('E:/workdata/YOUR_PROJECT/cleaned-data/parquet-registers/lmdb/**/*.parquet')"
)

# 3. Get a lazy dplyr table - nothing is loaded into R yet
lmdb <- tbl(con, "lmdb")

# Same three steps for LPR, with a different path and view name, e.g.:
# dbExecute(con, "CREATE VIEW lpr AS SELECT * FROM read_parquet('…/lpr/**/*.parquet')")
# lpr <- tbl(con, "lpr")

# When finished:
# dbDisconnect(con, shutdown = TRUE)
```

Do **not** feed medication extractors with `duckplyr::read_parquet_duckdb()` or
`fastreg::read_register()`. `fastreg` is still useful for SAS→parquet setup on
DST: [fastreg](https://dp-next.github.io/fastreg/). Tiny in-memory frames work
for toy runs (see the vignette).

## Recommended workflow

The usual study run is the five steps below (load lists → extract diagnoses →
extract medications with batch → merge → count). How many *conditions* you get
depends on your CSV code lists - not on the package API. The package itself
exposes a fixed set of **functions** (listed in the next section).

```r
library(regmorbidity)

# ---------------------------------------------------------------------------
# 0. Registers (you supply these)
#    On DST: open parquet via DuckDB as in the section above → objects `lpr`
#    and `lmdb`. Below we assume those already exist.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 1. Load code lists
#    The package already ships CSV code lists (one file per condition under
#    inst/extdata/codelists/). You do NOT pass file paths here unless you have
#    your own edited copies.
#
#    mb_codelist() with no arguments = load ALL bundled rows (ATC + ICD).
#    vocab = "ATC" / "ICD10" = keep only that vocabulary after load (a filter).
#    conditions = "..." = keep only named condition(s) (also a filter).
#
#    Medication extractors need ATC rows; diagnosis extractors need ICD rows.
#    So we load twice with different filters - still the same bundled CSVs.
# ---------------------------------------------------------------------------
rx_codes <- mb_codelist(vocab = "ATC")
# Example: only hypertension ICD rows for a small diagnosis demo.
# For a full study, drop `conditions=` to keep every bundled ICD condition:
#   dx_codes <- mb_codelist(vocab = "ICD10")
dx_codes <- mb_codelist(conditions = "hypertension", vocab = "ICD10")

# Own lists instead of bundled? Point at a folder of CSVs or one big CSV:
#   codes <- mb_codelist("path/to/my_codelists/")

# ---------------------------------------------------------------------------
# 2. Extract diagnoses from LPR
#    Writes one onset .rds per condition into outdir (here: data/dx/).
# ---------------------------------------------------------------------------
mb_extract_diagnosis(
  lpr,
  codes  = dx_codes,
  outdir = "data/dx",
  from   = as.Date("1995-01-01"),  # inclusive start; set to your study
  to     = as.Date("2018-12-31")   # optional inclusive end; not required - omit to leave open-ended
)

# ---------------------------------------------------------------------------
# 3. Extract medications from LMDB - prefer batch
#    Same idea: one onset .rds per ATC condition into data/rx/.
#    `from` is required on medication extractors (1997 recommended on DST).
#    `to` is optional (not required) - omit for no upper date bound.
# ---------------------------------------------------------------------------
mb_extract_medication_batch(
  lmdb,
  codes  = rx_codes,
  outdir = "data/rx",
  from   = as.Date("1997-01-01"),  # required on LMDB
  to     = as.Date("2018-12-31")   # optional; not required
)

# ---------------------------------------------------------------------------
# 4. Merge diagnosis + medication halves
#    Reads the two outdirs and returns one long data frame in memory
#    (person × condition × onset). No new files.
# ---------------------------------------------------------------------------
long <- mb_merge_all("data/dx", "data/rx")

# ---------------------------------------------------------------------------
# 5. Reshape and count
#    Wide = one column per condition (onset dates).
#    Count = how many conditions each person has as of a date (ever after onset).
# ---------------------------------------------------------------------------
wide <- mb_to_wide(long)
mb_count_conditions(wide, as_of = "2015-01-01")
```

For sequential extract, `keep_events`, prevalence lookback, exclusion timing,
and a worked miniature: `vignette("regmorbidity")`.

## All functions

Listed in the order most people meet them. One sentence each - **when to use**.

### Load list

- **`mb_codelist`** - load the bundled CSVs (or your own folder/CSV/data frame),
  then optionally filter with `conditions=` / `vocab=`; start here before any
  extract or QA. Does not invent lists - it reads them.

### Preflight QA

- **`mb_check_codes`** - for each code in your list, count how many register
  rows match (prefix); optional preflight before a long DST extract to catch
  typos / dead codes that would give silent zeros.
- **`mb_lookup`** - which conditions in the list claim this code?
- **`mb_overlap`** - which codes are shared between conditions?
- **`mb_inspect_code_lengths`** *(advanced)* - reports string lengths of code
  columns in the *register* sample (e.g. whether `atc2` is 3 characters);
  use before sequential medication extract / prefilter debugging - **not** for
  reviewing the code list.

### Extract

- **`mb_extract_diagnosis`** - LPR onset: one matching diagnosis is the
  condition; write one onset `.rds` per condition.
- **`mb_extract_medication_batch`** - **prefer** this for medication onset
  (one SQL pass, low RAM); does not keep raw events.
- **`mb_extract_medication`** - sequential backup; use when you need
  `keep_events = TRUE` (prevalence / debug) or one-condition disk checkpoints.

### Combine

- **`mb_merge_all`** - combine diagnosis + medication extract directories into
  one long table (the usual dx+rx merge).
- **`mb_merge_conditions`** - merge one condition’s two in-memory data frames
  (OR/AND / optional `logic` column).
- **`mb_load_conditions`** - load one extract outdir to long (skips
  `*_all_events.rds`); use for medication-only studies.

### Analyse

- **`mb_to_wide`** - long onset table → one row per person, one column per
  condition (onset dates).
- **`mb_count_conditions`** - ever-after condition counts as of a date (once
  met, always met).
- **`mb_prevalence`** - lookback prevalence on **raw events**, not onset
  `.rds`; needs events from sequential extract with `keep_events = TRUE`.

### Provisional

- **`mb_apply_exclusions`** - optional stage-2 Prior-style exclusions;
  incomplete vs Prior; skip until you need stage-2 exclusions (HTN still has C03 / HF / CKD gaps).

## More detail

`vignette("regmorbidity")` - worked miniature, batch vs sequential +
`keep_events`, `mb_check_codes` example, prevalence, merge logic, exclusion
timing, and why `mb_inspect_code_lengths` exists.

`?regmorbidity` - package help for every exported function.

Design notes (`ASSUMPTIONS_AND_LIMITATIONS.txt`, `DECISIONS.md`) live in the
source working tree if present; otherwise use the vignette and `?regmorbidity`.
Output is onset-only by default.
