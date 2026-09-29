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
# Install from GitHub (needs the remotes package once):
#   install.packages("remotes")
remotes::install_github("sara-schwartz/regmorbidity")
library(regmorbidity)
```

On Statistics Denmark (DST), the parquet + DuckDB path also needs **DBI**,
**duckdb**, and **dplyr** (listed under Suggests; install them if missing).

## What you supply

1. **LPR** - **one** combined diagnosis table with person id, ICD code,
   contact date (defaults: `pnr`, `C_DIAG`, `D_INDDTO`). The package does
   not merge LPR2 / LPR3 / psychiatric tables for you.
2. **LMDB** - person id, **full** ATC, dispensing date (defaults: `pnr`,
   `atc`, `eksd`). Medication extractors **require** `from` (a `Date`).
   There is no silent floor: you must pass it yourself.

   **Why we recommend `from = as.Date("1997-01-01")` on LMDB.** Until 1996,
   children's prescriptions were often recorded under the mother's CPR/PNR;
   from 1997 they are recorded under the child's own (Pottegard et al. 2016,
   NPR data-resource profile). Keeping pre-1997 rows can therefore attribute
   some child dispensings to the mother. Pass an earlier `from` only if your
   study needs 1995-96 *and* you have handled the mother-CPR issue another
   way. LMDB itself begins in 1995.

3. **Code lists** - call `mb_codelist()` (bundled), or pass your own folder /
   CSV / data frame.

### Bundled code lists (36 conditions)

The package ships **36 CSV files** under `inst/extdata/codelists/` (one file
per condition). Basenames:

alcohol, allergy, anemia, anorexia_bulimia, atrial_fibrillation, bipolar,
cancer, ckd, connective_tissue, dementia, diabetes, diverticular,
dyslipidemia, epilepsy, gout, hearing, heart_failure, hiv, hypertension,
ibd, ihd, liver, migraine, multiple_sclerosis, osteoporosis, pad, pain,
parkinson, prostate, pulmonary, schizophrenia, stroke, substance_abuse,
thyroid, ulcer, vision.

Typical columns: `condition`, `condition_label`, `category`, `vocab_id`
(`ATC` or `ICD10`), `code`, `exclude`, `min_prescriptions`, `window_days`,
`note`. Content is **first-pass ICD + ATC**; psychological distress is
archived out of the default set. You usually load everything once with
`codes <- mb_codelist()`; diagnosis and medication extractors then keep the
vocabulary they need (ICD vs ATC) internally. Optional filters:
`conditions = c("hypertension", "diabetes")` to keep named conditions only.

## On DST: parquet + DuckDB

On Statistics Denmark (DST), registers should be **parquet**, opened through
raw DuckDB + DBI + `dbplyr` - not loaded as SAS into R's memory. Medication
batch needs a `dbplyr` `tbl_lazy` from a plain
`DBI::dbConnect(duckdb::duckdb())` connection.

Replace the folder path below with **your** LMDB parquet folder on the project
drive (ask a colleague if unsure). The `**` is a DuckDB glob, not something you
fill in by hand: it means "this folder and every subfolder". `*.parquet` means
"every file ending in `.parquet`". Together, `'.../lmdb/**/*.parquet'` reads all
parquet files under that LMDB folder (typical when years are split into
subfolders). If everything sits in one flat folder with no subfolders, you can
use `'.../lmdb/*.parquet'` instead.

```r
library(DBI)
library(duckdb)
library(dplyr)

# 1. Open an empty DuckDB database in memory
con <- dbConnect(duckdb())

# 2. Point DuckDB at your LMDB parquet files (edit THIS path only)
#    **  = this folder and every subfolder
#    *.parquet = every file ending in .parquet
dbExecute(
  con,
  "CREATE VIEW lmdb AS
   SELECT * FROM read_parquet('E:/workdata/YOUR_PROJECT/cleaned-data/parquet-registers/lmdb/**/*.parquet')"
)

# 3. Get a lazy dplyr table - nothing is loaded into R yet
lmdb <- tbl(con, "lmdb")

# Same three steps for LPR, with a different path and view name, e.g.:
# dbExecute(con, "CREATE VIEW lpr AS SELECT * FROM read_parquet('.../lpr/**/*.parquet')")
# lpr <- tbl(con, "lpr")

# When finished:
# dbDisconnect(con, shutdown = TRUE)
```

Do **not** feed medication extractors with `duckplyr::read_parquet_duckdb()` or
`fastreg::read_register()`. `fastreg` is still useful for SAS->parquet setup on
DST: [fastreg](https://dp-next.github.io/fastreg/). Tiny in-memory frames work
for toy runs (see the vignette).

## Recommended workflow

The usual study run is: load lists -> extract diagnoses -> extract medications
(prefer batch) -> merge -> reshape -> count. How many *conditions* you get
depends on your CSV code lists, not on the package API.

Dates such as `1995-01-01`, `1997-01-01`, `2015-01-01`, and `2018-12-31` in
the examples below are **placeholders**. Replace them with your study window.

```r
library(regmorbidity)

# ---------------------------------------------------------------------------
# 0. Registers (you supply these)
#    On DST: open parquet via DuckDB as in the section above -> objects `lpr`
#    and `lmdb`. Below we assume those already exist.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 1. Load code lists
#    Bundled lists live under inst/extdata/codelists/ (36 CSVs). You do NOT
#    pass file paths here unless you have your own edited copies.
#
#    Default: load ALL bundled rows (ATC + first-pass ICD). Diagnosis and
#    medication extractors then keep the vocabulary they need themselves.
#
#    To run every bundled condition:
codes <- mb_codelist()
#
#    To run a named subset instead (same CSVs, filtered after load):
# codes <- mb_codelist(conditions = c("hypertension", "diabetes"))
#
#    Own lists instead of bundled? Point at a folder of CSVs or one big CSV:
# codes <- mb_codelist("path/to/my_codelists/")
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 2. Extract diagnoses from LPR
#    One matching diagnosis = the condition. Writes one onset .rds per
#    condition into outdir (here: data/dx/).
#    Tip: use a fresh outdir, or clear old .rds first. Resume skips files that
#    already exist, so a dirty folder of leftover results can silently leave
#    stale onsets in place.
# ---------------------------------------------------------------------------
mb_extract_diagnosis(
  lpr,
  codes  = codes,
  outdir = "data/dx",
  from   = as.Date("1995-01-01"),  # inclusive start; REPLACE with your study
  to     = as.Date("2018-12-31")   # optional inclusive end; omit to leave open
)

# ---------------------------------------------------------------------------
# 3. Extract medications from LMDB - prefer batch
#    One SQL pass, low RAM, one onset .rds per ATC condition into data/rx/.
#    `from` is required (1997 recommended on DST - mother-CPR; see above).
#    `to` is optional - omit for no upper date bound.
# ---------------------------------------------------------------------------
mb_extract_medication_batch(
  lmdb,
  codes  = codes,
  outdir = "data/rx",
  from   = as.Date("1997-01-01"),  # required on LMDB; REPLACE if your study differs
  to     = as.Date("2018-12-31")   # optional; not required
)

# ---------------------------------------------------------------------------
# 4. Merge diagnosis + medication halves
#    Reads the two outdirs and returns one long data frame in memory:
#    one row per person x condition with an onset_date. No new files.
#
#    Example shape of `long` (illustrative):
#
#    | pnr | condition     | onset_date |
#    |-----|---------------|------------|
#    | 1   | hypertension  | 2008-01-01 |
#    | 1   | diabetes      | 2010-03-15 |
#    | 2   | hypertension  | 2009-06-01 |
#
#    This table is what "which conditions" looks like. The count step below
#    does not re-list them; it only adds how many each person has.
# ---------------------------------------------------------------------------
long <- mb_merge_all("data/dx", "data/rx")

# ---------------------------------------------------------------------------
# 5. Reshape and count
#    mb_to_wide: one row per person, one column per condition (onset dates).
#    mb_count_conditions: same wide table, plus:
#      - n_conditions  = how many conditions are present as of `as_of`
#                        (ever after onset: once met, always met)
#      - multimorbid   = TRUE when n_conditions >= 2
#    It does NOT list which conditions those are; the onset-date columns from
#    mb_to_wide already carry that.
# ---------------------------------------------------------------------------
wide <- mb_to_wide(long)
mb_count_conditions(wide, as_of = "2015-01-01")  # REPLACE with your as-of date
```

Danish register ICD codes often carry a leading **D** (`DI10` in the register
vs `I10` in published lists). The package normalises both forms before
matching.

If column names differ from the defaults, pass `id_col`, `code_col`, and
`date_col` on the extract call. Codes that never match the register produce
empty extracts and silent zeros in downstream counts - see
`mb_check_codes` in the vignette before a long DST run.

For sequential extract, `keep_events`, prevalence lookback, exclusion timing,
helpers that test or compare code lists, and a worked miniature:
`vignette("regmorbidity")`.

## Batch vs sequential medication methods

Two medication extractors exist on purpose. Prefer **batch** for onset
studies; keep **sequential** when you need raw events or per-condition
checkpoints.

| | `mb_extract_medication_batch` (prefer) | `mb_extract_medication` (sequential) |
|---|---|---|
| How it runs | One SQL pass over the register | Loops conditions one by one |
| RAM | Low - DuckDB returns onsets only | Higher if you keep events |
| Onset `.rds` | Same shape | Same shape |
| Raw events | No | Yes, with `keep_events = TRUE` |
| Checkpoints | Per rule-batch (coarser) | Per condition (finer resume) |

**Why both exist.** Onset studies usually only need the date a condition was
first met. Batch asks DuckDB for that onset in one pass and never materialises
matching dispensings in R - that is the efficient default. Prevalence
("evidence within the last N days") and some debugging need the raw matching
dispensings; those come from sequential extract with `keep_events = TRUE`,
which writes `<condition>_all_events.rds` beside the onset file. Batch has no
equivalent side file, so prevalence cannot be fed from a batch-only outdir.
Both still write the same onset `.rds` shape, so merge / load / wide / count
stay unchanged either way.

## Core workflow functions

Only the functions you need for the recommended run. Helpers to test or
compare code lists (`mb_check_codes`, `mb_lookup`, `mb_overlap`,
`mb_inspect_code_lengths`), prevalence, exclusions, and sequential detail
live in `vignette("regmorbidity")`.

- **`mb_codelist`** - load the bundled CSVs (or your own folder / CSV /
  data frame). Start here. Default `mb_codelist()` loads all 36 conditions.
- **`mb_extract_diagnosis`** - LPR onset: one matching diagnosis is the
  condition; write one onset `.rds` per condition.
- **`mb_extract_medication_batch`** - **prefer** for medication onset (one
  SQL pass, low RAM); does not keep raw events.
- **`mb_extract_medication`** - sequential backup; use when you need
  `keep_events = TRUE` (prevalence / debug) or per-condition disk checkpoints.
- **`mb_merge_all`** - combine diagnosis + medication extract directories into
  one long table (person x condition x onset).
- **`mb_to_wide`** - long onset table -> one row per person, one onset-date
  column per condition.
- **`mb_count_conditions`** - returns the same wide table plus `n_conditions`
  and `multimorbid` (`>= 2`) as of a date. Does not list which conditions;
  those remain the onset-date columns from `mb_to_wide`.
- **`mb_load_conditions`** - load one extract outdir to long (skips
  `*_all_events.rds`); use for medication-only studies.

## More detail

`vignette("regmorbidity")` - deepening only (assumes you read this README
first): worked miniature, batch vs sequential + `keep_events`, medication
rules, prevalence, merge logic, exclusions, and code-list QA helpers.

`?regmorbidity` - package help for every exported function.

`ASSUMPTIONS_AND_LIMITATIONS.txt`, `DECISIONS.md`, and `METHODS.md` are
**local working-tree notes** (often gitignored). They are not part of a
normal GitHub install. If you installed from GitHub, use this README, the
vignette, and `?help`.

Output is onset-only by default.
