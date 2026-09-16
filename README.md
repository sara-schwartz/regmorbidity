# regmorbidity

Comorbidity indicators from Danish register data, defined in **editable CSV code
lists** rather than in code.

This is not an implementation of any published index. The extraction logic is
derived from Jie Zhang's `dgt_medication.R`, which is known to run on DST; its
condition definitions were in turn inspired by the code lists of Prior et al.
(2016). The bundled lists are a starting point to revise, not an instrument to
cite.

**Read [`ASSUMPTIONS_AND_LIMITATIONS.txt`](ASSUMPTIONS_AND_LIMITATIONS.txt)
before reporting any number from this package.**

Authors: Jie Zhang, Sara Schwartz (saras@clin.au.dk)

---

## Where the code lists are

`inst/extdata/codelists/` — **one CSV per condition**. All 15 are **ATC only**
(`vocab_id = ATC`); there are no ICD-10 rows yet.

```
allergy.csv       diabetes.csv      hypertension.csv  osteoporosis.csv  prostate.csv
bipolar.csv       distress.csv      ihd.csv           pain.csv          pulmonary.csv
dementia.csv      dyslipidemia.csv  migraine.csv      epilepsy.csv      thyroid.csv
```

Each file is transcribed from one filter block in `dgt_medication.R`, and the
`note` column records the source line. A whole condition fits on a screen, so
one file can be sent to one clinician to review. `pain.csv`:

```csv
condition,condition_label,category,vocab_id,code,exclude,min_prescriptions,window_days,note
pain,Painful condition,Musculoskeletal,ATC,N02A,FALSE,4,365,dgt_medication.R line 221
pain,Painful condition,Musculoskeletal,ATC,N02BA51,FALSE,4,365,dgt_medication.R line 221
pain,Painful condition,Musculoskeletal,ATC,N02BE,FALSE,4,365,dgt_medication.R line 221
pain,Painful condition,Musculoskeletal,ATC,M01A,FALSE,4,365,dgt_medication.R line 221
pain,Painful condition,Musculoskeletal,ATC,M02A,FALSE,4,365,dgt_medication.R line 221
```

| column | meaning |
|---|---|
| `condition` | the key. Will need to match the ICD-10 rows for the same disease. |
| `condition_label` | free text, for reading. |
| `category` | optional grouping. |
| `vocab_id` | `ATC` or `ICD10`. |
| `code` | a code **prefix**. `C09` matches C09AA05, C09DB01, … |
| `exclude` | `TRUE` removes matches of this code (e.g. one drug out of a class). |
| `min_prescriptions` | dispensings needed before the condition counts. |
| `window_days` | they must fall within this many days. Blank = any time apart. |
| `note` | free text. |

`min_prescriptions` and `window_days` describe the **condition**, so they must be
the same on every row of a file. Two different answers is rejected rather than
silently resolved.

### What is not here, and why

There are no ICD-10 rows yet, so every condition Prior defines by diagnosis is
currently missing or understated. That is the gap; see `TODO.txt` section 3.

What is *not* a gap: gout, alcohol problems, glaucoma and Parkinson's have no
ATC codes here, and it turns out they have none in Prior either. His composite
uses diagnosis only for all four — the medication half is commented out in his
Stata with the note *"nye medicindata"* — and glaucoma is not a separate
condition at all, but part of his vision problem (`H40`). Anemia nominally adds
`B03`, where he notes *"Ingen medicindata"*.

### Using your own lists

```r
codes <- mb_codelist("my_codelists/")        # a folder of one CSV per disease
codes <- mb_codelist("all_conditions.csv")   # or a single combined table
codes <- mb_codelist(my_data_frame)          # or an object built in R
```

A file with no `condition` column takes its condition from the file name, so
`gigt.csv` becomes condition `gigt`. Otherwise file names are free — the
`condition` column decides.

To start from the bundled lists and edit from there:

```r
mb_write_codelist(mb_codelist(), "my_codelists")   # 15 files to edit
codes <- mb_codelist("my_codelists")               # read them back
```

`mb_codelist()` validates on load and refuses lists that would fail quietly: a
condition whose rows disagree about the rule, a condition made only of
exclusions, duplicated codes, a code that another code already covers.

---

## What the package needs from you

Nothing but an object `dplyr` works on. The package never opens a file itself —
you hand it a table, and it only ever calls `filter()`, `select()` and
`collect()` on it.

That means all of these work the same way:

```r
lmdb <- fastreg::read_register("lmdb")          # lazy DuckDB table
lmdb <- arrow::open_dataset("lmdb/")            # arrow Dataset
lmdb <- readRDS("lmdb.rds")                     # a plain data frame
```

Neither `fastreg` nor Parquet is a requirement. They are how you would sensibly
read a register this size, not something the package depends on.

What it does require is three columns:

| | |
|---|---|
| person id | `id_col`, default `PNR` |
| the **full** ATC code | `code_col`, default `ATC` — not a truncated level |
| dispensing date | `date_col`, default `EKSD` |

`atc2` and `year` are optional. Without `atc2` the codes are matched against the
full code column directly and the package says so; without `year`, set
`year_min = NULL`.

**One honest caveat.** The test suite runs on plain data frames, because that is
what can be tested away from DST. Working on a lazy DuckDB or arrow table is the
design, not a verified fact. `mb_inspect_codes()` is a good first call on a new
setup for exactly that reason: it exercises the same verbs and costs nothing.

---

## Extracting

```r
library(regmorbidity)

lmdb <- read_parquet_duckdb(here(registers, "lmdb")) |>
  select(PNR = pnr, ATC = atc, atc2, EKSD = eksd, year)

mb_inspect_codes(lmdb)          # check the column lengths first

codes <- mb_codelist("my_codelists")

mb_extract_medication(
  lmdb,
  conditions = NULL,            # NULL = every condition in the list
  codes      = codes,
  outdir     = "data/conditions"
)
```

One `.rds` per condition is written **as soon as that condition finishes**, so a
run that is interrupted resumes rather than restarting:

```r
mb_extract_medication(lmdb, codes = codes, outdir = "data/conditions")  # skips what is done
```

A condition that fails is recorded in the returned summary and the loop carries
on — so always check `summary$status`, because a batch that finished may still
contain errors. Configuration mistakes (a column that does not exist) stop
before the loop instead.

Then:

```r
long <- mb_load_conditions("data/conditions")     # one row per person-condition
wide <- mb_to_wide(long)                          # one column per condition
cnt  <- mb_count_conditions(wide, as_of = "2015-01-01")
```

### The output

One row per person **who met the criterion**:

| column | |
|---|---|
| `PNR` | person id |
| `condition` | |
| `onset_date` | date of the qualifying dispensing — the 2nd, or the 4th for pain |
| `onset_code` | which ATC code it was |
| `source` | `medication`, `diagnosis` or `both` |
| `n_records` | matching records for that person |
| `first_date`, `last_date` | |

---

## Three things that change your numbers

**`year_min` defaults to 1997.** That is what `dgt_medication.R` uses throughout,
so it is the default here. LMDB covers 1995–2022, so the default **deliberately
discards 1995 and 1996**. The reason for 1997 is not documented in Jie's script;
it is carried forward because it is what the working pipeline does, not because
it has been checked. Studying the mid-1990s means `year_min = NULL`. A related
consequence: onset dates near the start of follow-up are left-censored and must
not be read as incidence.

**Same-day pickups are de-duplicated.** Two packages collected on one pharmacy
visit are two rows in LMDB. Counted naively they satisfy a two-prescription rule
on their own, making "two prescriptions" mean "one visit". Dispensings are
de-duplicated per person-date before counting, as the Stata original does
(`bys pid eksd: keep if _n == 1`). Jie's `flag_users()` did not, so counts here
run **lower** than from her script. `dedupe_same_day = FALSE` turns it off.

**`window_days = 365` follows the Stata original.** `multsygdom_udtrak.do` gives
every dispensing 365.25 days of coverage and requires two coverages to overlap
(four for pain) — which is the same rule as "2 dispensings within 365 days". The
codebook agrees: every ATC row is `time_frame = last_year`. Jie's `flag_users()`
had no window at all, so this package is **stricter** than `dgt_medication.R` and
closer to the Stata. Blank the column to reproduce her script instead.

---

## Checking the code lists themselves

Three functions answer questions about the definitions, without touching data.

**Do two conditions share codes?** With 39 conditions there are 741 pairs, which
is more than anyone checks by eye.

```r
mb_overlap(codes)
#>   ATC  C02C > C02CA   hypertension covers prostate
```

Prefix-aware, so it catches a condition whose code merely *contains* another's,
and ICD-normalised, so a list mixing `DJ45` and `J45` still shows the overlap.
An overlap is not automatically wrong — Prior's own lists overlap — but one
record creating two conditions should be a decision, not a discovery.

**What would this record be counted as?** The question you want answered when a
prevalence looks too high.

```r
mb_lookup("C02CA01", codes)
#>     query    condition vocab_id  code exclude
#> 1 C02CA01 hypertension      ATC  C02C   FALSE
#> 2 C02CA01     prostate      ATC C02CA   FALSE
```

**What changed since last time?** Both sides are normalised first, so reordered
rows, letter case and a blank-versus-`FALSE` exclude cell are not reported.
Rule changes are read out separately, because they move every person's result
where a code change usually moves a few.

```r
mb_compare("codelists_v1/", "codelists_v2/")
#> 3 change(s) across 3 condition(s). 15 -> 15 conditions, 61 -> 61 codes.
#> Rule changes — these change results for everyone:
#>   pain window_days: 365 -> 730
```

---

## The atc2 trap

`atc2` holds 3 characters. A 4-character pattern matched against it returns
nothing — no error, no warning, an empty result that looks like a finding.

This is why extraction is two-stage and the codes are applied twice: a cheap
`grepl` on `atc2` to narrow the register, then the real codes against the full
`ATC` column after `collect()`. `mb_inspect_codes()` prints the actual length
distribution of both columns; run it before trusting a new extract.
`mb_check_codes()` goes further and reports any code in your list that matches
nothing at all in the register — which is how you tell a true zero from a typo.

Test 14 asserts that truncating the recorded codes to 3 characters reproduces
Jie's literal `atc2` pattern for all 15 conditions. If you edit a list and that
test fails, do not ignore it: the pushed-down filter would silently drop rows the
codes should have caught.

---

## Working from a single condition

`flag_users()` keeps the calling convention from `dgt_medication.R`, so it drops
into existing scripts:

```r
aht_users <- flag_users(aht_filtered, ATC, min_prescriptions = 2)
aht_users <- flag_users(aht_filtered, ATC, min_prescriptions = 2,
                        window_days = 365)      # new
```

---

## Diagnoses

`mb_extract_diagnosis()` is the ICD-10 counterpart and takes the same arguments:

```r
mb_extract_diagnosis(lpr, codes = codes, outdir = "data/diagnoses")
```

It needs no two-stage filter — there is no `atc2` equivalent — and no
prescription rule: one recorded diagnosis is the condition, and onset is the
earliest one. Records are de-duplicated per person-day like the medication side,
so `n_records` means *days with a record* on both halves rather than rows on one
and days on the other; one contact coding both `I500` and `I509` is one day.
`dedupe_same_day = FALSE` turns it off.

Danish register codes carry a leading D (`DI50`), published code lists do not
(`I50`). Both are accepted. The discrimination is exact rather than a guess: a
WHO code is a letter followed by digits, so a Danish code is D followed by a
*letter*. `DD86` is the Danish form of `D86`; `D86` is already WHO form.

Combine the halves with:

```r
long <- mb_merge_all(diagnosis_dir  = "data/diagnoses",
                     medication_dir = "data/conditions",
                     codes = codes)
```

`OR` takes the earlier of the two dates, `AND` requires both and takes the
later. Which applies comes from an optional `logic` column in the code list,
defaulting to `OR`; Prior uses `AND` only for epilepsy.

The merged output carries `icd_date`, `rx_date`, `n_icd` and `n_rx` alongside
`onset_date`, so you can always see which half drove the result. The two counts
are over the whole period, not within the window — they answer "how much of this
did they have", not "how many met the rule".

**Not implemented:** Prior's exclusion rules — dyslipidaemia only if not IHD,
and so on — reference *other* conditions, so they cannot be applied one
condition at a time. See `TODO.txt` section 1(b).

Note that the bundled ATC files use Jie's condition names (`hypertension`,
`ihd`, …) while Sara's codebook uses `d_bt`, `d_ihd`. Whichever the ICD-10 lists
use, the two halves must agree — renaming the `condition` column is enough.

---

## Installing

```r
# from the folder containing regmorbidity/
install.packages("regmorbidity", repos = NULL, type = "source")
library(regmorbidity)
?regmorbidity          # grouped overview of every function
```

`NAMESPACE` and everything under `man/` are generated by roxygen2 from the
comments in `R/`. After changing an `@export` tag or any documentation, run
`roxygen2::roxygenise()` — do not edit those files by hand.

`R CMD check` passes clean on the built tarball. Check the tarball, not the
source directory: in directory mode the check does not derive `Author` and
`Maintainer` from `Authors@R` and reports them missing, which is an artefact
rather than a problem to fix.

```sh
R CMD build regmorbidity && R CMD check regmorbidity_0.1.0.tar.gz
```

---

## Tests

```sh
Rscript tests/test_regmorbidity.R
```

99 checks over synthetic data. The ones worth knowing about:

- the prescription rule, its window, and same-day de-duplication on both halves
- that stage 2 keeps `N02A` (pain) apart from `N02C` (migraine), and that
  anchoring rejects a mid-string match
- that the derived stage-1 prefixes still reproduce `dgt_medication.R` — if you
  edit a code list and this fails, the pushed-down filter would silently drop
  rows the codes should have caught
- the Danish D prefix, in both directions: `DD86` → `D86`, `D86` unchanged
- the OR / AND merge, including that `AND` with an empty half warns
- checkpoint, resume, and that one failing condition does not abort the batch
- that `mb_compare()` reports a real change and stays quiet about reordering

Nothing here has been run against real register data.

Nothing here has been run against real DST data.
