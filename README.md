# regmorbidity

Condition indicators, with dates, from Danish register data - defined in
**editable CSV code lists** rather than in code. It computes no score, and is
not an implementation of the Danish Multimorbidity Index or of any other
published index. The bundled lists take their starting point in Prior et al.
(2016) and are a starting point to revise, not an instrument to cite.

Authors: Jie Zhang, Sara Schwartz (saras@clin.au.dk)

## Installing

```r
# from the folder containing regmorbidity/
install.packages("regmorbidity", repos = NULL, type = "source")
library(regmorbidity)
```

## Quick start

```r
library(regmorbidity)

codes <- mb_codelist()          # the 15 bundled ATC code lists
mb_inspect_codes(lmdb)          # check the register's column lengths first

mb_extract_medication(
  lmdb,                          # any object dplyr can filter()/select()/collect()
  codes  = codes,
  outdir = "data/conditions"     # one .rds per condition, written as it finishes
)

long <- mb_load_conditions("data/conditions")   # one row per person-condition
wide <- mb_to_wide(long)                        # one column per condition
mb_count_conditions(wide, as_of = "2015-01-01")
```

`vignette("regmorbidity")` is the full walkthrough: how this compares to
Prior et al. (2016) - what is the same, what is deliberately different -
code lists and checking them, diagnoses and merging the two halves,
point-in-time prevalence with `mb_prevalence()`, what is still open and who
it is blocked on, and next steps.

`?regmorbidity` gives the grouped overview of all 19 functions.

## Before reporting a number

This package has never been run on real register data. The source repository
carries `ASSUMPTIONS_AND_LIMITATIONS.txt`, the short list of choices that move
a number, and `DECISIONS.md`, the full argument for each of them - read before
trusting output from a real run.
