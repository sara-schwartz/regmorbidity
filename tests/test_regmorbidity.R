# .............................................................................
# test_regmorbidity.R - the whole test suite
#
# PURPOSE
#   Prove that the rules the package claims to apply are the rules it actually
#   applies. Every check below corresponds to a decision documented in
#   ASSUMPTIONS_AND_LIMITATIONS.txt; if a check fails, either the code changed
#   or that document is now wrong.
#
# RUN     Rscript tests/test_regmorbidity.R
# INPUT   the bundled code lists; all patient data here is synthetic
# OUTPUT  PASS/FAIL lines on stdout. Nothing is written outside tempdir().
#
# WHAT THESE TESTS DO NOT COVER
#   Nothing here has touched real register data. Behaviour on the full LMDB -
#   especially memory during collect() for pain and hypertension - is unverified.
#   See TODO.txt section 5.
#
# WHY SO MUCH ATTENTION TO EMPTY RESULTS
#   The failure mode that matters in this package is silent: a wrong column, a
#   truncated code, an unanchored pattern all produce an empty result rather
#   than an error, and an empty result reads as a finding. Several tests below
#   exist only to prove that a specific wrong answer is NOT returned quietly.
#
# CONTENTS
#   1. Setup
#   1. Code lists
#   2. Validation
#   3. One CSV per disease
#   4. The prescription rule
#   5. Two-stage matching
#   6. Anchoring
#   7. Checkpoint and resume
#   8. Error handling
#   9. Combining output
#   10. flag_users()
#   11. Missing prefilter column
#   12. Short codes
#   13. mb_check_codes()
#   14. Provenance
#   15. from/to date window
# .............................................................................

suppressMessages({library(dplyr); library(data.table); library(rlang)})
# 1. Setup ----

# Works both from an installed package and straight from the source tree, so the
# suite can be run before the package is installable.
if (requireNamespace("regmorbidity", quietly = TRUE)) {
  library(regmorbidity)
  # Internals the suite reaches into. They are not public on purpose, so they
  # have to be fetched explicitly here - and R CMD check runs against the
  # INSTALLED package, which is the configuration that catches a missing one.
  for (nm in c("MB_CODELIST_COLS", "mb_run_conditions", "mb_restrict_ids",
               "mb_lookback_days", "mb_rule_batches", "mb_batch_query",
               "mb_is_lazy",
               # demoted from export (DECISIONS.md 1.3); suite still covers them.
               # mb_overlap / mb_lookup are public again (2026-09-28) — not listed here.
               "mb_validate_codelist", "mb_condition_logic", "mb_normalize_icd10",
               "mb_flag_users", "mb_prevalence_all", "mb_compare")) {
    assign(nm, get(nm, envir = asNamespace("regmorbidity")))
  }
  CODES <- system.file("extdata", "codelists", package = "regmorbidity")
} else {
  pkg <- if (dir.exists("../R")) ".." else "."
  for (f in list.files(file.path(pkg, "R"), full.names = TRUE)) source(f)
  CODES <- file.path(pkg, "inst", "extdata", "codelists")
}
stopifnot(dir.exists(CODES))

#' Report one check. Deliberately minimal - no testthat dependency, so the suite
#' runs anywhere R does, including a locked-down analysis server.
ok <- function(label, cond) {
  cat(if (isTRUE(cond)) "  PASS  " else "  FAIL  ", label, "\n")
}
# 1. Code lists ----

cat("\n=== 1. code list loads and validates ===\n")
codes <- mb_codelist(CODES)
ok("36 conditions (distress archived)", length(unique(codes$condition)) == 36)
ok("309 rows",      nrow(codes) == 309)
ok("pain needs 4",  unique(codes$min_prescriptions[codes$condition == "pain"]) == 4)
ok("window 365",    unique(codes$window_days[codes$condition == "hypertension"]) == 365)
ok("HTN excludes C02CA (prostate ATC)",
   any(codes$condition == "hypertension" & codes$code == "C02CA" & codes$exclude))
ok("distress absent from default bundled set (Sara 2026-09-28)",
   !"distress" %in% codes$condition)
ok("epilepsy ICD is G40-G41 only",
   identical(sort(unique(codes$code[codes$condition == "epilepsy" &
                                    codes$vocab_id == "ICD10"])),
             c("G40", "G41")))
ok("heart_failure is I50 only",
   identical(unique(codes$code[codes$condition == "heart_failure"]), "I50"))
ok("ckd first-pass is N18 only",
   identical(unique(codes$code[codes$condition == "ckd"]), "N18"))
ok("mb_codelist(conditions=) keeps named conditions",
   identical(sort(unique(mb_codelist(CODES, conditions = c("pain", "migraine"))$condition)),
             c("migraine", "pain")))
ok("mb_codelist(conditions=) refuses unknown names",
   inherits(try(mb_codelist(CODES, conditions = "not_a_condition"), silent = TRUE),
            "try-error"))

# C02CA is prostate; HTN must not count it even though C02C would match.
c02ca_lmdb <- data.frame(
  pnr = c("x","x","y","y"),
  atc = c("C02CA01","C02CA01","C09AA05","C09AA05"),
  atc2 = c("C02","C02","C09","C09"),
  eksd = as.Date(c("2010-01-01","2010-02-01","2010-01-01","2010-02-01")),
  stringsAsFactors = FALSE
)
c02ca_htn <- mb_extract_medication(c02ca_lmdb, "hypertension", codes = codes,
                                   from = as.Date("1990-01-01"),
                                   verbose = FALSE)[["hypertension"]]
c02ca_pros <- mb_extract_medication(c02ca_lmdb, "prostate", codes = codes,
                                    from = as.Date("1990-01-01"),
                                    verbose = FALSE)[["prostate"]]
ok("C02CA01 does not count as hypertension (exclude row)",
   !"x" %in% c02ca_htn$pnr)
ok("C09AA05 still counts as hypertension",
   "y" %in% c02ca_htn$pnr)
ok("C02CA01 still counts as prostate",
   "x" %in% c02ca_pros$pnr)
# 2. Validation ----

cat("\n=== 2. validation catches silent mistakes ===\n")
bad <- codes; bad$min_prescriptions[bad$condition == "hypertension"][1] <- 3L
ok("conflicting rule rejected",
   inherits(try(mb_validate_codelist(bad, quiet = TRUE), silent = TRUE), "try-error"))
bad2 <- codes; bad2$exclude[bad2$condition == "thyroid"] <- TRUE
ok("exclude-only condition rejected",
   inherits(try(mb_validate_codelist(bad2, quiet = TRUE), silent = TRUE), "try-error"))
# 3. One CSV per disease ----

cat("\n=== 3. one CSV per disease ===\n")
d <- file.path(tempdir(), "per_disease")
unlink(d, recursive = TRUE)
paths <- mb_write_codelist(codes, d)
ok("36 files written", length(paths) == 36)
back <- mb_codelist(d)
ok("round-trips identically",
   isTRUE(all.equal(codes[order(codes$condition, codes$code), MB_CODELIST_COLS],
                    back[order(back$condition, back$code), MB_CODELIST_COLS],
                    check.attributes = FALSE)))

# a user-added file with no `condition` column: name comes from the file name
write.csv(data.frame(vocab_id = "ATC", code = "M04", min_prescriptions = 2,
                     window_days = 365),
          file.path(d, "gigt.csv"), row.names = FALSE)
back2 <- mb_codelist(d)
ok("condition taken from file name", "gigt" %in% back2$condition)
# 4. The prescription rule ----

cat("\n=== 4. the 2-prescription rule ===\n")
# p1: two C09 a month apart          -> qualifies, onset = 2nd
# p2: two C09 three years apart      -> outside a 365-day window, no onset
# p3: two boxes on the SAME day      -> one dispensing after dedupe, no onset
# p4: one C09 only                   -> no onset
lmdb <- data.frame(
  pnr  = c("p1","p1", "p2","p2", "p3","p3", "p4"),
  atc  = c("C09AA05","C09CA01", "C09AA05","C09AA05", "C09AA05","C09CA01", "C09AA05"),
  atc2 = "C09",
  eksd = as.Date(c("2010-01-10","2010-02-10", "2010-01-10","2013-01-10",
                   "2010-05-01","2010-05-01", "2010-01-10")),
  year = c(2010,2010, 2010,2013, 2010,2010, 2010),
  stringsAsFactors = FALSE
)
res <- mb_extract_medication(lmdb, "hypertension", codes = codes, verbose = FALSE,
                           from = as.Date("1990-01-01"))[["hypertension"]]
ok("only p1 qualifies", identical(res$pnr, "p1"))
ok("onset = 2nd dispensing", res$onset_date == as.Date("2010-02-10"))
ok("p3 same-day pair does not qualify", !"p3" %in% res$pnr)

res_nodedupe <- mb_extract_medication(lmdb, "hypertension", codes = codes,
                                   dedupe_same_day = FALSE, verbose = FALSE,
                           from = as.Date("1990-01-01"))[["hypertension"]]
ok("without dedupe, p3 would qualify", "p3" %in% res_nodedupe$pnr)

res_nowin <- mb_extract_medication(lmdb, "hypertension",
                                codes = transform(codes, window_days = NA_integer_),
                                verbose = FALSE,
                           from = as.Date("1990-01-01"))[["hypertension"]]
ok("with no window, p2 qualifies too", all(c("p1","p2") %in% res_nowin$pnr))
# 5. Two-stage matching ----

cat("\n=== 5. stage 2 separates codes sharing a stage-1 prefix ===\n")
# N02A = pain, N02C = migraine. Both sit under atc2 = "N02".
lmdb2 <- data.frame(
  pnr  = c(rep("a", 4), rep("b", 2)),
  atc  = c("N02AA05","N02AA05","N02AA05","N02AA05", "N02CC01","N02CC01"),
  atc2 = "N02",
  eksd = as.Date(c("2010-01-01","2010-02-01","2010-03-01","2010-04-01",
                   "2010-01-01","2010-02-01")),
  year = 2010, stringsAsFactors = FALSE
)
r <- mb_extract_medication(lmdb2, c("pain","migraine"), codes = codes, verbose = FALSE,
                           from = as.Date("1990-01-01"))
ok("pain finds only a", identical(r$pain$pnr, "a"))
ok("migraine finds only b", identical(r$migraine$pnr, "b"))
ok("pain onset = 4th (min_prescriptions = 4)",
   r$pain$onset_date == as.Date("2010-04-01"))
# 6. Anchoring ----

cat("\n=== 6. anchoring ===\n")
# An unanchored 'C09' would also match a code that merely contains it.
lmdb3 <- data.frame(pnr = c("x","x"), atc = c("XC09AA","XC09AA"), atc2 = "C09",
                    eksd = as.Date(c("2010-01-01","2010-02-01")), year = 2010,
                    stringsAsFactors = FALSE)
r3 <- mb_extract_medication(lmdb3, "hypertension", codes = codes, verbose = FALSE,
                           from = as.Date("1990-01-01"))[["hypertension"]]
ok("mid-string match rejected", nrow(r3) == 0)
# 7. Checkpoint and resume ----

cat("\n=== 7. checkpoint and resume ===\n")
od <- file.path(tempdir(), "out"); unlink(od, recursive = TRUE)
s1 <- mb_extract_medication(lmdb, c("hypertension","dyslipidemia"), codes = codes,
                         outdir = od, verbose = FALSE,
                           from = as.Date("1990-01-01"))
ok("two files written", length(list.files(od, pattern = "\\.rds$")) == 2)
ok("both marked done", all(s1$status == "done"))
s2 <- mb_extract_medication(lmdb, c("hypertension","dyslipidemia"), codes = codes,
                         outdir = od, verbose = FALSE,
                           from = as.Date("1990-01-01"))
ok("second run skips both", all(s2$status == "skipped"))
s3 <- mb_extract_medication(lmdb, c("hypertension","dyslipidemia"), codes = codes,
                         outdir = od, resume = FALSE, verbose = FALSE,
                           from = as.Date("1990-01-01"))
ok("resume = FALSE re-runs", all(s3$status == "done"))
# 8. Error handling ----

cat("\n=== 8. a failing condition does not kill the batch ===\n")
lmdb_bad <- rbind(lmdb,
                  data.frame(pnr = "p9", atc = "C10AA01", atc2 = "C10",
                             eksd = as.Date("2010-01-01"), year = 2010))
lmdb_bad$eksd <- "not-a-date"   # both conditions now match, both fail to parse
sfail <- suppressWarnings(
  mb_extract_medication(lmdb_bad, c("hypertension","dyslipidemia"), codes = codes,
                     outdir = od, resume = FALSE, verbose = FALSE,
                           from = as.Date("1990-01-01")))
ok("errors recorded, loop continues", nrow(sfail) == 2 && all(sfail$status == "error"))

# A missing column is a configuration error and should stop before the loop.
ok("missing column fails fast",
   inherits(try(mb_extract_medication(lmdb, "hypertension", codes = codes,
                                   date_col = "NOPE", verbose = FALSE,
                           from = as.Date("1990-01-01")),
                silent = TRUE), "try-error"))
# 9. Combining output ----

cat("\n=== 9. load back, widen, count ===\n")
unlink(od, recursive = TRUE)
mb_extract_medication(lmdb2, c("pain","migraine"), codes = codes, outdir = od,
                   verbose = FALSE,
                           from = as.Date("1990-01-01"))
long <- mb_load_conditions(od)
ok("2 person-condition rows", nrow(long) == 2)
wide <- mb_to_wide(long, id_col = "pnr")
ok("one column per condition", all(c("pain","migraine") %in% names(wide)))
cnt <- mb_count_conditions(wide, id_col = "pnr")
ok("each person has 1 condition", all(cnt$n_conditions == 1))
ok("nobody multimorbid", !any(cnt$multimorbid))
cnt_early <- mb_count_conditions(wide, id_col = "pnr", as_of = "2010-01-15")
ok("as_of excludes later onsets", sum(cnt_early$n_conditions) == 0)
# 10. flag_users() ----

cat("\n=== 10. flag_users() keeps Jie's calling convention ===\n")
disp <- data.frame(pnr = c("q","q","q"), atc = "C09AA05",
                   eksd = as.Date(c("2010-01-01","2010-02-01","2010-03-01")),
                   stringsAsFactors = FALSE)
f1 <- mb_flag_users(disp, atc, min_prescriptions = 2)
f2 <- mb_flag_users(disp, "atc", min_prescriptions = 2)
ok("bare column name works", nrow(f1) == 1 && f1$onset_date == as.Date("2010-02-01"))
ok("string column name works", identical(f1, f2))
ok("n_records reported",f1$n_records == 3)
# 11. Missing prefilter column ----

cat("\n=== 11. missing prefilter column falls back ===\n")
lmdb4 <- lmdb[, c("pnr","atc","eksd","year")]
r4 <- mb_extract_medication(lmdb4, "hypertension", codes = codes, verbose = FALSE,
                           from = as.Date("1990-01-01"))[["hypertension"]]
ok("still finds p1 without atc2", identical(r4$pnr, "p1"))
# 12. Short codes ----

cat("\n=== 12. codes shorter than the prefilter column ===\n")
short <- rbind(codes[codes$condition == "hypertension", ][1, ])
short$code <- "C0"; short$condition <- "short_code"
lmdb5 <- data.frame(pnr = c("z","z"), atc = c("C09AA05","C01DA02"), atc2 = c("C09","C01"),
                    eksd = as.Date(c("2010-01-01","2010-02-01")), year = 2010,
                    stringsAsFactors = FALSE)
r5 <- mb_extract_medication(lmdb5, "short_code", codes = rbind(codes, short),
                         verbose = FALSE,
                           from = as.Date("1990-01-01"))[["short_code"]]
ok("2-char code matched via full atc column", nrow(r5) == 1 && r5$pnr == "z")
# 13. mb_check_codes() ----

cat("\n=== 13. mb_check_codes finds codes that match nothing ===\n")
chk <- mb_check_codes(lmdb, codes = codes[codes$condition %in% c("hypertension","thyroid"), ])
ok("C09A found",  chk$found[chk$code == "C09A"])
ok("H03A missing", !chk$found[chk$code == "H03A"])
# 14. Provenance ----

cat("\n=== 14. stage-1 prefixes still reproduce dgt_medication.R ===\n")
# If someone edits a code list so its 3-character truncation no longer covers
# Jie's atc2 pattern, the pushed-down filter silently drops rows the codes
# should have caught. Lock the mapping down.
jie_atc2 <- list(
  hypertension = "C02|C03|C04|C07|C08|C09", dyslipidemia = "C10",
  ihd = "C01", diabetes = "A10", thyroid = "H03", pulmonary = "R03",
  allergy = "R06|R01", prostate = "C02|G04", osteoporosis = "M05|G03|H05",
  pain = "N02|M01|M02", migraine = "N02", epilepsy = "N03",
  bipolar = "N05", dementia = "N06"
  # distress (N06A) archived out of default LTC set — Sara 2026-09-28
)
bad <- character(0)
for (cond in names(jie_atc2)) {
  # ICD rows now share files with ATC; provenance check is ATC-only
  cc <- codes$code[codes$condition == cond & !codes$exclude &
                   codes$vocab_id == "ATC"]
  derived <- sort(unique(substr(cc, 1L, 3L)))
  expected <- sort(unique(strsplit(jie_atc2[[cond]], "|", fixed = TRUE)[[1]]))
  if (!identical(derived, expected)) bad <- c(bad, cond)
}
ok("ATC conditions match Jie's atc2 patterns", length(bad) == 0)
if (length(bad)) cat("     mismatched:", paste(bad, collapse = ", "), "\n")
# 15. from/to date window ----

cat("\n=== 15. from/to date window (from required on LMDB) ===\n")
ok("medication from formal still defaults to NULL (runtime error if left)",
   is.null(formals(mb_extract_medication)$from))
ok("to defaults to NULL",   is.null(formals(mb_extract_medication)$to))
ok("batch from formal still defaults to NULL (runtime error if left)",
   is.null(formals(mb_extract_medication_batch)$from))
ok("diagnosis from defaults to NULL",
   is.null(formals(mb_extract_diagnosis)$from))

old <- data.frame(pnr = c("o","o"), atc = "C09AA05", atc2 = "C09",
                  eksd = as.Date(c("1995-01-01","1995-06-01")),
                  stringsAsFactors = FALSE)
ok("NULL from errors on medication extract",
   inherits(try(mb_extract_medication(old, "hypertension", codes = codes,
                           verbose = FALSE), silent = TRUE), "try-error"))
ok("NULL from errors on medication batch extract",
   inherits(try(mb_extract_medication_batch(old, codes = codes,
                           verbose = FALSE), silent = TRUE), "try-error"))
ok("early from keeps the entire register, including 1995",
   nrow(mb_extract_medication(old, "hypertension", codes = codes,
                           from = as.Date("1990-01-01"),
                           verbose = FALSE)[["hypertension"]]) == 1)
ok("from = 1997-01-01 drops 1995 dispensings",
   nrow(mb_extract_medication(old, "hypertension", codes = codes,
                           from = as.Date("1997-01-01"),
                           verbose = FALSE)[["hypertension"]]) == 0)

# Inclusive bounds on both ends, and a mid-window person who qualifies.
win <- data.frame(
  pnr = c("a","a", "b","b", "c","c"),
  atc = "C09AA05", atc2 = "C09",
  eksd = as.Date(c("1999-12-01","2000-01-15",
                   "2000-06-01","2000-07-01",
                   "2001-01-01","2001-02-01")),
  stringsAsFactors = FALSE)
win_res <- mb_extract_medication(win, "hypertension", codes = codes,
                                 from = as.Date("2000-01-01"),
                                 to   = as.Date("2000-12-31"),
                                 verbose = FALSE)[["hypertension"]]
ok("from/to inclusive: person fully inside the window qualifies",
   "b" %in% win_res$pnr)
ok("from/to inclusive: person with both dates after `to` is dropped",
   !"c" %in% win_res$pnr)
ok("from/to inclusive: person whose only in-window date cannot meet min=2 alone",
   !"a" %in% win_res$pnr)

ok("bad from errors clearly",
   inherits(try(mb_extract_medication(old, "hypertension", codes = codes,
                                   from = "not-a-date", verbose = FALSE),
                silent = TRUE), "try-error"))

cat("\nDone.\n")


# 16. ICD-10 code normalisation ----

cat("\n=== 16. the Danish D prefix ===\n")
# The discrimination that matters: DD86 is the Danish form of D86, but D86 is
# already WHO form. Stripping blindly would turn D86 into 86.
ok("DE10 -> E10", mb_normalize_icd10("DE10") == "E10")
ok("DD86 -> D86", mb_normalize_icd10("DD86") == "D86")
ok("D86 unchanged (second char is a digit)", mb_normalize_icd10("D86") == "D86")
ok("I50 unchanged", mb_normalize_icd10("I50") == "I50")
ok("dots and case handled", mb_normalize_icd10("dj30.1") == "J301")


# 17. mb_extract_diagnosis() ----

cat("\n=== 17. mb_extract_diagnosis ===\n")
dx_codes <- data.frame(
  condition = c("hf", "hf", "connective"),
  vocab_id  = "ICD10",
  code      = c("I50", "DI50", "D86"),   # same code both conventions
  stringsAsFactors = FALSE
)
lpr <- data.frame(
  pnr      = c("a", "a", "b", "c"),
  C_DIAG   = c("DI500", "DI509", "DD86", "DI10"),
  D_INDDTO = as.Date(c("2005-03-01", "2001-01-01", "2009-05-05", "2000-01-01")),
  stringsAsFactors = FALSE
)
dx <- mb_extract_diagnosis(lpr, codes = dx_codes, verbose = FALSE)
ok("heart failure finds a", identical(dx$hf$pnr, "a"))
ok("onset is the EARLIEST diagnosis", dx$hf$onset_date == as.Date("2001-01-01"))
ok("both records counted", dx$hf$n_records == 2)
ok("source labelled", dx$hf$source == "diagnosis")
ok("D86 written without the D prefix still matches DD86",
   identical(dx$connective$pnr, "b"))
ok("DI10 not picked up by either condition",
   !"c" %in% c(dx$hf$pnr, dx$connective$pnr))


# 18. Merging the two halves ----

cat("\n=== 18. OR / AND merge ===\n")
dxh <- data.frame(pnr = c("p","q"), condition = "x", source = "diagnosis",
                  onset_date = as.Date(c("2005-01-01","2007-01-01")),
                  stringsAsFactors = FALSE)
rxh <- data.frame(pnr = c("p","r"), condition = "x", source = "medication",
                  onset_date = as.Date(c("2003-01-01","2008-01-01")),
                  stringsAsFactors = FALSE)

m_or <- mb_merge_conditions(dxh, rxh, logic = "OR")
ok("OR keeps everyone", setequal(m_or$pnr, c("p","q","r")))
ok("OR takes the earlier date",
   m_or$onset_date[m_or$pnr == "p"] == as.Date("2003-01-01"))
ok("source = both where both halves present",
   m_or$source[m_or$pnr == "p"] == "both")
ok("source = diagnosis where only that half",
   m_or$source[m_or$pnr == "q"] == "diagnosis")

m_and <- mb_merge_conditions(dxh, rxh, logic = "AND")
ok("AND keeps only p", identical(m_and$pnr, "p"))
ok("AND takes the LATER date",
   m_and$onset_date == as.Date("2005-01-01"))

ok("AND with an empty half warns and returns nothing",
   nrow(suppressWarnings(
     mb_merge_conditions(dxh, NULL, logic = "AND"))) == 0)
ok("OR with one half still works",
   nrow(mb_merge_conditions(dxh, NULL, logic = "OR")) == 2)
ok("bad logic rejected",
   inherits(try(mb_merge_conditions(dxh, rxh, logic = "XOR"), silent = TRUE),
            "try-error"))


# 19. logic column and merging from disk ----

cat("\n=== 19. logic from the code list ===\n")
lg <- data.frame(condition = c("a","b"), vocab_id = "ICD10",
                 code = c("I50","G40"), logic = c(NA, "AND"),
                 stringsAsFactors = FALSE)
lg <- mb_codelist(lg, validate = FALSE)
ok("blank logic defaults to OR", mb_condition_logic(lg, "a") == "OR")
ok("AND is read from the column", mb_condition_logic(lg, "b") == "AND")

dxd <- file.path(tempdir(), "dx"); unlink(dxd, recursive = TRUE)
rxd <- file.path(tempdir(), "rx"); unlink(rxd, recursive = TRUE)
mb_extract_diagnosis(lpr, codes = dx_codes, outdir = dxd, verbose = FALSE)
all_merged <- mb_merge_all(diagnosis_dir = dxd, codes = dx_codes,
                           verbose = FALSE)
ok("merge_all works with only a diagnosis half", nrow(all_merged) == 2)
ok("both conditions present",
   setequal(all_merged$condition, c("hf", "connective")))


# 20. mb_compare() ----

cat("\n=== 20. code list diffing ===\n")
base_cl <- data.frame(
  condition = c("dm","dm","bt","gone"),
  vocab_id  = c("atc","ICD10","atc","atc"),
  code      = c("A10A","E10","C09","M04"),
  min_prescriptions = c(2,NA,2,2),
  window_days = c(365,NA,365,365),
  stringsAsFactors = FALSE)

ok("identical lists report nothing",
   nrow(mb_compare(base_cl, base_cl, verbose = FALSE)) == 0)

# formatting-only differences must NOT be reported
noise <- base_cl
noise$code <- tolower(noise$code)
noise <- noise[rev(seq_len(nrow(noise))), ]
ok("case and row order are not changes",
   nrow(mb_compare(base_cl, noise, verbose = FALSE)) == 0)

edited <- base_cl[base_cl$condition != "gone", ]
edited$code[edited$condition == "dm" & edited$vocab_id == "atc"] <- "A10B"
edited <- rbind(edited, data.frame(condition = "new", vocab_id = "ICD10",
                                   code = "I50", min_prescriptions = NA,
                                   window_days = NA, stringsAsFactors = FALSE))
edited$window_days[edited$condition == "bt"] <- 730

d <- mb_compare(base_cl, edited, verbose = FALSE)
ok("condition removed detected",
   any(d$change == "condition removed" & d$condition == "gone"))
ok("condition added detected",
   any(d$change == "condition added" & d$condition == "new"))
ok("code removed detected",
   any(d$change == "code removed" & d$item == "A10A"))
ok("code added detected",
   any(d$change == "code added" & d$item == "A10B"))
ok("rule change detected with old and new",
   any(d$change == "rule changed" & d$item == "window_days" &
       d$old == "365" & d$new == "730"))
ok("untouched vocab not reported",
   !any(d$condition == "dm" & d$item == "E10"))

# a code moving from include to exclude is one change, not two
flip <- base_cl
flip$exclude <- c(FALSE, FALSE, TRUE, FALSE)
f <- mb_compare(base_cl, flip, verbose = FALSE)
ok("include -> exclude is one 'exclude changed' row",
   nrow(f) == 1 && f$change == "exclude changed" && f$item == "C09")

# labels are cosmetic and off by default
lab <- base_cl; lab$condition_label <- "Renamed"
ok("label change hidden by default",
   nrow(mb_compare(base_cl, lab, verbose = FALSE)) == 0)
ok("label change shown on request",
   any(mb_compare(base_cl, lab, include_labels = TRUE,
                  verbose = FALSE)$change == "label changed"))


# 21. mb_overlap() and mb_lookup() ----

cat("\n=== 21. overlaps between conditions ===\n")
ov_codes <- rbind(
  data.frame(condition="hyp",  vocab_id="ATC",   code=c("C02C","C09A")),
  data.frame(condition="pros", vocab_id="ATC",   code="C02CA"),
  data.frame(condition="kol",  vocab_id="ICD10", code=c("DJ44","DJ45")),
  data.frame(condition="asth", vocab_id="ICD10", code="DJ45",
             stringsAsFactors = FALSE))

ov <- mb_overlap(ov_codes, verbose = FALSE)
ok("prefix overlap found (C02C covers C02CA)",
   any(ov$code_a == "C02C" & ov$code_b == "C02CA"))
ok("identical overlap found (DJ45 in two conditions)",
   any(ov$relation == "identical" & ov$code_a == "DJ45"))
ok("non-overlapping codes not reported", !any(ov$code_a == "C09A"))

ok("disjoint list reports nothing",
   nrow(mb_overlap(data.frame(condition=c("a","b"), vocab_id="ATC",
                              code=c("A10","C09"), stringsAsFactors=FALSE),
                   verbose = FALSE)) == 0)

# an exclusion narrows a condition; it is not an overlap with another
excl <- rbind(ov_codes[ov_codes$condition %in% c("hyp","pros"), ],
              data.frame(condition="pros", vocab_id="ATC", code="C02C",
                         stringsAsFactors = FALSE))
excl$exclude <- excl$code == "C02C" & excl$condition == "pros"
ok("exclude rows ignored",
   !any(mb_overlap(excl, verbose = FALSE)$code_b == "C02C"))

cat("\n=== 21b. what does one code become ===\n")
lk <- mb_lookup("C02CA01", ov_codes)
ok("one dispensing maps to two conditions",
   setequal(lk$condition, c("hyp","pros")))
ok("ICD lookup works across the D convention",
   nrow(mb_lookup("J450", ov_codes)) == 2)
ok("unknown code returns nothing", nrow(mb_lookup("ZZZZ", ov_codes)) == 0)


# 22. Diagnosis de-duplication and merge counts ----

cat("\n=== 22. counts on both halves ===\n")
lpr2 <- data.frame(
  pnr      = c("a","a","a","b"),
  C_DIAG   = c("DI500","DI509","DI50","DI50"),   # two codes, same day
  D_INDDTO = as.Date(c("2005-01-01","2005-01-01","2008-01-01","2001-01-01")),
  stringsAsFactors = FALSE)
hf <- data.frame(condition="hf", vocab_id="ICD10", code="I50",
                 stringsAsFactors = FALSE)

d_dedup <- mb_extract_diagnosis(lpr2, codes = hf, verbose = FALSE)$hf
ok("same-day records count once by default",
   d_dedup$n_records[d_dedup$pnr == "a"] == 2)
d_raw <- mb_extract_diagnosis(lpr2, codes = hf, dedupe_same_day = FALSE,
                           verbose = FALSE)$hf
ok("without dedup they count separately",
   d_raw$n_records[d_raw$pnr == "a"] == 3)
ok("onset unaffected by dedup",
   d_dedup$onset_date[d_dedup$pnr == "a"] == as.Date("2005-01-01"))

rxh2 <- data.frame(pnr = "a", condition = "hf", onset_date = as.Date("2004-01-01"),
                   n_records = 7L, stringsAsFactors = FALSE)
mg <- mb_merge_conditions(d_dedup, rxh2, logic = "OR")
ok("n_icd carried into the merge", mg$n_icd[mg$pnr == "a"] == 2)
ok("n_rx carried into the merge",  mg$n_rx[mg$pnr == "a"] == 7)
ok("n_rx is NA where there is no medication half",
   is.na(mg$n_rx[mg$pnr == "b"]))


# 23. NAMESPACE is generated, and exports only what it should ----

cat("\n=== 23. the public surface ===\n")
# NAMESPACE is generated by roxygen2 from the @export tags, so it cannot drift
# from them. What it CAN do is quietly grow: an @export tag added to a helper
# makes it public, and a public function is a promise. This checks the surface
# is the one we meant.
pkg_root <- if (dir.exists("../R")) ".." else if (dir.exists("R")) "." else NA

if (is.na(pkg_root)) {
  cat("  SKIP   source tree not found (running from an installed package)\n")
} else {
  ns_lines <- readLines(file.path(pkg_root, "NAMESPACE"), warn = FALSE)
  ok("NAMESPACE is generated, not hand-edited",
     any(grepl("^# Generated by roxygen2", ns_lines)))

  declared <- sort(unlist(lapply(parse(file.path(pkg_root, "NAMESPACE")),
    function(e) if (is.call(e) && identical(e[[1]], as.name("export")))
                  as.character(e[[2]]) else NULL)))

  expected <- sort(c(
    "mb_codelist", "mb_write_codelist",
    "mb_inspect_codes", "mb_check_codes",
    "mb_lookup", "mb_overlap",
    "mb_extract_medication", "mb_extract_medication_batch",
    "mb_extract_diagnosis",
    "mb_merge_conditions", "mb_merge_all",
    "mb_load_conditions", "mb_to_wide",
    "mb_count_conditions", "mb_prevalence",
    "mb_apply_exclusions"))

  ok("exactly the intended 16 functions are public",
     identical(declared, expected))
  if (!identical(declared, expected)) {
    cat("     unexpectedly public:", paste(setdiff(declared, expected), collapse = ", "), "\n")
    cat("     expected but absent:", paste(setdiff(expected, declared), collapse = ", "), "\n")
  }

  # the per-condition workers and the driver are implementation, not API
  ok("internal helpers stay internal",
     !any(c("mb_extract_one", "mb_diagnose_one", "mb_run_conditions",
            "mb_onset", "mb_pattern", "mb_half",
            "mb_validate_codelist", "mb_condition_logic", "mb_normalize_icd10",
            "mb_flag_users", "mb_prevalence_all",
            "mb_compare") %in% declared))
}

cat("\nDone.\n")


# 24. Regressions found in the 2026-09-05 review ----

cat("\n=== 24. review regressions ===\n")
rev_codes <- data.frame(condition = "dm", vocab_id = "ATC", code = "A10A",
                        min_prescriptions = 2, window_days = 365,
                        stringsAsFactors = FALSE)
rev_lmdb <- data.frame(pnr = c("a","a"), atc = "A10AB01", atc2 = "A10",
                       eksd = as.Date(c("2005-01-01","2005-02-01")), year = 2005,
                       stringsAsFactors = FALSE)

# side files must not be read back as conditions
sf <- file.path(tempdir(), "sidefiles"); unlink(sf, recursive = TRUE)
mb_extract_medication(rev_lmdb, codes = rev_codes, outdir = sf,
                   keep_events = TRUE, verbose = FALSE,
                           from = as.Date("1990-01-01"))
ok("a side file is written", file.exists(file.path(sf, "dm_all_events.rds")))
back <- mb_load_conditions(sf)
ok("mb_load_conditions ignores side files",
   nrow(back) == 1 && identical(unique(back$condition), "dm"))
ok("mb_merge_all ignores them too",
   nrow(mb_merge_all(medication_dir = sf, codes = rev_codes,
                     verbose = FALSE)) == 1)

# a worker field clashing with one the driver sets should not crash the run
clash <- mb_run_conditions("x", function(cond)
  list(data = data.frame(a = 1), condition = "CLASH", status = "CLASH",
       n_persons = 0, n_rows = 0), verbose = FALSE)
ok("worker cannot overwrite the driver's summary fields",
   clash$summary$condition == "x" && clash$summary$status == "done")

# duplicated person-condition rows are a warning, not a silent choice
dup <- data.frame(pnr = c("a","a"), condition = "c",
                  onset_date = as.Date(c("2005-01-01","2001-01-01")),
                  stringsAsFactors = FALSE)
ok("mb_to_wide warns on duplicate person rows",
   tryCatch({ mb_to_wide(dup); FALSE }, warning = function(w) TRUE))

# the package must stay ASCII-clean or R CMD check warns
src <- if (dir.exists("../R")) "../R" else "R"
non_ascii <- unlist(lapply(list.files(src, pattern = "\\.R$", full.names = TRUE),
  function(f) grep("[^\x01-\x7F]", readLines(f, warn = FALSE), value = TRUE)))
ok("no non-ASCII characters in R/", length(non_ascii) == 0)


# 25. mb_prevalence() ----

cat("\n=== 25. mb_prevalence() ===\n")

prev_codes <- data.frame(condition = "dm", vocab_id = "ATC", code = "A10A",
                         min_prescriptions = 2, window_days = 365,
                         stringsAsFactors = FALSE)
prev_lmdb <- data.frame(
  pnr  = c("a","a","a",  "b",         "c","c",         "d","d"),
  atc  = "A10AB01",
  atc2 = "A10",
  eksd = as.Date(c("2005-01-01","2005-02-01","2010-01-01",
                   "2005-01-01",
                   "2020-01-01","2020-02-01",
                   "2004-06-01","2005-03-01")),
  year = c(2005,2005,2010, 2005, 2020,2020, 2004,2005),
  stringsAsFactors = FALSE)

# a: qualifies by 2005-02-01 (2 within 365d), and has a later, irrelevant
#    dispensing in 2010
# b: only ever one dispensing - never qualifies
# c: qualifies, but only in 2020 - after as_of below
# d: qualifies over the full history (2004-06-01, 2005-03-01 are 273 days
#    apart), but the pair straddles the lookback boundary used below

pd <- file.path(tempdir(), "prevdisp"); unlink(pd, recursive = TRUE)
mb_extract_medication(prev_lmdb, codes = prev_codes, outdir = pd,
                   keep_events = TRUE, verbose = FALSE,
                           from = as.Date("1990-01-01"))
full <- readRDS(file.path(pd, "dm.rds"))
disp <- readRDS(file.path(pd, "dm_all_events.rds"))

as_of <- as.Date("2006-01-01")

## 25a. lookback = Inf reproduces mb_count_conditions() exactly ----
prev_inf <- mb_prevalence(disp, as_of = as_of, lookback = Inf,
                          min_prescriptions = 2, window_days = 365)

wide <- mb_to_wide(full, value = "onset_date")
counted <- mb_count_conditions(wide, as_of = as_of)

# Bug found while building this test: as.matrix() on a Date column formats it
# to a date STRING, so the as_of comparison in mb_count_conditions() was a
# string compared against a raw numeric day count - never true either way.
# Fixed in combine.R by converting each column explicitly instead.
ok("mb_count_conditions(as_of=) counts an onset that occurred before as_of",
   counted$n_conditions[counted$pnr == "a"] == 1)
ok("mb_count_conditions(as_of=) does not count an onset that occurred after as_of",
   counted$n_conditions[counted$pnr == "c"] == 0)

# one condition column, so n_conditions is 0/1 - "has dm as of as_of"
has_dm <- counted$pnr[counted$n_conditions >= 1]

ok("only people whose onset was by as_of are prevalent",
   setequal(prev_inf$pnr, has_dm))
ok("that is a and d, not b (never qualifies) or c (qualifies too late)",
   identical(sort(prev_inf$pnr), c("a", "d")))
ok("onset_date matches the full extraction, unaffected by the as_of cutoff",
   prev_inf$onset_date[prev_inf$pnr == "a"] ==
     full$onset_date[full$pnr == "a"])
ok("n_records is restricted to records at or before as_of, unlike the full run",
   prev_inf$n_records[prev_inf$pnr == "a"] == 2 &&
     full$n_records[full$pnr == "a"] == 3)

## 25b. a finite lookback excludes an onset that is too old ----
prev_short <- mb_prevalence(disp, as_of = as_of, lookback = 200,
                            min_prescriptions = 2, window_days = 365)
ok("a's dispensings (2005-01-01, 2005-02-01) are both outside a 200-day
    lookback from as_of (2006-01-01)",
   !"a" %in% prev_short$pnr)

## 25c. the window-boundary fix: a qualifying pair split across the naive
##      restriction boundary is still found, by widening the candidate set by
##      window_days and reporting the latest qualifying run instead of the
##      first (d's pair is 2004-06-01 and 2005-03-01, 273 days apart, and the
##      naive [as_of - lookback, as_of] restriction below would cut off
##      2004-06-01 alone)
prev_split <- mb_prevalence(disp, as_of = as.Date("2005-06-01"),
                            lookback = 300, id_col = pnr,
                            min_prescriptions = 2, window_days = 365)
ok("d's pair straddles the naive boundary but is still found",
   "d" %in% prev_split$pnr)
ok("the onset reported is the pair's later date, not the widened window's edge",
   prev_split$onset_date[prev_split$pnr == "d"] == as.Date("2005-03-01"))

## 25c2. the fix does not manufacture evidence that is not there: a pair that
##      completed genuinely before the true lookback boundary, even with the
##      widened candidate set, is still excluded
too_old <- mb_prevalence(disp, as_of = as.Date("2007-01-01"), lookback = 300,
                         min_prescriptions = 2, window_days = 365)
ok("d's pair (completing 2005-03-01) is really too old for this lookback",
   !"d" %in% too_old$pnr)

## 25c3. the fix only applies when there is a fixed span to widen by. With
##      window_days = NA (any distance apart) there is none, so the original
##      straddle limitation still applies - this is the one that remains in
##      ASSUMPTIONS_AND_LIMITATIONS.txt section 7
prev_na <- mb_prevalence(disp, as_of = as.Date("2005-06-01"), lookback = 300,
                         min_prescriptions = 2, window_days = NA)
ok("with window_days = NA there is no span to widen by, so d is still missed",
   !"d" %in% prev_na$pnr)

## 25d. min_prescriptions = 1 reduces to "any qualifying record in window",
##      the diagnosis-side use
diag_records <- data.frame(pnr = c("e","e","f"),
                           C_DIAG = c("DI500","DI509","DI500"),
                           D_INDDTO = as.Date(c("2004-01-01","2005-01-01",
                                                 "2020-01-01")),
                           stringsAsFactors = FALSE)
prev_diag <- mb_prevalence(diag_records, as_of = as_of, lookback = Inf,
                           code_col = C_DIAG, min_prescriptions = 1,
                           id_col = pnr, date_col = D_INDDTO)
ok("one qualifying diagnosis in the window is enough",
   setequal(prev_diag$pnr, "e"))

## 25e. as_of and lookback_days are carried onto the output for audit
ok("as_of is recorded on the result",
   all(prev_inf$as_of == as_of))
ok("lookback_days is NA when the lookback was Inf",
   all(is.na(prev_inf$lookback_days)))
ok("lookback_days is recorded when finite",
   all(prev_short$lookback_days == 200))

## 25f. bad arguments refuse rather than guess
ok("as_of must be length 1",
   inherits(tryCatch(mb_prevalence(disp, as_of = as.Date(c("2005-01-01",
                                                            "2006-01-01"))),
                     error = function(e) e), "error"))
ok("lookback must be positive",
   inherits(tryCatch(mb_prevalence(disp, as_of = as_of, lookback = 0),
                     error = function(e) e), "error"))


# 26. ids= cohort restriction (DECISIONS.md 7.5) ----

cat("\n=== 26. ids= restricts to a cohort before anything else runs ===\n")

ids_codes <- data.frame(condition = "dm", vocab_id = "ATC", code = "A10A",
                        min_prescriptions = 1, window_days = NA,
                        stringsAsFactors = FALSE)
ids_lmdb <- data.frame(pnr = c("a", "b", "c"), atc = "A10AB01",
                       atc2 = "A10",
                       eksd = as.Date("2010-01-01"), year = 2010,
                       stringsAsFactors = FALSE)

full_result <- mb_extract_medication(ids_lmdb, codes = ids_codes,
                                     outdir = NULL, verbose = FALSE,
                           from = as.Date("1990-01-01"))$dm
ok("ids = NULL (default) keeps everyone, unchanged",
   setequal(full_result$pnr, c("a", "b", "c")))

restricted <- mb_extract_medication(ids_lmdb, codes = ids_codes,
                                    outdir = NULL, verbose = FALSE,
                                    ids = c("a", "c"),
                           from = as.Date("1990-01-01"))$dm
ok("ids = restricts to exactly that cohort",
   setequal(restricted$pnr, c("a", "c")))
ok("a person not in ids never appears, even though they match the codes",
   !"b" %in% restricted$pnr)

ids_lpr <- data.frame(pnr = c("a", "b"), C_DIAG = "DI10",
                      D_INDDTO = as.Date("2010-01-01"),
                      stringsAsFactors = FALSE)
ids_dx_codes <- data.frame(condition = "hf", vocab_id = "ICD10", code = "I10",
                           exclude = FALSE, stringsAsFactors = FALSE)
restricted_dx <- mb_extract_diagnosis(ids_lpr, codes = ids_dx_codes,
                                      outdir = NULL, verbose = FALSE,
                                      ids = "a")$hf
ok("ids = works the same way on the diagnosis side",
   identical(restricted_dx$pnr, "a"))

## 26b. against a real lazy backend, not just a data frame - confirms the
##      restriction is an actual JOIN in the generated SQL, not a silent
##      fallback to collecting everyone and filtering in R afterward
if (requireNamespace("DBI", quietly = TRUE) &&
    requireNamespace("duckdb", quietly = TRUE) &&
    requireNamespace("dbplyr", quietly = TRUE)) {

  con <- DBI::dbConnect(duckdb::duckdb())
  DBI::dbWriteTable(con, "lmdb_lazy", ids_lmdb)
  lazy_lmdb <- dplyr::tbl(con, "lmdb_lazy")

  lazy_query <- mb_restrict_ids(lazy_lmdb, "pnr", c("a", "c"))
  sql_text <- toupper(as.character(dbplyr::remote_query(lazy_query)))
  # dbplyr translates semi_join() to a correlated "WHERE EXISTS" subquery
  # against a copied-in temp table, not the literal word JOIN - that is the
  # standard SQL form for a semi-join, and this is what proves the
  # restriction runs server-side rather than collecting everyone first.
  ok("mb_restrict_ids() pushes the restriction into the SQL (WHERE EXISTS)",
     grepl("WHERE", sql_text, fixed = TRUE) &&
       grepl("EXISTS", sql_text, fixed = TRUE))

  lazy_result <- mb_extract_medication(lazy_lmdb, codes = ids_codes,
                                       outdir = NULL, verbose = FALSE,
                                       ids = c("a", "c"),
                           from = as.Date("1990-01-01"))$dm
  ok("mb_extract_medication() against a lazy backend gives the same answer",
     setequal(lazy_result$pnr, c("a", "c")))

  DBI::dbDisconnect(con, shutdown = TRUE)
} else {
  cat("  SKIP   duckdb/dbplyr not installed - lazy-backend check skipped\n")
}


# 27. mb_prevalence_all() - every condition, its own lookback, one call ----

cat("\n=== 27. mb_prevalence_all() ===\n")

pa_codes <- data.frame(
  condition = c("dm", "htn"), vocab_id = "ATC", code = c("A10A", "C02"),
  min_prescriptions = 2, window_days = 365, stringsAsFactors = FALSE)

pa_lmdb <- data.frame(
  pnr  = c("a","a","a", "b",       "c","c",       "d","d",       "e","e"),
  atc  = c(rep("A10AB01", 3), "A10AB01", rep("C02CA01", 2),
           rep("A10AB01", 2), rep("C02CA01", 2)),
  atc2 = c(rep("A10", 4), rep("C02", 2), rep("A10", 2), rep("C02", 2)),
  eksd = as.Date(c("2005-01-01","2005-02-01","2010-01-01",
                   "2005-01-01",
                   "2020-01-01","2020-02-01",
                   "2004-06-01","2005-03-01",
                   "2005-01-01","2005-06-01")),
  year = c(2005,2005,2010, 2005, 2020,2020, 2004,2005, 2005,2005),
  stringsAsFactors = FALSE)
# e is htn's own qualifying person, onset well before as_of (2006-01-01) -
# c (also htn) qualifies too, but only in 2020, after as_of, so c alone
# would leave htn empty and silently absent from the combined result.

pa_dir <- file.path(tempdir(), "prevall"); unlink(pa_dir, recursive = TRUE)
mb_extract_medication(pa_lmdb, codes = pa_codes, outdir = pa_dir,
                      keep_events = TRUE, verbose = FALSE,
                           from = as.Date("1990-01-01"))

pa_as_of <- as.Date("2006-01-01")

## 27a. a single lookback value applies to every condition, and matches
##      calling mb_prevalence() by hand for the same condition
pa_single <- mb_prevalence_all(pa_dir, as_of = pa_as_of, lookback = Inf,
                               codes = pa_codes, verbose = FALSE)
disp_dm  <- readRDS(file.path(pa_dir, "dm_all_events.rds"))
solo_dm  <- mb_prevalence(disp_dm, as_of = pa_as_of, lookback = Inf,
                          min_prescriptions = 2, window_days = 365)
ok("mb_prevalence_all() with a single lookback matches a solo mb_prevalence() call",
   setequal(pa_single$pnr[pa_single$condition == "dm"], solo_dm$pnr))
ok("both conditions are present in the combined result",
   setequal(pa_single$condition, c("dm", "htn")))

## 27b. a per-condition lookback table - different conditions, different
##      values, exactly what B7/DECISIONS.md 4.1 asked for
pa_lb <- data.frame(condition = c("dm", "htn"),
                    lookback  = c("ever", "last_two_years"),
                    stringsAsFactors = FALSE)
pa_multi <- mb_prevalence_all(pa_dir, as_of = pa_as_of, lookback = pa_lb,
                              codes = pa_codes, verbose = FALSE)
ok("a lookback word ('ever') resolves to Inf",
   all(is.na(pa_multi$lookback_days[pa_multi$condition == "dm"])))
ok("a lookback word ('last_two_years') resolves to 730 days",
   all(pa_multi$lookback_days[pa_multi$condition == "htn"] == 730))

## 27c. a condition missing from the lookback table defaults to Inf, with a
##      message rather than an error
pa_lb_partial <- data.frame(condition = "dm", lookback = "ever",
                            stringsAsFactors = FALSE)
msg <- tryCatch({
  withCallingHandlers(
    mb_prevalence_all(pa_dir, as_of = pa_as_of, lookback = pa_lb_partial,
                      codes = pa_codes, verbose = FALSE),
    message = function(m) invokeRestart("muffleMessage")
  )
  "no error"
}, error = function(e) conditionMessage(e))
ok("a condition missing from the lookback table does not error",
   identical(msg, "no error"))

## 27d. conditions = restricts which files are read, and complains if one is
##      missing rather than silently skipping it
ok("an unknown condition name is refused, not silently dropped",
   inherits(tryCatch(mb_prevalence_all(pa_dir, as_of = pa_as_of,
                                       conditions = c("dm", "not_a_condition"),
                                       codes = pa_codes, verbose = FALSE),
                     error = function(e) e), "error"))

## 27e. an unrecognised lookback word refuses rather than silently becoming NA
ok("an unrecognised lookback word is refused",
   inherits(tryCatch(mb_lookback_days("last_year"), error = function(e) e),
            "error"))


# 28. mb_extract_medication_batch() ----

cat("\n=== 28. mb_extract_medication_batch() ===\n")

mb28_codes <- data.frame(
  condition = c("dm","dm", "htn","htn","htn", "pain","pain"),
  vocab_id  = "ATC",
  code      = c("A10A","A10B", "C02","C03","C09", "N02A","M01A"),
  min_prescriptions = c(2,2, 2,2,2, 4,4),
  window_days = 365,
  stringsAsFactors = FALSE)

mb28_lmdb <- data.frame(
  pnr  = c("a","a","a", "b",       "c","c",       "d","d",   "e","e","e","e"),
  atc  = c(rep("A10AB01", 3), "A10AB01", rep("C02CA01", 2),
          rep("A10AB01", 2), rep("N02AA01", 4)),
  eksd = as.Date(c("2005-01-01","2005-02-01","2010-01-01",
                   "2005-01-01",
                   "2020-01-01","2020-02-01",
                   "2004-06-01","2005-03-01",
                   "2018-01-01","2018-03-01","2018-05-01","2018-07-01")),
  year = c(2005,2005,2010, 2005, 2020,2020, 2004,2005, 2018,2018,2018,2018),
  stringsAsFactors = FALSE)
# a: qualifies for dm (2005-01-01, 2005-02-01 within 365d); 2010 dispensing
#    is irrelevant. b: only one dm dispensing - never qualifies. c: two htn
#    dispensings within 365d. d: dm pair straddling 2004/2005, 273 days
#    apart - qualifies. e: four pain dispensings within a year - qualifies
#    on the 4th, the only min_prescriptions = 4 rule.

## 28a. matches mb_extract_medication() exactly, on the same data
mb28_new <- mb_extract_medication_batch(mb28_lmdb, codes = mb28_codes,
                                        outdir = NULL, verbose = FALSE,
                           from = as.Date("1990-01-01"))
mb28_old <- mb_extract_medication(mb28_lmdb, codes = mb28_codes,
                                  outdir = NULL, verbose = FALSE,
                           from = as.Date("1990-01-01"))
for (cond in c("dm", "htn", "pain")) {
  a <- mb28_new[[cond]][order(mb28_new[[cond]]$pnr), ]
  b <- mb28_old[[cond]][order(mb28_old[[cond]]$pnr), ]
  ok(paste(cond, "- batch output columns match mb_extract_medication()'s"),
     identical(names(a), names(b)))
  ok(paste(cond, "- batch output agrees with mb_extract_medication()"),
     setequal(a$pnr, b$pnr) && isTRUE(all.equal(a$onset_date, b$onset_date)))
}

## 28b. exclusion codes (existing per-row exclude = TRUE mechanism)
mb28_codes_ex <- data.frame(
  condition = c("ex", "ex"), vocab_id = "ATC", code = c("A10", "A10B"),
  exclude = c(FALSE, TRUE), min_prescriptions = 2, window_days = 365,
  stringsAsFactors = FALSE)
mb28_lmdb_ex <- data.frame(
  pnr  = c("p1","p1", "p2","p2"),
  atc  = c("A10AB01","A10AB01", "A10BA01","A10BA01"),
  eksd = as.Date(c("2020-01-01","2020-02-01", "2020-01-01","2020-02-01")),
  year = 2020, stringsAsFactors = FALSE)
mb28_ex <- mb_extract_medication_batch(mb28_lmdb_ex, codes = mb28_codes_ex,
                                       outdir = NULL, verbose = FALSE,
                           from = as.Date("1990-01-01"))
ok("exclude = TRUE removes a code from within its own condition",
   "p1" %in% mb28_ex$ex$pnr && !("p2" %in% mb28_ex$ex$pnr))

## 28c. a condition matching nothing gets a correctly-shaped, empty result
mb28_codes_zero <- data.frame(
  condition = c("real", "nothing"), vocab_id = "ATC", code = c("A10", "Z99"),
  min_prescriptions = 2, window_days = 365, stringsAsFactors = FALSE)
mb28_zero <- mb_extract_medication_batch(mb28_lmdb, conditions = c("dm", "nothing"),
                                         codes = rbind(mb28_codes[mb28_codes$condition == "dm", ],
                                                      mb28_codes_zero["nothing" == mb28_codes_zero$condition, ]),
                                         outdir = NULL, verbose = FALSE,
                           from = as.Date("1990-01-01"))
ok("a condition matching nothing gets a 0-row result, not an error",
   nrow(mb28_zero$nothing) == 0)
ok("that empty result has the same 8 columns as a real one",
   identical(names(mb28_zero$nothing), names(mb28_zero$dm)))

## 28d. window_days = NA - any distance apart still qualifies
mb28_codes_na <- data.frame(condition = "any", vocab_id = "ATC", code = "A10",
                            min_prescriptions = 2, window_days = NA,
                            stringsAsFactors = FALSE)
mb28_lmdb_na <- data.frame(pnr = c("p1","p1"), atc = c("A10AB01","A10AB01"),
                          eksd = as.Date(c("2000-01-01","2020-01-01")),
                          year = c(2000,2020), stringsAsFactors = FALSE)
mb28_na <- mb_extract_medication_batch(mb28_lmdb_na, codes = mb28_codes_na,
                                       outdir = NULL,
                                       verbose = FALSE,
                           from = as.Date("1990-01-01"))
ok("window_days = NA - a 20-year gap still qualifies",
   "p1" %in% mb28_na$any$pnr)

## 28e. from is applied upstream, not as a post-hoc filter on onset -
##      the mechanism behind the bipolar anomaly found on real DST data,
##      2026-09-26 (TODO.txt section 5)
mb28_codes_ym <- data.frame(condition = "x", vocab_id = "ATC", code = "A10",
                           min_prescriptions = 2, window_days = 365,
                           stringsAsFactors = FALSE)
mb28_lmdb_ym <- data.frame(
  pnr = c("p1","p1"), atc = c("A10AB01","A10AB01"),
  eksd = as.Date(c("1996-12-01","1997-02-01")),
  stringsAsFactors = FALSE)
mb28_ym <- mb_extract_medication_batch(mb28_lmdb_ym, codes = mb28_codes_ym,
                                       outdir = NULL,
                                       from = as.Date("1997-01-01"),
                                       verbose = FALSE)
ok("from removes a pre-window anchor dispensing, not just a post-hoc onset filter",
   !("p1" %in% mb28_ym$x$pnr))

## 28f. resume is batch-granular: skip only when EVERY condition sharing a
##      rule already has its .rds; otherwise the whole batch reruns
mb28_dir <- file.path(tempdir(), "mb28resume"); unlink(mb28_dir, recursive = TRUE)
mb28_codes_r <- data.frame(condition = c("c1","c2"), vocab_id = "ATC",
                          code = c("A10","C02"), min_prescriptions = 2,
                          window_days = 365, stringsAsFactors = FALSE)
mb28_lmdb_r <- data.frame(pnr = "p1", atc = "A10AB01",
                         eksd = as.Date("2020-01-01"), year = 2020,
                         stringsAsFactors = FALSE)
s1 <- mb_extract_medication_batch(mb28_lmdb_r, codes = mb28_codes_r,
                                  outdir = mb28_dir, verbose = FALSE,
                           from = as.Date("1990-01-01"))
ok("first run: both conditions done", all(s1$status == "done"))

s2 <- mb_extract_medication_batch(mb28_lmdb_r, codes = mb28_codes_r,
                                  outdir = mb28_dir, resume = TRUE,
                                  verbose = FALSE,
                           from = as.Date("1990-01-01"))
ok("second run: both files exist, whole batch skipped",
   all(s2$status == "skipped"))

unlink(file.path(mb28_dir, "c1.rds"))
s3 <- mb_extract_medication_batch(mb28_lmdb_r, codes = mb28_codes_r,
                                  outdir = mb28_dir, resume = TRUE,
                                  verbose = FALSE,
                           from = as.Date("1990-01-01"))
ok("third run: one file missing, whole batch reruns (not per-condition)",
   all(s3$status == "done"))

## 28g. interop: mb_load_conditions()/mb_merge_all() work unchanged on a mix
##      of mb_extract_medication() and mb_extract_medication_batch() output
##      in the same outdir
mb28_dir2 <- file.path(tempdir(), "mb28interop"); unlink(mb28_dir2, recursive = TRUE)
mb28_codes_i <- data.frame(condition = c("old_fn","new_fn"), vocab_id = "ATC",
                          code = c("A10","C02"), min_prescriptions = 2,
                          window_days = 365, stringsAsFactors = FALSE)
mb28_lmdb_i <- data.frame(
  pnr = c("p1","p1","p2","p2"),
  atc = c("A10AB01","A10AB01","C02CA01","C02CA01"),
  eksd = as.Date(c("2020-01-01","2020-02-01","2020-01-01","2020-02-01")),
  year = 2020, stringsAsFactors = FALSE)
mb_extract_medication(mb28_lmdb_i, conditions = "old_fn", codes = mb28_codes_i,
                      outdir = mb28_dir2, verbose = FALSE,
                           from = as.Date("1990-01-01"))
mb_extract_medication_batch(mb28_lmdb_i, conditions = "new_fn", codes = mb28_codes_i,
                            outdir = mb28_dir2, verbose = FALSE,
                           from = as.Date("1990-01-01"))
mb28_combined <- mb_load_conditions(mb28_dir2)
ok("mb_load_conditions() reads output from both extraction functions together",
   setequal(unique(mb28_combined$condition), c("old_fn", "new_fn")))

## 28h. against a real lazy backend - confirms push-down (LAG/OVER/
##      PARTITION BY), not a silent fallback to collecting everyone first
if (requireNamespace("DBI", quietly = TRUE) &&
    requireNamespace("duckdb", quietly = TRUE) &&
    requireNamespace("dbplyr", quietly = TRUE)) {

  con <- DBI::dbConnect(duckdb::duckdb())
  DBI::dbWriteTable(con, "mb28_lazy", mb28_lmdb)
  mb28_lazy_tbl <- dplyr::tbl(con, "mb28_lazy")

  ok("mb_is_lazy() recognises a dbplyr tbl, not a plain data frame",
     mb_is_lazy(mb28_lazy_tbl) && !mb_is_lazy(mb28_lmdb))

  mb28_batch <- mb_rule_batches(mb_codelist(mb28_codes, vocab = "ATC"),
                                c("dm", "htn"))[["2|365"]]
  mb28_q <- mb_batch_query(mb28_lazy_tbl, mb_codelist(mb28_codes, vocab = "ATC"),
                           mb28_batch, "pnr", "atc", "eksd", NULL, NULL,
                           TRUE, NULL)
  mb28_sql <- toupper(as.character(dbplyr::remote_query(mb28_q)))
  ok("batch query pushes down: LAG/OVER/PARTITION BY all in the generated SQL",
     grepl("LAG(", mb28_sql, fixed = TRUE) &&
       grepl("OVER (", mb28_sql, fixed = TRUE) &&
       grepl("PARTITION BY", mb28_sql, fixed = TRUE))

  mb28_lazy_res <- mb_extract_medication_batch(mb28_lazy_tbl, codes = mb28_codes,
                                               outdir = NULL, verbose = FALSE,
                           from = as.Date("1990-01-01"))
  for (cond in c("dm", "htn", "pain")) {
    a <- mb28_lazy_res[[cond]][order(mb28_lazy_res[[cond]]$pnr), ]
    b <- mb28_new[[cond]][order(mb28_new[[cond]]$pnr), ]
    ok(paste(cond, "- lazy DuckDB backend agrees with the local data.frame path"),
       setequal(a$pnr, b$pnr) && isTRUE(all.equal(a$onset_date, b$onset_date)))
  }

  DBI::dbDisconnect(con, shutdown = TRUE)
} else {
  cat("  SKIP   duckdb/dbplyr not installed - lazy-backend check skipped\n")
}

cat("\nDone.\n")
