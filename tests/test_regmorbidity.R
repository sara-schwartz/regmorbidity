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
#   15. year_min default
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
  for (nm in c("MB_CODELIST_COLS", "mb_run_conditions")) {
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
ok("15 conditions", length(unique(codes$condition)) == 15)
ok("61 rows",       nrow(codes) == 61)
ok("pain needs 4",  unique(codes$min_prescriptions[codes$condition == "pain"]) == 4)
ok("window 365",    unique(codes$window_days[codes$condition == "hypertension"]) == 365)
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
ok("15 files written", length(paths) == 15)
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
  PNR  = c("p1","p1", "p2","p2", "p3","p3", "p4"),
  ATC  = c("C09AA05","C09CA01", "C09AA05","C09AA05", "C09AA05","C09CA01", "C09AA05"),
  atc2 = "C09",
  EKSD = as.Date(c("2010-01-10","2010-02-10", "2010-01-10","2013-01-10",
                   "2010-05-01","2010-05-01", "2010-01-10")),
  year = c(2010,2010, 2010,2013, 2010,2010, 2010),
  stringsAsFactors = FALSE
)
res <- mb_extract_medication(lmdb, "hypertension", codes = codes, verbose = FALSE)[["hypertension"]]
ok("only p1 qualifies", identical(res$PNR, "p1"))
ok("onset = 2nd dispensing", res$onset_date == as.Date("2010-02-10"))
ok("p3 same-day pair does not qualify", !"p3" %in% res$PNR)

res_nodedupe <- mb_extract_medication(lmdb, "hypertension", codes = codes,
                                   dedupe_same_day = FALSE, verbose = FALSE)[["hypertension"]]
ok("without dedupe, p3 would qualify", "p3" %in% res_nodedupe$PNR)

res_nowin <- mb_extract_medication(lmdb, "hypertension",
                                codes = transform(codes, window_days = NA_integer_),
                                verbose = FALSE)[["hypertension"]]
ok("with no window, p2 qualifies too", all(c("p1","p2") %in% res_nowin$PNR))
# 5. Two-stage matching ----

cat("\n=== 5. stage 2 separates codes sharing a stage-1 prefix ===\n")
# N02A = pain, N02C = migraine. Both sit under atc2 = "N02".
lmdb2 <- data.frame(
  PNR  = c(rep("a", 4), rep("b", 2)),
  ATC  = c("N02AA05","N02AA05","N02AA05","N02AA05", "N02CC01","N02CC01"),
  atc2 = "N02",
  EKSD = as.Date(c("2010-01-01","2010-02-01","2010-03-01","2010-04-01",
                   "2010-01-01","2010-02-01")),
  year = 2010, stringsAsFactors = FALSE
)
r <- mb_extract_medication(lmdb2, c("pain","migraine"), codes = codes, verbose = FALSE)
ok("pain finds only a", identical(r$pain$PNR, "a"))
ok("migraine finds only b", identical(r$migraine$PNR, "b"))
ok("pain onset = 4th (min_prescriptions = 4)",
   r$pain$onset_date == as.Date("2010-04-01"))
# 6. Anchoring ----

cat("\n=== 6. anchoring ===\n")
# An unanchored 'C09' would also match a code that merely contains it.
lmdb3 <- data.frame(PNR = c("x","x"), ATC = c("XC09AA","XC09AA"), atc2 = "C09",
                    EKSD = as.Date(c("2010-01-01","2010-02-01")), year = 2010,
                    stringsAsFactors = FALSE)
r3 <- mb_extract_medication(lmdb3, "hypertension", codes = codes, verbose = FALSE)[["hypertension"]]
ok("mid-string match rejected", nrow(r3) == 0)
# 7. Checkpoint and resume ----

cat("\n=== 7. checkpoint and resume ===\n")
od <- file.path(tempdir(), "out"); unlink(od, recursive = TRUE)
s1 <- mb_extract_medication(lmdb, c("hypertension","dyslipidemia"), codes = codes,
                         outdir = od, verbose = FALSE)
ok("two files written", length(list.files(od, pattern = "\\.rds$")) == 2)
ok("both marked done", all(s1$status == "done"))
s2 <- mb_extract_medication(lmdb, c("hypertension","dyslipidemia"), codes = codes,
                         outdir = od, verbose = FALSE)
ok("second run skips both", all(s2$status == "skipped"))
s3 <- mb_extract_medication(lmdb, c("hypertension","dyslipidemia"), codes = codes,
                         outdir = od, resume = FALSE, verbose = FALSE)
ok("resume = FALSE re-runs", all(s3$status == "done"))
# 8. Error handling ----

cat("\n=== 8. a failing condition does not kill the batch ===\n")
lmdb_bad <- rbind(lmdb,
                  data.frame(PNR = "p9", ATC = "C10AA01", atc2 = "C10",
                             EKSD = as.Date("2010-01-01"), year = 2010))
lmdb_bad$EKSD <- "not-a-date"   # both conditions now match, both fail to parse
sfail <- suppressWarnings(
  mb_extract_medication(lmdb_bad, c("hypertension","dyslipidemia"), codes = codes,
                     outdir = od, resume = FALSE, verbose = FALSE))
ok("errors recorded, loop continues", nrow(sfail) == 2 && all(sfail$status == "error"))

# A missing column is a configuration error and should stop before the loop.
ok("missing column fails fast",
   inherits(try(mb_extract_medication(lmdb, "hypertension", codes = codes,
                                   date_col = "NOPE", verbose = FALSE),
                silent = TRUE), "try-error"))
# 9. Combining output ----

cat("\n=== 9. load back, widen, count ===\n")
unlink(od, recursive = TRUE)
mb_extract_medication(lmdb2, c("pain","migraine"), codes = codes, outdir = od,
                   verbose = FALSE)
long <- mb_load_conditions(od)
ok("2 person-condition rows", nrow(long) == 2)
wide <- mb_to_wide(long, id_col = "PNR")
ok("one column per condition", all(c("pain","migraine") %in% names(wide)))
cnt <- mb_count_conditions(wide, id_col = "PNR")
ok("each person has 1 condition", all(cnt$n_conditions == 1))
ok("nobody multimorbid", !any(cnt$multimorbid))
cnt_early <- mb_count_conditions(wide, id_col = "PNR", as_of = "2010-01-15")
ok("as_of excludes later onsets", sum(cnt_early$n_conditions) == 0)
# 10. flag_users() ----

cat("\n=== 10. flag_users() keeps Jie's calling convention ===\n")
disp <- data.frame(PNR = c("q","q","q"), ATC = "C09AA05",
                   EKSD = as.Date(c("2010-01-01","2010-02-01","2010-03-01")),
                   stringsAsFactors = FALSE)
f1 <- mb_flag_users(disp, ATC, min_prescriptions = 2)
f2 <- mb_flag_users(disp, "ATC", min_prescriptions = 2)
ok("bare column name works", nrow(f1) == 1 && f1$onset_date == as.Date("2010-02-01"))
ok("string column name works", identical(f1, f2))
ok("n_records reported",f1$n_records == 3)
# 11. Missing prefilter column ----

cat("\n=== 11. missing prefilter column falls back ===\n")
lmdb4 <- lmdb[, c("PNR","ATC","EKSD","year")]
r4 <- mb_extract_medication(lmdb4, "hypertension", codes = codes, verbose = FALSE)[["hypertension"]]
ok("still finds p1 without atc2", identical(r4$PNR, "p1"))
# 12. Short codes ----

cat("\n=== 12. codes shorter than the prefilter column ===\n")
short <- rbind(codes[codes$condition == "hypertension", ][1, ])
short$code <- "C0"; short$condition <- "short_code"
lmdb5 <- data.frame(PNR = c("z","z"), ATC = c("C09AA05","C01DA02"), atc2 = c("C09","C01"),
                    EKSD = as.Date(c("2010-01-01","2010-02-01")), year = 2010,
                    stringsAsFactors = FALSE)
r5 <- mb_extract_medication(lmdb5, "short_code", codes = rbind(codes, short),
                         verbose = FALSE)[["short_code"]]
ok("2-char code matched via full ATC column", nrow(r5) == 1 && r5$PNR == "z")
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
  distress = "N06", bipolar = "N05", dementia = "N06"
)
bad <- character(0)
for (cond in names(jie_atc2)) {
  cc <- codes$code[codes$condition == cond & !codes$exclude]
  derived <- sort(unique(substr(cc, 1L, 3L)))
  expected <- sort(unique(strsplit(jie_atc2[[cond]], "|", fixed = TRUE)[[1]]))
  if (!identical(derived, expected)) bad <- c(bad, cond)
}
ok("all 15 conditions match Jie's atc2 patterns", length(bad) == 0)
if (length(bad)) cat("     mismatched:", paste(bad, collapse = ", "), "\n")
# 15. year_min default ----

cat("\n=== 15. year_min defaults to 1997 ===\n")
ok("default is 1997", formals(mb_extract_medication)$year_min == 1997)
old <- data.frame(PNR = c("o","o"), ATC = "C09AA05", atc2 = "C09",
                  EKSD = as.Date(c("1995-01-01","1995-06-01")), year = 1995,
                  stringsAsFactors = FALSE)
ok("1995 dropped by default",
   nrow(mb_extract_medication(old, "hypertension", codes = codes,
                           verbose = FALSE)[["hypertension"]]) == 0)
ok("year_min = NULL keeps 1995",
   nrow(mb_extract_medication(old, "hypertension", codes = codes, year_min = NULL,
                           verbose = FALSE)[["hypertension"]]) == 1)
ok("missing year column errors clearly",
   inherits(try(mb_extract_medication(old[, c("PNR","ATC","atc2","EKSD")],
                                   "hypertension", codes = codes,
                                   verbose = FALSE), silent = TRUE), "try-error"))

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
  PNR      = c("a", "a", "b", "c"),
  C_DIAG   = c("DI500", "DI509", "DD86", "DI10"),
  D_INDDTO = as.Date(c("2005-03-01", "2001-01-01", "2009-05-05", "2000-01-01")),
  stringsAsFactors = FALSE
)
dx <- mb_extract_diagnosis(lpr, codes = dx_codes, verbose = FALSE)
ok("heart failure finds a", identical(dx$hf$PNR, "a"))
ok("onset is the EARLIEST diagnosis", dx$hf$onset_date == as.Date("2001-01-01"))
ok("both records counted", dx$hf$n_records == 2)
ok("source labelled", dx$hf$source == "diagnosis")
ok("D86 written without the D prefix still matches DD86",
   identical(dx$connective$PNR, "b"))
ok("DI10 not picked up by either condition",
   !"c" %in% c(dx$hf$PNR, dx$connective$PNR))


# 18. Merging the two halves ----

cat("\n=== 18. OR / AND merge ===\n")
dxh <- data.frame(PNR = c("p","q"), condition = "x", source = "diagnosis",
                  onset_date = as.Date(c("2005-01-01","2007-01-01")),
                  stringsAsFactors = FALSE)
rxh <- data.frame(PNR = c("p","r"), condition = "x", source = "medication",
                  onset_date = as.Date(c("2003-01-01","2008-01-01")),
                  stringsAsFactors = FALSE)

m_or <- mb_merge_conditions(dxh, rxh, logic = "OR")
ok("OR keeps everyone", setequal(m_or$PNR, c("p","q","r")))
ok("OR takes the earlier date",
   m_or$onset_date[m_or$PNR == "p"] == as.Date("2003-01-01"))
ok("source = both where both halves present",
   m_or$source[m_or$PNR == "p"] == "both")
ok("source = diagnosis where only that half",
   m_or$source[m_or$PNR == "q"] == "diagnosis")

m_and <- mb_merge_conditions(dxh, rxh, logic = "AND")
ok("AND keeps only p", identical(m_and$PNR, "p"))
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
  vocab_id  = c("ATC","ICD10","ATC","ATC"),
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
edited$code[edited$condition == "dm" & edited$vocab_id == "ATC"] <- "A10B"
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
  PNR      = c("a","a","a","b"),
  C_DIAG   = c("DI500","DI509","DI50","DI50"),   # two codes, same day
  D_INDDTO = as.Date(c("2005-01-01","2005-01-01","2008-01-01","2001-01-01")),
  stringsAsFactors = FALSE)
hf <- data.frame(condition="hf", vocab_id="ICD10", code="I50",
                 stringsAsFactors = FALSE)

d_dedup <- mb_extract_diagnosis(lpr2, codes = hf, verbose = FALSE)$hf
ok("same-day records count once by default",
   d_dedup$n_records[d_dedup$PNR == "a"] == 2)
d_raw <- mb_extract_diagnosis(lpr2, codes = hf, dedupe_same_day = FALSE,
                           verbose = FALSE)$hf
ok("without dedup they count separately",
   d_raw$n_records[d_raw$PNR == "a"] == 3)
ok("onset unaffected by dedup",
   d_dedup$onset_date[d_dedup$PNR == "a"] == as.Date("2005-01-01"))

rxh2 <- data.frame(PNR = "a", condition = "hf", onset_date = as.Date("2004-01-01"),
                   n_records = 7L, stringsAsFactors = FALSE)
mg <- mb_merge_conditions(d_dedup, rxh2, logic = "OR")
ok("n_icd carried into the merge", mg$n_icd[mg$PNR == "a"] == 2)
ok("n_rx carried into the merge",  mg$n_rx[mg$PNR == "a"] == 7)
ok("n_rx is NA where there is no medication half",
   is.na(mg$n_rx[mg$PNR == "b"]))


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
    "mb_codelist", "mb_validate_codelist", "mb_write_codelist",
    "mb_inspect_codes", "mb_check_codes", "mb_overlap", "mb_lookup",
    "mb_compare", "mb_extract_medication", "mb_extract_diagnosis", "mb_flag_users",
    "mb_normalize_icd10", "mb_merge_conditions", "mb_merge_all",
    "mb_condition_logic", "mb_load_conditions", "mb_to_wide",
    "mb_count_conditions"))

  ok("exactly the intended 18 functions are public",
     identical(declared, expected))
  if (!identical(declared, expected)) {
    cat("     unexpectedly public:", paste(setdiff(declared, expected), collapse = ", "), "\n")
    cat("     expected but absent:", paste(setdiff(expected, declared), collapse = ", "), "\n")
  }

  # the per-condition workers and the driver are implementation, not API
  ok("internal helpers stay internal",
     !any(c("mb_extract_one", "mb_diagnose_one", "mb_run_conditions",
            "mb_onset", "mb_pattern", "mb_half") %in% declared))
}

cat("\nDone.\n")


# 24. Regressions found in the 2026-09-05 review ----

cat("\n=== 24. review regressions ===\n")
rev_codes <- data.frame(condition = "dm", vocab_id = "ATC", code = "A10A",
                        min_prescriptions = 2, window_days = 365,
                        stringsAsFactors = FALSE)
rev_lmdb <- data.frame(PNR = c("a","a"), ATC = "A10AB01", atc2 = "A10",
                       EKSD = as.Date(c("2005-01-01","2005-02-01")), year = 2005,
                       stringsAsFactors = FALSE)

# side files must not be read back as conditions
sf <- file.path(tempdir(), "sidefiles"); unlink(sf, recursive = TRUE)
mb_extract_medication(rev_lmdb, codes = rev_codes, outdir = sf,
                   save_dispensings = TRUE, verbose = FALSE)
ok("a side file is written", file.exists(file.path(sf, "dm_dispensings.rds")))
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
dup <- data.frame(PNR = c("a","a"), condition = "c",
                  onset_date = as.Date(c("2005-01-01","2001-01-01")),
                  stringsAsFactors = FALSE)
ok("mb_to_wide warns on duplicate person rows",
   tryCatch({ mb_to_wide(dup); FALSE }, warning = function(w) TRUE))

# the package must stay ASCII-clean or R CMD check warns
src <- if (dir.exists("../R")) "../R" else "R"
non_ascii <- unlist(lapply(list.files(src, pattern = "\\.R$", full.names = TRUE),
  function(f) grep("[^\x01-\x7F]", readLines(f, warn = FALSE), value = TRUE)))
ok("no non-ASCII characters in R/", length(non_ascii) == 0)
