# .............................................................................
# diagnosis.R - deriving conditions from hospital diagnoses
#
# PURPOSE
#   The ICD-10 counterpart to extract.R. One row per person per condition, with
#   the date of the first qualifying diagnosis.
#
# INPUT   lpr   - diagnosis records, lazy or in memory
#         codes - a code list from mb_codelist(), ICD10 rows
# OUTPUT  one .rds per condition in `outdir`, plus a summary data frame
#
# WHY THIS IS SIMPLER THAN THE MEDICATION SIDE
#   No two-stage filter: there is no atc2 equivalent for diagnoses, so the codes
#   go straight at the diagnosis column. No prescription rule either - one
#   recorded diagnosis is the condition. What replaces those in Prior is the
#   recovery window (Ever / last two years / last five years), which this
#   package does not implement; see ASSUMPTIONS_AND_LIMITATIONS section 7.
#
# THE LEADING D
#   Danish register ICD-10 codes carry a leading D: DE10, DI50, DD86. Published
#   code lists do not: E10, I50, D86. Both forms are accepted and normalised to
#   the WHO form before matching - see mb_normalize_icd10(), which is the one
#   piece of real subtlety in this file.
#
# CONTENTS
#   1. Normalising ICD-10 codes
#   2. Extracting many conditions
#        2.1 Filter to the condition
#        2.2 First diagnosis per person
# .............................................................................


# 1. Normalising ICD-10 codes ----

#' Strip the Danish D prefix, leaving a WHO ICD-10 code
#'
#' Danish register codes prepend D to the WHO code, so I50 is stored as DI50 and
#' D86 as DD86. A code list may be written either way - Jie's script uses DE10,
#' Prior's Web Table 1 uses E10 - so both sides are normalised before matching
#' rather than requiring the user to know which convention a file follows.
#'
#' The discrimination is exact, not a guess. A WHO code is a letter followed by
#' digits, so a Danish code is D followed by a LETTER. `DD86` is the Danish form
#' of `D86`; `D86` is already WHO form because its second character is a digit.
#' Matching on `^D[A-Z]` therefore never strips a D that belongs to the code.
#'
#' @param x Character vector of ICD-10 codes, in either convention.
#' @return The same codes in WHO form, upper-cased and stripped of whitespace.
#' @export
mb_normalize_icd10 <- function(x) {
  x <- toupper(gsub("[[:space:].]", "", as.character(x)))
  danish <- grepl("^D[A-Z]", x)
  x[danish] <- substring(x[danish], 2L)
  x
}


# 2. Extracting many conditions ----

#' Extract diagnosis-based conditions from a hospital register
#'
#' Runs one condition at a time and writes each result to its own `.rds`, so an
#' interrupted run resumes where it stopped. Same contract as
#' [mb_extract_medication()], so the two halves can be combined with
#' [mb_merge_conditions()].
#'
#' @param lpr The diagnosis register. A lazy table (arrow / duckplyr / dbplyr)
#'   or a plain data frame.
#' @param conditions Conditions to extract. `NULL` = every ICD-10 condition in
#'   `codes`.
#' @param codes A code list; see [mb_codelist()].
#' @param outdir Directory for one `.rds` per condition. `NULL` = return the
#'   results as a named list instead of writing anything.
#' @param resume Skip conditions whose `.rds` is already in `outdir`.
#' @param id_col,code_col,date_col Column names. `code_col` and `date_col`
#'   default to `C_DIAG`/`D_INDDTO`, the LPR2 names used in `streamlineMM.R`
#'   and DST's own field names, kept as-is. `id_col` defaults to lowercase
#'   `pnr`, matching this package's own convention rather than her renamed
#'   `PNR` - see DECISIONS.md.
#' @param year_col,year_min Optional year restriction. Unlike the medication
#'   side this defaults to `NULL`: Prior applies no year floor to diagnoses, and
#'   an earlier diagnosis is exactly what a lookback is for.
#' @param dedupe_same_day Count several diagnosis records on one date as one.
#'   On by default so `n_records` means the same thing on both halves - days
#'   with a record - rather than rows, which would make the two counts
#'   incomparable. One contact coding both I500 and I509 is one day either way.
#' @param keep_records Also write `<condition>_records.rds`, the matched
#'   diagnoses before the first-per-person step. Large.
#' @param ids Restrict to these person ids before anything else runs, via a
#'   join rather than a literal list - see [mb_restrict_ids()]. `NULL`
#'   (default) keeps everyone.
#' @param verbose Print progress.
#' @return With `outdir`: a summary data frame, invisibly. Without: a named list
#'   of per-condition data frames, with the summary attached as an attribute.
#'   **Always check `summary$status`.**
#' @export
mb_extract_diagnosis <- function(lpr,
                              conditions   = NULL,
                              codes        = mb_codelist(),
                              outdir       = NULL,
                              resume       = TRUE,
                              id_col       = "pnr",
                              code_col     = "C_DIAG",
                              date_col     = "D_INDDTO",
                              year_col     = "year",
                              year_min     = NULL,
                              dedupe_same_day = TRUE,
                              keep_records = FALSE,
                              ids          = NULL,
                              verbose      = TRUE) {

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.", call. = FALSE)
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.", call. = FALSE)
  }

  codes <- mb_codelist(codes, vocab = "ICD10", validate = TRUE)

  if (is.null(conditions)) {
    conditions <- unique(codes$condition)
  } else {
    unknown <- setdiff(conditions, unique(codes$condition))
    if (length(unknown)) {
      stop("Not in the code list: ", paste(unknown, collapse = ", "),
           "\nAvailable: ",
           paste(sort(unique(codes$condition)), collapse = ", "), call. = FALSE)
    }
  }

  # Setup mistakes fail identically for every condition, so catch them once.
  available_cols <- mb_cols(lpr)
  if (!is.null(available_cols)) {
    for (nm in c(id_col, code_col, date_col)) {
      if (!nm %in% available_cols) {
        stop("Column '", nm, "' is not in the data. Present: ",
             paste(utils::head(available_cols, 30), collapse = ", "),
             call. = FALSE)
      }
    }
    if (!is.null(year_min) && !year_col %in% available_cols) {
      stop("year_min = ", year_min, " but there is no '", year_col,
           "' column.", call. = FALSE)
    }
  }

  if (!is.null(outdir) && !dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  run <- mb_run_conditions(
    conditions = conditions,
    outdir     = outdir,
    resume     = resume,
    verbose    = verbose,
    unit       = "diagnoses",
    worker     = function(condition) {
      mb_diagnose_one(
        lpr = lpr, condition = condition, codes = codes,
        id_col = id_col, code_col = code_col, date_col = date_col,
        year_col = year_col, year_min = year_min,
        dedupe_same_day = dedupe_same_day, keep_records = keep_records,
        ids = ids
      )
    }
  )

  mb_finish_run(run, outdir, verbose)
}


#' Extract one diagnosis-based condition
#' @keywords internal
mb_diagnose_one <- function(lpr, condition, codes,
                            id_col, code_col, date_col,
                            year_col, year_min, dedupe_same_day = TRUE,
                            keep_records = FALSE, ids = NULL) {

  rules         <- codes[codes$condition == condition, , drop = FALSE]
  include_codes <- mb_normalize_icd10(unique(rules$code[!rules$exclude]))
  exclude_codes <- mb_normalize_icd10(unique(rules$code[rules$exclude]))

  if (!length(include_codes)) {
    stop("no include codes for condition '", condition, "'", call. = FALSE)
  }
  pattern <- mb_pattern(include_codes)

  ## 2.1 Filter to the condition ----
  # The pushed-down filter has to run on the raw column, which still carries the
  # Danish D. Matching "D?" plus the WHO pattern lets the database do the work
  # without a normalising pass over the whole register first.
  query <- mb_restrict_ids(lpr, id_col, ids)
  if (!is.null(year_min)) {
    query <- dplyr::filter(query, .data[[year_col]] >= !!year_min)
  }
  raw_pattern <- paste0("^D?(", substring(pattern, 3L))
  query <- dplyr::filter(query, grepl(!!raw_pattern, .data[[code_col]]))
  query <- dplyr::select(query, dplyr::all_of(c(id_col, code_col, date_col)))

  records <- dplyr::collect(query)

  if (!nrow(records)) {
    return(mb_null_diagnosis(id_col, condition, length(include_codes), pattern))
  }

  # Exact match, now that both sides are in WHO form. The pushed-down filter
  # above is deliberately loose - this is where the definition is applied.
  normalized <- mb_normalize_icd10(records[[code_col]])
  records <- records[grepl(pattern, normalized), , drop = FALSE]

  if (length(exclude_codes) && nrow(records)) {
    normalized <- mb_normalize_icd10(records[[code_col]])
    records <- records[!grepl(mb_pattern(exclude_codes), normalized), ,
                       drop = FALSE]
  }

  n_rows <- nrow(records)
  if (!n_rows) {
    return(mb_null_diagnosis(id_col, condition, length(include_codes), pattern))
  }

  ## 2.2 First diagnosis per person ----
  dt <- data.table::as.data.table(records[, c(id_col, code_col, date_col)])
  data.table::setnames(dt, c(id_col, code_col, date_col),
                       c("id", "code", "date"))
  dt$code <- mb_normalize_icd10(dt$code)
  dt$date <- mb_as_date(dt$date, date_col)
  dt <- dt[!is.na(dt$id) & !is.na(dt$date) & nzchar(dt$code), ]

  if (!nrow(dt)) {
    return(mb_null_diagnosis(id_col, condition, length(include_codes), pattern))
  }

  # Sorted before de-duplicating so which code is reported as onset_code is
  # deterministic, rather than whatever order the backend returned.
  data.table::setorderv(dt, c("id", "date", "code"))
  if (dedupe_same_day) dt <- unique(dt, by = c("id", "date"))

  people <- dt[, list(onset_date = date[1L],
                      onset_code = code[1L],
                      n_records  = .N,
                      first_date = date[1L],
                      last_date  = date[.N]), by = "id"]

  data.table::setnames(people, "id", id_col)
  people <- as.data.frame(people, stringsAsFactors = FALSE)
  people$condition <- condition
  people$source    <- "diagnosis"
  people <- people[, c(id_col, "condition", "source", "onset_date",
                       "onset_code", "n_records", "first_date", "last_date")]

  list(data          = people,
       extra         = list("_records" = if (keep_records) records else NULL),
       n_codes       = length(include_codes),
       pattern       = pattern,
       n_rows        = n_rows,
       n_persons_any = nrow(people),
       n_persons     = nrow(people))
}


#' @keywords internal
mb_null_diagnosis <- function(id_col, condition, n_codes, pattern) {
  out <- data.frame(
    id         = character(0),
    condition  = character(0),
    source     = character(0),
    onset_date = as.Date(character(0)),
    onset_code = character(0),
    n_records  = integer(0),
    first_date = as.Date(character(0)),
    last_date  = as.Date(character(0)),
    stringsAsFactors = FALSE
  )
  names(out)[1] <- id_col

  list(data = out, extra = list(), n_codes = n_codes, pattern = pattern,
       n_rows = 0L, n_persons_any = 0L, n_persons = 0L)
}
