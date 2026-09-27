# .............................................................................
# extract.R - deriving medication conditions from the dispensing register
#
# PURPOSE
#   Turn a code list plus LMDB into one row per person per condition, with the
#   date the person met the condition's prescription rule.
#
# INPUT   lmdb  - dispensing records, lazy (arrow/duckplyr/dbplyr) or in memory
#         codes - a code list from mb_codelist()
# OUTPUT  one .rds per condition in `outdir`, plus a summary data frame
#
# WHY ONE CONDITION AT A TIME
#   A single pass matching every condition at once does not finish on the full
#   register. Each condition is written to disk the moment it is done, so an
#   interrupted run resumes instead of starting over, and a condition that fails
#   does not discard the fifteen that already succeeded.
#
# WHY THE FILTER RUNS TWICE
#   The two-stage filter is Jie Zhang's, from dgt_medication.R, and it is the
#   part of this package that is known to run on DST:
#
#     lmdb |> filter(year >= 1997 & grepl('C01', atc2)) |> select(...) |> collect()
#     ... |> filter(grepl('C01DA', ATC))
#
#   Stage 1 is a cheap grepl on the 3-character level-2 column. It pushes down
#   into parquet/duckdb, so the register is narrowed before anything enters RAM.
#   Stage 2 applies the real codes to the full ATC column after collect().
#
#   Stage 1 on its own is not a definition and must never be used as one: atc2
#   holds 3 characters, so a 4-character code matches nothing there and returns
#   an empty result with no error. C01 is every cardiac drug; C01DA is nitrates.
#
# TWO DELIBERATE CHANGES TO dgt_medication.R
#   1. Patterns are anchored with "^", so a code can only match at the start of
#      an ATC code. This is Jie's own idiom - see the commented-out line 60,
#      ## atc_antihypertens <- "^C03|^C02CA|^G04CA|^C07|^C08|^C09"
#   2. Codes come from the code list rather than being written into each call,
#      so a condition is named in exactly one place. Both the duplicated
#      flag_users() bugs at lines 175 and 196 of that script were possible only
#      because the condition was named twice.
#
# SEE ALSO
#   ASSUMPTIONS_AND_LIMITATIONS.txt - what the rules mean and where they came
#   from. TODO.txt section 0 - the three rules still to be agreed with Jie.
#
# CONTENTS
#   1. Constants
#   2. Flagging users of one medication
#   3. Extracting many conditions
#        3.1 Fail on setup mistakes before the loop
#        3.2 One condition at a time
#        3.3 Stage 1: cheap filter, pushed down into the query
#        3.4 Stage 2: the real codes, against the full ATC code
#        3.5 The prescription rule
#   4. Small helpers
# .............................................................................
# 1. Constants ----

# The defaults for year_min, prefilter_col and prefilter_len are written as
# literals in the signature rather than as named constants. A named constant
# would be invisible: it is not exported, so args(mb_extract_medication) and the
# help page would show a symbol the reader cannot resolve. For a default
# argument the value IS the documentation, and the reasoning lives in @param.
#
#   year_min = 1997        what dgt_medication.R uses throughout. LMDB starts
#                          in 1995, so this drops two years on purpose. The
#                          reason for 1997 is not documented anywhere in that
#                          script - see TODO.txt section 6.
#   prefilter_col = "atc2" ATC level 2, exactly 3 characters (C09, A10, N06).
#   prefilter_len = 3      verify with mb_inspect_codes() on a new extract.
# 2. Flagging users of one medication ----

#' Flag people who meet a medication criterion
#'
#' Takes dispensing records already filtered to one condition and returns one
#' row per person who reached `min_prescriptions` dispensings, with the date of
#' the qualifying dispensing as `onset_date` - the 2nd, or the 4th for pain, not
#' the first.
#'
#' Replaces the `flag_users()` used in `dgt_medication.R` and takes the same
#' first two arguments, so `flag_users(aht_filtered, ATC, min_prescriptions = 2)`
#' still works unchanged. It adds `window_days` and same-day de-duplication,
#' both of which come from the Stata original.
#'
#' @param data Dispensing records for one condition, already collected.
#' @param code_col Column holding the ATC code. Bare name or string.
#' @param min_prescriptions Dispensings needed before the condition counts.
#' @param window_days They must fall within this many days of each other.
#'   `NA` = no window, i.e. any distance apart, which is what
#'   `dgt_medication.R` does.
#' @param id_col,date_col Person id and dispensing date columns. Bare or string.
#' @param dedupe_same_day Count several packages picked up on one day as one
#'   dispensing. The Stata original does this
#'   (`bys pid eksd: keep if _n == 1`). Leaving it off lets a single pharmacy
#'   visit satisfy a two-prescription rule on its own.
#' @param keep_all Also return people who matched the codes but never reached
#'   `min_prescriptions` (`onset_date` is `NA`). Use it to count how many people
#'   the rule excluded.
#' @return One row per person: id, `onset_date`, `onset_code`,
#'   `n_records`, `first_date`, `last_date`.
#' @export
mb_flag_users <- function(data,
                       code_col          = atc,
                       min_prescriptions = 2L,
                       window_days       = NA,
                       id_col            = pnr,
                       date_col          = eksd,
                       dedupe_same_day   = TRUE,
                       keep_all          = FALSE) {

  # ensym() accepts both a bare column name and a string, so an explicit call
  # naming the actual column - e.g. flag_users(aht_filtered, ATC, ...) from
  # dgt_medication.R, whatever case that column happens to be - still works
  # unchanged; only the DEFAULT below is lowercase.
  code_col <- rlang::as_name(rlang::ensym(code_col))
  id_col   <- rlang::as_name(rlang::ensym(id_col))
  date_col <- rlang::as_name(rlang::ensym(date_col))

  mb_flag_users_impl(data, code_col, min_prescriptions, window_days,
                id_col, date_col, dedupe_same_day, keep_all)
}


#' The implementation behind [mb_flag_users()], taking column names as strings
#'
#' Separate because quasiquotation does not survive a plain function call: an
#' internal caller already holds the column names as strings, and passing them
#' through the `ensym()` path would capture the variable name rather than what
#' it contains. `mb_flag_users()` is the thin wrapper that does the capturing;
#' everything the function actually does happens here.
#' @keywords internal
mb_flag_users_impl <- function(data, code_col, min_prescriptions, window_days,
                          id_col, date_col, dedupe_same_day, keep_all,
                          latest = FALSE) {

  for (nm in c(id_col, code_col, date_col)) {
    if (!nm %in% names(data)) {
      stop("Column '", nm, "' not found. Present: ",
           paste(names(data), collapse = ", "), call. = FALSE)
    }
  }

  dispensings <- data.table::as.data.table(data[, c(id_col, code_col, date_col)])
  data.table::setnames(dispensings, c(id_col, code_col, date_col),
                       c("id", "code", "date"))

  dispensings$code <- toupper(trimws(as.character(dispensings$code)))
  dispensings$date <- mb_as_date(dispensings$date, date_col)

  # A row with no date cannot take part in a window, and a row with no id cannot
  # be attributed to anyone. Dropping them here keeps them out of the counts.
  dispensings <- dispensings[!is.na(dispensings$id) &
                             !is.na(dispensings$date) &
                             nzchar(dispensings$code), ]

  if (dedupe_same_day) {
    # De-duplication is per person-DATE, not per drug: two different
    # antihypertensives collected on one visit are one dispensing, which is what
    # stops a single pharmacy visit satisfying a two-prescription rule. Same as
    # the Stata (bys pid eksd: keep if _n == 1).
    #
    # Sorted first because unique() keeps whichever row comes first, and without
    # a sort that is whatever order the backend returned - so onset_code could
    # differ between two runs of the same query.
    data.table::setorderv(dispensings, c("id", "date", "code"))
    dispensings <- unique(dispensings, by = c("id", "date"))
  }

  users <- mb_onset(dispensings,
                    min_prescriptions = min_prescriptions,
                    window_days       = window_days,
                    keep_all          = keep_all,
                    latest            = latest)

  data.table::setnames(users, "id", id_col)
  as.data.frame(users, stringsAsFactors = FALSE)
}


#' Onset date under the "n dispensings within a window" rule
#'
#' Onset is the date of the n-th dispensing in the first run of n that falls
#' inside `window_days`. With `n = 1` this is the first dispensing; with
#' `window_days = NA` the n dispensings may be any distance apart.
#'
#' This reproduces the Stata original, written differently.
#' `multsygdom_udtrak.do` lines 141-148 give every dispensing 365.25 days of
#' coverage (`expand 2` / `+365.25`) and require at least two coverages to
#' overlap (`k < -1`), four for pain (`k < -1 - 2*("l_pain")`). Two intervals of
#' equal length overlap exactly when their start dates are less than that length
#' apart, so "2 whose intervals overlap" and "2 within 365 days" are the same
#' rule.
#'
#' @param latest Report the LAST qualifying run instead of the first. Used by
#'   `mb_prevalence()` only: extraction always wants the true, first-ever
#'   onset, but a windowed prevalence check wants the most recent evidence,
#'   not whichever qualifying run happens to be earliest in a widened
#'   candidate set. Does not change what counts as qualifying, only which
#'   qualifying run is reported when there is more than one.
#' @keywords internal
mb_onset <- function(dt, min_prescriptions, window_days, keep_all = FALSE,
                     latest = FALSE) {

  n <- as.integer(min_prescriptions)
  if (is.na(n) || n < 1L) n <- 1L

  data.table::setorderv(dt, c("id", "date"))

  if (n > 1L) {
    # Sorted by (id, date), so the n-th dispensing of a run ending at row i has
    # its first at row i-n+1. Comparing the two dates is the whole rule.
    lag_date <- data.table::shift(dt$date, n - 1L, type = "lag")
    lag_id   <- data.table::shift(dt$id,   n - 1L, type = "lag")

    # Without the id check, the lag would reach back into the previous person's
    # rows and credit them to this one.
    same_person <- !is.na(lag_id) & lag_id == dt$id
    gap_days    <- as.numeric(dt$date - lag_date)

    qualifies <- same_person & !is.na(gap_days) &
                 (is.na(window_days) | gap_days <= as.numeric(window_days))
  } else {
    qualifies <- rep(TRUE, nrow(dt))
  }
  dt$mb_ok <- qualifies

  # Counted over every matching dispensing, not only the qualifying ones, so
  # n_records answers "how much of this drug did they collect" rather than
  # "how many rows passed the rule".
  counts <- dt[, list(n_records = .N,
                      first_date      = date[1L],
                      last_date       = date[.N]), by = "id"]

  qualifying <- dt[dt$mb_ok, ]
  if (nrow(qualifying)) {
    # Rows are already in date order, so the first qualifying row per person is
    # the earliest date the rule was met, and the last is the most recent.
    # length(x), not .N: .N is only meaningful written directly in a
    # data.table j-expression, not inside a helper function called from one.
    pick <- if (latest) function(x) x[length(x)] else function(x) x[1L]
    onset <- qualifying[, list(onset_date = pick(date),
                               onset_code = pick(code)), by = "id"]
    users <- merge(counts, onset, by = "id", all.x = TRUE)
  } else {
    users <- counts
    users$onset_date <- as.Date(NA)
    users$onset_code <- NA_character_
  }

  if (!keep_all) users <- users[!is.na(users$onset_date), ]

  users[, c("id", "onset_date", "onset_code", "n_records",
            "first_date", "last_date")]
}
# 3. Extracting many conditions ----

#' Extract medication-based conditions from a dispensing register
#'
#' Runs one condition at a time and writes each result to its own `.rds` as soon
#' as it is done, so a run that is interrupted picks up where it stopped.
#'
#' @param lmdb The dispensing register. A lazy table (arrow / duckplyr / dbplyr)
#'   or a plain data frame. Needs a person id, a full ATC code and a date column.
#' @param conditions Conditions to extract. `NULL` = every ATC condition in
#'   `codes`.
#' @param codes A code list; see [mb_codelist()].
#' @param outdir Directory for one `.rds` per condition. `NULL` = return the
#'   results as a named list instead of writing anything.
#' @param resume Skip conditions whose `.rds` is already in `outdir`.
#' @param id_col,code_col,date_col Column names. `code_col` must hold the
#'   *full* ATC code, not a truncated level.
#' @param prefilter_col Short ATC level column for the pushed-down first stage,
#'   or `NULL` to filter on `code_col` directly.
#' @param prefilter_len Characters held by `prefilter_col`. Check with
#'   [mb_inspect_codes()] before trusting it.
#' @param year_col,year_min Restrict to dispensings from this year onward.
#'   Defaults to 1997, as `dgt_medication.R` uses throughout, which deliberately
#'   drops 1995-96. Set `year_min = NULL` to keep every year, or if the data has
#'   no year column.
#' @param dedupe_same_day,keep_all_users Passed to [mb_flag_users()].
#' @param save_dispensings Also write `<condition>_dispensings.rds`, the matched
#'   records before the prescription rule is applied - the equivalent of the
#'   `*_filtered.rds` files in `dgt_medication.R`. Large.
#' @param ids Restrict to these person ids before anything else runs, via a
#'   join rather than a literal list - see [mb_restrict_ids()]. `NULL`
#'   (default) keeps everyone. Use this when the run is for a defined study
#'   cohort rather than the whole register; it is a second, independent lever
#'   on memory use alongside `prefilter_col`.
#' @param verbose Print progress.
#' @return With `outdir`: a summary data frame, invisibly. Without: a named list
#'   of per-condition data frames, with the summary attached as an attribute.
#'   **Always check `summary$status`** - a batch that finished may still contain
#'   errors.
#' @export
mb_extract_medication <- function(lmdb,
                               conditions       = NULL,
                               codes            = mb_codelist(),
                               outdir           = NULL,
                               resume           = TRUE,
                               id_col           = "pnr",
                               code_col         = "atc",
                               date_col         = "eksd",
                               prefilter_col    = "atc2",
                               prefilter_len    = 3L,
                               year_col         = "year",
                               year_min         = 1997,
                               dedupe_same_day  = TRUE,
                               keep_all_users   = FALSE,
                               save_dispensings = FALSE,
                               ids              = NULL,
                               verbose          = TRUE) {

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.", call. = FALSE)
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.", call. = FALSE)
  }

  codes <- mb_codelist(codes, vocab = "ATC", validate = TRUE)

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
  ## 3.1 Fail on setup mistakes before the loop ----
  # A wrong column name would fail identically for every condition. Catching it
  # here turns fifteen warnings into one message naming the columns that do
  # exist. Only a message improvement: the run fails either way.
  #
  # Skipped when the backend will not report its columns, so a backend this does
  # not understand can never block a run that would have worked.
  available_cols <- mb_cols(lmdb)
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
           "' column. Pass year_col = <name>, or year_min = NULL to keep ",
           "every year.", call. = FALSE)
    }
  }

  if (!is.null(outdir) && !dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  # A missing prefilter column is not an error: matching the full code column
  # directly gives the same answer, only slower.
  if (!is.null(prefilter_col) && !mb_has_col(lmdb, prefilter_col)) {
    if (verbose) {
      message("Column '", prefilter_col, "' not found - matching on '",
              code_col, "' directly.")
    }
    prefilter_col <- NULL
  }
  ## 3.2 One condition at a time ----
  run <- mb_run_conditions(
    conditions = conditions,
    outdir     = outdir,
    resume     = resume,
    verbose    = verbose,
    unit       = "dispensings",
    worker     = function(condition) {
      mb_extract_one(
        lmdb = lmdb, condition = condition, codes = codes,
        id_col = id_col, code_col = code_col, date_col = date_col,
        prefilter_col = prefilter_col, prefilter_len = prefilter_len,
        year_col = year_col, year_min = year_min,
        dedupe_same_day = dedupe_same_day, keep_all_users = keep_all_users,
        keep_dispensings = save_dispensings, ids = ids
      )
    }
  )

  mb_finish_run(run, outdir, verbose)
}


#' Extract one condition: two-stage filter, then the prescription rule
#' @keywords internal
mb_extract_one <- function(lmdb, condition, codes,
                           id_col, code_col, date_col,
                           prefilter_col, prefilter_len,
                           year_col, year_min,
                           dedupe_same_day, keep_all_users,
                           keep_dispensings = FALSE, ids = NULL) {

  rules         <- codes[codes$condition == condition, , drop = FALSE]
  include_codes <- unique(rules$code[!rules$exclude])
  exclude_codes <- unique(rules$code[rules$exclude])

  if (!length(include_codes)) {
    stop("no include codes for condition '", condition, "'", call. = FALSE)
  }

  # Validated as constant within a condition by mb_validate_codelist(), so the
  # first value is the value.
  min_prescriptions <- unique(rules$min_prescriptions)[1]
  window_days       <- unique(rules$window_days)[1]
  ## 3.3 Stage 1: cheap filter, pushed down into the query ----
  # Only usable when every code is at least as long as the short column: a
  # 2-character code cannot be expressed as a set of 3-character prefixes.
  use_prefilter <- !is.null(prefilter_col) &&
                   all(nchar(include_codes) >= prefilter_len)

  query <- mb_restrict_ids(lmdb, id_col, ids)

  if (!is.null(year_min)) {
    if (!mb_has_col(lmdb, year_col)) {
      stop("year_min = ", year_min, " but there is no '", year_col,
           "' column. Pass year_col = <name>, or year_min = NULL to keep ",
           "every year.", call. = FALSE)
    }
    query <- dplyr::filter(query, .data[[year_col]] >= !!year_min)
  }

  if (use_prefilter) {
    stage1_pattern <- mb_pattern(substr(include_codes, 1L, prefilter_len))
    query <- dplyr::filter(query, grepl(!!stage1_pattern,
                                        .data[[prefilter_col]]))
  } else {
    stage1_pattern <- mb_pattern(include_codes)
    query <- dplyr::filter(query, grepl(!!stage1_pattern, .data[[code_col]]))
  }

  # Only three columns are collected. Everything else stays in the file.
  query <- dplyr::select(query, dplyr::all_of(c(id_col, code_col, date_col)))
  dispensings <- dplyr::collect(query)

  stage2_pattern <- mb_pattern(include_codes)

  if (!nrow(dispensings)) {
    return(mb_null_result(id_col, condition, length(include_codes),
                          stage2_pattern, min_prescriptions, window_days))
  }
  ## 3.4 Stage 2: the real codes, against the full ATC code ----
  full_codes <- toupper(trimws(as.character(dispensings[[code_col]])))
  dispensings <- dispensings[grepl(stage2_pattern, full_codes), , drop = FALSE]

  if (length(exclude_codes)) {
    full_codes <- toupper(trimws(as.character(dispensings[[code_col]])))
    dispensings <- dispensings[!grepl(mb_pattern(exclude_codes), full_codes), ,
                               drop = FALSE]
  }

  n_rows <- nrow(dispensings)
  if (!n_rows) {
    return(mb_null_result(id_col, condition, length(include_codes),
                          stage2_pattern, min_prescriptions, window_days))
  }
  ## 3.5 The prescription rule ----
  users <- mb_flag_users_impl(dispensings,
                         code_col          = code_col,
                         min_prescriptions = min_prescriptions,
                         window_days       = window_days,
                         id_col            = id_col,
                         date_col          = date_col,
                         dedupe_same_day   = dedupe_same_day,
                         keep_all          = keep_all_users)

  users$condition <- condition
  users$source    <- "medication"
  users <- users[, c(id_col, "condition", "source", "onset_date", "onset_code",
                     "n_records", "first_date", "last_date")]

  list(data              = users,
       extra             = list("_dispensings" =
                                 if (keep_dispensings) dispensings else NULL),
       n_codes           = length(include_codes),
       pattern           = stage2_pattern,
       n_rows            = n_rows,
       # n_persons_any counts everyone who matched the codes; n_persons only
       # those who met the rule. The gap between them is what the rule excluded.
       n_persons_any     = length(unique(dispensings[[id_col]])),
       n_persons         = nrow(users),
       min_prescriptions = min_prescriptions,
       window_days       = window_days)
}
# 4. Small helpers ----

#' Anchored alternation, e.g. c("N02A","N02BE") -> "^(N02A|N02BE)"
#'
#' Anchored because an unanchored pattern also matches in the middle of a code.
#' One alternation rather than a loop because it is a single pass in the
#' database engine.
#' @keywords internal
mb_pattern <- function(codes) {
  codes <- unique(codes[!is.na(codes) & nzchar(codes)])
  paste0("^(", paste(codes, collapse = "|"), ")")
}


#' Restrict a query to a set of person ids, before anything else runs
#'
#' semi_join() rather than filter(id_col %in% ids): a large `ids` vector
#' passed as a literal IN-list can translate badly or not at all on some lazy
#' backends, where a JOIN against a copied-in table is the documented,
#' reliable way to filter by a large local vector. copy = TRUE lets this work
#' whether `query` is local or remote - dplyr copies `ids` across only when it
#' needs to, and does nothing extra for a plain data frame.
#'
#' Filtering by a known cohort BEFORE the code filter, rather than after, is
#' the pattern used by other DST projects reading the same registers (see
#' DECISIONS.md 7.5): it is a second, independent lever on the RAM problem,
#' unrelated to whether a rule's window function pushes down.
#' @keywords internal
mb_restrict_ids <- function(query, id_col, ids) {
  if (is.null(ids)) return(query)
  id_tbl <- data.frame(unique(ids), stringsAsFactors = FALSE)
  names(id_tbl) <- id_col
  dplyr::semi_join(query, id_tbl, by = id_col, copy = TRUE)
}


#' Coerce a dispensing date, refusing rather than guessing
#'
#' A column that parses to all-NA would otherwise produce zero qualifying users
#' and look like a real, empty finding.
#' @keywords internal
mb_as_date <- function(v, label = "date") {
  if (inherits(v, "Date")) return(v)
  if (inherits(v, "POSIXt")) return(as.Date(v))
  out <- suppressWarnings(as.Date(v))
  if (all(is.na(out)) && length(v)) {
    stop("Could not read column '", label, "' as a Date. Convert it before ",
         "extracting.", call. = FALSE)
  }
  out
}


#' Column names of a data frame or lazy table, or NULL if they cannot be read
#'
#' NULL means "unknown", not "none" - callers must not treat it as an empty set,
#' or an unfamiliar backend would look like a table with no columns.
#' @keywords internal
mb_cols <- function(x) {
  nms <- try(colnames(x), silent = TRUE)
  if (inherits(nms, "try-error") || is.null(nms)) {
    nms <- try(names(x), silent = TRUE)
  }
  if (inherits(nms, "try-error") || is.null(nms) || !length(nms)) return(NULL)
  as.character(nms)
}


#' @keywords internal
mb_has_col <- function(x, col) {
  if (is.null(col)) return(FALSE)
  nms <- mb_cols(x)
  !is.null(nms) && col %in% nms
}


#' The empty result, shaped exactly like a real one
#'
#' Returned when a condition matches nothing, so downstream rbind() and column
#' selection behave the same whether or not anyone had the condition.
#' @keywords internal
mb_null_result <- function(id_col, condition, n_codes, pattern,
                           min_prescriptions, window_days) {
  out <- data.frame(
    id              = character(0),
    condition       = character(0),
    onset_date      = as.Date(character(0)),
    onset_code      = character(0),
    n_records = integer(0),
    first_date      = as.Date(character(0)),
    last_date       = as.Date(character(0)),
    stringsAsFactors = FALSE
  )
  names(out)[1] <- id_col

  list(data = out, extra = list(), n_codes = n_codes, pattern = pattern,
       n_rows = 0L, n_persons_any = 0L, n_persons = 0L,
       min_prescriptions = min_prescriptions, window_days = window_days)
}


#' One row of the run summary
#'
#' Records the rule that was applied and the pattern that was used, not just the
#' counts, so a result can be audited months later without rerunning it.
#' @keywords internal
mb_summary_row <- function(condition, status, file = NA_character_,
                           seconds = NA_real_, n_codes = NA_integer_,
                           pattern = NA_character_, n_rows = NA_integer_,
                           n_persons_any = NA_integer_, n_persons = NA_integer_,
                           min_prescriptions = NA_integer_,
                           window_days = NA_integer_, message = NA_character_) {
  data.frame(
    condition         = condition,
    status            = status,
    n_codes           = n_codes,
    pattern           = pattern,
    n_rows            = n_rows,
    n_persons_any     = n_persons_any,
    n_persons         = n_persons,
    min_prescriptions = min_prescriptions,
    window_days       = window_days,
    seconds           = round(seconds, 1),
    file              = file,
    message           = message,
    stringsAsFactors  = FALSE
  )
}
