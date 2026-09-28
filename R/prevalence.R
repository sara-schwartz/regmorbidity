# .............................................................................
# prevalence.R - condition status at a point in time, not just onset
#
# PURPOSE
#   mb_count_conditions(as_of=) answers "has this person ever met the rule by
#   T", which is the same as Prior's estimand only for conditions that never
#   resolve. This answers a different question: is there a qualifying record
#   close enough to T, using the SAME rule extraction already applies, re-run
#   against a record set restricted to a lookback window.
#
# INPUT   one condition's raw records - mb_extract_medication(keep_events=)
#         or mb_extract_diagnosis(keep_events=) - one row per matching record
# OUTPUT  one row per person prevalent at `as_of` under `lookback`
#
# WHY THIS DOES NOT DUPLICATE THE WINDOW RULE
#   mb_flag_users_impl() already IS the rule, factored out in extract.R. This
#   restricts the record set to a window and hands it to that same function
#   unchanged, so the rule cannot drift between the two paths.
#
# THE BOUNDARY PROBLEM, AND WHY THIS WIDENS THE CANDIDATE WINDOW
#   A qualifying run can start up to window_days before it completes (see
#   mb_onset()). Restricting candidates to exactly [as_of - lookback, as_of]
#   can therefore cut a qualifying run in half - one record just inside, its
#   partner just outside - and silently miss someone who is plausibly still
#   under active treatment. Fixed by pulling in an extra window_days on the
#   left, running the SAME rule against that wider set, but reporting the
#   LATEST qualifying run rather than the first (see mb_onset(latest=)) - the
#   first qualifying run in a widened set is not necessarily the one that
#   makes someone currently prevalent. Anyone whose latest qualifying run
#   still completes before the true lookback boundary is dropped afterwards.
#
#   This only applies when there is a fixed span to widen by: min_prescriptions
#   > 1 (a single record cannot be "split"), window_days is not NA (no fixed
#   span to buffer by when the two records may be any distance apart), and
#   lookback is finite (nothing to widen against otherwise). Outside those
#   conditions this reduces to the plain restrict-then-reapply version, which
#   is exact: lookback = Inf reproduces mb_count_conditions() exactly, because
#   there is no boundary to straddle in the first place.
#
# SEE ALSO
#   DECISIONS.md 4.1, ASSUMPTIONS_AND_LIMITATIONS.txt section 7.
#
# CONTENTS
#   1. Prevalence at a point in time
#   2. Prevalence across every condition, in one call
# .............................................................................
# 1. Prevalence at a point in time ----

#' Is a condition still recorded, as of a date
#'
#' `mb_count_conditions(as_of = )` counts a condition from its onset date
#' onward forever - once met, always met. This answers a different question:
#' does a qualifying record fall in `[as_of - lookback, as_of]`. With
#' `lookback = Inf` it reproduces `mb_count_conditions(as_of = )` exactly,
#' because both then ask only whether onset happened by `as_of`.
#'
#' Applies the SAME rule as extraction - `min_prescriptions` within
#' `window_days` of each other, with same-day de-duplication - by delegating
#' to the internal function extraction itself uses. Works unchanged for
#' diagnosis records: pass `min_prescriptions = 1`, for which `window_days`
#' has no effect, and any qualifying record in the window counts.
#'
#' When `min_prescriptions > 1`, `window_days` is finite and `lookback` is
#' finite, the candidate record set is widened by `window_days` on the left
#' before the rule is applied, and the MOST RECENT qualifying run is reported
#' rather than the first - otherwise a qualifying run split across the
#' `as_of - lookback` boundary would be missed even though the person may
#' still be under active treatment. See `ASSUMPTIONS_AND_LIMITATIONS.txt`
#' section 7 for what this does and does not fix.
#'
#' @param data One condition's raw matching records - what
#'   `mb_extract_medication(keep_events = TRUE)` or
#'   `mb_extract_diagnosis(keep_events = TRUE)` write to
#'   `<condition>_all_events.rds`. Not the
#'   extraction result itself, which has already applied the rule and
#'   discarded the records.
#' @param as_of Single date. The index date to evaluate prevalence at.
#' @param lookback How many days back from `as_of` a qualifying record may
#'   fall. `Inf` (default) means no lower bound, equivalent to
#'   `mb_count_conditions(as_of = )`. Anders's sheet gives per-condition
#'   values ("Diagnosis time frame": Ever / last two years / last five years
#'   for cancer); the medication side uses one year. There is no code list
#'   column for this yet - the code lists are frozen - so it is an argument
#'   here, one call per condition, until that is unfrozen. See DECISIONS.md
#'   4.1.
#' @param code_col Column holding the code. Bare name or string.
#' @param min_prescriptions Records needed within `window_days` of each other.
#'   Use `1` for diagnosis records, where any one qualifying record counts.
#' @param window_days They must fall within this many days of each other.
#'   `NA` = any distance apart. Has no effect when `min_prescriptions = 1`.
#' @param id_col,date_col Person id and record date columns. Bare or string.
#' @param dedupe_same_day Count several same-day records as one. See
#'   [mb_flag_users()].
#' @return One row per person prevalent at `as_of`: id, `onset_date` (the
#'   qualifying date within the window, not necessarily the person's
#'   first-ever record), `onset_code`, `n_records`, `first_date`, `last_date`,
#'   `as_of`, `lookback_days`. People with no qualifying record in the window
#'   are absent, the same convention [mb_extract_medication()] uses.
#' @export
mb_prevalence <- function(data, as_of, lookback = Inf,
                          code_col          = atc,
                          min_prescriptions = 2L,
                          window_days       = NA,
                          id_col            = pnr,
                          date_col          = eksd,
                          dedupe_same_day   = TRUE) {

  code_col <- rlang::as_name(rlang::ensym(code_col))
  id_col   <- rlang::as_name(rlang::ensym(id_col))
  date_col <- rlang::as_name(rlang::ensym(date_col))

  as_of <- mb_as_date(as_of, "as_of")
  if (length(as_of) != 1L || is.na(as_of)) {
    stop("`as_of` must be a single, non-missing date.", call. = FALSE)
  }
  if (length(lookback) != 1L || is.na(lookback) || lookback <= 0) {
    stop("`lookback` must be a single positive number of days, or Inf.",
         call. = FALSE)
  }

  if (!date_col %in% names(data)) {
    stop("Column '", date_col, "' not found. Present: ",
         paste(names(data), collapse = ", "), call. = FALSE)
  }

  n <- suppressWarnings(as.integer(min_prescriptions))
  widen <- is.finite(lookback) && !is.na(n) && n > 1L && is.finite(window_days)
  buffer <- if (widen) window_days else 0

  dates <- mb_as_date(data[[date_col]], date_col)
  keep  <- !is.na(dates) & dates <= as_of
  if (is.finite(lookback)) keep <- keep & dates >= (as_of - lookback - buffer)

  restricted <- data[keep, , drop = FALSE]

  users <- mb_flag_users_impl(restricted, code_col, min_prescriptions,
                              window_days, id_col, date_col,
                              dedupe_same_day, keep_all = FALSE,
                              latest = widen)

  if (widen) {
    # The wider candidate set can surface a qualifying run that completes
    # before the true lookback boundary - reaching back this far was only
    # ever meant to keep a run from being cut in half, not to find older
    # evidence. Drop it now that the rule has had the full, untruncated run
    # to work with.
    users <- users[users$onset_date >= (as_of - lookback), , drop = FALSE]
  }

  # A bare scalar assignment errors against a zero-row data frame rather than
  # recycling to zero rows, so rep() it explicitly.
  users$as_of         <- rep(as_of, nrow(users))
  users$lookback_days <- rep(if (is.finite(lookback)) lookback else NA_real_,
                             nrow(users))
  users
}
# 2. Prevalence across every condition, in one call ----

# Anders's sheet states lookback as a word, not a number of days. Named here,
# the same way MB_TIME_FRAME_DAYS in codelist.R names window_days's own
# word-to-number rule, so both translations are visible and editable in one
# place rather than copied wherever someone needs to type a lookback.
MB_LOOKBACK_DAYS <- c(ever = Inf, last_two_years = 730, last_five_years = 1825)

#' Turn a word like "ever" or "last_two_years" into a number of days
#'
#' Accepts a plain number already, so a `lookback` column can freely mix
#' words and numbers. Case-insensitive, spaces and hyphens folded to `_` so
#' "Last two years" and "last-two-years" both match.
#' @keywords internal
mb_lookback_days <- function(x) {
  is_word <- is.na(suppressWarnings(as.numeric(x)))
  words <- tolower(trimws(gsub("[-[:space:]]+", "_", as.character(x[is_word]))))
  unknown <- setdiff(words, names(MB_LOOKBACK_DAYS))
  if (length(unknown)) {
    stop("Not a recognised lookback word: ", paste(unique(unknown), collapse = ", "),
         "\nUse a number of days, Inf, or one of: ",
         paste(names(MB_LOOKBACK_DAYS), collapse = ", "), call. = FALSE)
  }
  out <- suppressWarnings(as.numeric(x))
  out[is_word] <- unname(MB_LOOKBACK_DAYS[words])
  out
}

#' Prevalence across every condition, with its own lookback, in one call
#'
#' [mb_prevalence()] answers this for one condition at a time, with
#' `lookback` typed in by hand at every call. This is the same question
#' asked across a whole set of conditions at once, each using its own
#' `min_prescriptions` / `window_days` (from `codes`, the same list
#' extraction used) and its own `lookback` (from Anders's sheet, or
#' wherever - the code lists are frozen, so this stays a separate argument
#' rather than a code list column; see DECISIONS.md 4.1).
#'
#' @param dir Directory holding `<condition>_all_events.rds` (written by
#'   `mb_extract_medication(keep_events = TRUE)` or
#'   `mb_extract_diagnosis(keep_events = TRUE)`) - the raw records
#'   `mb_prevalence()` needs, not the extraction result itself.
#' @param as_of Single date, the index date - the same for every condition.
#' @param lookback A single value (number of days, a word - `"ever"` /
#'   `"last_two_years"` / `"last_five_years"` - or `Inf`), applied to every
#'   condition; or a data
#'   frame with columns `condition` and `lookback` giving a different value
#'   per condition. A condition present in `dir` but missing from a named
#'   `lookback` table defaults to `Inf` (ever/never), with a message.
#' @param codes A code list; see [mb_codelist()]. Supplies
#'   `min_prescriptions` and `window_days` per condition - the same rule
#'   extraction already applied, so it cannot drift between the two.
#' @param conditions Conditions to include. `NULL` = every condition with
#'   both a raw-records file in `dir` and a row in `codes`.
#' @param id_col,code_col,date_col,dedupe_same_day Passed through to
#'   [mb_prevalence()] for every condition.
#' @param verbose Print progress.
#' @return One combined data frame, one row per person-condition prevalent
#'   at `as_of`, `condition` column added - the same columns
#'   [mb_prevalence()] returns for one condition.
#' @keywords internal
mb_prevalence_all <- function(dir, as_of, lookback = Inf,
                              codes           = mb_codelist(),
                              conditions      = NULL,
                              id_col          = pnr,
                              code_col        = atc,
                              date_col        = eksd,
                              dedupe_same_day = TRUE,
                              verbose         = TRUE) {

  id_col   <- rlang::as_name(rlang::ensym(id_col))
  code_col <- rlang::as_name(rlang::ensym(code_col))
  date_col <- rlang::as_name(rlang::ensym(date_col))

  if (!dir.exists(dir)) stop("No such directory: ", dir, call. = FALSE)

  # Raw-records files, not extraction results - the opposite exclusion to
  # mb_condition_files(), which is for the results and skips these.
  files <- list.files(dir, pattern = "_all_events\\.rds$")
  file_condition <- sub("_all_events\\.rds$", "", files)

  if (!is.null(conditions)) {
    missing_files <- setdiff(conditions, file_condition)
    if (length(missing_files)) {
      stop("No raw-records file for: ", paste(missing_files, collapse = ", "),
           "\n(Looked for <condition>_all_events.rds in ",
           dir, ")", call. = FALSE)
    }
    keep <- file_condition %in% conditions
    files <- files[keep]
    file_condition <- file_condition[keep]
  }
  if (!length(files)) {
    stop("No <condition>_all_events.rds files in ", dir,
         call. = FALSE)
  }

  rules <- unique(codes[, c("condition", "min_prescriptions", "window_days")])
  missing_rules <- setdiff(file_condition, rules$condition)
  if (length(missing_rules)) {
    stop("No entry in `codes` for: ", paste(missing_rules, collapse = ", "),
         call. = FALSE)
  }

  # One lookback value per condition, defaulting to Inf (ever/never) for any
  # condition not named - the same default mb_prevalence() itself has.
  if (is.data.frame(lookback)) {
    lb <- stats::setNames(mb_lookback_days(lookback$lookback), lookback$condition)
    missing_lb <- setdiff(file_condition, names(lb))
    if (length(missing_lb) && verbose) {
      message("No lookback given for: ", paste(missing_lb, collapse = ", "),
              " - defaulting to Inf (ever/never).")
    }
    lookback_for <- function(cond) {
      if (cond %in% names(lb)) lb[[cond]] else Inf
    }
  } else {
    single <- mb_lookback_days(lookback)
    lookback_for <- function(cond) single
  }

  parts <- lapply(seq_along(files), function(i) {
    cond <- file_condition[i]
    if (verbose) message("Prevalence: ", cond)

    data <- readRDS(file.path(dir, files[i]))
    rule <- rules[rules$condition == cond, , drop = FALSE][1L, ]

    # inject() + !! splices the already-resolved strings in as literals
    # before the call happens. Without it, mb_prevalence()'s own
    # rlang::ensym(code_col) would capture the SYMBOL NAME "code_col" from
    # this call site - not its value "atc" - and silently look for a column
    # literally called "code_col".
    result <- rlang::inject(mb_prevalence(
      data, as_of = as_of, lookback = lookback_for(cond),
      code_col = !!code_col, id_col = !!id_col, date_col = !!date_col,
      min_prescriptions = rule$min_prescriptions,
      window_days       = rule$window_days,
      dedupe_same_day   = dedupe_same_day
    ))
    if (!nrow(result)) return(NULL)
    result$condition <- cond
    result
  })

  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (!length(parts)) {
    warning("Every condition was empty - no one prevalent at as_of under ",
            "these lookbacks.", call. = FALSE)
    return(NULL)
  }

  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}
