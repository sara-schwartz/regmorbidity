# .............................................................................
# extract_batch.R - single-pass, batched medication extraction
#
# PURPOSE
#   mb_extract_medication() runs one condition at a time: for every condition
#   it filters, collect()s the matching raw dispensing rows into R, then
#   applies the "N dispensings within window_days" rule locally via
#   data.table (mb_flag_users_impl() / mb_onset() in extract.R). That works,
#   but does one full register scan and one collect() per condition.
#
#   mb_extract_medication_batch() does the same extraction, but batches every
#   condition sharing a prescription rule (min_prescriptions, window_days)
#   into ONE query: a UNION ALL of independently-filtered branches, followed
#   by ONE shared window-function pass computing "the previous qualifying
#   dispensing" per (person, condition). Only the final onset rows are ever
#   collect()ed - never the raw dispensings. Confirmed correct against Jie's
#   flag_users() for all 15 real conditions on the real DST register,
#   2026-09-26 - see TODO.txt section 5 and section 8, and DECISIONS.md.
#
# WHAT THIS DOES NOT DO (v1, deliberately)
#   No prefilter_col/prefilter_len - there is nothing to shrink, since raw
#   rows are never collected here at all. No keep_all_users - would change
#   the .rds shape (NA-onset rows) in a way not yet decided for the merge
#   functions. No keep_events - reintroduces exactly the raw-collect
#   cost this file exists to avoid. For prevalence inputs, use sequential
#   mb_extract_medication(conditions = "x", keep_events = TRUE, from/to/ids
#   as needed) instead.
#
# WHY A SEPARATE FILE, NOT A CHANGE TO mb_extract_medication()
#   Genuinely different code path (batched SQL vs. per-condition R loop),
#   not a variant of the existing function - kept apart on purpose so the
#   proven, tested one-condition-at-a-time path is untouched.
#
# CONTENTS
#   1. Grouping conditions into rule-batches
#   2. The batched query
#   3. The batch driver (checkpoint/resume at batch granularity)
#   4. mb_extract_medication_batch()
# .............................................................................


# 1. Grouping conditions into rule-batches ----

#' Split a code list into batches sharing one prescription rule
#'
#' Grouped by (min_prescriptions, window_days), not min_prescriptions alone.
#' mb_validate_codelist() already guarantees both are constant within a
#' condition, so this is always a safe, lossless split - and it matters
#' because the window comparison inside mb_batch_query() becomes a single
#' literal per batch. All 15 bundled ATC conditions currently fall into two
#' batches, (2, 365) and (4, 365); a future condition with a different
#' window becomes a third batch automatically, not a code change.
#' @keywords internal
mb_rule_batches <- function(codes, conditions) {
  codes <- codes[codes$condition %in% conditions, , drop = FALSE]
  rules <- unique(codes[, c("condition", "min_prescriptions", "window_days")])

  key <- paste(rules$min_prescriptions, rules$window_days, sep = "|")
  lapply(split(seq_len(nrow(rules)), key), function(idx) {
    list(conditions        = rules$condition[idx],
        min_prescriptions = rules$min_prescriptions[idx[1]],
        window_days       = rules$window_days[idx[1]])
  })
}


# 2. The batched query ----

#' Is this a lazy (remote) table, or a plain local data frame?
#'
#' window_order() only works on a dbplyr tbl_lazy - it errors outright on a
#' local data frame ("Did you mean arrange()?"), and a local group_by() %>%
#' lag() without an explicit arrange() first is silently positional, not
#' date-ordered (confirmed directly, not assumed). So the window-function
#' step needs two paths; this is the gate between them.
#' @keywords internal
mb_is_lazy <- function(x) inherits(x, "tbl_lazy")


#' One condition's filtered, deduplicated, labelled branch
#'
#' Builds one arm of the UNION ALL in mb_batch_query(). Same-day
#' de-duplication happens here, per branch, as a pushed-down GROUP BY +
#' MIN() - not after collect(), because there is no raw collect() in this
#' design at all. The final select() runs on BOTH paths (deduped or not) so
#' every branch has identical column order before the union - UNION ALL
#' aligns by position, not by name, so this is a correctness requirement,
#' not tidiness.
#' @keywords internal
mb_batch_branch <- function(base, condition, include_codes, exclude_codes,
                           id_col, code_col, date_col, dedupe_same_day) {

  include_patt <- mb_pattern(include_codes)
  branch <- dplyr::filter(base, grepl(!!include_patt, .data[[code_col]]))

  if (length(exclude_codes)) {
    exclude_patt <- mb_pattern(exclude_codes)
    branch <- dplyr::filter(branch, !grepl(!!exclude_patt, .data[[code_col]]))
  }

  if (dedupe_same_day) {
    # One row per person-date, keeping the alphabetically-first code on a
    # tie - the same tie-break mb_flag_users_impl()'s sort-then-unique
    # gives, done here as a pushed-down GROUP BY/MIN instead of a local sort.
    #
    # suppressWarnings(): on a LOCAL data frame with zero rows (a condition
    # matching nothing), dplyr type-checks this expression against an empty
    # prototype and min() warns "no non-missing arguments" even though the
    # actual result is correctly zero rows - confirmed harmless, not
    # masking a real problem, by testing the zero-row case directly.
    branch <- suppressWarnings(
      branch %>%
        dplyr::group_by(.data[[id_col]], .data[[date_col]]) %>%
        dplyr::summarise(!!code_col := min(.data[[code_col]]), .groups = "drop")
    )
  }

  branch <- dplyr::select(branch, dplyr::all_of(c(id_col, code_col, date_col)))
  dplyr::mutate(branch, condition = !!condition)
}


#' Build the full lazy query for one rule-batch
#'
#' Returns an UNCOLLECTED table: one row per person per condition that
#' matched anything, with n_records/first_date/last_date (every matching
#' dispensing) and onset_date/onset_code (the qualifying pair's date and
#' the code that triggered it, NA if the person never qualified) computed
#' alongside each other - all inside the one collect() the caller does once.
#'
#' from/to are applied here, upstream of every branch, before the window
#' function runs - not as a filter on the computed onset afterward. A
#' pre-window dispensing can be the earlier half of a qualifying pair whose
#' LATER half falls inside the window; filtering post hoc would keep that
#' pair wrongly (2026-09-26 DST session, the bipolar anomaly - TODO.txt
#' section 5).
#' @keywords internal
mb_batch_query <- function(lmdb, codes, batch, id_col, code_col, date_col,
                          from, to, dedupe_same_day, ids) {

  base <- mb_restrict_ids(lmdb, id_col, ids)
  base <- mb_restrict_dates(base, date_col, from, to)

  branches <- lapply(batch$conditions, function(condition) {
    rules         <- codes[codes$condition == condition, , drop = FALSE]
    include_codes <- unique(rules$code[!rules$exclude])
    exclude_codes <- unique(rules$code[rules$exclude])
    mb_batch_branch(base, condition, include_codes, exclude_codes,
                    id_col, code_col, date_col, dedupe_same_day)
  })
  matched <- Reduce(dplyr::union_all, branches)
  lazy    <- mb_is_lazy(lmdb)

  # Every matching (deduped) dispensing, independent of whether the rule
  # was ever met - the same meaning n_records/first_date/last_date already
  # have in mb_extract_one() (R/extract.R).
  counts <- matched %>%
    dplyr::group_by(.data[[id_col]], .data$condition) %>%
    dplyr::summarise(n_records  = dplyr::n(),
                     first_date = min(.data[[date_col]]),
                     last_date  = max(.data[[date_col]]),
                     .groups    = "drop")

  offset <- as.integer(batch$min_prescriptions) - 1L

  windowed <- matched %>% dplyr::group_by(.data[[id_col]], .data$condition)
  windowed <- if (lazy) {
    windowed %>%
      dbplyr::window_order(.data[[date_col]]) %>%
      dplyr::mutate(prev = dplyr::lag(.data[[date_col]], n = !!offset))
  } else {
    windowed %>%
      dplyr::arrange(.data[[id_col]], .data$condition, .data[[date_col]]) %>%
      dplyr::mutate(prev = dplyr::lag(.data[[date_col]], n = !!offset))
  }
  windowed <- dplyr::ungroup(windowed)
  windowed <- if (is.na(batch$window_days)) {
    dplyr::filter(windowed, !is.na(.data$prev))
  } else {
    dplyr::filter(windowed, !is.na(.data$prev),
                 .data[[date_col]] - .data$prev <= !!batch$window_days)
  }

  # The qualifying row itself - earliest date, and on a same-day tie the
  # alphabetically-first code, matching the tie-break used everywhere else.
  onset_pick <- windowed %>% dplyr::group_by(.data[[id_col]], .data$condition)
  onset_pick <- if (lazy) {
    onset_pick %>% dbplyr::window_order(.data[[date_col]], .data[[code_col]])
  } else {
    dplyr::arrange(onset_pick, .data[[id_col]], .data$condition,
                   .data[[date_col]], .data[[code_col]])
  }
  onset_pick <- onset_pick %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(!!id_col := .data[[id_col]], condition = .data$condition,
                     onset_date = .data[[date_col]],
                     onset_code = .data[[code_col]])

  # Left, not inner: keeps everyone who matched anything even without a
  # qualifying pair (onset columns NA) - the caller decides what to do with
  # those rows (v1 drops them, matching mb_extract_medication()'s default).
  dplyr::left_join(counts, onset_pick, by = c(id_col, "condition"))
}


#' The empty result for one condition, shaped exactly like a real one
#'
#' Empty batch result shaped like a real onset table (includes `source`).
#' Matches [mb_null_result()] in extract.R.
#' @keywords internal
mb_null_result_batch <- function(id_col, condition) {
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
  out
}


# 3. The batch driver (checkpoint/resume at batch granularity) ----

#' Run a worker over rule-batches, checkpointing per condition but resuming
#' per batch
#'
#' A single query computes every condition in one batch at once and cannot
#' be resumed mid-query the way mb_run_conditions()'s per-condition loop
#' can (R/run.R) - so resume here works at the coarser granularity a
#' batched query actually allows: skip a whole batch only when EVERY
#' member condition's `.rds` already exists; otherwise rerun the whole
#' batch, overwriting any of its conditions already done.
#'
#' A real, stated consequence: one condition's data causing the batch's
#' query to fail now fails every condition sharing its rule too - true
#' per-condition isolation only existed in the old driver.
#'
#' `worker(batch)` must return a named list, one entry per
#' `batch$conditions`, each `list(data = <data frame>, n_rows = ,
#' n_persons_any = , n_persons = )`.
#' @keywords internal
mb_run_batches <- function(batches, worker, outdir, resume, verbose) {

  summaries <- list()
  results   <- list()

  for (batch in batches) {
    conds <- batch$conditions
    files <- if (is.null(outdir)) rep(NA_character_, length(conds)) else
      file.path(outdir, paste0(conds, ".rds"))

    if (!is.null(outdir) && resume && all(file.exists(files))) {
      if (verbose) {
        message("[rule min_prescriptions=", batch$min_prescriptions,
                " window_days=", batch$window_days, "] ", length(conds),
                " condition(s) - already done, skipping")
      }
      for (j in seq_along(conds)) {
        summaries[[length(summaries) + 1]] <- mb_summary_row(
          conds[j], status = "skipped", file = files[j],
          min_prescriptions = batch$min_prescriptions,
          window_days = batch$window_days)
      }
      next
    }

    if (verbose) {
      message("[rule min_prescriptions=", batch$min_prescriptions,
              " window_days=", batch$window_days, "] ", length(conds),
              " condition(s) ...")
    }
    started_at <- Sys.time()

    # try() rather than letting it propagate: a batch failing must not
    # throw away every OTHER batch's results, even though it does take
    # down every condition sharing its own rule (see @keywords doc above).
    outcome <- try(worker(batch), silent = TRUE)
    elapsed_secs <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))

    if (inherits(outcome, "try-error")) {
      msg <- conditionMessage(attr(outcome, "condition"))
      warning("Rule batch (min_prescriptions=", batch$min_prescriptions,
             ", window_days=", batch$window_days, ") failed: ", msg,
             call. = FALSE)
      for (j in seq_along(conds)) {
        summaries[[length(summaries) + 1]] <- mb_summary_row(
          conds[j], status = "error", file = files[j], seconds = elapsed_secs,
          message = msg, min_prescriptions = batch$min_prescriptions,
          window_days = batch$window_days)
      }
      next
    }

    for (j in seq_along(conds)) {
      condition <- conds[j]
      one <- outcome[[condition]]

      if (!is.null(outdir)) {
        saveRDS(one$data, files[j])
      } else {
        results[[condition]] <- one$data
      }

      summaries[[length(summaries) + 1]] <- mb_summary_row(
        condition, status = "done", file = files[j], seconds = elapsed_secs,
        n_rows = one$n_rows, n_persons_any = one$n_persons_any,
        n_persons = one$n_persons, min_prescriptions = batch$min_prescriptions,
        window_days = batch$window_days)

      if (verbose) {
        message("    ", condition, ": ",
               format(one$n_persons, big.mark = ","), " persons from ",
               format(one$n_rows, big.mark = ","), " dispensings")
      }
    }
  }

  summary <- do.call(rbind, summaries)
  rownames(summary) <- NULL
  list(results = results, summary = summary)
}


# 4. mb_extract_medication_batch() ----

#' Extract every medication condition sharing a prescription rule in one query
#'
#' [mb_extract_medication()] runs one condition at a time. This runs every
#' condition sharing the same `min_prescriptions`/`window_days` rule in ONE
#' query instead - a UNION of independently-filtered branches, one shared
#' window-function pass - so 15 conditions become 2 queries, not 15, and no
#' raw dispensing row is ever collected into R. Confirmed correct against
#' [mb_flag_users()] for all 15 real conditions on the real DST register,
#' 2026-09-26 (TODO.txt section 5).
#'
#' An opt-in alternative to [mb_extract_medication()], not a replacement -
#' see `@details` for what it does not (yet) do.
#'
#' @param lmdb The dispensing register. A lazy table (dbplyr `tbl_lazy`,
#'   e.g. from `DBI::dbConnect()` + `dplyr::tbl()`) or a plain data frame.
#'   Needs a person id, a full ATC code and a date column.
#' @param conditions Conditions to extract. `NULL` = every ATC condition in
#'   `codes`.
#' @param codes A code list; see [mb_codelist()].
#' @param outdir Directory for one `.rds` per condition. `NULL` = return
#'   the results as a named list instead of writing anything.
#' @param resume Skip a whole rule-batch when every one of its conditions'
#'   `.rds` already exists. Coarser than [mb_extract_medication()]'s
#'   per-condition resume - see `@details`.
#' @param id_col,code_col,date_col Column names. `code_col` must hold the
#'   *full* ATC code, not a truncated level.
#' @param from Inclusive lower date bound on `date_col` (**required** `Date`;
#'   error if `NULL`). Applied early, upstream of every branch and the
#'   window-function pass. For LMDB pass `from = as.Date("1997-01-01")` -
#'   see [mb_extract_medication()] and DECISIONS.md.
#' @param to Inclusive upper date bound (`Date` or `NULL`). `NULL` (default)
#'   = no upper bound.
#' @param dedupe_same_day Count several same-day dispensings as one. See
#'   [mb_flag_users()].
#' @param ids Restrict to these person ids before anything else runs, via
#'   [mb_restrict_ids()]. `NULL` (default) keeps everyone.
#' @param verbose Print progress.
#'
#' @details
#' Not implemented in this function, on purpose:
#' * `prefilter_col`/`prefilter_len` - nothing to shrink, since no raw row
#'   is ever collected here.
#' * `keep_all_users` - would change the `.rds` shape (NA-onset rows kept)
#'   in a way not yet decided for the merge functions.
#' * `keep_events` - reintroduces exactly the raw-collect memory cost
#'   this function exists to avoid. For prevalence inputs use sequential
#'   `mb_extract_medication(conditions = "x", keep_events = TRUE)`
#'   (plus `from`/`to`/`ids` as needed) instead.
#'
#' Resume is coarser than [mb_extract_medication()]'s: a batch (every
#' condition sharing one prescription rule) is skipped only when ALL of
#' its conditions already have a `.rds`, since a single query cannot be
#' resumed partway through the way a per-condition loop can. One
#' consequence worth knowing before trusting a resumed run: a failure in
#' one condition's data fails every condition sharing its rule, not just
#' itself - true per-condition isolation only exists in
#' [mb_extract_medication()].
#'
#' @return With `outdir`: a summary data frame, invisibly. Without: a named
#'   list of per-condition data frames, with the summary attached as an
#'   attribute. Output columns match [mb_extract_medication()]'s exactly,
#'   so [mb_merge_conditions()]/[mb_merge_all()]/[mb_load_conditions()]
#'   work unchanged on either function's output.
#' @export
mb_extract_medication_batch <- function(lmdb,
                                        conditions      = NULL,
                                        codes           = mb_codelist(),
                                        outdir          = NULL,
                                        resume          = TRUE,
                                        id_col          = "pnr",
                                        code_col        = "atc",
                                        date_col        = "eksd",
                                        from            = NULL,
                                        to              = NULL,
                                        dedupe_same_day = TRUE,
                                        ids             = NULL,
                                        verbose         = TRUE) {

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.", call. = FALSE)
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

  available_cols <- mb_cols(lmdb)
  if (!is.null(available_cols)) {
    for (nm in c(id_col, code_col, date_col)) {
      if (!nm %in% available_cols) {
        stop("Column '", nm, "' is not in the data. Present: ",
             paste(utils::head(available_cols, 30), collapse = ", "),
             call. = FALSE)
      }
    }
  }

  from <- mb_require_lmdb_from(from)
  to   <- mb_as_bound(to, "to")

  if (!is.null(outdir) && !dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  batches <- mb_rule_batches(codes, conditions)

  worker <- function(batch) {
    result    <- mb_batch_query(lmdb, codes, batch, id_col, code_col,
                                date_col, from, to,
                                dedupe_same_day, ids)
    collected <- dplyr::collect(result)

    stats::setNames(lapply(batch$conditions, function(condition) {
      rows <- collected[collected$condition == condition, , drop = FALSE]

      if (!nrow(rows)) {
        return(list(data = mb_null_result_batch(id_col, condition),
                    n_rows = 0L, n_persons_any = 0L, n_persons = 0L))
      }

      n_persons_any <- nrow(rows)

      out <- data.frame(
        id         = rows[[id_col]],
        condition  = condition,
        source     = "medication",
        onset_date = rows$onset_date,
        onset_code = rows$onset_code,
        n_records  = rows$n_records,
        first_date = rows$first_date,
        last_date  = rows$last_date,
        stringsAsFactors = FALSE
      )
      names(out)[1] <- id_col
      # v1 drops the NA-onset rows, matching mb_extract_medication()'s
      # default (keep_all_users = FALSE, deferred here - see @details).
      out <- out[!is.na(out$onset_date), , drop = FALSE]

      list(data = out, n_rows = sum(rows$n_records),
          n_persons_any = n_persons_any, n_persons = nrow(out))
    }), batch$conditions)
  }

  run <- mb_run_batches(batches, worker, outdir, resume, verbose)
  mb_finish_run(run, outdir, verbose)
}
