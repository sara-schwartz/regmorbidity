# .............................................................................
# exclusions.R - Prior's cross-condition exclusion rules (optional stage 2)
#
# PURPOSE
#   After every condition has its own onset row, Prior's composites drop some
#   of them when another condition is also present (lipid unless IHD,
#   hypertension medication unless IHD/HF/(kidney for diuretics), distress
#   unless other mental disorders). That cannot live in a per-condition CSV
#   and cannot run inside extract. This is the optional second pass over the
#   assembled long or wide table.
#
# INPUT   long (pnr, condition, onset_date, ...) or wide (one date column per
#         condition), from mb_merge_all() / mb_load_conditions() / mb_to_wide()
# OUTPUT  the same shape, with excluded person-condition rows (long) or dates
#         set to NA (wide)
#
# WHAT THIS IS NOT
#   Not baked into mb_extract_*(). Not a substitute for Prior's time-varying
#   spells. Content of the rules is still a joint call with Jie (DECISIONS.md
#   5.1); this is a provisional implementation so the shape exists.
#
# CONDITION NAMES
#   Mapped to the bundled CSV keys (dyslipidemia, ihd, hypertension, distress,
#   bipolar, dementia). Conditions Prior names that we do not yet have
#   (heart_failure, kidney, emotion, schizophrenia, eating_disorder) cause that
#   part of a rule to be skipped with a message.
#
# CONTENTS
#   1. Public entry point
#   2. Rule table and helpers
# .............................................................................


# 1. Public entry point ----

#' Apply Prior's cross-condition exclusion rules (optional stage 2)
#'
#' Runs after extraction and merge. Operates on **condition presence** using
#' onset dates already in a long table (`pnr`, `condition`, `onset_date`) or a
#' wide table from [mb_to_wide()]. Does **not** re-query the register and is
#' not called from [mb_extract_medication()] / [mb_extract_medication_batch()] /
#' [mb_extract_diagnosis()].
#'
#' Default timing is `"as_of"` when `as_of` is given, otherwise `"ever_ever"`.
#' If you pass `timing = "as_of"` with `as_of = NULL`, the function falls back
#' to `"ever_ever"` with a message.
#'
#' @param long_or_wide Long result (`condition` + `onset_date` columns) or wide
#'   result (one date column per condition).
#' @param as_of Optional index date. With `timing = "as_of"`, a condition is
#'   treated as present only when its onset is on or before this date.
#' @param timing One of `"as_of"`, `"onset_order"`, `"ever_ever"`. See Details.
#' @param id_col Person id column name.
#' @param quiet Suppress informational messages about skipped / approximate
#'   rules.
#' @return The same shape as `long_or_wide`. Long: excluded person-condition
#'   rows are dropped. Wide: excluded cells are set to `NA`.
#'
#' @details
#' **Timing**
#' \describe{
#'   \item{`ever_ever`}{Drop the target if the person ever has any excluding
#'     condition (any onset date).}
#'   \item{`onset_order`}{Drop the target only if an excluding condition's
#'     onset is on or before the target's onset.}
#'   \item{`as_of`}{Drop the target if, at `as_of`, both the target and an
#'     excluding condition are present (onset \eqn{\le} `as_of`).}
#' }
#'
#' **Default Prior rules (package condition keys)**
#' \itemize{
#'   \item `dyslipidemia` excluded if `ihd` (Prior `r_lipid`).
#'   \item `hypertension` excluded if `ihd` - **approximate**. Prior's `r_bt`
#'     also needs heart failure, and a separate diuretic arm that further
#'     excludes kidney disease. Those conditions are not in the bundled ATC
#'     lists, and C03 is folded into the hypertension CSV, so the Prior
#'     diuretic arm cannot be expressed until lists are split (clinician + Jie).
#'   \item `distress` excluded if `bipolar` or `dementia`. Prior also names
#'     emotion, schizophrenia, and eating disorder, which are not yet
#'     extractable here; those limbs are skipped with a message.
#' }
#'
#' A rule whose **target** is absent from the data is skipped. An excluding
#' condition absent from the data is skipped for that rule (provisional).
#'
#' @seealso [mb_merge_all()], [mb_to_wide()], [mb_count_conditions()],
#'   `DECISIONS.md` section 5.1.
#' @export
mb_apply_exclusions <- function(long_or_wide,
                                as_of = NULL,
                                timing = c("as_of", "onset_order", "ever_ever"),
                                id_col = "pnr",
                                quiet = FALSE) {

  stopifnot(is.data.frame(long_or_wide))
  if (!id_col %in% names(long_or_wide)) {
    stop("Column '", id_col, "' not found.", call. = FALSE)
  }

  timing_arg_missing <- missing(timing)
  timing <- if (timing_arg_missing) {
    if (!is.null(as_of)) "as_of" else "ever_ever"
  } else {
    match.arg(timing, c("as_of", "onset_order", "ever_ever"))
  }

  if (identical(timing, "as_of") && is.null(as_of)) {
    if (!quiet) {
      message("mb_apply_exclusions: timing = \"as_of\" but as_of is NULL; ",
              "falling back to ever_ever.")
    }
    timing <- "ever_ever"
  }

  if (!is.null(as_of)) {
    as_of <- as.Date(as_of)
    if (length(as_of) != 1L || is.na(as_of)) {
      stop("'as_of' must be a single Date (or NULL).", call. = FALSE)
    }
  }

  is_long <- all(c("condition", "onset_date") %in% names(long_or_wide))
  long <- if (is_long) {
    mb_exclusions_as_long(long_or_wide, id_col)
  } else {
    mb_exclusions_wide_to_long(long_or_wide, id_col)
  }

  present <- unique(as.character(long$condition))
  rules <- mb_exclusion_rules()
  drop_keys <- character(0)

  for (rule in rules) {
    msgs <- mb_exclusion_apply_one(long, rule, present, timing, as_of, id_col)
    drop_keys <- c(drop_keys, msgs$drop_keys)
    if (!quiet) {
      for (m in msgs$messages) message("mb_apply_exclusions: ", m)
    }
  }

  if (length(drop_keys)) {
    row_key <- paste(long[[id_col]], long$condition, sep = "\r")
    long <- long[!row_key %in% unique(drop_keys), , drop = FALSE]
    rownames(long) <- NULL
  }

  if (is_long) {
    # Preserve any extra columns from the original long table by re-filtering
    # the caller's rows rather than returning the stripped long.
    orig_key <- paste(long_or_wide[[id_col]], long_or_wide$condition, sep = "\r")
    keep_key <- paste(long[[id_col]], long$condition, sep = "\r")
    out <- long_or_wide[orig_key %in% keep_key, , drop = FALSE]
    rownames(out) <- NULL
    out
  } else {
    mb_exclusions_long_to_wide(long, long_or_wide, id_col)
  }
}


# 2. Rule table and helpers ----

#' Prior's three rules, mapped to bundled condition keys
#' @keywords internal
mb_exclusion_rules <- function() {
  list(
    list(
      name = "dyslipidemia",
      target = "dyslipidemia",
      exclude_if = "ihd",
      also_needed = character(0),
      approximate = FALSE,
      approx_note = NULL
    ),
    list(
      name = "hypertension",
      target = "hypertension",
      exclude_if = "ihd",
      also_needed = c("heart_failure", "kidney"),
      approximate = TRUE,
      approx_note = paste0(
        "hypertension rule is approximate: Prior also excludes on heart ",
        "failure, and diuretics (C03) further exclude on kidney disease. ",
        "Those conditions are not in the bundled lists, and C03 is folded ",
        "into hypertension.csv, so the Prior diuretic arm cannot be ",
        "expressed until lists are split (clinician + Jie)."
      )
    ),
    list(
      name = "distress",
      target = "distress",
      exclude_if = c("bipolar", "dementia"),
      also_needed = c("emotion", "schizophrenia", "eating_disorder"),
      approximate = TRUE,
      approx_note = paste0(
        "distress rule is partial: Prior also excludes on emotion, ",
        "schizophrenia, and eating disorder, which are not yet extractable ",
        "from the bundled ATC-only lists."
      )
    )
  )
}


#' @keywords internal
mb_exclusions_as_long <- function(df, id_col) {
  out <- data.frame(
    id = as.character(df[[id_col]]),
    condition = as.character(df$condition),
    onset_date = as.Date(df$onset_date),
    stringsAsFactors = FALSE
  )
  names(out)[1] <- id_col
  if (anyNA(out$condition) || any(!nzchar(out$condition))) {
    stop("Long input has missing `condition` values.", call. = FALSE)
  }
  out
}


#' @keywords internal
mb_exclusions_wide_to_long <- function(wide, id_col) {
  cond_cols <- setdiff(names(wide), c(id_col, "n_conditions", "multimorbid"))
  if (!length(cond_cols)) {
    stop("Wide input has no condition columns.", call. = FALSE)
  }
  parts <- lapply(cond_cols, function(cond) {
    dates <- as.Date(wide[[cond]])
    keep <- !is.na(dates)
    if (!any(keep)) {
      return(data.frame(
        character(0), character(0), as.Date(character(0)),
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      as.character(wide[[id_col]][keep]),
      cond,
      dates[keep],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, parts)
  if (is.null(out) || !nrow(out)) {
    out <- data.frame(
      character(0), character(0), as.Date(character(0)),
      stringsAsFactors = FALSE
    )
  }
  names(out) <- c(id_col, "condition", "onset_date")
  rownames(out) <- NULL
  out
}


#' Rebuild wide from filtered long, preserving id order and non-condition cols
#' @keywords internal
mb_exclusions_long_to_wide <- function(long, wide_template, id_col) {
  cond_cols <- setdiff(names(wide_template),
                       c(id_col, "n_conditions", "multimorbid"))
  out <- wide_template
  for (cond in cond_cols) {
    out[[cond]] <- as.Date(NA)
  }
  if (nrow(long)) {
    for (i in seq_len(nrow(long))) {
      id <- long[[id_col]][i]
      cond <- long$condition[i]
      if (!cond %in% cond_cols) next
      row <- match(id, out[[id_col]])
      if (!is.na(row)) out[[cond]][row] <- long$onset_date[i]
    }
  }
  # Drop derived count columns; caller can recompute after exclusions.
  out$n_conditions <- NULL
  out$multimorbid <- NULL
  out
}


#' Apply one rule; return drop keys and messages
#' @keywords internal
mb_exclusion_apply_one <- function(long, rule, present, timing, as_of, id_col) {
  messages <- character(0)
  drop_keys <- character(0)

  if (!rule$target %in% present) {
    messages <- c(messages,
                  paste0("rule '", rule$name, "': target '", rule$target,
                         "' not in data - skipped."))
    return(list(drop_keys = drop_keys, messages = messages))
  }

  if (isTRUE(rule$approximate) && !is.null(rule$approx_note)) {
    messages <- c(messages, rule$approx_note)
  }

  missing_needed <- setdiff(rule$also_needed, present)
  if (length(missing_needed)) {
    messages <- c(messages,
                  paste0("rule '", rule$name, "': also needs ",
                         paste(missing_needed, collapse = ", "),
                         " (not in data) - those limbs skipped."))
  }

  excluders <- rule$exclude_if[rule$exclude_if %in% present]
  missing_ex <- setdiff(rule$exclude_if, present)
  if (length(missing_ex)) {
    messages <- c(messages,
                  paste0("rule '", rule$name, "': excluding condition(s) ",
                         paste(missing_ex, collapse = ", "),
                         " not in data - skipped for this rule."))
  }
  if (!length(excluders)) {
    messages <- c(messages,
                  paste0("rule '", rule$name, "': no usable excluding ",
                         "conditions in data - rule not applied."))
    return(list(drop_keys = drop_keys, messages = messages))
  }

  targets <- long[long$condition == rule$target, , drop = FALSE]
  if (!nrow(targets)) {
    return(list(drop_keys = drop_keys, messages = messages))
  }

  for (ex in excluders) {
    excl <- long[long$condition == ex, , drop = FALSE]
    if (!nrow(excl)) next
    # Earliest onset per person if duplicates somehow exist
    excl <- excl[order(excl[[id_col]], excl$onset_date), , drop = FALSE]
    excl <- excl[!duplicated(excl[[id_col]]), , drop = FALSE]
    idx <- match(targets[[id_col]], excl[[id_col]])
    hit <- !is.na(idx)
    if (!any(hit)) next

    t_onset <- targets$onset_date[hit]
    e_onset <- excl$onset_date[idx[hit]]
    pids <- targets[[id_col]][hit]

    drop <- switch(timing,
      ever_ever = rep(TRUE, sum(hit)),
      onset_order = !is.na(e_onset) & !is.na(t_onset) & e_onset <= t_onset,
      as_of = !is.na(e_onset) & !is.na(t_onset) &
        e_onset <= as_of & t_onset <= as_of,
      rep(FALSE, sum(hit))
    )
    if (any(drop)) {
      drop_keys <- c(drop_keys, paste(pids[drop], rule$target, sep = "\r"))
    }
  }

  n_drop <- length(unique(drop_keys))
  if (n_drop) {
    messages <- c(messages,
                  paste0("rule '", rule$name, "': dropped ", n_drop,
                         " person-row(s) of '", rule$target, "' (",
                         timing, ")."))
  }

  list(drop_keys = unique(drop_keys), messages = messages)
}

