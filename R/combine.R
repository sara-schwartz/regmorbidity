# .............................................................................
# combine.R - putting the per-condition files back together
#
# PURPOSE
#   mb_extract_medication() writes one .rds per condition, which is what makes a
#   long run resumable but is not what anyone analyses. These functions read
#   them back into the two shapes an analysis actually uses: long (one row per
#   person-condition) and wide (one column per condition).
#
# INPUT   the outdir used for mb_extract_medication()
# OUTPUT  data frames; nothing is written
#
# WHY THE SHAPES ARE SEPARATE
#   Long keeps every detail - which code triggered onset, how many dispensings,
#   first and last date - and is what you check before trusting anything. Wide
#   throws all of that away to keep one date per condition, and is what goes
#   into a regression. Going long -> wide is lossy on purpose.
#
# CONTENTS
#   1. Reading the extraction output
#   2. Reshaping
#   3. Counting
# .............................................................................
# 1. Reading the extraction output ----

#' The .rds files in a directory that are conditions
#'
#' `save_dispensings` and `keep_records` write `<condition>_dispensings.rds` and
#' `<condition>_records.rds` alongside the results. Those are raw records with a
#' different shape, so reading them back as conditions fails - which it did,
#' until this was factored out and used by both readers.
#' @keywords internal
mb_condition_files <- function(dir, full.names = TRUE) {
  files <- list.files(dir, pattern = "\\.rds$")
  files <- files[!grepl("_(records|dispensings)\\.rds$", files)]
  if (full.names) file.path(dir, files) else files
}


#' Load the per-condition .rds files written by mb_extract_medication()
#'
#' @param dir The `outdir` used for [mb_extract_medication()].
#' @param conditions Conditions to load. `NULL` = every `.rds` in `dir`.
#' @return One long data frame, one row per person-condition. `NULL` if every
#'   file was empty.
#' @export
mb_load_conditions <- function(dir, conditions = NULL) {

  if (!dir.exists(dir)) stop("No such directory: ", dir, call. = FALSE)

  files <- mb_condition_files(dir)
  if (!length(files)) stop("No condition .rds files in ", dir, call. = FALSE)

  if (!is.null(conditions)) {
    wanted  <- file.path(dir, paste0(conditions, ".rds"))
    # Silently returning fewer conditions than asked for would quietly
    # undercount multimorbidity, so a missing file is an error.
    not_yet <- conditions[!file.exists(wanted)]
    if (length(not_yet)) {
      stop("Not extracted yet: ", paste(not_yet, collapse = ", "),
           call. = FALSE)
    }
    files <- wanted
  }

  parts <- lapply(files, function(f) {
    one_condition <- readRDS(f)
    if (!nrow(one_condition)) return(NULL)
    # The file name is authoritative. A stale `condition` column can survive a
    # renamed code list; the file name cannot.
    one_condition$condition <- tools::file_path_sans_ext(basename(f))
    one_condition
  })

  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (!length(parts)) {
    warning("Every file was empty - no one met the criteria.", call. = FALSE)
    return(NULL)
  }

  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}
# 2. Reshaping ----

#' Turn the long result into one column per condition
#'
#' @param long Output of [mb_load_conditions()].
#' @param id_col Person id column name.
#' @param value Which column to spread: `"onset_date"` (default) or
#'   `"n_dispensations"`.
#' @return A wide data frame: one row per person, one column per condition.
#'   `NA` means the person never met that condition's criteria.
#' @export
mb_to_wide <- function(long, id_col = "pnr", value = "onset_date") {

  stopifnot(is.data.frame(long))
  for (nm in c(id_col, "condition", value)) {
    if (!nm %in% names(long)) {
      stop("Column '", nm, "' not found in `long`.", call. = FALSE)
    }
  }

  # Only people who met at least one condition appear here, because that is all
  # the extraction wrote. Join onto the study population to get the healthy
  # ones back - they are not absent from the data, they are absent from these
  # files by construction.
  ids        <- sort(unique(long[[id_col]]))
  conditions <- sort(unique(long$condition))

  out <- data.frame(ids, stringsAsFactors = FALSE)
  names(out) <- id_col

  for (condition in conditions) {
    one <- long[long$condition == condition, , drop = FALSE]
    # match() silently takes the first row. One person should appear once per
    # condition; more than that means two runs were stacked, and quietly
    # keeping whichever came first would be an arbitrary choice of onset date.
    if (anyDuplicated(one[[id_col]])) {
      warning("Condition '", condition, "': ",
              sum(duplicated(one[[id_col]])), " duplicated person row(s). ",
              "Keeping the first of each - check the input.", call. = FALSE)
    }
    out[[condition]] <- one[[value]][match(ids, one[[id_col]])]
  }
  out
}
# 3. Counting ----

#' Count how many of the conditions each person has
#'
#' @param wide Output of [mb_to_wide()].
#' @param id_col Person id column name.
#' @param as_of Optional date. Only conditions with an onset on or before this
#'   date are counted, so the index can be read at a point in time.
#' @return `wide` with `n_conditions` and `multimorbid` added.
#'
#' @section What this does not do:
#' Conditions are counted independently. There is no composite logic - no
#' exclusion of dyslipidemia when IHD is present, no requirement of both a
#' diagnosis and a prescription for epilepsy - so a count from here is not
#' comparable to a published multimorbidity count. See
#' `ASSUMPTIONS_AND_LIMITATIONS.txt` section 7 and `TODO.txt` section 0(c).
#'
#' `as_of` counts onsets up to a date but cannot express recovery: once a
#' condition starts here, it never ends. The Stata original splits follow-up
#' into spells that can close again. That is a different estimand, not an
#' approximation of this one.
#' @export
mb_count_conditions <- function(wide, id_col = "pnr", as_of = NULL) {

  # Recomputing over an already-counted table would otherwise treat
  # n_conditions as a condition.
  condition_cols <- setdiff(names(wide),
                            c(id_col, "n_conditions", "multimorbid"))

  # as.matrix() on a data frame of Date columns silently formats them to
  # date STRINGS, not the numeric day count `as_of` is compared against below
  # - so every as_of comparison was a string against a number and effectively
  # always false. Converting column by column keeps the numeric day count
  # instead. vapply() drops the matrix dimension when wide has exactly one
  # row, so that case is reshaped back explicitly.
  onset_dates <- vapply(wide[, condition_cols, drop = FALSE],
                        function(x) as.numeric(as.Date(x)),
                        numeric(nrow(wide)))
  if (is.null(dim(onset_dates))) {
    onset_dates <- matrix(onset_dates, nrow = nrow(wide),
                          dimnames = list(NULL, condition_cols))
  }

  has_condition <- if (is.null(as_of)) {
    !is.na(onset_dates)
  } else {
    !is.na(onset_dates) & onset_dates <= as.numeric(as.Date(as_of))
  }

  out <- wide
  out$n_conditions <- rowSums(has_condition)
  out$multimorbid  <- out$n_conditions >= 2L
  out
}
