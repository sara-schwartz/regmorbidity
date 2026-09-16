# .............................................................................
# inspect.R - checks to run BEFORE an extraction
#
# PURPOSE
#   Both functions here exist to separate a real zero from a broken query. That
#   distinction is invisible in the output of an extraction: a condition that
#   nobody has and a condition whose codes match nothing both come back as an
#   empty file, with no error and no warning.
#
# INPUT   the dispensing register, and (for mb_check_codes) a code list
# OUTPUT  printed diagnostics and a small data frame; nothing is written
#
# WHEN TO RUN THEM
#   mb_inspect_codes()  once per new extract, before anything else
#   mb_check_codes()    whenever a code list is edited, and before interpreting
#                       any zero
#
# CONTENTS
#   1. Constants
#   2. Column shape
#   3. Code list against the register
# .............................................................................
# 1. Constants ----

# Enough rows to see every code length that occurs, small enough to collect from
# a lazy table in seconds. This is a sample, not a census - it can miss a rare
# malformed value.
MB_INSPECT_SAMPLE_ROWS <- 1e5

# Columns shorter than this cannot hold a full ATC code, so matching a real code
# against them is guaranteed to fail.
MB_SHORT_COLUMN_NCHAR <- 3L
# 2. Column shape ----

#' Report how long the code columns actually are
#'
#' The failure this catches is silent. `atc2` holds 3 characters, so
#' `grepl("N02A", atc2)` matches nothing, raises no error, and returns an empty
#' result that reads as a finding. Run this once on any new extract and confirm
#' the lengths are what the code lists assume.
#'
#' @param lmdb The dispensing register (lazy table or data frame).
#' @param cols Columns to inspect. Defaults to the full ATC column and the short
#'   level column.
#' @param n_sample Rows to pull for the check.
#' @return A data frame of column, code length and frequency, invisibly. Also
#'   printed, with example values.
#' @export
mb_inspect_codes <- function(lmdb,
                             cols = c("ATC", "atc2"),
                             n_sample = MB_INSPECT_SAMPLE_ROWS) {

  cols <- cols[vapply(cols, function(cc) mb_has_col(lmdb, cc), logical(1))]
  if (!length(cols)) {
    stop("None of those columns are in the data. Available: ",
         paste(utils::head(mb_cols(lmdb), 40), collapse = ", "), call. = FALSE)
  }

  sampled <- dplyr::collect(
    utils::head(dplyr::select(lmdb, dplyr::all_of(cols)), n_sample))

  per_column <- list()

  for (cc in cols) {
    values <- toupper(trimws(as.character(sampled[[cc]])))
    values <- values[!is.na(values) & nzchar(values)]
    lengths_seen <- table(nchar(values))

    cat("\n", cc, " - ", format(length(values), big.mark = ","),
        " non-missing values in the sample\n", sep = "")

    if (!length(lengths_seen)) {
      cat("  (all missing)\n")
      next
    }

    # Examples matter more than counts here: they show whether a level is stored
    # cumulatively ("C09") or in isolation ("09C"), which no count would reveal.
    for (k in names(lengths_seen)) {
      examples <- utils::head(unique(values[nchar(values) == as.integer(k)]), 4)
      cat("  ", k, " chars: ",
          format(as.integer(lengths_seen[[k]]), big.mark = ","),
          "  e.g. ", paste(examples, collapse = ", "), "\n", sep = "")
    }

    per_column[[cc]] <- data.frame(column = cc,
                                   nchar  = as.integer(names(lengths_seen)),
                                   n      = as.integer(lengths_seen),
                                   stringsAsFactors = FALSE)
  }

  out <- do.call(rbind, per_column)
  rownames(out) <- NULL

  if (any(out$nchar <= MB_SHORT_COLUMN_NCHAR)) {
    cat("\nMatch codes longer than", MB_SHORT_COLUMN_NCHAR, "characters against",
        "a column that holds\nfull-length codes. On a short column they match",
        "nothing, and you get an\nempty result rather than an error.\n")
  }

  invisible(out)
}
# 3. Code list against the register ----

#' Check a code list against the codes actually present in the register
#'
#' Finds codes that match nothing at all: a typo, a code retired before the
#' study period, a drug never marketed in Denmark, or a whole ATC chapter
#' missing from the extract. The Stata original notes that chapters S and V were
#' absent from its data, which is why glaucoma came back empty there - a fact
#' invisible from the result alone.
#'
#' @param lmdb The dispensing register.
#' @param codes A code list.
#' @param code_col Full ATC column name.
#' @param conditions Conditions to check. `NULL` = all.
#' @return A data frame with one row per code and the number of dispensings it
#'   matches.
#' @export
mb_check_codes <- function(lmdb, codes = mb_codelist(), code_col = "ATC",
                           conditions = NULL) {

  codes <- mb_codelist(codes, vocab = "ATC", validate = FALSE)
  if (!is.null(conditions)) {
    codes <- codes[codes$condition %in% conditions, , drop = FALSE]
  }

  # Counting distinct codes first turns a scan of the whole register into a join
  # against a few thousand rows, which is what makes this cheap enough to run
  # every time a code list changes.
  present <- dplyr::collect(dplyr::count(lmdb, .data[[code_col]], name = "n"))
  names(present)[names(present) == code_col] <- "code_full"
  present$code_full <- toupper(trimws(as.character(present$code_full)))

  codes$n_dispensings <- vapply(seq_len(nrow(codes)), function(i) {
    matches <- grepl(mb_pattern(codes$code[i]), present$code_full)
    sum(present$n[matches])
  }, numeric(1))

  codes$found <- codes$n_dispensings > 0

  n_missing <- sum(!codes$found)
  if (n_missing) {
    message(n_missing, " code(s) match nothing in the register: ",
            paste(codes$code[!codes$found], collapse = ", "))
  } else {
    message("All ", nrow(codes), " codes match at least one dispensing.")
  }

  codes[, c("condition", "code", "exclude", "n_dispensings", "found")]
}
