# .............................................................................
# codelist.R - reading, checking and writing condition code lists
#
# PURPOSE
#   A code list says which register codes define which condition, and under what
#   rule. Keeping that in editable CSV rather than in R code is the whole point
#   of this package: a clinician can review one file per disease without reading
#   any code, and a reviewer can diff the definitions between two studies.
#
# INPUT   .csv files (one combined table, or one file per condition), or a
#         data frame the caller built themselves.
# OUTPUT  A validated code list data frame, used by mb_extract_medication().
#
# THE FORMAT
#   One row = one code prefix belonging to one condition.
#
#     condition          key, e.g. "hypertension". Must be identical across
#                        vocabularies, because it is what joins the ATC and
#                        ICD-10 halves of the same disease.
#     condition_label    human-readable name. Defaults to `condition`.
#     category           optional grouping, e.g. "Circulatory".
#     vocab_id           "ATC" or "ICD10".
#     code               code PREFIX. "C09" matches C09AA05, C09DB01, ...
#     exclude            TRUE = remove matches of this code. Default FALSE.
#     min_prescriptions  dispensings required before the condition counts.
#                        ATC rows only; meaningless for a diagnosis.
#     window_days        they must fall within this many days. NA = no window.
#     logic              how the diagnosis and medication halves combine:
#                        "OR" (either is enough) or "AND" (both required).
#                        Blank means OR. Prior uses AND only for epilepsy.
#     note               free text.
#
#   ATC and ICD-10 rows sit side by side in the same table so the diagnosis half
#   can be added later without changing this file or the format collaborators
#   have already learned.
#
# WHY min_prescriptions AND window_days LIVE ON EVERY ROW
#   They describe the CONDITION, not the individual code, so they are repeated
#   and must agree. Storing them per row keeps one condition in one file with no
#   second lookup table to maintain - at the cost of a consistency check, which
#   mb_validate_codelist() performs.
#
# SEE ALSO
#   ASSUMPTIONS_AND_LIMITATIONS.txt sections 1 and 3 for where the bundled
#   values come from and which of them are decisions rather than transcriptions.
#
# CONTENTS
#   1. Constants
#   2. Loading
#   3. Normalising
#   4. Validating
#        4.1 Vocabulary
#        4.2 ATC code shape
#        4.3 One rule per condition
#        4.4 Duplicates
#        4.5 Redundancy (warning, not an error)
#   5. Writing
# .............................................................................
# 1. Constants ----

#' Columns a code list is written out with, in order.
#' @keywords internal
MB_CODELIST_COLS <- c("condition", "condition_label", "category", "vocab_id",
                      "code", "exclude", "min_prescriptions", "window_days",
                      "logic", "note")

# The codebook writes its rule as a word ("last_year") rather than a number.
# Named here so the translation is visible and editable in one place.
MB_TIME_FRAME_DAYS <- c(last_year = 365L, "1y" = 365L, "365" = 365L)

# A condition whose file states no rule still has to require something, or every
# person who ever touched the drug once would count.
MB_DEFAULT_MIN_PRESCRIPTIONS <- 1L

# An ATC code is 1-7 characters (WHO ATC levels 1-5). Anything outside that is a
# typo, not a code. See ASSUMPTIONS_AND_LIMITATIONS section 6.
MB_ATC_MIN_NCHAR <- 1L
MB_ATC_MAX_NCHAR <- 7L
# 2. Loading ----

#' Load a code list
#'
#' @param x One of:
#'   * `NULL` - the ATC code lists bundled with the package. A starting point,
#'     not an authority: revise them for your own study.
#'   * a data frame - used as-is.
#'   * a path to a `.csv` file - one table holding every condition.
#'   * a path to a directory - every `.csv` in it is read and stacked. Use this
#'     to keep one file per disease. If a file has no `condition` column, the
#'     file name is used as the condition (`hypertension.csv` -> "hypertension").
#' @param vocab Optionally keep only one vocabulary, e.g. `"ATC"`.
#' @param validate Run [mb_validate_codelist()] before returning. Leave on
#'   unless you are deliberately inspecting a broken list.
#' @return A data frame with the columns described at the top of this file.
#' @export
mb_codelist <- function(x = NULL, vocab = NULL, validate = TRUE) {

  if (is.null(x)) {
    x <- system.file("extdata", "codelists", package = "regmorbidity")
    if (!nzchar(x) || !dir.exists(x)) {
      stop("The bundled code lists were not found. Pass a file or directory ",
           "to `mb_codelist()` explicitly, e.g. mb_codelist(\"my_codelists/\").",
           call. = FALSE)
    }
  }

  # Directory before file: a folder of per-disease CSVs is the format we expect
  # collaborators to hand back, so it should be the cheapest thing to pass in.
  if (is.data.frame(x)) {
    out <- mb_normalize_codelist(x, source_label = "<data frame>")
  } else if (length(x) == 1L && dir.exists(x)) {
    out <- mb_read_codelist_dir(x)
  } else if (length(x) == 1L && file.exists(x)) {
    out <- mb_normalize_codelist(mb_read_csv(x), source_label = basename(x))
  } else {
    stop("`x` is neither a data frame nor an existing file or directory: ",
         paste(x, collapse = ", "), call. = FALSE)
  }

  if (!is.null(vocab)) {
    keep <- toupper(out$vocab_id) %in% toupper(vocab)
    # Silently returning zero rows here would look like "no conditions defined"
    # much later, in the extraction, with nothing pointing back to this line.
    if (!any(keep)) {
      stop("No rows left after filtering to vocab_id = ",
           paste(vocab, collapse = "/"), ". Present: ",
           paste(sort(unique(out$vocab_id)), collapse = ", "), call. = FALSE)
    }
    out <- out[keep, , drop = FALSE]
  }

  if (validate) mb_validate_codelist(out)
  out
}


#' Read every .csv in a directory and stack them into one code list
#'
#' One file per disease is the reviewable unit: small enough to send to a
#' clinician, and it makes a change to one condition a one-file diff.
#' @keywords internal
mb_read_codelist_dir <- function(dir) {

  files <- list.files(dir, pattern = "\\.csv$", full.names = TRUE,
                      ignore.case = TRUE)
  if (!length(files)) {
    stop("No .csv files found in directory: ", dir, call. = FALSE)
  }

  parts <- lapply(files, function(f) {
    df <- mb_read_csv(f)
    if (!nrow(df)) return(NULL)

    # A per-disease file may leave the condition implicit in its name, so
    # somebody can add a condition by dropping in "gigt.csv" with two columns
    # and no boilerplate.
    condition_missing <-
      !"condition" %in% names(df) ||
      all(is.na(df$condition)) ||
      all(!nzchar(trimws(as.character(df$condition))))

    if (condition_missing) {
      df$condition <- tools::file_path_sans_ext(basename(f))
    }
    mb_normalize_codelist(df, source_label = basename(f))
  })

  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (!length(parts)) stop("Every .csv in ", dir, " was empty.", call. = FALSE)

  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}


#' Read one CSV with everything as character
#'
#' Codes are read as text on purpose. Left to guess, R turns a column of ATC
#' codes into something else the moment one of them looks numeric, and leading
#' zeros in ICD-10 codes are lost silently.
#' @keywords internal
mb_read_csv <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, colClasses = "character",
                  na.strings = c("NA", ""), check.names = TRUE)
}
# 3. Normalising ----

#' Fill in defaults, accept older column names, and coerce types
#'
#' Everything that makes a hand-edited file usable happens here, so the rest of
#' the package can assume one exact shape.
#' @keywords internal
mb_normalize_codelist <- function(df, source_label = "<code list>") {

  df <- as.data.frame(df, stringsAsFactors = FALSE)
  names(df) <- tolower(names(df))

  # Accept `n_prescriptions` as a synonym for min_prescriptions.
  if (!"min_prescriptions" %in% names(df) && "n_prescriptions" %in% names(df)) {
    df$min_prescriptions <- df$n_prescriptions
  }

  # Accept the `time_frame` shorthand from codes_multimorbidity_dk_ss.csv, so
  # the existing codebook can be loaded without being rewritten first.
  if (!"window_days" %in% names(df) && "time_frame" %in% names(df)) {
    time_frame <- tolower(trimws(as.character(df$time_frame)))
    df$window_days <- unname(MB_TIME_FRAME_DAYS[time_frame])
  }

  missing_cols <- setdiff(c("vocab_id", "code"), names(df))
  if (length(missing_cols)) {
    stop(source_label, ": missing required column(s): ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if (!"condition" %in% names(df)) {
    stop(source_label, ": no `condition` column, and the condition could not ",
         "be taken from the file name.", call. = FALSE)
  }

  as_trimmed_chr <- function(v) trimws(as.character(v))

  has <- function(col) col %in% names(df)

  out <- data.frame(
    condition       = as_trimmed_chr(df$condition),
    condition_label = if (has("condition_label")) as_trimmed_chr(df$condition_label) else NA_character_,
    category        = if (has("category")) as_trimmed_chr(df$category) else NA_character_,
    vocab_id        = toupper(as_trimmed_chr(df$vocab_id)),
    # Whitespace inside a code is always a typo, and one that would otherwise
    # make the code match nothing at all without any error.
    code            = toupper(gsub("[[:space:]]", "", as_trimmed_chr(df$code))),
    exclude         = if (has("exclude")) mb_as_logical(df$exclude) else FALSE,
    min_prescriptions = if (has("min_prescriptions")) suppressWarnings(as.integer(df$min_prescriptions)) else NA_integer_,
    window_days     = if (has("window_days")) suppressWarnings(as.integer(df$window_days)) else NA_integer_,
    logic           = if (has("logic")) toupper(as_trimmed_chr(df$logic)) else NA_character_,
    note            = if (has("note")) as_trimmed_chr(df$note)
                      else if (has("coding_definition")) as_trimmed_chr(df$coding_definition)
                      else NA_character_,
    stringsAsFactors = FALSE
  )

  label_blank <- is.na(out$condition_label) | !nzchar(out$condition_label)
  out$condition_label[label_blank] <- out$condition[label_blank]

  out$exclude[is.na(out$exclude)] <- FALSE
  out$min_prescriptions[is.na(out$min_prescriptions)] <- MB_DEFAULT_MIN_PRESCRIPTIONS

  # Blank rows are what a spreadsheet leaves behind after someone deletes a
  # code. Dropping them here keeps them out of every downstream count.
  out[out$vocab_id != "" & !is.na(out$code) & nzchar(out$code), , drop = FALSE]
}


#' Read the many ways a person writes TRUE in a spreadsheet
#'
#' Blank means FALSE: an empty `exclude` cell is the normal case, not a mistake.
#' @keywords internal
mb_as_logical <- function(v) {
  if (is.logical(v)) return(v)
  v <- tolower(trimws(as.character(v)))
  out <- rep(NA, length(v))
  out[v %in% c("true", "t", "yes", "y", "1")] <- TRUE
  out[v %in% c("false", "f", "no", "n", "0")] <- FALSE
  out[is.na(v) | v == ""] <- FALSE
  as.logical(out)
}
# 4. Validating ----

#' Check a code list for the mistakes that fail silently
#'
#' Deliberately narrow: it looks for errors that produce a plausible but wrong
#' result rather than a crash. A condition whose own rows disagree about the
#' rule, a condition made only of exclusions, duplicated codes, a code another
#' code already covers. Every one of those yields a number that looks fine.
#'
#' It does NOT check that the codes are clinically right. Use `mb_check_codes()`
#' against the register for whether a code matches anything at all, and a human
#' for whether it should.
#'
#' @param codes A code list.
#' @param quiet Suppress the "looks fine" message.
#' @return `codes`, invisibly. Stops on any problem found.
#' @keywords internal
mb_validate_codelist <- function(codes, quiet = FALSE) {

  stopifnot(is.data.frame(codes))

  missing_cols <- setdiff(c("condition", "vocab_id", "code", "exclude",
                            "min_prescriptions", "window_days"), names(codes))
  if (length(missing_cols)) {
    stop("Code list is missing column(s): ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  # Collected rather than thrown one at a time, so a person fixing a hand-edited
  # file sees every problem in one pass instead of rerunning after each fix.
  problems <- character(0)
  ## 4.1 Vocabulary ----
  bad_vocab <- setdiff(unique(codes$vocab_id), c("ATC", "ICD10"))
  if (length(bad_vocab)) {
    problems <- c(problems, paste0(
      "Unknown vocab_id: ", paste(bad_vocab, collapse = ", "),
      " (expected ATC or ICD10)"))
  }
  ## 4.2 ATC code shape ----
  atc <- codes[codes$vocab_id == "ATC", , drop = FALSE]
  if (nrow(atc)) {
    bad_len <- unique(atc$code[nchar(atc$code) < MB_ATC_MIN_NCHAR |
                               nchar(atc$code) > MB_ATC_MAX_NCHAR])
    if (length(bad_len)) {
      problems <- c(problems, paste0(
        "ATC codes must be ", MB_ATC_MIN_NCHAR, "-", MB_ATC_MAX_NCHAR,
        " characters: ", paste(bad_len, collapse = ", ")))
    }
    bad_chr <- unique(atc$code[grepl("[^A-Z0-9]", atc$code)])
    if (length(bad_chr)) {
      problems <- c(problems, paste0(
        "ATC codes contain unexpected characters: ",
        paste(bad_chr, collapse = ", ")))
    }
  }
  ## 4.3 One rule per condition ----
  # The rule belongs to the condition, so two different answers means one of the
  # rows is wrong. Picking either one silently would be a coin flip that changes
  # prevalence, so this is an error rather than a warning.
  by_condition <- split(codes, codes$condition)

  for (nm in names(by_condition)) {
    d <- by_condition[[nm]]

    if (length(unique(d$min_prescriptions)) > 1L) {
      problems <- c(problems, paste0(
        "Condition '", nm, "' has conflicting min_prescriptions: ",
        paste(sort(unique(d$min_prescriptions)), collapse = ", ")))
    }

    # Mixing NA with a number is also a conflict: "no window" and "365 days" are
    # different rules, not a missing value to be filled in.
    window_conflict <-
      length(unique(d$window_days[!is.na(d$window_days)])) > 1L ||
      (any(is.na(d$window_days)) && any(!is.na(d$window_days)))
    if (window_conflict) {
      problems <- c(problems, paste0(
        "Condition '", nm, "' has conflicting window_days: ",
        paste(unique(d$window_days), collapse = ", ")))
    }

    stated_logic <- unique(d$logic[!is.na(d$logic) & nzchar(d$logic)])
    if (length(stated_logic) > 1L) {
      problems <- c(problems, paste0(
        "Condition '", nm, "' has conflicting logic: ",
        paste(stated_logic, collapse = ", ")))
    }
    if (length(stated_logic) && !all(stated_logic %in% c("OR", "AND"))) {
      problems <- c(problems, paste0(
        "Condition '", nm, "' has logic '", stated_logic[1],
        "' (expected OR or AND)"))
    }

    if (nrow(d) && all(d$exclude)) {
      problems <- c(problems, paste0(
        "Condition '", nm, "' consists only of exclude = TRUE rows, so it can ",
        "never match anything."))
    }
  }
  ## 4.4 Duplicates ----
  dup <- duplicated(codes[, c("condition", "vocab_id", "code", "exclude")])
  if (any(dup)) {
    problems <- c(problems, paste0(
      sum(dup), " duplicated condition/code row(s), e.g. ",
      paste(utils::head(paste0(codes$condition[dup], ":", codes$code[dup]), 3),
            collapse = ", ")))
  }

  if (length(problems)) {
    stop("Code list problems:\n  - ", paste(problems, collapse = "\n  - "),
         call. = FALSE)
  }
  ## 4.5 Redundancy (warning, not an error) ----
  # "C09" already covers "C09AA05", so listing both changes nothing. Harmless in
  # itself, but usually the visible half of a paste error, so it is worth
  # saying out loud without refusing to run.
  for (nm in names(by_condition)) {
    include_codes <- unique(by_condition[[nm]]$code[!by_condition[[nm]]$exclude])
    for (code in include_codes) {
      shadowed <- setdiff(
        include_codes[startsWith(include_codes, code) &
                      nchar(include_codes) > nchar(code)],
        code)
      if (length(shadowed)) {
        warning("Condition '", nm, "': code '", code, "' already covers ",
                paste(shadowed, collapse = ", "),
                " - the longer code(s) add nothing.", call. = FALSE)
      }
    }
  }

  if (!quiet) {
    message("Code list OK: ", nrow(codes), " rows, ",
            length(unique(codes$condition)), " conditions, vocab ",
            paste(sort(unique(codes$vocab_id)), collapse = "/"), ".")
  }
  invisible(codes)
}
# 5. Writing ----

#' Write a code list out as one .csv per condition
#'
#' The point is review: one small file per disease can be sent to one clinician,
#' and a change to a definition shows up as a one-file diff rather than a moved
#' row in a 40-condition spreadsheet. Read the folder back with
#' `mb_codelist("<dir>")`.
#'
#' @param codes A code list.
#' @param dir Directory to write into. Created if needed.
#' @param overwrite Overwrite existing files of the same name.
#' @return The file paths written, invisibly.
#' @export
mb_write_codelist <- function(codes, dir, overwrite = FALSE) {

  codes <- mb_codelist(codes, validate = FALSE)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

  conditions <- unique(codes$condition)
  paths <- file.path(dir, paste0(conditions, ".csv"))

  # Refusing by default matters here: these files are the ones a collaborator
  # may have spent an afternoon editing, and overwriting them is unrecoverable.
  already_there <- file.exists(paths)
  if (any(already_there) && !overwrite) {
    stop("These files already exist (use overwrite = TRUE): ",
         paste(basename(paths[already_there]), collapse = ", "), call. = FALSE)
  }

  for (i in seq_along(conditions)) {
    one_condition <- codes[codes$condition == conditions[i],
                           MB_CODELIST_COLS, drop = FALSE]
    utils::write.csv(one_condition, paths[i], row.names = FALSE, na = "")
  }

  message("Wrote ", length(paths), " condition file(s) to ", dir)
  invisible(paths)
}
