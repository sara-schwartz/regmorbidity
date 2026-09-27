# .............................................................................
# merge.R - combining the diagnosis and medication halves of a condition
#
# PURPOSE
#   Most conditions are defined as "diagnosis AND/OR prescription". Once both
#   halves have been extracted, they have to be reduced to one date per person.
#
# INPUT   two per-condition tables, from mb_extract_diagnosis() and
#         mb_extract_medication()
# OUTPUT  one row per person with the combined onset date and where it came from
#
# THE TWO RULES
#   OR   the condition starts at whichever came first. Prior uses this for every
#        condition except one.
#   AND  both are required, and the condition starts at the later of the two -
#        the point at which both criteria are finally met. Prior uses this only
#        for epilepsy: "Diagnosis AND prescription of anti-epileptics".
#
#   This follows merge_condition() in Jie Zhang's streamlineMM.R, which is the
#   reference implementation.
#
# WHAT THIS DOES NOT DO
#   Prior's exclusion rules - dyslipidaemia only if not IHD, hypertension only
#   if not IHD or heart failure, distress only if no other mental disorder -
#   reference OTHER conditions and so cannot be applied here, one condition at a
#   time. They need a pass over the assembled table. Not implemented; see
#   TODO.txt section 1(b).
#
# CONTENTS
#   1. Merging one condition
#   2. Merging everything on disk
# .............................................................................


# 1. Merging one condition ----

#' Combine the diagnosis and medication halves of one condition
#'
#' @param diagnosis Output of [mb_extract_diagnosis()] for one condition, or `NULL`
#'   if the condition has no diagnosis half.
#' @param medication Output of [mb_extract_medication()] for one condition, or
#'   `NULL` if the condition is diagnosis-only.
#' @param logic `"OR"` (either half is enough, earliest date wins) or `"AND"`
#'   (both required, later date wins).
#' @param id_col Person id column name.
#' @param condition Condition name for the output. Taken from the inputs when
#'   not given.
#' @return One row per person: id, `condition`, `source`, `onset_date`,
#'   `icd_date`, `rx_date`, `n_icd`, `n_rx`. `source` is `"diagnosis"`,
#'   `"medication"` or `"both"`. People who do not meet the rule are dropped.
#'
#'   `n_icd` and `n_rx` count matching records over the whole period, not within
#'   the window - they answer "how much of this did they have", not "how many
#'   met the rule". With de-duplication on, both count days rather than rows.
#' @export
mb_merge_conditions <- function(diagnosis, medication, logic = "OR",
                                id_col = "pnr", condition = NULL) {

  logic <- toupper(logic)
  if (!logic %in% c("OR", "AND")) {
    stop("`logic` must be \"OR\" or \"AND\", not \"", logic, "\".", call. = FALSE)
  }

  if (is.null(diagnosis) && is.null(medication)) {
    stop("Both halves are NULL - nothing to merge.", call. = FALSE)
  }
  if (is.null(condition)) {
    condition <- unique(c(diagnosis$condition, medication$condition))[1]
  }

  # AND with a missing half can never be satisfied. Returning an empty table
  # rather than erroring keeps a batch running, but it is worth saying out loud
  # because it almost always means a half failed to extract.
  if (logic == "AND" && (is.null(diagnosis) || is.null(medication) ||
                         !nrow(as.data.frame(diagnosis)) ||
                         !nrow(as.data.frame(medication)))) {
    warning("Condition '", condition, "' uses AND but one half is empty - ",
            "no one can meet it.", call. = FALSE)
    return(mb_empty_merge(id_col))
  }

  icd <- mb_half(diagnosis,  id_col, "icd_date")
  rx  <- mb_half(medication, id_col, "rx_date")

  ids <- union(icd[[id_col]], rx[[id_col]])
  if (!length(ids)) return(mb_empty_merge(id_col))

  icd_i <- match(ids, icd[[id_col]])
  rx_i  <- match(ids, rx[[id_col]])
  icd_date <- icd$icd_date[icd_i]
  rx_date  <- rx$rx_date[rx_i]
  n_icd    <- icd$n[icd_i]
  n_rx     <- rx$n[rx_i]

  onset <- if (logic == "OR") {
    # pmin drops NA only when told to; with both missing it must stay missing.
    suppressWarnings(pmin(icd_date, rx_date, na.rm = TRUE))
  } else {
    ifelse(!is.na(icd_date) & !is.na(rx_date),
           pmax(icd_date, rx_date), NA)
  }
  onset <- as.Date(onset, origin = "1970-01-01")

  source <- ifelse(!is.na(icd_date) & !is.na(rx_date), "both",
            ifelse(!is.na(icd_date), "diagnosis", "medication"))

  out <- data.frame(ids, condition, source, onset, icd_date, rx_date,
                    n_icd, n_rx, stringsAsFactors = FALSE)
  names(out) <- c(id_col, "condition", "source", "onset_date",
                  "icd_date", "rx_date", "n_icd", "n_rx")

  out <- out[!is.na(out$onset_date), , drop = FALSE]
  rownames(out) <- NULL
  out
}


#' One half, reduced to id + date, tolerating NULL
#' @keywords internal
mb_half <- function(x, id_col, date_name) {
  empty <- data.frame(character(0), as.Date(character(0)), integer(0),
                      stringsAsFactors = FALSE)
  names(empty) <- c(id_col, date_name, "n")
  if (is.null(x)) return(empty)

  x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (!nrow(x)) return(empty)
  if (!id_col %in% names(x)) {
    stop("Column '", id_col, "' not found in one of the halves.", call. = FALSE)
  }
  if (!"onset_date" %in% names(x)) {
    stop("Column 'onset_date' not found in one of the halves.", call. = FALSE)
  }

  out <- x[, c(id_col, "onset_date")]
  names(out)[2] <- date_name
  out$n <- if ("n_records" %in% names(x)) x$n_records else NA_integer_
  out
}


#' @keywords internal
mb_empty_merge <- function(id_col) {
  out <- data.frame(character(0), character(0), character(0),
                    as.Date(character(0)), as.Date(character(0)),
                    as.Date(character(0)), integer(0), integer(0),
                    stringsAsFactors = FALSE)
  names(out) <- c(id_col, "condition", "source", "onset_date",
                  "icd_date", "rx_date", "n_icd", "n_rx")
  out
}


#' Which rule applies to a condition
#'
#' Read from an optional `logic` column in the code list, defaulting to `"OR"`.
#' Keeping it in the code list means the one condition that behaves differently
#' says so in its own file, rather than in a special case in the code.
#'
#' @param codes A code list.
#' @param condition Condition name.
#' @return `"OR"` or `"AND"`.
#' @export
mb_condition_logic <- function(codes, condition) {
  rows <- codes[codes$condition == condition, , drop = FALSE]
  if (!nrow(rows)) {
    stop("Condition '", condition, "' is not in the code list.", call. = FALSE)
  }
  if (!"logic" %in% names(rows)) return("OR")

  logic <- toupper(unique(rows$logic[!is.na(rows$logic) & nzchar(rows$logic)]))
  if (!length(logic)) return("OR")
  if (length(logic) > 1L) {
    stop("Condition '", condition, "' has conflicting logic: ",
         paste(logic, collapse = ", "), call. = FALSE)
  }
  logic
}


# 2. Merging everything on disk ----

#' Merge every condition from the two extraction directories
#'
#' @param diagnosis_dir,medication_dir The `outdir`s used for
#'   [mb_extract_diagnosis()] and [mb_extract_medication()]. Either may be `NULL`.
#' @param codes A code list, used for the per-condition `logic` column.
#' @param conditions Conditions to merge. `NULL` = every condition present in
#'   either directory.
#' @param id_col Person id column name.
#' @param verbose Print one line per condition.
#' @return One long data frame, one row per person-condition.
#' @export
mb_merge_all <- function(diagnosis_dir = NULL, medication_dir = NULL,
                         codes = mb_codelist(), conditions = NULL,
                         id_col = "pnr", verbose = TRUE) {

  if (is.null(diagnosis_dir) && is.null(medication_dir)) {
    stop("Give at least one of diagnosis_dir and medication_dir.", call. = FALSE)
  }

  present <- function(dir) {
    if (is.null(dir) || !dir.exists(dir)) return(character(0))
    tools::file_path_sans_ext(mb_condition_files(dir, full.names = FALSE))
  }

  dx_have <- present(diagnosis_dir)
  rx_have <- present(medication_dir)

  if (is.null(conditions)) conditions <- sort(union(dx_have, rx_have))
  if (!length(conditions)) {
    stop("No extracted conditions found in either directory.", call. = FALSE)
  }

  read_one <- function(dir, condition, have) {
    if (!condition %in% have) return(NULL)
    readRDS(file.path(dir, paste0(condition, ".rds")))
  }

  parts <- lapply(conditions, function(condition) {
    dx <- read_one(diagnosis_dir,  condition, dx_have)
    rx <- read_one(medication_dir, condition, rx_have)
    if (is.null(dx) && is.null(rx)) return(NULL)

    logic <- tryCatch(mb_condition_logic(codes, condition),
                      error = function(e) "OR")

    merged <- mb_merge_conditions(dx, rx, logic = logic, id_col = id_col,
                                  condition = condition)
    if (verbose) {
      message(sprintf("  %-22s %-4s %s persons", condition, logic,
                      format(nrow(merged), big.mark = ",")))
    }
    if (!nrow(merged)) NULL else merged
  })

  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (!length(parts)) {
    warning("No one met any condition.", call. = FALSE)
    return(NULL)
  }

  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}
