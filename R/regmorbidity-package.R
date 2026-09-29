# .............................................................................
# regmorbidity-package.R - the package-level help page, and the imports
#
# PURPOSE
#   Two jobs. It carries the overview a person sees when they type ?regmorbidity,
#   which is where a grouped list of the functions belongs - visible to users
#   rather than buried in NAMESPACE. And it holds the @importFrom tags, so
#   NAMESPACE can be generated in full and never edited by hand.
#
# CONTENTS
#   1. Package documentation
# .............................................................................


# 1. Package documentation ----

#' regmorbidity: comorbidity from Danish register data
#'
#' @description
#' Conditions are defined in editable CSV code lists - one file per condition -
#' rather than in code. The package reads those lists, extracts the conditions
#' from a dispensing register and a hospital register, and combines the two
#' halves.
#'
#' This is not an implementation of any published index. The bundled lists
#' (ATC + first-pass ICD; distress archived out of the default set) take their
#' starting point in those of Prior et al. (2016) and are meant to be revised,
#' not cited. The number of conditions comes from your CSV lists; the package API is a fixed set of exported functions. Design notes
#' (`ASSUMPTIONS_AND_LIMITATIONS.txt`, `DECISIONS.md`) live in the source
#' working tree if present; otherwise use `vignette("regmorbidity")` and this
#' help page.
#'
#' @section Recommended workflow:
#' \describe{
#'   \item{[mb_codelist()]}{load/filter lists (`conditions=`, `vocab=`)}
#'   \item{[mb_extract_diagnosis()]}{LPR onset}
#'   \item{[mb_extract_medication_batch()]}{preferred medication onset (SQL,
#'     low RAM)}
#'   \item{[mb_merge_all()]}{combine dx+rx extract directories -> long}
#'   \item{[mb_to_wide()] / [mb_count_conditions()]}{ever-after counts as of
#'     a date}
#' }
#'
#' @section Backup / special:
#' \describe{
#'   \item{[mb_extract_medication()]}{sequential; use when you need
#'     `keep_events=TRUE` or one-condition debug}
#'   \item{[mb_merge_conditions()]}{merge one condition's two data frames
#'     (in memory)}
#'   \item{[mb_load_conditions()]}{load one extract outdir to long (skips
#'     `*_all_events.rds`); medication-only studies}
#'   \item{[mb_prevalence()]}{lookback prevalence on **raw events**, not
#'     onset `.rds`}
#' }
#'
#' @section Preflight QA (optional, before long DST runs):
#' \describe{
#'   \item{[mb_check_codes()]}{for each list code, count register rows that
#'     match (prefix); catch typos / dead codes before a long extract}
#'   \item{[mb_lookup()]}{which conditions claim this code?}
#'   \item{[mb_overlap()]}{codes shared between conditions}
#' }
#'
#' @section Advanced:
#' \describe{
#'   \item{[mb_inspect_code_lengths()]}{reports string lengths of code columns
#'     in the *register* sample (e.g. is `atc2` 3 characters?). Use before
#'     sequential medication extract / prefilter debugging. Not for reviewing
#'     the code list.}
#' }
#'
#' @section Provisional:
#' \describe{
#'   \item{[mb_apply_exclusions()]}{optional stage-2; incomplete vs Prior;
#'     not part of the recommended workflow. HTN still has C03 / HF / CKD gaps.}
#' }
#'
#' @section Authors:
#' Jie Zhang and Sara Schwartz (saras@@clin.au.dk).
#'
#' @importFrom rlang .data
#' @importFrom rlang :=
#' @importFrom data.table .N
#' @importFrom dplyr %>%
#' @keywords internal
"_PACKAGE"
