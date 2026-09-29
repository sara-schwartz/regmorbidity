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
#' This is not an implementation of any published index. The bundled lists take
#' their starting point in those of Prior et al. (2016) and are meant to be
#' revised, not cited. Read `ASSUMPTIONS_AND_LIMITATIONS.txt` before reporting
#' any number from it.
#'
#' @section Code lists:
#' \describe{
#'   \item{[mb_codelist()]}{load from a folder of per-condition CSVs, one file,
#'     or a data frame (validates on load)}
#' }
#'
#' @section Checks to run first:
#' \describe{
#'   \item{[mb_inspect_codes()]}{how long are the code columns really - the
#'     check that catches a 4-character pattern aimed at a 3-character column}
#'   \item{[mb_check_codes()]}{do these codes match anything at all, which is
#'     how a true zero is told from a typo}
#'   \item{[mb_lookup()]}{which conditions a register code would count as}
#'   \item{[mb_overlap()]}{codes shared between conditions (QA before counting)}
#' }
#'
#' @section Extraction:
#' \describe{
#'   \item{[mb_extract_medication_batch()]}{recommended ATC path: every
#'     condition sharing one prescription rule in one query, collecting only
#'     final onset rows}
#'   \item{[mb_extract_medication()]}{sequential backup - one condition at a
#'     time, checkpointed per condition; use when you need
#'     `keep_events` or finer resume}
#'   \item{[mb_extract_diagnosis()]}{the ICD-10 half, same contract as the
#'     sequential medication extractor}
#' }
#'
#' @section Downstream (one path):
#' \describe{
#'   \item{[mb_merge_conditions()]}{combine the two halves of one condition
#'     (OR or AND)}
#'   \item{[mb_merge_all()]}{every condition, from the two output folders}
#'   \item{[mb_load_conditions()]}{read the per-condition `.rds` files back}
#'   \item{[mb_to_wide()]}{one column per condition}
#'   \item{[mb_count_conditions()]}{`n_conditions` and `multimorbid`}
#'   \item{[mb_prevalence()]}{condition status at a point in time, with a
#'     lookback window, rather than ever-after-onset}
#'   \item{[mb_apply_exclusions()]}{optional stage-2 Prior exclusion rules on
#'     an assembled long or wide table (provisional; not baked into extract)}
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
