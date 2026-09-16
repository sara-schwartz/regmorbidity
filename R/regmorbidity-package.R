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
#'     or a data frame}
#'   \item{[mb_validate_codelist()]}{refuse lists that would fail silently}
#'   \item{[mb_write_codelist()]}{split a list into one CSV per condition, for
#'     review by someone who does not read R}
#' }
#'
#' @section Checks to run first:
#' \describe{
#'   \item{[mb_inspect_codes()]}{how long are the code columns really - the
#'     check that catches a 4-character pattern aimed at a 3-character column}
#'   \item{[mb_check_codes()]}{do these codes match anything at all, which is
#'     how a true zero is told from a typo}
#'   \item{[mb_overlap()]}{which conditions share codes}
#'   \item{[mb_lookup()]}{what would this one record be counted as}
#'   \item{[mb_compare()]}{what changed between two versions of a list}
#' }
#'
#' @section Extraction:
#' \describe{
#'   \item{[mb_extract_medication()]}{the ATC half, one condition at a time,
#'     checkpointed so an interrupted run resumes}
#'   \item{[mb_extract_diagnosis()]}{the ICD-10 half, same contract}
#'   \item{[mb_flag_users()]}{the prescription rule on its own, for records you
#'     have already filtered}
#'   \item{[mb_normalize_icd10()]}{`DI50` and `I50` are the same code}
#' }
#'
#' @section Putting it together:
#' \describe{
#'   \item{[mb_merge_conditions()]}{combine the two halves of one condition,
#'     OR or AND}
#'   \item{[mb_merge_all()]}{every condition, from the two output folders}
#'   \item{[mb_condition_logic()]}{which rule a condition uses}
#'   \item{[mb_load_conditions()]}{read the per-condition files back}
#'   \item{[mb_to_wide()]}{one column per condition}
#'   \item{[mb_count_conditions()]}{`n_conditions` and `multimorbid`}
#' }
#'
#' @section Authors:
#' Jie Zhang and Sara Schwartz (saras@@clin.au.dk).
#'
#' @importFrom rlang .data
#' @importFrom data.table .N
#' @keywords internal
"_PACKAGE"
