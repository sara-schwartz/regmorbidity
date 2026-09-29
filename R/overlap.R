# .............................................................................
# overlap.R - which conditions share codes, and what a code belongs to
#
# PURPOSE
#   A person can legitimately have two conditions at once. What is usually NOT
#   intended is for one dispensing or one diagnosis to create both of them.
#   These two functions make that visible before it reaches a result table.
#
# INPUT   a code list
# OUTPUT  a table of overlaps, or the conditions a single code belongs to
#
# WHY THIS EXISTS
#   Both overlaps we know about in this project were found by accident:
#   C02CA is half the prostate definition and sits inside hypertension's C02C,
#   and J45/J46 are the whole asthma definition and sit inside a J40-J47 COPD
#   range. Neither produces an error. Both inflate a multimorbidity count. A
#   thirty-nine condition list has more pairs than anyone will check by eye.
#
# WHAT AN OVERLAP IS NOT
#   Not automatically a bug. Prior's own lists overlap - his hypertension uses
#   all of C02, which contains the C02CA of his prostate definition. The point
#   is that it should be a decision, not a discovery.
#
# CONTENTS
#   1. Overlaps between conditions
#   2. What does this code belong to
# .............................................................................


# 1. Overlaps between conditions ----

#' Find codes shared between conditions
#'
#' Prefix matching means overlap is not only exact equality: a condition listing
#' `C02C` captures everything a condition listing `C02CA` captures. Both are
#' reported.
#'
#' @param codes A code list; see [mb_codelist()].
#' @param vocab Restrict to one vocabulary, e.g. `"ATC"`. `NULL` = both.
#' @param verbose Print a summary.
#' @return A data frame, one row per overlapping pair of codes: `vocab_id`,
#'   `condition_a`, `code_a`, `condition_b`, `code_b`, `relation`. Empty if the
#'   conditions are disjoint.
#' @keywords internal
#' @export
mb_overlap <- function(codes = mb_codelist(), vocab = NULL, verbose = TRUE) {

  codes <- mb_codelist(codes, vocab = vocab, validate = FALSE)
  # Exclusions narrow a condition rather than defining it, so an exclude row
  # sharing a code with another condition is not an overlap.
  codes <- codes[!codes$exclude, , drop = FALSE]

  found <- list()

  for (vocabulary in unique(codes$vocab_id)) {
    v <- codes[codes$vocab_id == vocabulary, , drop = FALSE]
    conditions <- sort(unique(v$condition))
    if (length(conditions) < 2L) next

    # ICD-10 is compared in WHO form, so a list that writes DJ45 in one
    # condition and J45 in another still shows the overlap. The reported codes
    # stay as written, so they can be found in the file.
    key <- if (vocabulary == "ICD10") mb_normalize_icd10(v$code) else v$code

    for (i in seq_len(length(conditions) - 1L)) {
      for (j in seq(i + 1L, length(conditions))) {
        a_rows <- which(v$condition == conditions[i])
        b_rows <- which(v$condition == conditions[j])

        for (ai in a_rows) for (bi in b_rows) {
          a <- key[ai]; b <- key[bi]
          relation <-
            if (identical(a, b)) "identical"
            else if (startsWith(b, a)) paste0(conditions[i], " covers ", conditions[j])
            else if (startsWith(a, b)) paste0(conditions[j], " covers ", conditions[i])
            else next

          found[[length(found) + 1L]] <- data.frame(
            vocab_id = vocabulary,
            condition_a = conditions[i], code_a = v$code[ai],
            condition_b = conditions[j], code_b = v$code[bi],
            relation = relation, stringsAsFactors = FALSE)
        }
      }
    }
  }

  out <- if (length(found)) do.call(rbind, found) else
    data.frame(vocab_id = character(0), condition_a = character(0),
               code_a = character(0), condition_b = character(0),
               code_b = character(0), relation = character(0),
               stringsAsFactors = FALSE)
  rownames(out) <- NULL

  if (verbose) {
    if (!nrow(out)) {
      message("No overlapping codes: every condition is disjoint.")
    } else {
      pairs <- unique(paste(out$condition_a, out$condition_b))
      message(nrow(out), " overlapping code(s) across ", length(pairs),
              " condition pair(s):")
      for (i in seq_len(nrow(out))) {
        message("  ", out$vocab_id[i], "  ", out$code_a[i],
                if (out$relation[i] == "identical") " = " else " > ",
                out$code_b[i], "   ", out$relation[i])
      }
      message("\nNot necessarily wrong - Prior's own lists overlap. ",
              "But one record should create two conditions only on purpose.")
    }
  }
  out
}


# 2. What does this code belong to ----

#' Which conditions would this code be counted as?
#'
#' The question you want answered when a prevalence looks too high. Takes a
#' code as it appears in the register - full length, and with the Danish D on
#' an ICD-10 code if that is how your data holds it.
#'
#' @param code One or more codes, e.g. `"C02CA01"` or `"DJ450"`.
#' @param codes A code list; see [mb_codelist()].
#' @return A data frame of the matching rows: `query`, `condition`, `vocab_id`,
#'   `code`, `exclude`. Empty if the code belongs to nothing.
#' @export
mb_lookup <- function(code, codes = mb_codelist()) {

  codes <- mb_codelist(codes, validate = FALSE)
  query <- toupper(gsub("[[:space:].]", "", as.character(code)))

  hits <- lapply(query, function(q) {
    # ICD-10 is compared in WHO form on both sides so DJ450 finds a list that
    # says J45, and a list that says DJ45 finds it too.
    atc <- codes$vocab_id == "ATC" & startsWith(q, codes$code)
    icd <- codes$vocab_id == "ICD10" &
           startsWith(mb_normalize_icd10(q), mb_normalize_icd10(codes$code))

    matched <- codes[atc | icd, c("condition", "vocab_id", "code", "exclude"),
                     drop = FALSE]
    if (!nrow(matched)) return(NULL)
    cbind(query = q, matched, stringsAsFactors = FALSE)
  })

  hits <- hits[!vapply(hits, is.null, logical(1))]
  out <- if (length(hits)) do.call(rbind, hits) else
    data.frame(query = character(0), condition = character(0),
               vocab_id = character(0), code = character(0),
               exclude = logical(0), stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}
