# .............................................................................
# compare.R - what changed between two versions of a code list
#
# PURPOSE
#   Code lists are the study definition. When one changes, the change IS the
#   methods paragraph, and it has to be visible. A CSV diff shows moved rows and
#   reformatted numbers; this shows what actually changed about the conditions.
#
# INPUT   two code lists, in any form mb_codelist() accepts
# OUTPUT  one row per change, plus a printed summary
#
# WHY NOT JUST USE git diff
#   Because both sides are normalised first. "c09" and "C09 " and "C09" are the
#   same code, an empty exclude cell and FALSE are the same rule, and a condition
#   whose rows were reordered has not changed at all. A text diff reports all
#   three as changes and buries the one that matters.
#
# CONTENTS
#   1. Comparing
#        1.1 Conditions added and removed
#        1.2 Codes added and removed
#        1.3 Rules
#   2. Summarising
# .............................................................................


# 1. Comparing ----

# Condition-level fields that change what the extraction DOES.
MB_RULE_FIELDS <- c("min_prescriptions", "window_days", "logic")

# Condition-level fields that only change how it reads.
MB_LABEL_FIELDS <- c("condition_label", "category")


#' What changed between two code lists
#'
#' Both sides are normalised through [mb_codelist()] first, so reordering,
#' whitespace, letter case and a blank-versus-`FALSE` exclude cell are not
#' reported as changes. What is reported is a change to the definition.
#'
#' @param old,new Two code lists. Either may be a directory of per-condition
#'   CSVs, a single CSV, or a data frame - anything [mb_codelist()] accepts.
#' @param include_labels Also report changes to `condition_label` and
#'   `category`, which do not affect what is extracted. Off by default so the
#'   result is a list of things that change results.
#' @param verbose Print a summary.
#' @return A data frame with one row per change: `condition`, `vocab_id`,
#'   `change`, `item`, `old`, `new`. Zero rows means the two are equivalent.
#' @keywords internal
mb_compare <- function(old, new, include_labels = FALSE, verbose = TRUE) {

  old <- mb_codelist(old, validate = FALSE)
  new <- mb_codelist(new, validate = FALSE)

  changes <- list()
  add <- function(...) changes[[length(changes) + 1L]] <<- mb_change_row(...)

  old_conditions <- unique(old$condition)
  new_conditions <- unique(new$condition)

  ## 1.1 Conditions added and removed ----
  for (cond in setdiff(new_conditions, old_conditions)) {
    d <- new[new$condition == cond, , drop = FALSE]
    add(cond, paste(sort(unique(d$vocab_id)), collapse = "/"),
        "condition added", "", "",
        paste(d$code, collapse = " "))
  }
  for (cond in setdiff(old_conditions, new_conditions)) {
    d <- old[old$condition == cond, , drop = FALSE]
    add(cond, paste(sort(unique(d$vocab_id)), collapse = "/"),
        "condition removed", "",
        paste(d$code, collapse = " "), "")
  }

  ## 1.2 Codes added and removed ----
  # Keyed on condition + vocabulary + code, so a code that flips between
  # include and exclude is reported as a changed rule rather than as a removal
  # and an unrelated addition.
  for (cond in intersect(old_conditions, new_conditions)) {
    o <- old[old$condition == cond, , drop = FALSE]
    n <- new[new$condition == cond, , drop = FALSE]

    for (vocab in union(o$vocab_id, n$vocab_id)) {
      ov <- o[o$vocab_id == vocab, , drop = FALSE]
      nv <- n[n$vocab_id == vocab, , drop = FALSE]

      for (code in setdiff(nv$code, ov$code)) {
        add(cond, vocab, "code added", code, "", code)
      }
      for (code in setdiff(ov$code, nv$code)) {
        add(cond, vocab, "code removed", code, code, "")
      }
      for (code in intersect(ov$code, nv$code)) {
        o_excl <- ov$exclude[match(code, ov$code)]
        n_excl <- nv$exclude[match(code, nv$code)]
        if (!identical(o_excl, n_excl)) {
          add(cond, vocab, "exclude changed", code,
              tolower(as.character(o_excl)), tolower(as.character(n_excl)))
        }
      }
    }

    ## 1.3 Rules ----
    fields <- if (include_labels) c(MB_RULE_FIELDS, MB_LABEL_FIELDS)
              else MB_RULE_FIELDS

    for (field in fields) {
      o_val <- mb_single(o[[field]])
      n_val <- mb_single(n[[field]])
      if (!identical(o_val, n_val)) {
        change <- if (field %in% MB_RULE_FIELDS) "rule changed" else "label changed"
        add(cond, paste(sort(unique(n$vocab_id)), collapse = "/"),
            change, field, o_val, n_val)
      }
    }
  }

  out <- if (length(changes)) do.call(rbind, changes) else mb_change_row()[0, ]
  rownames(out) <- NULL
  out <- out[order(out$condition, out$change, out$item), , drop = FALSE]
  rownames(out) <- NULL

  if (verbose) mb_report_compare(out, old, new)
  out
}


#' One row of the change table
#' @keywords internal
mb_change_row <- function(condition = character(0), vocab_id = character(0),
                          change = character(0), item = character(0),
                          old = character(0), new = character(0)) {
  data.frame(condition = condition, vocab_id = vocab_id, change = change,
             item = item, old = old, new = new, stringsAsFactors = FALSE)
}


#' The one value a per-condition field should have, as text
#'
#' Returns "" for a field that is unset, and joins with "/" if a condition
#' disagrees with itself - which mb_validate_codelist() would reject, but
#' mb_compare() must survive, since comparing a broken list against a fixed one
#' is exactly when you need it.
#' @keywords internal
mb_single <- function(x) {
  if (is.null(x)) return("")
  x <- x[!is.na(x)]
  if (is.character(x)) x <- x[nzchar(x)]
  if (!length(x)) return("")
  paste(sort(unique(as.character(x))), collapse = "/")
}


# 2. Summarising ----

#' @keywords internal
mb_report_compare <- function(changes, old, new) {

  if (!nrow(changes)) {
    message("No differences: the two code lists are equivalent.")
    return(invisible(NULL))
  }

  n_cond <- length(unique(changes$condition))
  message(nrow(changes), " change(s) across ", n_cond, " condition(s). ",
          length(unique(old$condition)), " -> ",
          length(unique(new$condition)), " conditions, ",
          nrow(old), " -> ", nrow(new), " codes.")

  counts <- table(changes$change)
  for (nm in names(counts)) {
    message("  ", format(nm, width = 18), counts[[nm]])
  }

  # Rule changes deserve to be read out: they alter every person's result,
  # whereas a code change usually alters a few.
  rules <- changes[changes$change == "rule changed", , drop = FALSE]
  if (nrow(rules)) {
    message("\nRule changes - these change results for everyone:")
    for (i in seq_len(nrow(rules))) {
      message("  ", rules$condition[i], " ", rules$item[i], ": ",
              if (nzchar(rules$old[i])) rules$old[i] else "(unset)", " -> ",
              if (nzchar(rules$new[i])) rules$new[i] else "(unset)")
    }
  }
  invisible(NULL)
}
