# .............................................................................
# run.R - the shared per-condition driver
#
# PURPOSE
#   Both mb_extract_medication() and mb_extract_diagnosis() do the same thing around
#   the edges: loop over conditions, skip what is already on disk, time each
#   one, catch a failure without losing the rest, and build a summary row. Only
#   the work in the middle differs. That loop lives here once.
#
# INPUT   a vector of conditions and a worker function
# OUTPUT  a named list of results and a summary data frame
#
# WHY A DRIVER RATHER THAN TWO LOOPS
#   The resume behaviour is the part most likely to be got subtly wrong, and it
#   is the part whose failure is least visible - a run that silently skips a
#   condition looks exactly like a run that finished. One implementation, one
#   set of tests.
#
# CONTENTS
#   1. The driver
# .............................................................................


# 1. The driver ----

#' Run a worker over many conditions, checkpointing each to disk
#'
#' @param conditions Character vector of condition names.
#' @param worker A function of one argument (the condition name) returning a
#'   list with at least `data`, plus any fields to record in the summary. It may
#'   also return `extra`, a named list of side objects; each is written as
#'   `<condition><name>.rds`.
#' @param outdir Directory for one `.rds` per condition. `NULL` = keep results
#'   in memory and write nothing.
#' @param resume Skip conditions whose `.rds` already exists.
#' @param verbose Print progress.
#' @param unit Word for the progress line, e.g. "dispensings" or "diagnoses".
#' @return A list with `results` (named, `NULL` entries where written to disk)
#'   and `summary` (one row per condition).
#' @keywords internal
mb_run_conditions <- function(conditions, worker, outdir = NULL, resume = TRUE,
                              verbose = TRUE, unit = "records") {

  summaries <- vector("list", length(conditions))
  results   <- vector("list", length(conditions))
  names(results) <- conditions

  for (i in seq_along(conditions)) {
    condition <- conditions[i]
    out_file  <- if (is.null(outdir)) NA_character_ else
      file.path(outdir, paste0(condition, ".rds"))

    if (!is.null(outdir) && resume && file.exists(out_file)) {
      if (verbose) {
        message("[", i, "/", length(conditions), "] ", condition,
                " - already done, skipping")
      }
      summaries[[i]] <- mb_summary_row(condition, status = "skipped",
                                       file = out_file)
      next
    }

    if (verbose) message("[", i, "/", length(conditions), "] ", condition, " ...")
    started_at <- Sys.time()

    # try() rather than letting it propagate: on a run of forty conditions
    # taking hours, one bad code list must not throw away the rest.
    outcome <- try(worker(condition), silent = TRUE)

    elapsed_secs <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))

    if (inherits(outcome, "try-error")) {
      msg <- conditionMessage(attr(outcome, "condition"))
      warning("Condition '", condition, "' failed: ", msg, call. = FALSE)
      summaries[[i]] <- mb_summary_row(condition, status = "error",
                                       file = out_file, seconds = elapsed_secs,
                                       message = msg)
      next
    }

    if (!is.null(outdir)) {
      saveRDS(outcome$data, out_file)
      for (nm in names(outcome$extra)) {
        if (!is.null(outcome$extra[[nm]])) {
          saveRDS(outcome$extra[[nm]],
                  file.path(outdir, paste0(condition, nm, ".rds")))
        }
      }
    } else {
      results[[i]] <- outcome$data
    }

    # Only the fields the worker actually supplied; the rest stay NA. The four
    # the driver sets itself are excluded, because a worker returning one of
    # them would otherwise fail with "matched by multiple actual arguments",
    # which says nothing about what is wrong.
    ours <- c("condition", "status", "file", "seconds")
    summary_fields <- outcome[setdiff(intersect(names(outcome),
                                                names(formals(mb_summary_row))),
                                      ours)]
    summaries[[i]] <- do.call(mb_summary_row, c(
      list(condition = condition, status = "done", file = out_file,
           seconds = elapsed_secs),
      summary_fields))

    if (verbose) {
      message("    ", format(outcome$n_persons, big.mark = ","),
              " persons from ", format(outcome$n_rows, big.mark = ","),
              " ", unit, " (", round(elapsed_secs), "s)")
    }
  }

  summary <- do.call(rbind, summaries)
  rownames(summary) <- NULL

  list(results = results, summary = summary)
}


#' Shape the return value the same way for every extract_* function
#' @keywords internal
mb_finish_run <- function(run, outdir, verbose) {
  if (is.null(outdir)) {
    out <- run$results
    attr(out, "summary") <- run$summary
    return(out)
  }
  if (verbose) {
    n_errors <- sum(run$summary$status == "error")
    if (n_errors) {
      message("Finished with ", n_errors, " error(s) - see the summary.")
    }
  }
  invisible(run$summary)
}
