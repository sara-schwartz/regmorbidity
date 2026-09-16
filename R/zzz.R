# .............................................................................
# zzz.R - package-level housekeeping
#
# PURPOSE
#   Nothing here affects behaviour. It exists so that `R CMD check` stays quiet
#   about names that only look undefined.
#
# WHY THIS IS NEEDED
#   data.table refers to columns as bare names inside `dt[...]`, and rlang lets
#   mb_flag_users() take a bare column name. Static analysis cannot tell those from
#   forgotten variables, so each is declared here. A name in this list is NOT a
#   global variable - it is a column name evaluated inside a data.table or a
#   symbol captured by ensym().
#
#   Add to this list only when R CMD check complains, and only for names that
#   are genuinely columns or captured symbols.
# .............................................................................

utils::globalVariables(c(
  # column names used inside data.table expressions (extract.R)
  "id", "date", "code", "mb_ok", ".N",
  # bare column names accepted as defaults by mb_flag_users() (extract.R)
  "ATC", "PNR", "EKSD"
))
