#' print decision table from a "dec.table" object.
#' @description \code{print} method for class "\code{dec.table}"
#' @param x an object of class \code{"dec.table"}, a result of a call to \code{dec.table}.
#' @param ... Not used argument.
#' @details \code{print.dec.table} prints the decision table with legend keys.
#' @export
print.dec.table <- function(x, ...) {
  print(x$table)
  cat("\n")
  cat("Column : number of DLTs ;", "\n")
  cat("   Row : number of patients at current stage ;", "\n")
  cat("     E : Escalate to the next higher dose ;", "\n")
  cat("     S : Stay at the same dose ;", "\n")
  cat("     D : De-escalate to the previous lower dose ;", "\n")
  cat("    DU : De-escalate and never use this dose again.", "\n")
}
