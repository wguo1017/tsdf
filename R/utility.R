#' plot decision table from a "dec.table" object.
#' @description \code{plot} method for class "\code{dec.table}"
#' @param x an object of class \code{"dec.table"}, a result of a call to \code{dec.table}.
#' @param ... Not used argument.
#' @details \code{plot.dec.table} prints the decision boundarys.
#' @import graphics
#' @export
plot.dec.table <- function(x, ...) {
  n <- x$n
  nc <- cumsum(n)
  r <- x$E
  s <- x$D
  su <- x$DU
  col <- c(1:3, "cornsilk")
  plot(nc, r, xaxt = "n", ylab = "Boundary", xlab = 'Sample Size', ylim = c(0, max(su)+2), type = "o", col=col[1], main = "Decision Plot")
  points(nc, s + 1, type = "o", col = col[2])
  points(nc, su + 1, type = "o", col = col[3])
  axis(1, at = nc, labels = nc)
  polygon(c(nc, rev(nc)), c(r + 0.05, rev(s + 1 - 0.05)), col = "cornsilk", border = NA)
  legend("topleft", c("E (<=)", "D (>=)", "DU (>=)", "S"), fill = col)
}

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
  cat("    DU : De-escalate and never use this dose again","\n")
}

