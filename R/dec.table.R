#' generate three-stage dose-finding decision table
#' @description Generate three stage dose finding decision table
#' @param alpha.l Left-side overall type 1 error. Control the upper bound of dose escalation.
#' @param alpha.r Right-side overall type 1 error. Control the lower bound of dose de-escalatition.
#' @param alpha.u Right-side overall type 1 error. This also controls the lower bound of dose de-escalatition, but it is used to find lower bound for "DU".
#' @param pc A numeric vector of target toxicity. Should be a vector with 1 or 2(when the target is an interval).
#' @param pc.u A numeric vector of target toxicity which is used to obtain "DU" in the decision table.
#' @param n A vector of sample size at each stage. \code{sum(n)} is the total sample size.
#' @param sf The alpha-spending function to be used. \code{sf="OF"} or "\code{sf="Pocock"} uses spending function in R package \code{\link{gsDesign}}; or a user supplied spending function.
#' @param ... Not used argument.
#' @return An object of class "dec.table" is a list containing:
#'  \item{table}{the generated decision table.}
#'  \item{alpha.two}{a vector of true type 1 error for two-tailed test.}
#   \item{alpha.one}{a vector of true type 1 error for right-tailed test.}
#'  \item{E}{a vector of "E" bound.}
#'  \item{D}{a vector of "D" bound.}
#'  \item{DU}{a vector of "DU" bound.}
#'  \item{pc}{input; a vector of target toxicity}
#'  \item{pc.u}{input; a vector of target toxicity }
#'  \item{n}{input; a vector with sample size at each stage.}
#'  \item{sf}{input; the alpha-spending function used.}
#' @author Wenchuan Guo <wguo007@ucr.edu>
#' @import gsDesign
#' @import stats
#' @export
#' @examples
#' n <- rep(3, 3)
#' alpha.l <- 0.6
#' alpha.r <- 0.4
#' alpha.u <- 0.3
#' pc <- c(0.29, 0.31)
#' pc.u <- 0.3
#' # print out decision table
#' dec.table(alpha.l, alpha.r, alpha.u, pc, pc.u, n)$table

dec.table <- function(alpha.l, alpha.r, alpha.u, pc, pc.u, n, sf  = "Pocock") {
  # check
  err <- c(alpha.l, alpha.r, alpha.u)
  k <- length(n)
  if(sum(err < 0 | err > 1) != 0) {
    stop("'alpha' should between 0 and 1.")
  }
  if(length(pc) > 2) {
    stop("'pc''s length should less than 2 (two-tailed test).")
  }
  if(length(pc) == 1) {
    pc <- rep(pc, 2)
  }
  if(length(pc.u) != 1) {
    stop("'pc.u''s length should be 1 (right-tailed test).")
  }
  if(k != 3 & k != 2) {
    stop("This function only find two-stage/three-stage optimal design")
  }
  if(sf != "Pocock" & sf != "OF" & !is.function(sf)){
    stop("'sf' should be either 'OF', 'Pocock' or a user specified spending function" )
  }
  if(k == 3) {
    out.two <- three.opt(alpha.l, alpha.r, pc, n, sf)
    out.one <- right.three.opt(alpha.u, pc.u, n, sf)
  } else {
    out.two <- two.opt(alpha.l, alpha.r, pc, n, sf)
    out.one <- right.two.opt(alpha.u, pc.u, n, sf)
  }
  
  des <- list(E = out.two$out[1, 1:k], D = out.two$out[1, (k+1):(2*k)], DU = out.one$out[1, 1:k],  n = n, pc = pc, pc.u = pc.u, sf = sf, alpha.two = out.two$out[1, (2*k+1):(4*k)], alpha.one = out.one$out[1, (k+1):(2*k)])
  r <- des$E
  s <- des$D
  su <- des$DU
  ns <- length(n)
  nt <- sum(n)
  nc <- cumsum(n)
  ans <- matrix(0, nt+1, ns)
  rownames(ans) <- 0:nt
  colnames(ans) <- nc
  for(j in 1:length(r)){
    ans[ ,j][1:(r[j]+1)] <- "E"
    ans[ ,j][(s[j]+2):(nc[j]+1)] <- "D"
    ans[ ,j][(su[j]+2):(nc[j]+1)] <- "DU"
    ans[ ,j][which(ans[, j] ==0)] <- "S"
    ind.s <- (1:(nc[k]+1)) > nc[j]+1
    ans[ , ][ind.s] <- rep("0", sum(ind.s))
  }
  out <- c(des, list(table=as.table(ans)))
  class(out) <- "dec.table"
  return(out)
}


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
  col <- c("green", "gold2", "red", "beige")
  plot(nc, r, xaxt = "n", ylab = "Boundary", xlab = 'Sample Size', ylim = c(0, max(su)+2), type = "o", col=col[1], main = "Decision Plot", panel.first = grid())
  points(nc, s + 1, type = "o", col = col[2])
  points(nc, su + 1, type = "o", col = col[3])
  axis(1, at = nc, labels = nc)
  polygon(c(nc, rev(nc)), c(r + 0.05, rev(s + 1 - 0.05)), col = col[4], border = NA)
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
  cat("	  Row : Number of DLTs", "\n")
  cat("Column : Number of subjects", "\n")
  cat("     E : Escalate to the next higher dose", "\n")
  cat("     S : Stay at the same dose", "\n")
  cat("     D : De-escalate to the previous lower dose", "\n")
  cat("    DU : De-escalate and never use this dose again", "\n")
}