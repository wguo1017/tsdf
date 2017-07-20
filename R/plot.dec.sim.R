#' plot simulation results from a dec.sim object
#' @description Three plots are currently available: a plot of true toxicity at each dose level (\code{type = "s"}); a bar plot of the probability of selecting as the MTD for each dose level (\code{type = "prob"}); a bar plot of the average number of patients treated at each dose level (\code{type = "np"}); a bar plot of the average number of patients experienced DLT at each dose level (\code{type = "dlt"}) and \code{type = "all"} generates all above plots.
#' @param x an object of class \code{"dec.sim"} or \code{"sl.sim"}, a result of a call to \code{dec.sim} or \code{sl.sim}.
#' @param pt a vector with target toxicity for each scenario.
#' @param s scenario to be plotted. Defaults to 1.
#' @param type plot type. See descriptions above.
#' @param label a logical value indicating if values are shown on plot.
#' @param col graphical parameter \code{col}; see details \code{\link{par}}.
#' @param text.col plotting color of text shown.
#' @param cex graphical parameter \code{col}; see details \code{\link{par}}.
#' @param ... arguments to be passed to \code{\link{plot}} methods.
#' @import graphics
#' @export
#' @examples 
#' # generate decision table
#' dt <- dec.table(0.6,0.4,0.2,0.3,0.3,c(3,3,3))
#' # simulate trials from test data 
#' test.file <- system.file("extdata", "testS.csv", package = "tsdf")
#' out <- sl.sim(dt$table, test.file)
#' plot(out, pt=rep(0.3,2), s=1, type="all")
#' plot(out, pt=rep(0.3,2), s=2, type="prob")
#' plot(out, pt=rep(0.3,2), s=1, type="np")
#' plot(out, pt=rep(0.3,2), s=2, type="dlt")

plot.dec.sim <- function(x, pt, s = 1, type = c("all", "s", "prob", "np", "dlt"), label = TRUE, col = "cornflowerblue", text.col = "darkblue", cex = 1, ...) {
  type <- match.arg(type)
  if (class(x)[1] == "sl.sim") {
    obj <- x[[s]]
    ns <- length(x)
  }
  else {
    obj <- x
    ns <- 1
  }
  truep <- obj$truep
  n_dose <- length(truep)
  names.arg <- 1:n_dose
  scen <- paste("( Scenario", s, ")")
  if (missing(pt)) {
    print("No target toxicity/mtd shown")
  }
  else if (length(pt) != ns) {
    stop("true toxicity is a vector with length = # of cenarios")
  }
  else {
    mtd <- max(which(truep <= pt[s]))
    names.arg[mtd] <- paste(mtd, "(MTD)")
  }
  if (type == "all") {
    par(mfrow = c(2, 2), cex = cex)
    plot(x, pt, s, "s", label, col, text.col, cex, ...)
    plot(x, pt, s, "prob", label, col, text.col, cex, ...)
    plot(x, pt, s, "dlt", label, col, text.col, cex, ...)
    plot(x, pt, s, "np", label, col, text.col, cex, ...)
  }
  if (type == "prob") {
    ans <- obj$mtd.prob
    out <- barplot(rep(NA, n_dose), xlab = "Dose level", names.arg = names.arg, panel.first = box(), ylim = c(0, max(ans) * 1.1))
    grid()
    barplot(ans, add = TRUE, col = col)
    title(paste("Probability of selection", scen))
    if (label) {
      text(out, ans, labels = ans, pos = 3, cex = cex, col = text.col)
    }
  }
  if (type == "dlt") {
    ans <- obj$dlt
    out <- barplot(rep(NA, n_dose), xlab = "Dose level", names.arg = names.arg, panel.first = box(), ylim = c(0, max(ans) * 1.1))
    grid()
    barplot(ans, add = TRUE, col = col)
    title(paste("Average number of DLTs", scen))
    if (label) {
      text(out, ans, labels = ans, pos = 3, cex = cex, col = text.col)
    }
  }
  if (type == "np") {
    ans <- obj$n.patients
    out <- barplot(rep(NA, n_dose), xlab = "Dose level", names.arg = names.arg, panel.first = box(), ylim = c(0, max(ans) * 1.1))
    grid()
    barplot(ans, add = TRUE, col = col)
    title(paste("Average number of patients", scen))
    if (label) {
      text(out, ans, labels = ans, pos = 3, cex = cex, col = text.col)
    }
  }
  if (type == "s") {
    out <- plot(1:n_dose, truep, xlab = "Dose level", ylab = "", col = col, type = "b", ylim = c(min(truep) * 0.9, max(truep) * 1.1), xaxt = "n", pch = 19, panel.first = c(box(), grid()))
    axis(1, at = 1:n_dose, labels = names.arg)
    title(paste("True toxicity", scen))
    if (!missing(pt)) {
      abline(pt[s], 0, lty = 5)
      text(n_dose/2, pt[s], labels = paste("Target =", pt[s]), pos = 3, cex = cex, col = text.col)
    }
    if (label) {
      text(1:n_dose, truep, labels = truep, pos = 3, cex = cex, col = text.col)
    }
  }
}
