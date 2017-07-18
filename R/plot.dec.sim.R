#' plot simulation results from a dec.sim object
#' @description Three plots are currently available: a plot of true toxicity at each dose level (\code{type = "s"}); a bar plot of the probability of selecting as the MTD for each dose level (\code{type = "prob"}); a bar plot of the average number of patients treated at each dose level (\code{type = "np"}); a bar plot of the average number of patients experienced DLT at each dose level (\code{type = "dlt"}) and \code{type = "all"} generates all above plots.
#' @param x an object of class \code{"dec.sim"} or \code{"sl.sim"}, a result of a call to \code{dec.sim} or \code{sl.sim}.
#' @param pt a vector with target toxicity for each scenario.
#' @param s scenario to be plotted. Defaults to 1.
#' @param type plot type. See descriptions above.
#' @param show a logical value indicating if values are shown on plot.
#' @param col graphical parameter \code{col}; see details \code{\link{par}}.
#' @param text.col plotting color of text shown.
#' @param cex graphical parameter \code{col}; see details \code{\link{par}}.
#' @param ... arguments to be passed to \code{\link{plot}} methods.
#' @import graphics
#' @export

plot.dec.sim <- function(x, pt, s = 1, type = c("all", "s", "prob", "np", "dlt"), show = TRUE, col = "cornflowerblue", text.col = "darkblue", cex = 1, ...) {
  # check
  type <- match.arg(type)
  if(class(x)[1] == "sl.sim") {
    obj <- x[[s]]
    ns <- length(x)
  } else {
    obj <- x
    ns <- 1
  }
  if(length(pt) != ns) {
    stop("true toxicity is a vector with length = # of cenarios")
  }
  # initialization
  truep <- obj$truep
  n_dose <- length(truep)
  names.arg <- 1:n_dose
  scen <- paste("( Scenario", s, ")")
  if(missing(pt)) {
    warning("true toxicity is missing")
  } else {
    mtd <- max(which(truep <= pt[s]))
    names.arg[mtd] <- paste(mtd, "(MTD)")
  }
  # plot as specified "type"
  if(type == "all") {
    par(mfrow = c(2,2), cex = cex)
    plot(x, pt, s, "s", show, col, text.col, cex, ...)
    plot(x, pt, s, "prob", show, col, text.col, cex, ...)
    plot(x, pt, s, "dlt", show, col, text.col, cex, ...)
    plot(x, pt, s, "np", show, col, text.col, cex, ...)
  }
  if(type == "prob") {
    ans <- obj$mtd.prob
    out <- barplot(rep(NA, n_dose), xlab = "Dose level", names.arg = names.arg,  panel.first = box(), ylim = c(0, max(ans) * 1.1))
    grid()
    barplot(ans, add = TRUE, col = col)
    title(paste("Probability of selection", scen))
    if(show){
      text(out, ans, labels = ans, pos = 3, cex = cex, col = text.col)
    }
  }
  if(type == "dlt") {
    ans <- obj$dlt
    out <- barplot(rep(NA, n_dose), xlab = "Dose level", names.arg = names.arg,  panel.first = box(), ylim = c(0, max(ans) * 1.1))
    grid()
    barplot(ans, add = TRUE, col = col)
    title(paste("Average number of DLTs", scen))
    if(show){
      text(out, ans, labels = ans, pos = 3, cex = cex, col = text.col)
    }
  }
  if(type == "np") {
    ans <- obj$n.patients
    out <- barplot(rep(NA, n_dose), xlab = "Dose level", names.arg = names.arg,  panel.first = box(), ylim = c(0, max(ans) * 1.1))
    grid()
    barplot(ans, add = TRUE, col = col)
    title(paste("Average number of patients", scen))
    if(show){
      text(out, ans, labels = ans, pos = 3, cex = cex, col = text.col)
    }
  }
  if(type == "s") {
    out <- plot(1:n_dose, truep, xlab = "Dose level", ylab="", col = col, type = "b", ylim = c(min(truep) * 0.9, max(truep) * 1.1), xaxt="n", pch = 19, panel.first = c(box(), grid()))
    axis(1, at = 1:n_dose, labels = names.arg)
    title(paste("True toxicity", scen))
    if(!missing(pt)) {
      abline(pt[s], 0, lty = 5)
      text(n_dose/2, pt[s], labels = paste("Target =", pt[s]), pos = 3, cex = cex, col = text.col)
    }
    if(show){
      text(1:n_dose, truep, labels = truep, pos = 3, cex = cex, col = text.col)
    }
  }
}

