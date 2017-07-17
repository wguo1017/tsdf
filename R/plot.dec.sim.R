#' plot simulation results from a dec.sim object
#' @description Three plots are currently available: a plot of true toxicity at each dose level (\code{type="s"}); a bar plot of the probability of selecting as the MTD for each dose level (\code{type="prob"}) and a bar plot of the average number of patients treated at each dose level (\code{type="np"}).
#' @param x an object of class \code{"dec.sim"}, a result of a call to \code{dec.sim} or \code{sl.sim}.
#' @param pt target toxicity for each scenario.
#' @param s scenario to be plotted. Default is 1.
#' @param type plot type. See descriptions above.
#' @param show a logical value indicating if values are shown on plot.
#' @param ... Not used argument.
#' @import graphics
#' @export

plot.dec.sim <- function(x, pt, s = 1, type = c("s", "prob", "np"), show = TRUE, ...) {
  # check
  type <- match.arg(type)
  if(class(x)[1] == "sl.sim") {
    obj <- x[[s]]
  } else {
    obj <- x
  }
  # initialization
  truep <- obj$truep
  n_dose <- length(truep)
  names.arg <- 1:n_dose
  if(missing(pt)) {
    warning("true toxicity is missing")
  } else {
    mtd <- max(which(truep <= pt[s]))
    names.arg[mtd] <- paste(mtd, "(MTD)")
  }
  # plot as specified "type"
  if(type == "prob") {
    mtd.prob <- sapply(1:n_dose, function(ii) mean(obj$MTD==ii))
    out <- barplot(mtd.prob, xlab = "Dose level", col = "cornflowerblue", main = paste("Scenario", s, ":", "prob. of selected as the MTD "), names.arg = names.arg, panel.first = grid(), ylim = c(0, max(mtd.prob) + 0.1))
    if(show){
      text(out, mtd.prob, labels = mtd.prob, pos = 3, cex = 1, col = "darkblue")
    }
  }
  if(type == "np") {
    np <- obj$n.patients
    out <- barplot(np, xlab = "Dose level", col = "cornflowerblue", names.arg = names.arg,  panel.first = grid(), main = paste("Scenario", s, ":", "avg. number of patients treated"), ylim = c(0, max(np) + 0.8))
    if(show){
      text(out, np, labels = np, pos = 3, cex = 1, col = "darkblue")
    }
  }
  if(type == "s") {
    out <- plot(1:n_dose, truep, xlab = "Dose level", ylab="", col = "cornflowerblue", type = "o", main = paste("Scenario", s, ":", "true toxicity"), ylim = c(min(truep) - 0.03, max(truep) + 0.03), xaxt="n", panel.first = grid())
    axis(1, at = 1:n_dose, labels = names.arg)
    if(!missing(pt)) {
      abline(pt[s], 0, lty = 5)
      text(n_dose/2, pt[s], labels = paste("Target = ", pt[s]), pos = 3, cex = 1, col = "darkblue")
    }
    if(show){
      text(1:n_dose, truep, labels = truep, pos = 3, cex = 1, col = "darkblue")
    }
  }
}




