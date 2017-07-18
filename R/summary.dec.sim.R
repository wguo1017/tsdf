#' Summarizing simulation results from a dec.sim object
#' @description \code{summary} method for class \code{"dec.sim"}.
#' @param object an object of class \code{"dec.sim"}, a result of a call to \code{dec.sim} or \code{sl.sim}.
#' @param pt target toxicity for each scenario.
#' @param ... Not used argument.
#' @import stats
#' @export

summary.dec.sim <- function(object, pt, ...) {
  if(missing(pt)) {
    warning("Missing true toxicity; set to be 1.")
    pt <- 1
  }
  if(class(object)[1] == "sl.sim" & length(object) != length(pt)) {
    warning("pt length not equal to number of scenarios; only returns the first scenario stats")
  }
  cat("How to read the table :", "\n\n")
  cat("Level : Dose level", "\n")
  cat("Truth : True toxicity probability", "\n")
  cat("MTD : The probability of selecting as the MTD", "\n")
  cat("DLT : The average number of patients experienced DLT", "\n")
  cat("NP : The average number of patients treated", "\n\n")
  ns <- length(pt)
  avg.np <- rep(0, ns)
  avg.dose <- rep(0, ns)
  avg.prob <- rep(0, ns)
  out <- vector("list", ns)
  for(i in 1:ns) {
    if(class(object)[1] == "sl.sim"){
      ans <- object[[i]]
    } else {
      ans <- object
    }
    truep <- ans$truep
    n_dose <- length(truep)
    res <- matrix(0, n_dose, 5)
    colnames(res) <- c("Level", "Truth", "MTD", "DLT", "NP")
    res[, 1] <- 1:n_dose
    res[, 2] <- truep
    res[, 3] <- ans$mtd.prob
    res[, 4] <- ans$dlt
    res[, 5] <- ans$n.patients
    # calculate stats
    avg.np[i] <- sum(res[, 5])
    avg.dose[i] <- res[, 3][max(which(truep <= pt[i]))]
    avg.prob[i] <-sum(res[, 5][truep <= pt[i]])/sum(res[, 5])
    # table legends
    cat("Scenario ", i, "\n\n")
    cat("Simulation results are based on", ans$nsim, "simulated trials","\n")
    cat("Starting dose level is", ans$start.level, "\n")
    cat("Target toxicity probability =", pt[i], "\n")
    cat("Average number of subjects =", avg.np[i], "\n")
    cat("Probability of selecting the true MTD =", avg.dose[i], "\n")
    cat("Probability of patients treated at or below the true MTD =", avg.prob[i], "\n\n")
    print(as.data.frame(res))
    cat("\n")
  }
}
