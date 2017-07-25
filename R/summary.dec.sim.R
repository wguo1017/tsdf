#' Summarizing simulation results from a dec.sim object
#' @description \code{summary} method for class \code{"dec.sim"}.
#' @details \code{summary} is used for formating important statistics for dose-finding simulation. Giving the target toxicity, it returns the probability of selecting current dose level as the MTD, probability of selecting the true MTD, probability of subjects treated at or below the true MTD, etc. The MTD is defined as the highest dose level such that the toxicity probability is less than target toxicity probability, if target is less than the smallest probability, then the lowest dose level is set as MTD. For example, if target is 0.3 and true toxicity for five doses are 0.1, 0.25, 0.35, 0.40, then MTD is dose 2.
#' @param object an object of class \code{"dec.sim"}, a result of a call to \code{dec.sim} or \code{sl.sim}.
#' @param pt target toxicity for each scenario.
#' @param ... Not used argument.
#' @examples 
#' test.file <- system.file("extdata", "testS.csv", package = "tsdf")
#' dt <- dec.table(0.6,0.4,0.2,0.3,0.3,c(3,3,3))
#' out <- sl.sim(dt$table, test.file)
#' pt <- c(0.3, 0.4)
#' summary(out, pt)
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
  cat("What does each column represent ?", "\n\n")
  cat("Level : Dose level", "\n")
  cat("Truth : True toxicity probability", "\n")
  cat("MTD : The probability of selecting current dose level as the MTD", "\n")
  cat("DLT : The average number of subjects experienced DLT at current dose level", "\n")
  cat("NP : The average number of subjects treated at current dose level", "\n\n")
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
    if(sum(truep <= pt[i]) == 0) {
      mtd.dose <- 1
      avg.dose[i] <- res[, 3][mtd.dose]
      avg.prob[i] <- res[,5][1] / sum(res[, 5])
    } else {
      mtd.dose <- max(which(truep <= pt[i]))
      avg.dose[i] <- res[, 3][mtd.dose]
      avg.prob[i] <-sum(res[, 5][truep <= pt[i]])/sum(res[, 5])
    }
    # table legends
    cat("Scenario ", i, "\n\n")
    cat("Simulation results are based on", ans$nsim, "simulated trials","\n")
    cat("Starting dose level is", ans$start.level, "; MTD is dose", mtd.dose, "\n")
    cat("Target toxicity probability =", pt[i], "\n")
    cat("Average number of subjects =", avg.np[i], "\n")
    cat("Probability of selecting the true MTD =", avg.dose[i], "\n")
    cat("Probability of no selection =", mean(is.na(ans$mtd)), "\n")
    cat("Probability of subjects treated at or below the true MTD =", avg.prob[i], "\n\n")
    print(as.data.frame(res))
    cat("\n")
  }
}