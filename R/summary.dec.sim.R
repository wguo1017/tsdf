#' Summarizing simulation results from a dec.sim object
#' @description \code{summary} method for class \code{"dec.sim"}.
#' @param object an object of class \code{"dec.sim"}, a result of a call to \code{dec.sim} or \code{sl.sim}.
#' @param pt target toxicity for each scenario.
#' @param ... Not used argument.
#' @import stats
#' @export

summary.dec.sim <- function(object, pt, ...) {
  if(missing(pt)) {
    stop("Missing true toxicity.")
  }
  if(class(object)[1] == "sl.sim" & length(object) != length(pt)) {
    warning("only returns the first scenario stats")
  }
  cat("Truth: True Toxicity Probability;", "\n")
  cat("Prob: The probability of selecting each dose level as the MTD;", "\n")
  cat("NP:  The average number of patients treated at each dose level.", "\n\n")
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
    res <- matrix(0, n_dose, 4)
    colnames(res) <- c("Level", "Truth", "Prob.", "NP")
    res[, 1] <- 1:n_dose
    res[, 2] <- truep
    res[, 3] <- sapply(1:n_dose, function(ii) mean(ans$MTD==ii))
    res[, 4] <- ans$n.patients
    # calculate stats
    avg.np[i] <- sum(res[, 4])
    avg.dose[i] <- res[, 3][max(which(truep <= pt[i]))]
    avg.prob[i] <-sum(res[, 4][truep <= pt[i]])/sum(res[, 4])
    # table legends
    cat("Scenario ", i, "\n\n")
    cat("The average number of subjects:", avg.np[i], "\n")
    cat("The probability of selecting the true MTD:", avg.dose[i], "\n")
    cat("The probability of patients treated at or below the true MTD:", avg.prob[i], "\n\n")
    print(as.data.frame(res))
    cat("\n")
  }
}
