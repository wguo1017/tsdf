#' run dose-finding simulations
#' @description Run dose-finding simulations based on a customized decision table.
#' @param truep A vector of length k (the number of doses being considered in the trial), with values equal to the true probabilities of toxicity at the dose levels.
#' @param decTable A customized decision table. (same format as output of \code{\link{dec.table}})
#' @param start.level Starting dose level. Defaults to 1, i.e. the lowest dose level.
#' @param nsim The number of simulation trials. Defaults to 1000.
#' @return The functions \code{\link{summary.dec.sim}} is used to obtain and print a summary table of the results (recommended). An object of class \code{"dec.sim"} is a list containing:
#'  \item{mtd}{a vector of dose levels giving the recommended maximum tolerated dose (MTD) at the end of the trial.}
#'  \item{mtd.prob}{a vector of length \code{k} giving the average proportions of selected as MTD at each dose level}
#'  \item{n.patients}{the average number of patients dosed at each level.}
#'  \item{dlt}{the average number of DLTs experienced at each dose level}
#'  \item{truep}{input; true probabilities of toxicity.}
#'  \item{start.level}{input; starting dose level.}
#'  \item{nsim}{input; number of simulated trails.}
#' @author Wenchuan Guo <wguo007@ucr.edu>
#' @import stats
#' @export
#' @examples
#' truep <- c(0.3, 0.45, 0.5, 0.6)
#' res <- dec.table(0.6,0.4,0.2,0.3,0.3,c(3,3,3))
#' out <- dec.sim(truep, res$table, start.level=2, nsim=1000)
#' summary(out, pt=0.3)

dec.sim  <- function(truep, decTable, start.level = 1, nsim = 1000) {
  # initialization
  n_dose <- length(truep)
  input.sample <- as.numeric(colnames(decTable))
  maxn <- nrow(decTable) - 1
  nc <- ncol(decTable)
  if(input.sample[length(input.sample)] != maxn) {
    stop("Please check decision table format")
  }
  add <- c(input.sample[1], diff(input.sample))
  mtd <- rep(0, nsim)
  np <- matrix(0, nrow = nsim, ncol = n_dose)
  dlt <- matrix(0, nsim, n_dose)
  for(i in 1:nsim){
    # dose need to be removed
    rm.dose <- NULL
    # current dose level
    dose <- start.level
    # stage current dose at
    sample.stage <- rep(1, n_dose)
    while(mtd[i] == 0 & !is.na(mtd[i])){
      add.sample <- add[sample.stage[dose]]
      np[i, dose] <- np[i, dose] + add.sample
      dlt[i, dose] <- dlt[i, dose] + sum(rbinom(add.sample, 1, truep[dose]))
      des <- decTable[dlt[i, dose]+1, sample.stage[dose]]
      sample.stage[dose] <- sample.stage[dose] + 1
      # case: stay (S)
      if(des == "S"){
        if(np[i, dose] != maxn) {
          dose <- dose
        } 
        else {
          mtd[i] <- dose
        } 
      }
      # case : escalate (E)
      if(des == "E") {
        # check if next dose level is available
        if((dose+1) %in% rm.dose){
          mtd[i] <- dose
        } else {
          # check if can escalate
          if(dose != n_dose){
            if(np[i, dose + 1] != maxn) {
              dose <- dose + 1
            } else {
              next_des <- decTable[dlt[i, dose + 1] + 1, nc]
              if(next_des == "D" | next_des == "DU"){
                mtd[i] <- dose
              } else {
                mtd[i] <- dose + 1
              }
            }
          } else {
            mtd[i] <- dose
          }	
        }
      }
      # case : de-escalate (D)
      if(des == "D") {
        if(dose != 1) {
          if(np[i, dose-1] != maxn) {
            dose <- dose - 1
          } else {
            mtd[i] <- dose -1 
          }
        } else {
          mtd[i] <- NA
        }
      }
      # case : de-escalate (DU)
      if(des == "DU") {
        if(dose != 1) {
          # can not go back to this dose again
          rm.dose <- c(rm.dose, dose)
          if(np[i, dose-1] != maxn) {
            dose <- dose - 1
          } else {
            mtd[i] <- dose -1 
          }
        } else {
          mtd[i] <- NA
        }
      }
    }		
  }
  mtd.prob <- sapply(1:n_dose, function(ii) sum(mtd == ii, na.rm = TRUE)/nsim)
  out <- list(mtd = mtd, mtd.prob = mtd.prob, dlt = colMeans(dlt), n.patients = colMeans(np), truep = truep, start.level = start.level, nsim = nsim)
  class(out) <- "dec.sim"
  return(out)
}


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
  } else {
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
  } else {
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
    out <- plot(1:n_dose, truep, xlab = "Dose level", ylab = "Toxicity", col = col, type = "b", ylim = c(min(truep) * 0.9, max(truep) * 1.1), xaxt = "n", pch = 19, panel.first = c(box(), grid()))
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