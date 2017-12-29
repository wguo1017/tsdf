#' run dose-finding simulations
#' @description Run dose-finding simulations based on a customized decision table.
#' @details Assume there are $d$ dose levels to be studied. Denote the cumulative number of patients treated and cumulative number of DLTs at the current dose level are $n_i$ and $m_i$, respectively. $n_{max}$ is the maximum number of patients permitted to be treated at each dose level. The procedure is as follows
#' \itemize{
#' \item{Step 1 : }{Update cumulative number of DLTs $m_i$ and total number of patients $n_i$ treated at the current dose and use the decision table to make a decision: if decision is ``S" --> step 2; if decision is ``D" or ``DU'' --> step 3; if decision is ``E'' --> step 4}
#' \item{Step 2 : }{If $n_i = n_{max}$, declare dose i as the MTD; otherwise, update $m_i$ and $n_i$ with additional cohort of patients and go to Step 1.}
#' \item{Step 3 : }{If the current dose level is the highest dose level, then: if n_i < n_max, update $m_i$ and $n_i$ with additional cohort of patients and go to Step 1; otherwise, stop the trial and declare that the MTD is higher than the highest dose level (inconclusive); If the current dose is not the lowest dose, then: if $n_{i-1} < n_{max}$, update $m_{i-1}$ and $n_{i-1}$ with additional cohort of patients and set the current dose level to be the next lower dose level, and go to Step 1; otherwise, stop the trial and declare the next lower dose level is the MTD; Additionally, if the decision is ``DU'', record this dose level as DU and never treat additional patients at the current dose level again.}
#' \item{Step 4 : }{If the current dose level is the highest dose level, then: if $n_i < n_{max}$, update $m_i$ and $n_i$ with additional cohort of patients and go to Step 1; otherwise, stop the trial and declare that the MTD is higher than the highest dose level (inconclusive); If the next higher dose level is of status DU, then: if $n_i < n_{max}$, update $m_i$ and $n_i$ with additional cohort of patients and go to step 1; otherwise stop, the current dose level is MTD; Otherwise: if $n_{i+1} < n_{max}$, update $m_{i+1}$ and  $n_{i+1}$ with additional cohort of patients, set the current dose level to be next higher dose level, and go to step 1; else, the current dose level is the MTD. }
#' }
#' @param truep a vector of length k (the number of doses being considered in the trial), with values equal to the true probabilities of toxicity at the dose levels.
#' @param decTable a customized decision table. (same format as output of \code{\link{dec.table}})
#' @param start.level starting dose level. Defaults to 1, i.e. the lowest dose level.
#' @param nsim the number of simulation trials. Defaults to 1000.
#' @return the functions \code{\link{summary.dec.sim}} is used to obtain and print a summary table of the results (recommended). An object of class \code{"dec.sim"} is a list containing:
#'  \item{mtd}{a vector of dose levels giving the recommended maximum tolerated dose (MTD) at the end of the trial.}
#'  \item{mtd.prob}{a vector of length \code{k} giving the average proportions of selected as MTD at each dose level.}
#'  \item{over.prob}{a vector of length \code{k} giving the average proportions of selected as over the MTD at each dose level.}
#'  \item{n.patients}{the average number of patients dosed at each level.}
#'  \item{dlt}{the average number of DLTs experienced at each dose level.}
#'  \item{truep}{input; true probabilities of toxicity.}
#'  \item{start.level}{input; starting dose level.}
#'  \item{nsim}{input; number of simulated trails.}
#' @author Wenchuan Guo <wguo007@ucr.edu>
#' @import stats
#' @export
#' @examples
#' truep <- c(0.3, 0.45, 0.5, 0.6)
#' res <- dec.table(0.6,0.4,0.2,0.3,c(3,3,3))
#' out <- dec.sim(truep, res$table, start.level = 2, nsim=1000)
#' summary(out, pt = 0.3)

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
  np <- dlt <- des <- over <- matrix(0, nsim, n_dose)
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
      des[i, dose] <- as.character(decTable[dlt[i, dose]+1, sample.stage[dose]])
      sample.stage[dose] <- sample.stage[dose] + 1
      # case: stay (S)
      if(des[i, dose] == "S"){
        if(np[i, dose] != maxn) {
          dose <- dose
        } 
        else {
          mtd[i] <- dose
        } 
      }
      # case : escalate (E)
      if(des[i, dose] == "E") {
        if(dose != n_dose) {
          # check if next dose level is available
          if((dose+1) %in% rm.dose) {
            if(np[i, dose] != maxn) {
              dose <- dose
            } else {
              mtd[i] <- dose
            }
          } else {
            if(np[i, dose + 1] != maxn) {
              dose <- dose + 1
            } else {
              mtd[i] <- dose
            }
          }
        } else {
          mtd[i] <- "U"
        }
      }  
      # case : de-escalate (D)
      if(des[i, dose] == "D") {
        if(dose != 1) {
          if(np[i, dose-1] != maxn) {
            dose <- dose - 1
          } else {
            mtd[i] <- dose -1 
          }
        } else {
          mtd[i] <- "L"
        }
      }
      # case : de-escalate (DU)
      if(des[i, dose] == "DU") {
        if(dose != 1) {
          # can not go back to this dose again
          rm.dose <- c(rm.dose, dose)
          if(np[i, dose-1] != maxn) {
            dose <- dose - 1
          } else {
            mtd[i] <- dose -1 
          }
        } else {
          mtd[i] <- "L"
        }
      }
    }		
  }
  # calculate mtd probabilty 
  mtd.prob <- sapply(c(1:n_dose, "L", "U"), function(ii) mean(mtd == ii))
  # calculate power
  over <- matrix(0, nsim, n_dose)
  over[des == "D" | des == "DU"] <- 1
  over[des == "S" | des == "E"] <- 0.5
  max.over <- apply(over, 1, max)
  for(i in 1 : nsim) {
    if(max.over[i] == 1) {
      iid <- which(over[i, ] == 1)[1]
    } else {
      iid <- min(max(which(over[i, ] == 0.5)) + 1, n_dose)
    }
    over[i, iid : n_dose] <- 1
  }
  over[over == 0.5] <- 0
  over.prob <- colMeans(over)
  out <- list(mtd = mtd, mtd.prob = mtd.prob, dlt = colMeans(dlt), n.patients = colMeans(np), over.prob = over.prob, truep = truep, start.level = start.level, nsim = nsim)
  class(out) <- "dec.sim"
  return(out)
}
