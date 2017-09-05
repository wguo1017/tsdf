#' run dose-finding simulations
#' @description Run dose-finding simulations based on a customized decision table.
#' @details Denote the number of patients treated by current dose level as m_i and the maximum number of objects treated at each dose level as m_max. The procedure is as follows
#' \itemize{
#' \item{Step 1 : }{Denote the dose level being used to treat patients as the current dose level(dose i). Accrue and treat k patients at the current dose level.}
#' \item{Step 2 : }{Record number of DLTs and use the decision table to make a decision: if decision is ``S" --> step 3; if decision is ``D" or ``DU'' --> step 4; if decision is ``E" --> step 4}
#' \item{Step 3 : }{If m_i = m_max, declare dose i as the MTD; otherwise, go to step 1. Need to de-escalate dose since the MTD has been exceeded. If current dose is the lowest dose, then stop the trail and declare the MTD is lower than the lowest dose level (inconclusive); If the next-lower level has m_max objects, then stop the trial and declare the MTD is the next lower dose level, otherwise, set the current dose level to be the next-lower dose level and go to step 2; If the decision is ``DU'', never use current dose level again.}
#' \item{Step 4 : }{Need to escalate dose level. If the current dose level is the highest dose level, then stop the trial and declare that the MTD is higher than the highest dose level (inconclusive); otherwise, set the current dose level to be next dose level and go to step 2.}
#' }
#' @param truep a vector of length k (the number of doses being considered in the trial), with values equal to the true probabilities of toxicity at the dose levels.
#' @param decTable a customized decision table. (same format as output of \code{\link{dec.table}})
#' @param start.level starting dose level. Defaults to 1, i.e. the lowest dose level.
#' @param nsim the number of simulation trials. Defaults to 1000.
#' @return the functions \code{\link{summary.dec.sim}} is used to obtain and print a summary table of the results (recommended). An object of class \code{"dec.sim"} is a list containing:
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
            mtd[i] <- "U"
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
          mtd[i] <- "L"
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
          mtd[i] <- "L"
        }
      }
    }		
  }
  mtd.prob <- sapply(c(1:n_dose, "L", "U"), function(ii) mean(mtd == ii))
  out <- list(mtd = mtd, mtd.prob = mtd.prob, dlt = colMeans(dlt), n.patients = colMeans(np), truep = truep, start.level = start.level, nsim = nsim)
  class(out) <- "dec.sim"
  return(out)
}
