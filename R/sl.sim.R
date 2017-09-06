#' Dose-finding simulations for a list of scenarios
#' @description Run dose-finding simulations based on a customized decision table for a list of scenarios.
#' @param decTable A customized decision table. (same format as output of \code{dec.table})
#' @param file The name of the file which the data are to be read from. See details in \code{\link{read.table}}.
#' @param header A logical value indicating whether the file contains the names of the variables as its first line. Default is \code{FALSE}. See details in \code{\link{read.table}}.
#' @param sep The field separator character. Default is \code{","}. See details in \code{\link{read.table}}.
#' @param ... arguments to be passed to \code{\link{read.table}} methods.
#' @details In each line of the input file, the parameters must be ordered in accordance as follows: \code{pt}, \code{start.level}, \code{nsim}, \code{truep}. See details in \code{\link{read.table}}. The algorithm for dose-finding is described in \code{\link{dec.sim}}.
#' @return The function \code{\link{summary}} is used to obtain and print a summary table of the results. An object of class \code{"dec.sim"} (1 scenario) or \code{"sl.sim"} (more than 1 scenarios)is a list containing:
#'  \item{MTD}{A vector of dose levels giving the recommended maximum tolerated dose (MTD) at the end of the trial.}
#'  \item{n.patients}{The average number of patients dosed at each level.}
#   \item{alpha.one}{a vector of true type 1 error for right-tailed test.}
#'  \item{truep}{input; true probabilities of toxicity.}
#'  \item{start.level}{input; starting dose level.}
#'  \item{nsim}{input; number of simulated trails.}
#' @author Wenchuan Guo <wguo007@ucr.edu>
#' @import stats
#' @import utils
#' @export
#' @examples
#' dt <- dec.table(0.6,0.4,0.2,0.3,c(3,3,3))
#' test.file <- system.file("extdata", "testS.csv", package = "tsdf")
#' # use a customized decision table
#' table.file <- system.file("extdata", "decTable.csv", package = "tsdf")
#' dec <- read.table(table.file, sep=",", col.names=c(3,4,8,10), row.names = 1, check.names = FALSE)
#' out1 <- sl.sim(dt$table, test.file)
#' out2 <- sl.sim(dec, test.file)

sl.sim <- function(decTable, file, header = TRUE, sep = ",", ...) {
  sl <- read.table(file, header, sep, ...)
  out <- lapply(1:nrow(sl), function(i){
    info <- sl[i, ][!is.na(sl[i, ])]
    n.info <- length(info)
    dec.sim(info[3:n.info], decTable, info[1], info[2])
  })
  if(nrow(sl) > 1) {
    class(out) <- c("sl.sim", "dec.sim")
  } else {
    class(out) <- "dec.sim"
  }
  return(out)
}

