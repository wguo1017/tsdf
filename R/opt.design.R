#' Zhong's 2-/3- stage Phase II design
#' @description calculate optimal 2-/3-stage design given by Bob Zhong
#' @details The two-stage design setup is: n1 patients are treated in the first stage. At the end of the first stage, either the trial continues to the second stage or inefficacy is concluded and the trial is stopped (early termination), depending on the number of responses observed at the first stage. If the trial does continue to the second stage, additional n2 patients are treated. Three-stage design is an extension of two-stage design where one stage is added between Stage 1 and 2. The left-side rejection region is response <= r_i for i = 1, 2, or 3 and right-side rejection region is  response > s. Alpha-spending method is added to two-/three-stage designs. \code{opt.design} supports Hwang-Shih-DeCani spending function. You can change the definition of \code{HSD} function to use a different spending function. 
#' @param alpha1 left-side overall type I error. 
#' @param alpha2 right-side overall type I error.
#' @param beta type II error
#' @param pc a numeric vector of response rate. should be a vector with length 1 or 2.
#' @param pe alternative hypothesis.
#' @param stage 2 or 3. default to 2 (2-stage design).
#' @param stop.eff logical flag. default to \code{FALSE}. if \code{stop.eff = TRUE}, the trial may stop for efficacy at interim.
#' @param frac_n1 proportion of n1. for 2-stage design, default to \code{c(0.3, 0.6)}, i.e. the range of \code{n1} is \code{0.2*n} to \code{0.5*n}. for 3-stage design, default to \code{c(0.2, 0.3)}, i.e. the range of \code{n1} is \code{0.2*n} to \code{0.3*n}
#' @param frac_n2 proportion of n2. Used for 3-stage design. default to \code{c(0.2, 0.4)}.
#' @param sf.param  a single real value specifying the gamma parameter for which Hwang-Shih-DeCani spending is to be computed; allowable range is [-40, 40]. Increasing this parameter implies that more error is spent early stage and less is available in late stage. For two-stage designs, default to \code{NULL}(alpha-spending is not used); for three-stage designs, default to 4.
#' @param show logical. If \code{TRUE}, current sample size is shown as total sample size increase. 
#' @param nmax maximum sample size. default to 100.
#' @param n.choice stop criterion for the search of feasible designs. stop if number of designs is more than \code{n.choice}
#' @param ... not used argument.
#' @return An object of class "opt.design" is a list containing:
#'  \item{bdry}{rejection regions}
#'  \item{error}{true type 1/2 errors}
#'  \item{n}{sample size at each stage} 
#'  \item{complete}{complete list of feasible designs}
#'  \item{alpha1}{input; left-side type 1 error}
#'  \item{alpha2}{input; right-side type 1 error}
#'  \item{beta}{input; type 2 error}
#'  \item{pc}{input; a vector of response rate.}
#'  \item{pe}{input; a vector of alternative response rate}
#'  \item{sf}{input; the alpha-spending function used}
#'  \item{stage}{input; two- or three- stage design is used}
#' @author Wenchuan Guo <wguo1017@gmail.com>, Jianan Hui <jiananhuistat@gmail.com>
#' @references
#' Zhong. (2012) Single-arm Phase IIA clinical trials with go/no-go decisions. \emph{Contemporary Clinical Trials}, \bold{33}, 1272--1279.
#' @import stats
#' @export
#' @examples 
#'  alpha1 <- 0.15
#'  alpha2 <- 0.10
#'  beta <- 0.15
#'  pc <- 0.25
#'  pe <- pc + 0.20
#'  # calculate optimal two-stage design without using alpha-spending
#'  opt.design(alpha1, alpha2, beta, pc, pe, stage=2)
#'  \dontrun{
#'  # calculate optimal two-stage design with Pocock-like spending function 
#'  opt.design(alpha1, alpha2, beta, pc, pt, stage = 2, sf.param = 1)
#'  
#'  # calculate optimal three-stage design with =O’Brien-Fleming like spending function
#'  opt.design(alpha1, alpha2, beta, pc, pt, stage = 3, sf.param = -4)
#'  }

opt.design <- function(alpha1, alpha2, beta, pc, pe, stage = 2, stop.eff = FALSE, frac_n1 = NULL, frac_n2 = NULL, sf.param = NULL, show = FALSE, nmax = 100, n.choice = 1, ...) {
  if(stage !=2 & stage !=3){
    stop("only support two and three stage designs")
  }
  if(alpha1 > 1 | alpha1 < 0){
    stop("'alpha1' should between 0 and 1")
  }
  if(alpha2 > 1 | alpha2 < 0){
    stop("'alpha2' should between 0 and 1")
  }
  if(beta > 1 | beta < 0){
    stop("'beta' should between 0 and 1")
  }
  if(stage == 2) {
    if(is.null(frac_n1)) frac_n1 <- c(0.3, 0.6)
    out <- zhong.two(alpha1, alpha2, beta, pc, pe, stop.eff, sf.param, show, nmax, n.choice, frac_n1)
  } 
  if(stage == 3) {
    if(is.null(sf.param)) sf.param <- 4
    if(is.null(frac_n1)) frac_n1 <- c(0.2, 0.3)
    if(is.null(frac_n2)) frac_n2 <- c(0.2, 0.4)
    out <- zhong.three(alpha1, alpha2, beta, pc, pe, frac_n1, frac_n2, sf.param, stop.eff, show, nmax)
  }
  input <- list(alpha1 = alpha1, alpha2 = alpha2, beta = beta, pc = pc, pe = pe, sf.param = sf.param, stage = stage, frac_n1 = frac_n1, frac_n2 = frac_n2)
  out <- c(out, input)
  class(out) <- "opt.design"
  return(out)
}

