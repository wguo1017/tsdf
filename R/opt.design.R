#' Zhong's 2-/3- stage Phase II design
#' @description calculate optimal 2-/3-stage design given by Bob Zhong
#' @param alpha1 left-side overall type I error. 
#' @param alpha2 right-side overall type I error.
#' @param beta type II error
#' @param pc a numeric vector of response rate. Should be a vector with 1 or 2.
#' @param pt alternative hypothesis
#' @param stage 2 or 3. default to 2 (2-stage design).
#' @param frac_n1 proportion of n1. Used for 3-stage design. default to \code{c(0.2, 0.3)}, i.e. the range of \code{n1} is \code{0.2*n} to \code{0.3*n}
#' @param frac_n2 proportion of n2. Used for 3-stage design. default to \code{c(0.2, 0.4)}.
#' @param sf.param  A single real value specifying the gamma parameter for which Hwang-Shih-DeCani spending is to be computed; allowable range is [-40, 40]. Details in\code{\link{gsDesign}}. for two-stage designs, default to \code{NULL}, alpha-spending is not used; for three-stage designs, default to 4.
#' @param show logical. If \code{TRUE}, current sample size is shown as total sample size increase. 
#' @param nmax maximum sample size. default to 100.
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
#'  \item{pt}{input; a vector of alternative response rate}
#'  \item{sf}{input; the alpha-spending function used}
#'  \item{stage}{input; two- or three- stage design is used}
#' @author Wenchuan Guo <wguo007@ucr.edu>
#' @references
#' Zhong. (2012) Single-arm Phase IIA clinical trials with go/no-go decisions. \emph{Contemporary Clinical Trials}, \bold{33}, 1272--1279.
#' @import stats
#' @export
#' @examples 
#' \dontrun{
#'  alpha1 <- 0.15
#'  alpha2 <- 0.10
#'  beta <- 0.15
#'  pc <- 0.25
#'  pt <- pc + 0.20
#'  opt.design(alpha1, alpha2, beta, pc, pt, stage=2, show=1)
#'  opt.design(alpha1, alpha2, beta, pc, pt, stage=3, show=1)
#' }

opt.design <- function(alpha1, alpha2, beta, pc, pt, stage = 2, frac_n1 = c(0.2, 0.3), frac_n2 = c(0.2,0.4), sf.param = NULL, show = FALSE, nmax = 100, ...) {
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
  if(frac_n1[2] + frac_n2[2] > 1){
    warning(" n1+n2 can not greater than total sample size, set to default fractions")
    frac_n1 <- c(0.2, 0.3)
    frac_n2 <- c(0.2,0.4)
  }
  if(stage == 2) {
    out <- zhong.two(alpha1, alpha2, beta, pc, pt, sf.param, show, nmax = 100)
  } else {
    if(is.null(sf.param)) sf.param <- 4
    out <- zhong.three(alpha1, alpha2, beta, pc, pt, frac_n1, frac_n2, sf.param, show, nmax)
  }
  input <- list(alpha1 = alpha1, alpha2 = alpha2, beta = beta, pc = pc, pt = pt, sf.param = sf.param, stage = stage)
  out <- c(out, input)
  class(out) <- "opt.design"
  return(out)
}

