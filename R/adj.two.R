#' Zhong's 2-/3- stage Phase II design
#' @description adjust Zhong's 2-/3-stage design for over-/under-running
#' @details To be added
#' @param n1 sample size at stage 1.
#' @param r1 inefficacy boundary at stage 1.
#' @param s1 efficacy boundary at stage 1. if no early stopping for efficacy, \code{s1} should equal to \code{n1}.
#' @param n2 sample size at stage 2.
#' @param alpha1 left-side overall type I error. 
#' @param alpha2 right-side overall type I error.
#' @param beta type II error.
#' @param pc a numeric vector of response rate. should be a vector with length 1 or 2.
#' @param pe alternative hypothesis.
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
#' @import stats
#' @export
#' @examples 
#'  n1 <- 22
#'  r1 <- 6
#'  s1 <- 22
#'  n2 <- 24
#'  pc <- 0.4
#'  pe <- pc + 0.15
#'  alpha1 <- 0.3
#'  alpha2 <- 0.1
#'  beta <- 0.2
#'  out <- adj.two(n1, r1, s1, n2, alpha1, alpha2, beta, pc, pe)


adj.two <- function(n1, r1, s1, n2, alpha1, alpha2, beta, pc, pe, ...){
  if(length(pc) == 1) {
    pc <- rep(pc, 2)
  }
  comb <- NULL
  err <- NULL
  L1 <- pbinom(r1, n1, pc[1])
  R1 <- 1 - pbinom(s1, n1, pc[2])
  t1 <- (r1+1) : s1
  for(r2 in (s1 + n2) : (r1+1)){
    L2 <- L1 + sum(dbinom(t1, n1, pc[1]) * pbinom(r2 - t1, n2, pc[1]))
    if(L2 <= alpha1) {
      s2_l <- ifelse(s1!=n1, max(r2, s1), r2)
      for(s2 in s2_l : (s1+n2)) {
        R2 <- R1 + sum(dbinom(t1, n1, pc[2]) * (1 - pbinom(s2 - t1, n2, pc[2])))
        if(R2 <= alpha2) {
          Rpe <- sum(dbinom(t1, n1, pe) * (1 - pbinom(s2 - t1, n2, pe))) + 1 - pbinom(s1, n1, pe)
          comb <- rbind(c(r1, r2, s1, s2, n1, n2), comb)
          err <- rbind(c(L1, L2, R1, R2, 1- Rpe), err)
        } 
        else next
      }
    } else next 
  }
  out <- round(cbind(err, comb), 4)
  # out <- out[out[, 5] < (beta + 0.2), ]
  out <- out[order(-out[, 2], -out[, 4], out[, 5]), ]
  out <- as.matrix(do.call(rbind, by(out, out[, 7], FUN=function(x) head(x, 1))))
  colnames(out) <- c("alpha11", "alpha12", "alpha21", "alpha22", "beta", "r1", "r2", "s1", "s2", "n1", "n2")
  out <- out[order(-out[, 2], -out[, 4], out[, 5]), ]
  opt <- out[1, ]
  res <- list(bdry = opt[6:9], error = opt[1:5], n = opt[10:11], complete = out)
  input <- list(alpha1 = alpha1, alpha2 = alpha2, beta = beta, pc = pc, pe = pe, stage = 2)
  res <- c(res, input)
  class(res) <- "opt.design"
  return(res)
}