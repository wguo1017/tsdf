#' @keywords internal
right.opt <- function(alpha, pc, n, sf = "Pocock", ...){
  # initialization
  nc <- cumsum(n)
  nt <- nc[3]

  if(sf == "Pocock" | sf == "OF") {
    as_right <- cumsum(gsDesign(test.type = 1, alpha = alpha, timing = nc/nt, sfu = sf)$upper$spend)
  } else {
    as_right <- sf(alpha, timing = nc/nt)
  }
  comb <- NULL
  err <- NULL
  beta <- NULL
  # boundary of s1 (0, n1-1]
  s1_bdry <- 0:(n[1]-1)
  out_s1 <- 1-pbinom(s1_bdry, n[1], pc)
  ind_s1 <- which(out_s1 <= as_right[1])
  # check if condition holds
  if(length(ind_s1)==0) {
    return(print("No optimal design (right side)"))
  } else {
    s1 <- s1_bdry[ind_s1]
    out_s1 <- out_s1[ind_s1]
  }
  # loop over all s1
  for(s1i in 1:length(s1)){
    t1 <- 0:s1[s1i]
    # boundary of s2
    s2_bdry <- s1[s1i]:(nc[2]-1)
    out_s2 <- sapply(s2_bdry, function(s2) sum(dbinom(t1, n[1], pc) * (1 - pbinom(s2-t1, n[2], pc)))) + out_s1[s1i]
    ind_s2 <- which(out_s2 <= as_right[2])
    # check if conditions hold
    if(length(ind_s2) == 0) {
      next
    } else {
      s2 <- s2_bdry[ind_s2]
      out_s2 <- out_s2[ind_s2]
    }
    # loop over all s2 given s1
    for(s2i in 1:length(s2)){
      # boundary of s3
      s3_bdry <- s2[s2i]:(nc[3]-1)
      out_s3 <- sapply(t1, function(tt){
        t2 <- 0:(s2[s2i]-tt)
        return(sapply(s3_bdry, function(s3) dbinom(tt, n[1], pc)  *sum(dbinom(t2, n[2], pc)*(1 - pbinom(s3-tt-t2, n[3], pc)))))
      })
      # check if nrow(out_s3)==1
      if(is.vector(out_s3)) {
        out_s3 <- sum(out_s3) + out_s2[s2i]
      }  else {
        out_s3 <- rowSums(out_s3) + out_s2[s2i]
      }
      ind_s3 <-  which(out_s3 <= as_right[3])
      # check if conditions hold
      if(length(ind_s3) == 0) {
        next
      } else {
        s3 <- s3_bdry[ind_s3]
        out_s3 <- out_s3[ind_s3]
      }
      # calculate type 2 error
      beta <- sapply(t1, function(tt){
        t2 <- 0:(s2[s2i]-tt)
        return(sapply(s3, function(ss3) dbinom(tt, n[1], pc+0.2) * sum(dbinom(t2, n[2], pc+0.2) * pbinom(ss3-tt-t2, n[3], pc+0.2))))
      })
      if(is.vector(beta)){
        beta <- sum(beta)
      } else {
        beta <- rowSums(beta)
      }
      # save feasible designs & errors
      comb <- rbind(comb, cbind(s1[s1i], s2[s2i], s3))
      err <- rbind(err, cbind(out_s1[s1i], out_s2[s2i], out_s3, beta))
    }
  }
  # merge results
  if(is.null(comb)){
    return(print("No optimal design (right side)"))
  } else {
    rs <- c("s1", "s2", "s3")
    err1 <- c("alpha11", "alpha12", "alpha13")
    err2 <- c("beta")
    out <- cbind(comb, err)
    out <- out[order(-out[,4], -out[,4], -out[,6], out[,7]), ]
    colnames(out) <- c(rs, err1, err2)
  }
  out <- list(out = out, pc = pc, n = n, sf = sf, alpha = alpha)
  class(out) <- "1opt"
  return(out)
}
