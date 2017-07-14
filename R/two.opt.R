#' @keywords internal
two.opt <- function(alpha1, alpha2, pc, n, sf = "Pocock", ...){
  # initialization
  nc <- cumsum(n)
  nt <- nc[3]
  if(sf == "Pocock" | sf == "OF") {
    as_left <- cumsum(gsDesign(test.type = 1, alpha = alpha1, timing = nc/nt, sfu = sf)$upper$spend)
    as_right <- cumsum(gsDesign(test.type = 1, alpha = alpha2, timing = nc/nt, sfu = sf)$upper$spend)
  } else {
    as_left <- sf(alpha1, timing = nc/nt)
    as_right <- sf(alpha2, timing = nc/nt)
  }
  comb <- NULL
  err <- NULL
  beta1 <- beta2 <- NULL
  # boundary of r1: [0, n1]
  r1_bdry <- 0:nc[1]
  out_r1 <- pbinom(r1_bdry, n[1], pc[1])
  ind_r1 <- which(out_r1 <= as_left[1])
  # check if conditions hold
  if(length(ind_r1)==0) {
    stop("No optimal design")
  } else {
    r1 <- r1_bdry[ind_r1]
    out_r1 <- out_r1[ind_r1]
  }
  # loop over all r1
  for(r1i in 1:length(r1)){
    # boundary of s1 (r1, n1-1]
    s1_bdry <- r1[r1i]:(n[1]-1)
    out_s1 <- 1 - pbinom(s1_bdry, n[1], pc[2])
    ind_s1 <- which(out_s1 <= as_right[1])
    # check if condition holds
    if(length(ind_s1) == 0) {
      next
    } else {
      s1 <- s1_bdry[ind_s1]
      out_s1 <- out_s1[ind_s1]
    }
    # loop over all s1 given r1
    for(s1i in 1:length(s1)){
      # boundary of r1: [r1, n1+n2-1]
      r2_bdry <- r1[r1i]:(nc[2]-1)
      t1 <- (r1[r1i]+1):s1[s1i]
      out_r2 <- sapply(r2_bdry, function(r2) sum(dbinom(t1, n[1], pc[1]) * pbinom(r2-t1, n[2], pc[1]))) + pbinom(r1[r1i], n[1], pc[1])
      ind_r2 <- which(out_r2 <= as_left[2])
      # check if conditions hold
      if(length(ind_r2) == 0) {
        next
      } else {
        r2 <- r2_bdry[ind_r2]
        out_r2 <- out_r2[ind_r2]
      }
      # loop over all r2 given r1, s1
      for(r2i in 1:length(r2)){
        # boundary of s2
        s2_bdry <- r2[r2i]:nc[2]
        out_s2 <- sapply(s2_bdry, function(s2) sum(dbinom(t1, n[1], pc[2]) * (1-pbinom(s2-t1, n[2], pc[2])))) + 1-pbinom(s1[s1i], n[1], pc[2])
        ind_s2 <- which(out_s2 <= as_right[2])
        # check if conditions hold
        if(length(ind_s2) == 0) {
          next
        } else {
          s2 <- s2_bdry[ind_s2]
          out_s2 <- out_s2[ind_s2]
        }
        # loop over all s2 give r2, r1, s1
        for(s2i in 1:length(s2)){
          # boundary of r3: [r2, n1+n2+n3=n]
          r3_bdry <- r2[r2i]:nc[3]
          out_r3 <- sapply(t1, function(tt){
            t2 <- (r2[r2i]-tt+1):(s2[s2i]-tt)
            return(sapply(r3_bdry, function(r3) dbinom(tt, n[1], pc[1]) * sum(dbinom(t2, n[2], pc[1]) * pbinom(r3-tt-t2, n[3], pc[1]))))
          })
          # check if nrow(out_r3)==1
          if(is.vector(out_r3)) {
            out_r3 <- sum(out_r3) + out_r2[r2i]
          }  else {
            out_r3 <- rowSums(out_r3) + out_r2[r2i]
          }
          ind_r3 <- which(out_r3 <= as_left[3])
          # check if conditions hold
          if(length(ind_r3) == 0) {
            next
          } else {
            r3 <- r3_bdry[ind_r3]
            out_r3 <- out_r3[ind_r3]
          }
          # loop over all r3 given r1, r2, s1, s2
          for(r3i in 1:length(r3)){
            # boundary of s3
            s3_bdry <- r3[r3i]:(nc[3]-1)
            out_s3 <- sapply(t1, function(tt){
              t2 <- (r2[r2i]-tt+1):(s2[s2i]-tt)
              return(sapply(s3_bdry, function(s3) dbinom(tt, n[1], pc[2]) * sum(dbinom(t2, n[2], pc[2]) * (1-pbinom(s3-tt-t2, n[3], pc[2])))))
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
            beta1 <- sapply(t1, function(tt){
              t2 <- (r2[r2i]-tt+1):(s2[s2i]-tt)
              return(sapply(s3, function(ss3) dbinom(tt, n[1], pc[1]-0.2) * sum(dbinom(t2, n[2], pc[1]-0.2) * (1-pbinom(ss3-tt-t2, n[3], pc[1]-0.2)))))
            })
            beta2 <- sapply(t1, function(tt){
              t2 <- (r2[r2i]-tt+1):(s2[s2i]-tt)
              return(sapply(s3, function(ss3) dbinom(tt, n[1], pc[2]+0.2) * sum(dbinom(t2, n[2], pc[2]+0.2) * pbinom(ss3-tt-t2, n[3], pc[2]+0.2))))
            })
            if(is.vector(beta1)){
              beta1 <- sum(beta1)
              beta2 <- sum(beta2)
            } else {
              beta1 <- rowSums(beta1)
              beta2 <- rowSums(beta2)
            }
            # save feasible designs & errors
            comb <- rbind(comb, cbind(r1[r1i], r2[r2i], r3[r3i], s1[s1i], s2[s2i], s3))
            err <- rbind(err, cbind(out_r1[r1i], out_r2[r2i], out_r3[r3i], out_s1[s1i], out_s2[s2i], out_s3, beta1, beta2))
          }
        }
      }
    }
  }
  # merge results
  if(is.null(comb)){
    stop("No optimal design")
  } else {
    rs <- c("r1", "r2", "r3", "s1", "s2", "s3")
    err1 <- c("alpha11", "alpha12", "alpha13", "alpha21", "alpha22", "alpha23")
    err2 <- c("beta1", "beta2")
    out <- cbind(comb, err)
    out <- out[order(-out[,7], -out[,10], -out[,8], -out[,11], -out[,9], -out[,12]), ]
    colnames(out) <- c(rs, err1, err2)
  }
  out <- list(out = out, pc = pc, n = n, alpha = c(alpha1, alpha2), sf = sf)
  class(out) <- "2opt"
  return(out)
}
