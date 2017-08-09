#' @keywords internal
three.opt <- function(alpha1, alpha2, pc, n, sf = "Pocock", ...){
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
  # boundary of r1: [0, n1]
  r1_bdry <- 0:nc[1]
  out_r1 <- unique(pbinom(r1_bdry, n[1], pc[1]))
  ind_r1 <- which(out_r1 <= as_left[1])
  # check if conditions hold
  if(length(ind_r1)==0) stop("No optimal design (left side)")
  r1 <- r1_bdry[max(ind_r1)]
  out_r1 <- out_r1[max(ind_r1)]
  # boundary of s1 (r1, n1-1]
  s1_bdry <- r1:(n[1]-1)
  out_s1 <- unique(1 - pbinom(s1_bdry, n[1], pc[2]))
  ind_s1 <- which(out_s1 <= as_right[1])
  # check if condition holds 
  if(length(ind_s1) == 0) stop("No optimal design (right side)")
  s1 <- s1_bdry[min(ind_s1)]
  out_s1 <- out_s1[min(ind_s1)]
  # boundary of r1: [r1, n1+n2-1]
  r2_bdry <- r1:(nc[2]-1)
  t1 <- (r1+1):s1
  out_r2 <- sapply(r2_bdry, function(r2) sum(dbinom(t1, n[1], pc[1]) * pbinom(r2-t1, n[2], pc[1]))) + out_r1
  out_r2 <- unique(out_r2)
  ind_r2 <- which(out_r2 <= as_left[2])
  # check if conditions hold
  if(length(ind_r2) == 0) stop("No optimal design (left side)")
  r2 <- r2_bdry[max(ind_r2)]
  out_r2 <- out_r2[max(ind_r2)]
  # boundary of s2
  s2_bdry <- r2:nc[2]
  out_s2 <- sapply(s2_bdry, function(s2) sum(dbinom(t1, n[1], pc[2]) * (1-pbinom(s2-t1, n[2], pc[2])))) + out_s1
  out_s2 <- unique(out_s2)
  ind_s2 <- which(out_s2 <= as_right[2])
  # check if conditions hold
  if(length(ind_s2) == 0) stop("No optimal design (right side)")
  s2 <- s2_bdry[min(ind_s2)]
  out_s2 <- out_s2[min(ind_s2)]
  # boundary of r3: [r2, n1+n2+n3=n]
  r3_bdry <- r2:nc[3]
  out_r3 <- sapply(t1, function(tt){
    t2 <- (r2-tt+1):(s2-tt)
    return(sapply(r3_bdry, function(r3) dbinom(tt, n[1], pc[1]) * sum(dbinom(t2, n[2], pc[1]) * pbinom(r3-tt-t2, n[3], pc[1]))))
  })
  # check if nrow(out_r3)==1
  if(is.vector(out_r3)) {
    out_r3 <- sum(out_r3) + out_r2
  }  else {
    out_r3 <- rowSums(out_r3) + out_r2
  } 
  out_r3 <-unique(out_r3)
  ind_r3 <- which(out_r3 <= as_left[3])
  # check if conditions hold
  if(length(ind_r3) == 0) stop("No optimal design (left side)")
  r3 <- r3_bdry[max(ind_r3)]
  out_r3 <- out_r3[max(ind_r3)]
  # boundary of s3
  s3_bdry <- r3:(nc[3]-1)
  out_s3 <- sapply(t1, function(tt){
    t2 <- (r2-tt+1):(s2-tt)
    return(sapply(s3_bdry, function(s3) dbinom(tt, n[1], pc[2]) * sum(dbinom(t2, n[2], pc[2]) * (1-pbinom(s3-tt-t2, n[3], pc[2])))))
  })
  # check if nrow(out_s3)==1
  if(is.vector(out_s3)) {
    out_s3 <- sum(out_s3) + out_s2
  }  else {
    out_s3 <- rowSums(out_s3) + out_s2
  } 
  out_s3 <- unique(out_s3)
  ind_s3 <-  which(out_s3 <= as_right[3])
  # check if conditions hold
  if(length(ind_s3) == 0) stop("No optimal design (right side)")
  s3 <- s3_bdry[min(ind_s3)]
  out_s3 <- out_s3[min(ind_s3)]
  # save feasible designs & errors
  bdry <- c(r1, r2, r3, s1, s2, s3)
  err <- c(out_r1, out_r2, out_r3, out_s1, out_s2, out_s3)
  # merge results
  names(bdry) <- c("r1", "r2", "r3", "s1", "s2", "s3")
  names(err) <- c("alpha11", "alpha12", "alpha13", "alpha21", "alpha22", "alpha23")
  out <- list(bdry = bdry, error = err, pc = pc, n = n, alpha = c(alpha1, alpha2), sf = sf)
  class(out) <- "2opt"
  return(out)
}

#' @keywords internal
two.opt <- function(alpha1, alpha2, pc, n, sf = "Pocock", ...){
  # initialization
  nc <- cumsum(n)
  nt <- nc[2]
  if(sf == "Pocock" | sf == "OF") {
    as_left <- cumsum(gsDesign(k = 2, test.type = 1, alpha = alpha1, timing = nc/nt, sfu = sf)$upper$spend)
    as_right <- cumsum(gsDesign(k = 2, test.type = 1, alpha = alpha2, timing = nc/nt, sfu = sf)$upper$spend)
  } else {
    as_left <- sf(alpha1, timing = nc/nt)
    as_right <- sf(alpha2, timing = nc/nt)
  }
  comb <- NULL
  err <- NULL
  # boundary of r1: [0, n1]
  r1_bdry <- 0:nc[1]
  out_r1 <- unique(pbinom(r1_bdry, n[1], pc[1]))
  ind_r1 <- which(out_r1 <= as_left[1])
  # check if conditions hold
  if(length(ind_r1)==0) stop("No optimal design (left side)")
  r1 <- r1_bdry[max(ind_r1)]
  out_r1 <- out_r1[max(ind_r1)]
  # boundary of s1 (r1, n1-1]
  s1_bdry <- r1:(n[1]-1)
  out_s1 <- unique(1 - pbinom(s1_bdry, n[1], pc[2]))
  ind_s1 <- which(out_s1 <= as_right[1])
  # check if condition holds 
  if(length(ind_s1) == 0) stop("No optimal design (right side)")
  s1 <- s1_bdry[min(ind_s1)]
  out_s1 <- out_s1[min(ind_s1)]
  # boundary of r1: [r1, n1+n2-1]
  r2_bdry <- r1:(nc[2]-1)
  t1 <- (r1+1):s1
  out_r2 <- sapply(r2_bdry, function(r2) sum(dbinom(t1, n[1], pc[1]) * pbinom(r2-t1, n[2], pc[1]))) + out_r1
  out_r2 <- unique(out_r2)
  ind_r2 <- which(out_r2 <= as_left[2])
  # check if conditions hold
  if(length(ind_r2) == 0) stop("No optimal design (left side)")
  r2 <- r2_bdry[max(ind_r2)]
  out_r2 <- out_r2[max(ind_r2)]
  # boundary of s2
  s2_bdry <- r2:nc[2]
  out_s2 <- sapply(s2_bdry, function(s2) sum(dbinom(t1, n[1], pc[2]) * (1-pbinom(s2-t1, n[2], pc[2])))) + out_s1
  out_s2 <- unique(out_s2)
  ind_s2 <- which(out_s2 <= as_right[2])
  # check if conditions hold
  if(length(ind_s2) == 0) stop("No optimal design (right side)")
  s2 <- s2_bdry[min(ind_s2)]
  out_s2 <- out_s2[min(ind_s2)]
  bdry <- c(r1, r2, s1, s2)
  err <- c(out_r1, out_r2, out_s1, out_s2)
  # merge results
  names(bdry) <- c("r1", "r2", "s1", "s2")
  names(err) <- c("alpha11", "alpha12", "alpha21", "alpha22")
  out <- list(bdry= bdry, error = err, pc = pc, n = n, alpha = c(alpha1, alpha2), sf = sf)
  class(out) <- "2opt"
  return(out)
}

#' @keywords internal
right.two.opt <- function(alpha, pc, n, sf = "Pocock", ...){
  # initialization
  nc <- cumsum(n)
  nt <- nc[2]
  
  if(sf == "Pocock" | sf == "OF") {
    as_right <- cumsum(gsDesign(k = 2, test.type = 1, alpha = alpha, timing = nc/nt, sfu = sf)$upper$spend)
  } else {
    as_right <- sf(alpha, timing = nc/nt)
  }
  # boundary of s1 (0, n1-1]
  s1_bdry <- 0:(n[1]-1)
  out_s1 <- 1-pbinom(s1_bdry, n[1], pc)
  ind_s1 <- which(out_s1 <= as_right[1])
  # check if condition holds 
  if(length(ind_s1)==0) stop("No optimal design (right side)")
  s1 <- s1_bdry[min(ind_s1)]
  out_s1 <- out_s1[min(ind_s1)]
  t1 <- 0:s1
  # boundary of s2
  s2_bdry <- s1:(nc[2]-1)
  out_s2 <- sapply(s2_bdry, function(s2) sum(dbinom(t1, n[1], pc) * (1 - pbinom(s2-t1, n[2], pc)))) + out_s1
  out_s2 <- unique(out_s2)
  ind_s2 <- which(out_s2 <= as_right[2])
  # check if conditions hold
  if(length(ind_s2) == 0) stop("No optimal design (right side)")
  s2 <- s2_bdry[min(ind_s2)]
  out_s2 <- out_s2[min(ind_s2)]
  bdry <- c(s1, s2)
  err <- c(out_s1, out_s2)
  # merge results
  names(bdry) <- c("s1", "s2")
  names(err) <- c("alpha11", "alpha12")
  out <- list(bdry = bdry, error = err, pc = pc, n = n, sf = sf, alpha = alpha)
  class(out) <- "1opt"
  return(out)
}

#' @keywords internal
right.three.opt <- function(alpha, pc, n, sf = "Pocock", ...){
  # initialization
  nc <- cumsum(n)
  nt <- nc[3]
  if(sf == "Pocock" | sf == "OF") {
    as_right <- cumsum(gsDesign(test.type = 1, alpha = alpha, timing = nc/nt, sfu = sf)$upper$spend)
  } else {
    as_right <- sf(alpha, timing = nc/nt)
  }
  # boundary of s1 (0, n1-1]
  s1_bdry <- 0:(n[1]-1)
  out_s1 <- 1-pbinom(s1_bdry, n[1], pc)
  ind_s1 <- which(out_s1 <= as_right[1])
  # check if condition holds 
  if(length(ind_s1)==0) stop("No optimal design (right side)")
  s1 <- s1_bdry[min(ind_s1)]
  out_s1 <- out_s1[min(ind_s1)]
  # loop over all s1 
  t1 <- 0:s1
  # boundary of s2
  s2_bdry <- s1:(nc[2]-1)
  out_s2 <- sapply(s2_bdry, function(s2) sum(dbinom(t1, n[1], pc) * (1 - pbinom(s2-t1, n[2], pc)))) + out_s1
  out_s2 <- unique(out_s2)
  ind_s2 <- which(out_s2 <= as_right[2])
  # check if conditions hold
  if(length(ind_s2) == 0) stop("No optimal design (right side)")
  s2 <- s2_bdry[min(ind_s2)]
  out_s2 <- out_s2[min(ind_s2)]
  # boundary of s3
  s3_bdry <- s2:(nc[3]-1)
  out_s3 <- sapply(t1, function(tt){
    t2 <- 0:(s2-tt)
    return(sapply(s3_bdry, function(s3) dbinom(tt, n[1], pc)  *sum(dbinom(t2, n[2], pc)*(1 - pbinom(s3-tt-t2, n[3], pc)))))
  })
  # check if nrow(out_s3)==1
  if(is.vector(out_s3)) {
    out_s3 <- sum(out_s3) + out_s2
  } else {
    out_s3 <- rowSums(out_s3) + out_s2
  } 
  out_s3 <- unique(out_s3)
  ind_s3 <-  which(out_s3 <= as_right[3])
  # check if conditions hold
  if(length(ind_s3) == 0) stop("No optimal design (right side)")
  s3 <- s3_bdry[min(ind_s3)]
  out_s3 <- out_s3[min(ind_s3)]
  # save feasible designs & errors
  bdry <- c(s1, s2, s3)
  err <- c(out_s1, out_s2, out_s3)
  # merge results
  names(bdry) <- c("s1", "s2", "s3")
  names(err) <- c("alpha11", "alpha12", "alpha13")
  out <- list(bdry = bdry, error = err, pc = pc, n = n, sf = sf, alpha = alpha)
  class(out) <- "1opt"
  return(out)
}