#' @keywords internal
three.opt <- function(alpha1, alpha2, pt, n, sf.param, pe.par, ...){
  # initialization
  nc <- cumsum(n)
  nt <- nc[3]
  as_left <- HSD(alpha1, nc/nt, sf.param)
  as_right <- HSD(alpha2, nc/nt, sf.param)
  # boundary of r1: [0, n1]
  r1_bdry <- 0:nc[1]
  out_r1 <- unique(pbinom(r1_bdry, n[1], pt[1]))
  ind_r1 <- which(out_r1 <= as_left[1])
  # check if conditions hold
  if(length(ind_r1)==0) stop("No optimal design (left side)")
  r1 <- r1_bdry[max(ind_r1)]
  out_r1 <- out_r1[max(ind_r1)]
  # boundary of s1 (r1, n1-1]
  s1_bdry <- r1:(n[1]-1)
  out_s1 <- unique(1 - pbinom(s1_bdry, n[1], pt[2]))
  ind_s1 <- which(out_s1 <= as_right[1])
  # check if condition holds 
  if(length(ind_s1) == 0) stop("No optimal design (right side)")
  s1 <- s1_bdry[min(ind_s1)]
  out_s1 <- out_s1[min(ind_s1)]
  # boundary of r1: [r1, n1+n2-1]
  r2_bdry <- r1:(nc[2]-1)
  t1 <- (r1+1):s1
  out_r2 <- sapply(r2_bdry, function(r2) sum(dbinom(t1, n[1], pt[1]) * pbinom(r2-t1, n[2], pt[1]))) + out_r1
  out_r2 <- unique(out_r2)
  ind_r2 <- which(out_r2 <= as_left[2])
  # check if conditions hold
  if(length(ind_r2) == 0) stop("No optimal design (left side)")
  r2 <- r2_bdry[max(ind_r2)]
  out_r2 <- out_r2[max(ind_r2)]
  # boundary of s2
  s2_bdry <- r2:nc[2]
  out_s2 <- sapply(s2_bdry, function(s2) sum(dbinom(t1, n[1], pt[2]) * (1-pbinom(s2-t1, n[2], pt[2])))) + out_s1
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
    return(sapply(r3_bdry, function(r3) dbinom(tt, n[1], pt[1]) * sum(dbinom(t2, n[2], pt[1]) * pbinom(r3-tt-t2, n[3], pt[1]))))
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
    return(sapply(s3_bdry, function(s3) dbinom(tt, n[1], pt[2]) * sum(dbinom(t2, n[2], pt[2]) * (1-pbinom(s3-tt-t2, n[3], pt[2])))))
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
  pe <- pt[2] + pe.par
  emp_power <- 1 - pbinom(s1, n[1], pe) + sum(dbinom(t1, n[1], pe) * (1-pbinom(s2-t1, n[2], pe))) + sapply(t1, function(tt){
    t2 <- (r2-tt+1):(s2-tt)
    return(dbinom(tt, n[1], pe) * sum(dbinom(t2, n[2], pe) * (1-pbinom(s3-tt-t2, n[3], pe))))
  })
  # merge results
  names(bdry) <- c("r1", "r2", "r3", "s1", "s2", "s3")
  names(err) <- c("alpha11", "alpha12", "alpha13", "alpha21", "alpha22", "alpha23")
  out <- list(bdry = bdry, error = err, pt = pt, n = n, alpha = c(alpha1, alpha2), beta = 1 - emp_power, sf.param = sf.param)
  class(out) <- "2opt"
  return(out)
}

#' @keywords internal
two.opt <- function(alpha1, alpha2, pt, n, sf.param, pe.par, ...){
  # initialization
  nc <- cumsum(n)
  nt <- nc[2]
  as_left <- HSD(alpha1, nc/nt, sf.param)
  as_right <- HSD(alpha2, nc/nt, sf.param)
  comb <- NULL
  err <- NULL
  # boundary of r1: [0, n1]
  r1_bdry <- 0:nc[1]
  out_r1 <- unique(pbinom(r1_bdry, n[1], pt[1]))
  ind_r1 <- which(out_r1 <= as_left[1])
  # check if conditions hold
  if(length(ind_r1)==0) stop("No optimal design (left side)")
  r1 <- r1_bdry[max(ind_r1)]
  out_r1 <- out_r1[max(ind_r1)]
  # boundary of s1 (r1, n1-1]
  s1_bdry <- r1:(n[1]-1)
  out_s1 <- unique(1 - pbinom(s1_bdry, n[1], pt[2]))
  ind_s1 <- which(out_s1 <= as_right[1])
  # check if condition holds 
  if(length(ind_s1) == 0) stop("No optimal design (right side)")
  s1 <- s1_bdry[min(ind_s1)]
  out_s1 <- out_s1[min(ind_s1)]
  # boundary of r1: [r1, n1+n2-1]
  r2_bdry <- r1:(nc[2]-1)
  t1 <- (r1+1):s1
  out_r2 <- sapply(r2_bdry, function(r2) sum(dbinom(t1, n[1], pt[1]) * pbinom(r2-t1, n[2], pt[1]))) + out_r1
  out_r2 <- unique(out_r2)
  ind_r2 <- which(out_r2 <= as_left[2])
  # check if conditions hold
  if(length(ind_r2) == 0) stop("No optimal design (left side)")
  r2 <- r2_bdry[max(ind_r2)]
  out_r2 <- out_r2[max(ind_r2)]
  # boundary of s2
  s2_bdry <- r2:nc[2]
  out_s2 <- sapply(s2_bdry, function(s2) sum(dbinom(t1, n[1], pt[2]) * (1-pbinom(s2-t1, n[2], pt[2])))) + out_s1
  out_s2 <- unique(out_s2)
  ind_s2 <- which(out_s2 <= as_right[2])
  # check if conditions hold
  if(length(ind_s2) == 0) stop("No optimal design (right side)")
  s2 <- s2_bdry[min(ind_s2)]
  out_s2 <- out_s2[min(ind_s2)]
  bdry <- c(r1, r2, s1, s2)
  # calculate type-2 error with pt = pt + 0.2
  pe <- pt[2] + pe.par
  err <- c(out_r1, out_r2, out_s1, out_s2)
  emp_power <- 1 - pbinom(s1, n[1], pe) + sum(dbinom(t1, n[1], pe) * (1-pbinom(s2-t1, n[2], pe)))
  # merge results
  names(bdry) <- c("r1", "r2", "s1", "s2")
  names(err) <- c("alpha11", "alpha12", "alpha21", "alpha22")
  out <- list(bdry= bdry, error = err, pt = pt, n = n, alpha = c(alpha1, alpha2), beta = 1 - emp_power, sf.param = sf.param)
  class(out) <- "2opt"
  return(out)
}

#' @keywords internal
right.two.opt <- function(alpha, pt, n, sf.param, ...){
  # initialization
  nc <- cumsum(n)
  nt <- nc[2]
  as_right <- HSD(alpha, nc/nt, sf.param)
  # boundary of s1 (0, n1-1]
  s1_bdry <- 0:(n[1]-1)
  out_s1 <- 1-pbinom(s1_bdry, n[1], pt)
  ind_s1 <- which(out_s1 <= as_right[1])
  # check if condition holds 
  if(length(ind_s1)==0) stop("No optimal design (right side)")
  s1 <- s1_bdry[min(ind_s1)]
  out_s1 <- out_s1[min(ind_s1)]
  t1 <- 0:s1
  # boundary of s2
  s2_bdry <- s1:(nc[2]-1)
  out_s2 <- sapply(s2_bdry, function(s2) sum(dbinom(t1, n[1], pt) * (1 - pbinom(s2-t1, n[2], pt)))) + out_s1
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
  out <- list(bdry = bdry, error = err, pt = pt, n = n, sf.param = sf.param, alpha = alpha)
  class(out) <- "1opt"
  return(out)
}

#' @keywords internal
right.three.opt <- function(alpha, pt, n, sf.param, ...){
  # initialization
  nc <- cumsum(n)
  nt <- nc[3]
  as_right <- HSD(alpha, nc/nt, sf.param)
  # boundary of s1 (0, n1-1]
  s1_bdry <- 0:(n[1]-1)
  out_s1 <- 1-pbinom(s1_bdry, n[1], pt)
  ind_s1 <- which(out_s1 <= as_right[1])
  # check if condition holds 
  if(length(ind_s1)==0) stop("No optimal design (right side)")
  s1 <- s1_bdry[min(ind_s1)]
  out_s1 <- out_s1[min(ind_s1)]
  # loop over all s1 
  t1 <- 0:s1
  # boundary of s2
  s2_bdry <- s1:(nc[2]-1)
  out_s2 <- sapply(s2_bdry, function(s2) sum(dbinom(t1, n[1], pt) * (1 - pbinom(s2-t1, n[2], pt)))) + out_s1
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
    return(sapply(s3_bdry, function(s3) dbinom(tt, n[1], pt)  *sum(dbinom(t2, n[2], pt)*(1 - pbinom(s3-tt-t2, n[3], pt)))))
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
  out <- list(bdry = bdry, error = err, pt = pt, n = n, sf.param = sf.param, alpha = alpha)
  class(out) <- "1opt"
  return(out)
}

#' @keywords internal
three.cal <- function(alpha1, alpha2, beta, pc, pt, n, sf.param = 1, ...){
  # initialization
  nc <- cumsum(n)
  nt <- nc[3]
  as_left <- HSD(alpha1, nc/nt, sf.param)
  comb <- NULL
  err <- NULL
  # boundary of r1: [0, n1]
  r1_bdry <- 0 : nc[1]
  out_r1 <- pbinom(r1_bdry, n[1], pc[1])
  ind_r1 <- which(out_r1 <= as_left[1])
  # check if conditions hold
  if(length(ind_r1)==0) return(NULL)
  r1 <- r1_bdry[max(ind_r1)]
  out_r1 <- out_r1[max(ind_r1)]
  # boundary of r1: [r1, n1+n2-1]
  r2_bdry <- r1 : (nc[2]-1)
  t1 <- (r1+1) : n[1]
  out_r2 <- sapply(r2_bdry, function(r2) sum(dbinom(t1, n[1], pc[1]) * pbinom(r2-t1, n[2], pc[1]))) + out_r1
  ind_r2 <- which(out_r2 <= as_left[2])
  # check if conditions hold
  if(length(ind_r2) == 0) return(NULL)
  r2 <- r2_bdry[max(ind_r2)]
  out_r2 <- out_r2[max(ind_r2)]
  # boundary of r3: [r2, n1+n2+n3=n]
  r3_bdry <- r2 : nc[3]
  out_r3 <- sapply(t1, function(tt){
    t2 <- (r2-tt+1) : n[2]
    return(sapply(r3_bdry, function(r3) dbinom(tt, n[1], pc[1]) * sum(dbinom(t2, n[2], pc[1]) * pbinom(r3-tt-t2, n[3], pc[1]))))
  })
  # check if nrow(out_r3)==1
  if(is.vector(out_r3)) out_r3 <- sum(out_r3) + out_r2
  else out_r3 <- rowSums(out_r3) + out_r2
  ind_r3 <- which(out_r3 <= as_left[3])
  # check if conditions hold
  if(length(ind_r3) == 0) return(NULL)
  r3 <- r3_bdry[max(ind_r3)]
  out_r3 <- out_r3[max(ind_r3)]
  # boundary of s
  s_bdry <- (r3+1) : nc[3]
  out_s <- sapply(t1, function(tt){
    t2 <- (r2-tt+1) : n[2]
    return(sapply(s_bdry, function(s) dbinom(tt, n[1], pc[2]) * sum(dbinom(t2, n[2], pc[2]) * (1-pbinom(s-tt-t2, n[3], pc[2])))))
  })
  # check if nrow(out_s3)==1
  if(is.vector(out_s)) out_s <- sum(out_s)
  else out_s <- rowSums(out_s)
  ind_s <-  which(out_s <= alpha2)
  # check if conditions hold
  if(length(ind_s) == 0) return(NULL)
  s <- s_bdry[min(ind_s)]
  out_s <- out_s[min(ind_s)]
  # calculate type 2 error 
  beta_t <- sapply(t1, function(tt){
    t2 <- (r2-tt+1) : n[2]
    return(dbinom(tt, n[1], pt) * sum(dbinom(t2, n[2], pt) * pbinom(s-tt-t2, n[3], pt)))
  })					
  beta_t <- sum(beta_t)
  ind_beta <- which(beta_t <= beta)
  if(length(ind_beta) == 0) return(NULL)
  # save feasible designs & errors
  comb <- c(r1, r2, r3, s)
  err <- c(out_r1, out_r2, out_r3, out_s, beta_t)
  # merge results
  rs <- c("r1", "r2", "r3", "s")
  err_t <- c("alpha11", "alpha12", "alpha13", "alpha2", "beta")
  out <- c(comb, err)
  names(out) <- c(rs, err_t)
  return(out)
}

#' @keywords internal
zhong.three <- function(alpha1, alpha2, beta, pc, pt, frac_n1 = c(0.28, 0.32), frac_n2 = c(0.28,0.32), sf.param = 1, show = TRUE, nmax = 100, ...) {
  if(length(pc) == 1) {
    pc <- rep(pc, 2)
  }
  nt <- 3
  out <- NULL
  comb.n <- NULL
  n_count <- 0
  while(n_count <= nt & nt <= nmax) {
    n1_bdry <- floor(nt*frac_n1[1]):ceiling(nt*frac_n1[2])
    n2_bdry <- floor(nt*frac_n2[1]):ceiling(nt*frac_n2[2])
    
    # remove n1=0, n2=0, from boundary 
    n1_bdry <- n1_bdry[n1_bdry!=0]
    n2_bdry <- n2_bdry[n2_bdry!=0]
    
    # find all possible outcomes of n1 n2 n3
    n1_n2 <- expand.grid(n1_bdry, n2_bdry)
    
    # remove n3=0
    sum_n1_n2 <- rowSums(n1_n2)
    id_rm <- sum_n1_n2 < nt 
    n3_bdry <- nt-sum_n1_n2[id_rm]
    n <- as.matrix(cbind(n1_n2[id_rm,], n3_bdry))
    colnames(n) <- c("n1", "n2", "n3")
    for(i in 1:nrow(n)){
      ans <- three.cal(alpha1, alpha2, beta, pc, pt, n[i, ], sf.param)
      if(!is.null(ans)) comb.n <- rbind(comb.n, n[i, ]) 
      out <- rbind(out, ans)
    }
    if(!is.null(out)) n_count <- 2 * nrow(out)
    nt <- nt+1
    if(show) print(paste("current sample size is", nt))
  }
  out <- cbind(out, comb.n)
  out <- out[order(-out[,5], -out[,6], -out[,7], -out[,8], out[,9]), ]
  bdry <- out[1, ][1:4]
  error <- out[1, ][5:9]
  nn <- out[1, ][10:12]
  names(bdry) <- c("r1", "r2", "r3", "s")
  names(nn) <- c("n1", "n2", "n3")
  names(error) <- c("alpha11", "alpha12", "alpha13", "alpha2", "beta")
  return(list(bdry = bdry, error = error, n = nn, complete = out))
}

#' @keywords internal
zhong.two <- function(alpha1, alpha2, beta, pc, pt, two.sided = TRUE, sf.param = NULL, show = TRUE, nmax = 100, n.ratio=1/2, ...){
  if(length(pc) == 1) {
    pc <- rep(pc, 2)
  }
  if(length(alpha1) == 1) {
    alpha1_ini <- c(1, alpha1)
    alpha2_ini <- c(1, alpha2)
  } else {
    stop("alpha1 and alpha2 should be a single value")
  }
  
  for(n in 2 : nmax) {
    comb <- NULL
    err <- NULL
    n_count <- 0
    for(n1 in 1:(n-1)) {
      n2 <- n - n1
      if(!is.null(sf.param)) {
        alpha1 <- HSD(alpha1_ini[2], c(n1, n)/n, sf.param)
        alpha2 <- HSD(alpha2_ini[2], c(n1, n)/n, sf.param)
      } else {
        alpha1 <- alpha1_ini
        alpha2 <- alpha2_ini
      }
      for(r1 in 0:(n1-1)) {
        if(two.sided) s1_lb <- r1 + 1
        else  s1_lb <- n1
        # check L1 
        L1 <- pbinom(r1, n1, pc[1])
        if(L1 <=  alpha1[1] & L1 <= alpha1[2]){
          for(s1 in s1_lb: n1){
            R1 <- 1 - pbinom(s1, n1, pc[2])
            if(R1 <= alpha2[1] & R1 <= alpha2[2]) {
              t1 <- (r1+1) : s1
              b <- dbinom(t1, n1, pc[1])
              for(r2 in r1 : (s1 + n2)){
                L2 <- L1 + sum(b * pbinom(r2 - t1, n2, pc[1]))
                if(L2 <= alpha1[2]) {
                  for(s2 in max(r2, s1) : (s1+n2)) {
                    R2 <- R1 + sum(b * (1 - pbinom(s2 - t1, n2, pc[2])))
                    # the beta is incorrect, need to be fixed
                    if(R2 <= alpha2[2]) {
                      Rpe <- sum(dbinom(t1, n1, pt) * (1 - pbinom(s2 - t1, n2, pt))) + 1 - pbinom(s1, n1, pt)
                    }
                    else next
                    if(Rpe >= 1 - beta){
                      comb <- rbind(c(r1, r2, s1, s2, n1, n2), comb)
                      err <- rbind(c(L1, L2, R1, R2, 1- Rpe), err)
                    } else next
                  }
                } else next	
              }
            } else next
          }
        } else next
      }
    }
    if(!is.null(comb)) n_count <- nrow(comb)
    if(n_count > n.ratio * n) break
    if(show) print(paste("current sample size is", n))
  }
  out <- cbind(err, comb)
  out <- out[order(-out[, 2], -out[, 4] , out[,5]),]
  out <- cbind(err, comb)
  colnames(out) <- c("alpha11", "alpha12", "alpha13", "alpha21", "alpha22", "alpha23", "beta", "r1", "r2", "r3", "s1", "s2", "s3", "n1", "n2", "n3")
  opt <- out[1, ]
  return(list(bdry = opt[5:7], error = opt[1:4], n = opt[8:11], complete = out))
}

#' @keywords internal
HSD <- function (alpha, t, param) {
  t[t > 1] <- 1
  spend <- if (param == 0) t * alpha else alpha * (1 - exp(-t * param))/(1 - exp(-param))
  return(spend)
}

