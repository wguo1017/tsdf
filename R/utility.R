#' plot decision table from a "dec.table" object.
#' @description \code{plot} method for class "\code{dec.table}"
#' @param x an object of class \code{"dec.table"}, a result of a call to \code{dec.table}.
#' @param ... Not used argument.
#' @details \code{plot.dec.table} prints the decision boundarys.
#' @import graphics
#' @export
#' @examples
#' truep <- c(0.3, 0.45, 0.5, 0.6)
#' out <- dec.table(0.6,0.4,0.2,0.3,c(3,3,3))
#' plot(out)
plot.dec.table <- function(x, ...) {
  n <- x$n
  nc <- cumsum(n)
  r <- x$E
  s <- x$D
  su <- x$DU
  col <- c("cadetblue1", "darkorchid2", "indianred2", "seagreen2")
  plot(nc, r, xaxt = "n", ylab = "Boundary", xlab = 'Sample Size', ylim = c(0, max(su)+2.3), type = "o", col=col[1], main = "Decision Plot", panel.first = grid())
  points(nc, s + 1, type = "o", col = col[2])
  points(nc, su + 1, type = "o", col = col[3])
  axis(1, at = nc, labels = nc)
  polygon(c(nc, rev(nc)), c(r + 0.05, rev(s + 1 - 0.05)), col = col[4], border = NA)
  polygon(c(nc, rev(nc)), c(rep(0,length(nc)), rev(r - 0.05)), col = col[1], border = NA)
  polygon(c(nc, rev(nc)), c(s + 1 + 0.05, rev(su + 1 - 0.05)), col = col[2], border = NA)
  polygon(c(nc, rev(nc)), c(su + 1 + 0.05, rev(rep(max(su)+2, length(nc)))), col = col[3], border = NA)
  legend("top", c("E (<=)", "D (>=)", "DU (>=)", "S"), horiz = TRUE, fill = col)
}


#' print decision table from a "dec.table" object.
#' @description \code{print} method for class "\code{dec.table}"
#' @param x an object of class \code{"dec.table"}, a result of a call to \code{dec.table}.
#' @param ... Not used argument.
#' @details \code{print.dec.table} prints the decision table with legend keys.
#' @export
#' @examples
#' print(dec.table(0.6,0.4,0.2,0.3,c(3,3,3)))
print.dec.table <- function(x, ...) {
  print(x$table)
  cat("\n")
  cat("	  Row : Number of DLTs", "\n")
  cat("Column : Number of subjects", "\n")
  cat("     E : Escalate to the next higher dose", "\n")
  cat("     S : Stay at the same dose", "\n")
  cat("     D : De-escalate to the previous lower dose", "\n")
  cat("    DU : De-escalate and never use this dose again", "\n")
  cat("\n")
  cat("Type 1 error : ", "\n")
  print(x$alpha.two)
  cat("Type 2 error : ", x$beta,"\n")
  cat("Right-side type 1 error : ", "\n")
  print(x$alpha.one)
}

#' plot simulation results from a dec.sim object
#' @description Three plots are currently available: a plot of true toxicity at each dose level (\code{type = "s"}); a bar plot of the probability of selecting as the MTD for each dose level (\code{type = "prob"}); a bar plot of the average number of patients treated at each dose level (\code{type = "np"}); a bar plot of the average number of patients experienced DLT at each dose level (\code{type = "dlt"}) and \code{type = "all"} generates all above plots.
#' @param x an object of class \code{"dec.sim"} or \code{"sl.sim"}, a result of a call to \code{dec.sim} or \code{sl.sim}.
#' @param pt a vector with target toxicity for each scenario.
#' @param s scenario to be plotted. Defaults to 1.
#' @param type plot type. See descriptions above.
#' @param label a logical value indicating if values are shown on plot.
#' @param col graphical parameter \code{col}; see details \code{\link{par}}.
#' @param text.col plotting color of text shown.
#' @param cex graphical parameter \code{col}; see details \code{\link{par}}.
#' @param ... arguments to be passed to \code{\link{plot}} methods.
#' @import graphics
#' @export
#' @examples 
#' # generate decision table
#' dt <- dec.table(0.6,0.4,0.2,0.3,c(3,3,3))
#' # simulate trials from test data 
#' test.file <- system.file("extdata", "testS.csv", package = "tsdf")
#' out <- sl.sim(dt$table, test.file)
#' plot(out, pt=rep(0.3,2), s=1, type="all")
#' plot(out, pt=rep(0.3,2), s=2, type="prob")
#' plot(out, pt=rep(0.3,2), s=1, type="np")
#' plot(out, pt=rep(0.3,2), s=2, type="dlt")
plot.dec.sim <- function(x, pt, s = 1, type = c("all", "s", "prob", "np", "dlt"), label = TRUE, col = "cornflowerblue", text.col = "darkblue", cex = 1, ...) {
  type <- match.arg(type)
  if (class(x)[1] == "sl.sim") {
    obj <- x[[s]]
    ns <- length(x)
  } else {
    obj <- x
    ns <- 1
  }
  truep <- obj$truep
  n_dose <- length(truep)
  names.arg <- 1:n_dose
  scen <- paste("( Scenario", s, ")")
  if (missing(pt)) {
    print("No target toxicity/mtd shown")
  }
  else if (length(pt) != ns) {
    stop("true toxicity is a vector with length = # of cenarios")
  } else {
    if(max(truep) < pt[s]) {
      print("MTD is higher than highest dose")
    } else if(min(truep) > pt[s]) {
      print("MTD is lower than lowest dose")
    } 
    else {
      mtd <- max(which(truep <= pt[s]))
      names.arg[mtd] <- paste(mtd, "(MTD)")
    }
  }
  if (type == "all") {
    par(mfrow = c(2, 2), cex = cex)
    plot(x, pt, s, "s", label, col, text.col, cex, ...)
    plot(x, pt, s, "prob", label, col, text.col, cex, ...)
    plot(x, pt, s, "dlt", label, col, text.col, cex, ...)
    plot(x, pt, s, "np", label, col, text.col, cex, ...)		
  }
  if (type == "prob") {
    ans <- obj$mtd.prob[1:n_dose]
    names(ans) <- NULL
    out <- barplot(rep(NA, n_dose), xlab = "Dose level", names.arg = names.arg, panel.first = box(), ylim = c(0, max(ans) * 1.1))
    grid()
    barplot(ans, add = TRUE, col = col)
    title(paste("Probability of selection", scen))
    if (label) {
      text(out, ans, labels = ans, pos = 3, cex = cex, col = text.col)
    }
  }
  if (type == "dlt") {
    ans <- obj$dlt
    out <- barplot(rep(NA, n_dose), xlab = "Dose level", names.arg = names.arg, panel.first = box(), ylim = c(0, max(ans) * 1.1))
    grid()
    barplot(ans, add = TRUE, col = col)
    title(paste("Average number of DLTs", scen))
    if (label) {
      text(out, ans, labels = ans, pos = 3, cex = cex, col = text.col)
    }
  }
  if (type == "np") {
    ans <- obj$n.patients
    out <- barplot(rep(NA, n_dose), xlab = "Dose level", names.arg = names.arg, panel.first = box(), ylim = c(0, max(ans) * 1.1))
    grid()
    barplot(ans, add = TRUE, col = col)
    title(paste("Average number of patients", scen))
    if (label) {
      text(out, ans, labels = ans, pos = 3, cex = cex, col = text.col)
    }
  }
  if (type == "s") {
    out <- plot(1:n_dose, truep, xlab = "Dose level", ylab = "Toxicity", col = col, type = "b", ylim = c(min(truep) * 0.9, max(truep) * 1.1), xaxt = "n", pch = 19, panel.first = c(box(), grid()))
    axis(1, at = 1:n_dose, labels = names.arg)
    title(paste("True toxicity", scen))
    if (!missing(pt)) {
      abline(pt[s], 0, lty = 5)
      text(n_dose/2, pt[s], labels = paste("Target =", pt[s]), pos = 3, cex = cex, col = text.col)
    }
    if (label) {
      text(1:n_dose, truep, labels = truep, pos = 3, cex = cex, col = text.col)
    }
  }
}


#' Summarizing simulation results from a dec.sim object
#' @description \code{summary} method for class \code{"dec.sim"}.
#' @details \code{summary} is used for formating important statistics for dose-finding simulation. Giving the target toxicity, it returns the probability of selecting current dose level as the MTD and over the MTD, probability of selecting the true MTD, probability of subjects treated at or below the true MTD, etc. The MTD is defined as the highest dose level such that the toxicity probability is less than target toxicity probability, if target is less than the smallest probability, then the lowest dose level is set as MTD. For example, if target is 0.3 and true toxicity for five doses are 0.1, 0.25, 0.35, 0.40, then MTD is dose 2.
#' @param object an object of class \code{"dec.sim"}, a result of a call to \code{dec.sim} or \code{sl.sim}.
#' @param pt target toxicity for each scenario.
#' @param ... Not used argument.
#' @examples 
#' test.file <- system.file("extdata", "testS.csv", package = "tsdf")
#' dt <- dec.table(0.6,0.4,0.2,0.3,c(3,3,3))
#' out <- sl.sim(dt$table, test.file)
#' pt <- c(0.3, 0.4)
#' summary(out, pt)
#' @import stats
#' @export
summary.dec.sim <- function(object, pt, ...) {
  if(missing(pt)) {
    warning("Missing true toxicity")
  }
  if(class(object)[1] == "sl.sim" & length(object) != length(pt)) {
    warning("pt length not equal to number of scenarios; only returns the first scenario stats")
  }
  ns <- length(pt)
  avg.np <- rep(0, ns)
  avg.dose <- rep(0, ns)
  avg.prob <- rep(0, ns)
  out <- vector("list", ns)
  for(i in 1:ns) {
    if(class(object)[1] == "sl.sim"){
      ans <- object[[i]]
    } else {
      ans <- object
    }
    truep <- ans$truep
    n_dose <- length(truep)
    res <- matrix(0, n_dose, 6)
    colnames(res) <- c("Level", "Truth", "MTD", "Over","DLT", "NP")
    res[, 1] <- 1:n_dose 
    res[, 2] <- truep
    res[, 3] <- ans$mtd.prob[1:n_dose]
    res[, 4] <- ans$over.prob[1:n_dose]
    res[, 5] <- ans$dlt
    res[, 6] <- ans$n.patients
    # calculate stats
    avg.np[i] <- sum(res[, 6])
    if(max(truep) < pt[i]) {
      mtd.dose <- "higher than highest dose"
      avg.dose[i] <- ans$mtd.prob[n_dose+2]
      avg.prob[i] <- 1
    } else if(min(truep) > pt[i]) {
      mtd.dose <- "lower than lowest dose"
      avg.dose[i] <- ans$mtd.prob[n_dose+1]
      avg.prob[i] <- 0
    } else {
      mtd.dose <- max(which(truep <= pt[i]))
      avg.dose[i] <- res[, 3][mtd.dose]
      avg.prob[i] <-sum(res[, 6][truep <= pt[i]])/sum(res[, 6])
    }
    out[[i]] <- list(dose.stats = res, prob.select = avg.dose[i], at.below.mtd = avg.prob[i], mtd = mtd.dose, nsim = ans$nsim, pt = pt[i], avg.np = avg.np[i], start.level = ans$start.level)
  }
  class(out) <- "summary.dec.sim"
  out
}

#' @export
print.summary.dec.sim <- function(x, ...) {
  cat("What does each column represent ?", "\n\n")
  cat("Level : Dose level", "\n")
  cat("Truth : True toxicity probability", "\n")
  cat("MTD   : The probability of selecting current dose level as the MTD", "\n")
  cat("Over  : The probability of selecting current dose level as over the MTD", "\n")
  cat("DLT   : The average number of subjects experienced DLT at current dose level", "\n")
  cat("NP    : The average number of subjects treated at current dose level", "\n\n")
  
  for(i in 1:length(x)){
    # table legends
    cat("Scenario ", i, "\n\n")
    cat("Simulation results are based on", x[[i]]$nsim, "simulated trials","\n")
    cat("Starting dose level is", x[[i]]$start.level, "; MTD is dose", x[[i]]$mtd, "\n")
    cat("Target toxicity probability =", x[[i]]$pt, "\n")
    cat("Average number of subjects =", x[[i]]$avg.np, "\n")
    cat("Probability of selecting the true MTD =", x[[i]]$prob.select, "\n")
    cat("Probability of subjects treated at or below the true MTD =", x[[i]]$at.below.mtd, "\n\n")
    print(as.data.frame(x[[i]]$dose.stats))
    cat("\n")
  }
}


#' print Zhong's design from a "opt.design" object.
#' @description \code{print} method for class "\code{opt.design}"
#' @param x an object of class \code{"opt.design"}, a result of a call to \code{opt.design}.
#' @param ... not used argument.
#' @export
#' @examples
#' alpha1 <- 0.20
#' alpha2 <- 0.1
#' beta <- 0.20
#' pc <- 0.5
#' pt <- pc + 0.2
#' out <- opt.design(alpha1, alpha2, beta, pc, pt, stage = 2, sf.param = 1)
#' print(out)
print.opt.design <- function(x, ...) {
  cat("\n Zhong's", x$stage,"stage Phase II design \n\n")
  cat("Minimal response rate: ", unique(x$pc), "\n")
  cat("Postulate response rate: ", x$pt, "\n")
  cat("Left-side type 1 error: ",x$alpha1, "\n")
  cat("Right-side type 1 error: ",x$alpha2, "\n")
  cat("Type 2 error: ",x$beta, "\n\n")
  cat("Notation:", "\n")
  cat("Left-side rejection region at stage i is response <= ri", "\n")
  cat("Right-side rejection region at stage", x$stage, "is response > s", "\n")
  cat("Sample size used at stage i is ni", "\n\n")
  cat("Optimal design : ", "\n")
  print(c(x$bdry, x$n))
  cat("True errors : ", "\n")
  print(x$error)
  cat("\n")
}
