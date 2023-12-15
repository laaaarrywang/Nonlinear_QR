coefL <- function (object, burn.perc = 0.5, nmc = 200, plot = TRUE, show.intercept = TRUE, 
                  reduce = TRUE, sparse = TRUE, ...) 
{
  nsamp <- object$dim[10]
  pars <- matrix(object$parsamp, ncol = nsamp)
  ss <- unique(round(nsamp * seq(burn.perc, 1, len = nmc + 
                                   1)[-1]))
  n <- object$dim[1]
  p <- object$dim[2]
  L <- object$dim[3]
  mid <- object$dim[4] + 1
  nknots <- object$dim[5]
  ngrid <- object$dim[6]
  a.sig <- object$hyper[1:2]
  a.kap <- matrix(object$hyper[-c(1:2)], nrow = 3)
  tau.g <- object$tau.g
  reg.ix <- object$reg.ix
  x.ce <- outer(rep(1, L), attr(object$x, "scaled:center"))
  x.sc <- outer(rep(1, L), attr(object$x, "scaled:scale"))
  base.bundle <- list()
  if (object$fbase.choice == 1) {
    base.bundle$q0 <- function(u, nu = Inf) return(1/(dt(qt(qrjoint:::unitFn(u), 
                                                            df = nu), df = nu) * qt(0.9, df = nu)))
    base.bundle$Q0 <- function(u, nu = Inf) return(qt(qrjoint:::unitFn(u), 
                                                      df = nu)/qt(0.9, df = nu))
    base.bundle$F0 <- function(x, nu = Inf) return(pt(x * 
                                                        qt(0.9, df = nu), df = nu))
  } else if (object$fbase.choice == 2) {
    base.bundle$q0 <- function(u, nu = Inf) return(1/dlogis(qlogis(qrjoint:::unitFn(u))))
    base.bundle$Q0 <- function(u, nu = Inf) return(qlogis(qrjoint:::unitFn(u)))
    base.bundle$F0 <- function(x, nu = Inf) return(plogis(x))
  } else {
    base.bundle$q0 <- function(u, nu = Inf) return(1/(dunif(qunif(u, 
                                                                  -1, 1), -1, 1)))
    base.bundle$Q0 <- function(u, nu = Inf) return(qunif(u, 
                                                         -1, 1))
    base.bundle$F0 <- function(x, nu = Inf) return(punif(x, 
                                                         -1, 1))
  }
  beta.samp <- sapply(ss, function(p1) estFnL(pars[, p1], object$x, 
                                             object$y, object$gridmats, L, mid, nknots, ngrid, a.kap, 
                                             a.sig, tau.g, reg.ix, reduce, x.ce, x.sc, base.bundle,sparse), 
                      simplify = "array")
  if (reduce) 
    tau.g <- tau.g[reg.ix]
  L <- length(tau.g)
  if (plot) {
    nr <- ceiling(sqrt(p + show.intercept))
    nc <- ceiling((p + show.intercept)/nr)
    cur.par <- par(no.readonly = TRUE)
    par(mfrow = c(nr, nc))
  }
  beta.hat <- array(0, c(L, p + 1, 3))
  plot.titles <- c("Intercept", object$xnames)
  beta.hat[, 1, ] <- qrjoint:::getBands(beta.samp[, 1, ], plot = (plot & 
                                                                    show.intercept), add = FALSE, x = tau.g, xlab = "tau", 
                                        ylab = "Coefficient", bty = "n")
  if (plot & show.intercept) 
    title(main = plot.titles[1])
  for (j in 2:(p + 1)) {
    beta.hat[, j, ] <- qrjoint:::getBands(beta.samp[, j, ], plot = plot, 
                                          add = FALSE, x = tau.g, xlab = "tau", ylab = "Coefficient", 
                                          bty = "n")
    if (plot) {
      title(main = plot.titles[j])
      abline(h = 0, lty = 2, col = 4)
    }
  }
  if (plot) 
    suppressWarnings(par(cur.par, no.readonly = TRUE))
  dimnames(beta.hat) <- list(tau = tau.g, beta = plot.titles, 
                             summary = c("b.lo", "b.med", "b.hi"))
  dimnames(beta.samp) <- list(tau = tau.g, beta = plot.titles, 
                              iter = 1:length(ss))
  invisible(list(beta.samp = beta.samp, beta.est = beta.hat))
  mid.red <- which(tau.g == object$tau.g[mid])
  parametric.list <- rbind(beta.samp[mid.red, , , drop = TRUE], 
                           sigma = qrjoint:::sigFn(pars[nknots * (p + 1) + (p + 1) + 1, ss], 
                                                   a.sig), nu = qrjoint:::nuFn(pars[(nknots + 1) * (p + 1) + 2, 
                                                                                    ss]))
  dimnames(parametric.list)[[1]][1 + 0:p] <- c("Intercept", 
                                               object$xnames)
  gamsignu <- t(apply(parametric.list, 1, quantile, pr = c(0.5, 
                                                           0.025, 0.975)))
  dimnames(gamsignu)[[2]] <- c("Estimate", "Lo95%", "Up95%")
  invisible(list(beta.samp = beta.samp, beta.est = beta.hat, 
                 parametric = gamsignu))
}
