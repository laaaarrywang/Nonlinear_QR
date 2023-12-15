summaryL = function (object, ntrace = 1000, burn.perc = 0.5, plot.dev = TRUE, 
                     more.details = FALSE, ...) 
{
  thin <- object$dim[9]
  nsamp <- object$dim[10]
  pars <- matrix(object$parsamp, ncol = nsamp)
  ss <- unique(pmax(1, round(nsamp * (1:ntrace/ntrace))))
  post.burn <- (ss > nsamp * burn.perc)
  dimpars <- object$dim
  dimpars[8] <- length(ss)
  n <- object$dim[1]
  p <- object$dim[2]
  ngrid <- object$dim[6]
  sm <- .C("DEVL", pars = as.double(pars[, ss]), x = as.double(object$x), 
           y = as.double(object$y), cens = as.integer(object$cens), 
           wt = as.double(object$wt), shrink = as.integer(object$shrink), 
           hyper = as.double(object$hyper), dim = as.integer(dimpars), 
           gridmats = as.double(object$gridmats), tau.g = as.double(object$tau.g), 
           devsamp = double(length(ss)), llsamp = double(length(ss) * 
                                                           n), pgsamp = double(length(ss) * ngrid * (p + 1)), 
           qlsamp = double(length(ss) * n), fbase.choice = as.integer(object$fbase.choice))
  deviance <- sm$devsamp
  ll <- matrix(sm$llsamp, ncol = length(ss))
  ql <- matrix(sm$qlsamp, ncol = length(ss))
  fit.waic <- qrjoint::waic(ll[, post.burn])
  pg <- matrix(sm$pgsamp, ncol = length(ss))
  prox.samp <- matrix(NA, p + 1, length(ss))
  for (i in 1:(p + 1)) {
    prox.samp[i, ] <- object$prox[apply(pg[(i - 1) * ngrid + 
                                             1:ngrid, ], 2, function(pr) sample(length(pr), 1, 
                                                                                prob = pr))]
  }
  if (more.details) {
    cur.par <- par(no.readonly = TRUE)
    par(mfrow = c(2, 2), mar = c(5, 4, 3, 2) + 0.1)
  }
  if (plot.dev) {
    plot(thin * ss, deviance, ty = "l", xlab = "Markov chain iteration", 
         ylab = "Deviance", bty = "n", main = "Fit trace plot", 
         )
    grid(col = "gray")
  }
  if (more.details) {
    ngrid <- length(object$prox)
    prior.grid <- exp(object$gridmats[nrow(object$gridmats), 
    ])
    lam.priorsamp <- qrjoint:::lamFn(sample(object$prox, ntrace, replace = TRUE, 
                                            prob = prior.grid))
    lam.prior.q <- quantile(lam.priorsamp, pr = c(0.025, 
                                                  0.5, 0.975))
    lam.samp <- qrjoint:::lamFn(prox.samp)
    a <- min(qrjoint:::lamFn(object$prox))
    b <- diff(range(qrjoint:::lamFn(object$prox))) * 1.2
    plot(thin * ss, lam.samp[1, ], ty = "n", ylim = a + c(0, 
                                                          b * (p + 1)), bty = "n", ann = FALSE, axes = FALSE)
    axis(1)
    for (i in 1:(p + 1)) {
      abline(h = b * (i - 1) + qrjoint:::lamFn(object$prox), col = "gray")
      abline(h = b * (i - 1) + lam.prior.q, col = "red", 
             lty = c(2, 1, 2))
      lines(thin * ss, b * (i - 1) + lam.samp[i, ], lwd = 1, 
            col = 4)
      if (i%%2) 
        axis(2, at = b * (i - 1) + qrjoint:::lamFn(object$prox[c(1, 
                                                                 ngrid)]), labels = round(object$prox[c(1, ngrid)], 
                                                                                          2), las = 1, cex.axis = 0.6)
      mtext(substitute(beta[index], list(index = i - 1)), 
            side = 4, line = 0.5, at = a + b * (i - 1) + 
              0.4 * b, las = 1)
    }
    title(xlab = "Markov chain iteration", ylab = "Proxmity posterior", 
          main = "Mixing over GP scaling")
    theta <- coda::as.mcmc(t(matrix(object$parsamp, ncol = nsamp)[, 
                                                            ss[post.burn]]))
    gg <- coda::geweke.diag(theta, 0.1, 0.5)
    zvals <- gg$z
    pp <- 2 * (1 - pnorm(abs(zvals)))
    plot(sort(pp), ylab = "Geweke p-values", xlab = "Parameter index (reordered)", 
         main = "Convergence diagnosis", ty = "h", col = 4, 
         ylim = c(0, 0.3), lwd = 2)
    abline(h = 0.05, col = 2, lty = 2)
    abline(a = 0, b = 0.1/length(pp), col = 2, lty = 2)
    mtext(c("BH-10%", "5%"), side = 4, at = c(0.1, 0.05), 
          line = 0.1, las = 0, cex = 0.6)
    npar <- length(object$par)
    suppressWarnings(image(1:npar, 1:npar, cor(theta), xlab = "Parameter index", 
                           ylab = "Parameter index", main = "Parameter correlation"))
    suppressWarnings(par(cur.par, no.readonly = TRUE))
  }
  invisible(list(deviance = deviance, pg = pg, prox = prox.samp, 
                 ll = ll, ql = ql, waic = fit.waic))
}
