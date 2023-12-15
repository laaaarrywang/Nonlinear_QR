estFnL <- function (par, x, y, gridmats, L, mid, nknots, ngrid, a.kap, 
                    a.sig, tau.g, reg.ix, reduce = TRUE, x.ce = 0, x.sc = 1, 
                    base.bundle,sparse = TRUE) {
  n <- length(y) # number of observations
  p <- ncol(x) # number of explanatory variables
  wKnot <- matrix(par[1:(nknots * (p + 1))], nrow = nknots) # initial value of w_j at m=6 supplementary knots.
  w0PP <- qrjoint:::ppFn0(wKnot[, 1], gridmats, L, nknots, ngrid) # wKnot[,1] is the initial value for w_0 at different knots
  w0 <- w0PP$w # \tilde{w_0}(\tau), L-dimensional
  #wPP <- apply(wKnot[, -1, drop = FALSE], 2, qrjoint:::ppFn, gridmats = gridmats, 
  #             L = L, nknots = nknots, ngrid = ngrid, a.kap = a.kap) # \tilde{w_j}(\tau), L-dimensional
  #wMat <- matrix(sapply(wPP, qrjoint:::extract, vn = "w"), ncol = p) # 121*13
  if (sparse){
    delta <- matrix(c(0,par[(nknots + 1)*(p+1) + 2 + 1:p]),ncol = p + 1)
    wKnot.delta = rbind(wKnot,delta)
    wPP <- apply(wKnot.delta[, -1, drop = FALSE], 2, ppFnL, gridmats = gridmats, 
                 L = L, nknots = nknots, ngrid = ngrid, a.kap = a.kap) # \tilde{w_j}(\tau), L-dimensional
    wMat <- matrix(sapply(wPP, qrjoint:::extract, vn = "w"), ncol = p) # 121*13
    
    indi <- sweep(abs(wMat), 2, delta[-1], FUN=">")
    wMat <- wMat*indi
  } else{
    wPP <- apply(wKnot[, -1, drop = FALSE], 2, qrjoint:::ppFn, gridmats = gridmats, 
                 L = L, nknots = nknots, ngrid = ngrid, a.kap = a.kap) # \tilde{w_j}(\tau), L-dimensional
    wMat <- matrix(sapply(wPP, qrjoint:::extract, vn = "w"), ncol = p) # 121*13
  }
  zeta0.dot <- exp(qrjoint:::shrinkFn(p) * (w0 - max(w0))) # proportional to derivative of zeta. One question: why delete max(w0)?
  zeta0 <- qrjoint:::trape(zeta0.dot[-c(1, L)], tau.g[-c(1, L)], L - 
                             2) # throw away the first and last grid (0 and 1), 
  zeta0.tot <- zeta0[L - 2]
  zeta0 <- c(0, tau.g[2] + (tau.g[L - 1] - tau.g[2]) * zeta0/zeta0.tot, 
             1) # zeta0.tot is the total area between tau.g[2] and tau.g[L-1]
  zeta0.dot <- (tau.g[L - 1] - tau.g[2]) * zeta0.dot/zeta0.tot
  zeta0.dot[c(1, L)] <- 0
  zeta0.ticks <- pmin(L - 1, pmax(1, sapply(zeta0, function(u) sum(tau.g <= 
                                                                     u)))) # for interpolation
  zeta0.dists <- (zeta0 - tau.g[zeta0.ticks])/(tau.g[zeta0.ticks + 
                                                       1] - tau.g[zeta0.ticks]) # 分子是zeta0中的点与紧邻的上一个tau.g中点的距离，分母是所在线段的长度
  vMat <- apply(wMat, 2, qrjoint:::transform.grid, ticks = zeta0.ticks, 
                dists = zeta0.dists) # 计算w(zeta(tau))
  reach <- nknots * (p + 1)
  gam0 <- par[reach + 1]
  reach <- reach + 1
  gam <- par[reach + 1:p]
  reach <- reach + p
  sigma <- qrjoint:::sigFn(par[reach + 1], a.sig)
  reach <- reach + 1
  nu <- qrjoint:::nuFn(par[reach + 1])
  b0dot <- sigma * base.bundle$q0(zeta0, nu) * zeta0.dot
  beta0.hat <- rep(NA, L)
  beta0.hat[mid:L] <- gam0 + qrjoint:::trape(b0dot[mid:L], tau.g[mid:L], 
                                             L - mid + 1)
  beta0.hat[mid:1] <- gam0 + qrjoint:::trape(b0dot[mid:1], tau.g[mid:1], 
                                             mid)
  vNorm <- sqrt(rowSums(vMat^2))
  a <- tcrossprod(vMat, x)
  aX <- apply(-a, 1, max)/vNorm
  aX[is.nan(aX)] <- Inf
  aTilde <- vMat/(aX * sqrt(1 + vNorm^2))
  ab0 <- b0dot * aTilde
  beta.hat <- kronecker(rep(1, L), t(gam))
  beta.hat[mid:L, ] <- beta.hat[mid:L, ] + apply(ab0[mid:L, 
                                                     , drop = FALSE], 2, qrjoint:::trape, h = tau.g[mid:L], len = L - 
                                                   mid + 1) # why mid?
  beta.hat[mid:1, ] <- beta.hat[mid:1, ] + apply(ab0[mid:1, 
                                                     , drop = FALSE], 2, qrjoint:::trape, h = tau.g[mid:1], len = mid) # why mid?
  beta.hat <- beta.hat/x.sc
  beta0.hat <- beta0.hat - rowSums(beta.hat * x.ce)
  betas <- cbind(beta0.hat, beta.hat)
  if (reduce) {
    betas <- betas[reg.ix, , drop = FALSE]
  }
  return(betas)
}
