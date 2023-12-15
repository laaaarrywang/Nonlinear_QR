ppFnL <- function (w.knot.delta, gridmats, L, nknots, ngrid, a.kap) 
{
  w.grid <- matrix(NA, L, ngrid)
  lpost.grid <- rep(NA, ngrid)
  w.knot <- w.knot.delta[-length(w.knot.delta)]
  delta <- w.knot.delta[length(w.knot.delta)]
  for (i in 1:ngrid) {
    A <- matrix(gridmats[1:(L * nknots), i], nrow = nknots)
    R <- matrix(gridmats[L * nknots + 1:(nknots * nknots), 
                         i], nrow = nknots)
    
    r <- qrjoint:::sum.sq(backsolve(R, w.knot, transpose = TRUE))
    w.grid[, i] <- colSums(A * w.knot)
    lpost.grid[i] <- (qrjoint:::logsum(-(nknots/2 + a.kap[1, ] + 1) * log1p(0.5 * r/a.kap[2, ]) 
                                       + a.kap[3, ] 
                                       + lgamma(a.kap[1, ] + nknots/2 + 1) 
                                       - lgamma(a.kap[1, ]) 
                                       - 0.5 * nknots * log(a.kap[2,]) 
                                       + log(pgamma(3/delta,nknots/2 + a.kap[1, ] + 1,scale = 1/(0.5 * r/a.kap[2, ]))) # sparsity is considered) 
                      - gridmats[nknots * (L + nknots) + 1, i] 
                      + gridmats[nknots * (L + nknots) + 2, i]))
  }
  lpost.sum <- qrjoint:::logsum(lpost.grid)
  post.grid <- exp(lpost.grid - lpost.sum)
  w <- c(w.grid %*% post.grid)
  return(list(w = w, lpost.sum = lpost.sum))
}