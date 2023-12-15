predL <- function (object, newdata = NULL, summarize = TRUE, burn.perc = 0.5, 
                   nmc = 200, sparse = TRUE, ...) 
{
  p <- object$dim[2]
  betas <- coefL(object, burn.perc = burn.perc, nmc = nmc, plot = FALSE,sparse)
  nsamp <- dim(betas$beta.samp)[3]
  L <- dim(betas$beta.samp)[1]
  if (is.null(newdata)) {
    Xpred <- cbind(1, sapply(1:p, function(r) object$x[, 
                                                       r] * attr(object$x, "scaled:scale")[r] + attr(object$x, 
                                                                                                     "scaled:center")[r]))
  }
  else {
    Xpred <- model.matrix(object$terms, data = newdata)
  }
  npred <- dim(Xpred)[1]
  pred <- array(NA, c(npred, L, nsamp))
  for (i in 1:nsamp) {
    pred[, , i] <- tcrossprod(Xpred, betas$beta.samp[, , 
                                                     i])
  }
  dimnames(pred) <- list(obs = rownames(Xpred), tau = round(object$tau.g[object$reg.ix], 
                                                            4), samp = 1:nsamp)
  if (summarize) {
    pred <- apply(pred, c(1, 2), quantile, probs = 0.5)
  }
  return(pred)
}