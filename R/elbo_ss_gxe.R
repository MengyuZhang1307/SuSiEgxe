# @title Get objective function from data and susie fit object.
# @param XtX a 2p by 2p matrix, K'K
# @param Xty a 2p vector, K'y,
# @param s a susie fit object
# @param yty a scaler, y'y, where y is centered to have mean 0
# @param n sample size
get_objective_ss_gxe = function (XtX, Xty, s, yty, n)
  Eloglik_ss_gxe(XtX,Xty,s,yty,n) - sum(s$KL)

# Expected loglikelihood for a susie fit.
Eloglik_ss_gxe = function (XtX, Xty, s, yty, n)
  -n/2*log(2*pi*s$sigma2) - 1/(2*s$sigma2) * get_ER2_ss_gxe(XtX,Xty,s,yty)

# Expected squared residuals.
get_ER2_ss_gxe = function (XtX, Xty, s, yty) {
  B = rep(s$alpha, 2) * s$mu # bl_bar
  XB2 = sum((B %*% XtX) * B)
  betabar = colSums(B)
  d = attr(XtX,"dKtK")
  d_postb2 = lapply(1:nrow(s$alpha), function(x) {s$mu2[[x]] * d %*% Diagonal(x = rep(s$alpha[x, ],2))}) # Posterior second moment.
  return(yty - 2*sum(betabar * Xty) + sum(betabar * (XtX %*% betabar)) -
         XB2 + sum(sapply(d_postb2, sum)))
}

# @title posterior expected loglikelihood for a single effect regression
# @param dXtX a 2p vector of diagonal elements of KtK
# @param Xty a 2p vector
# @param s2 the residual variance
# @param Eb the posterior mean of zeta (2p vector) (alpha * mu)
# @param Eb2 the posterior second moment of b (2p * 2p matrix) (alpha * mu2)
SER_posterior_e_loglik_ss_gxe = function (dXtX, Xty, s2, Eb, Eb2) {
  -0.5/s2 * (-2*sum(Eb*Xty) + sum(dXtX * Eb2))
}
 