# @title update each effect once
# @param XtX a 2p by 2p matrix
# @param XtX_inv a return object of function S_inverse_crossprod
# @param Xty a 2p vector
# @param s_init a list with elements sigma2, V, alpha, mu, Xr
# @param estimate_prior_variance boolean indicating whether to
#   estimate prior variance
# @param estimate_prior_method The method used for estimating prior
#   variance, 'optim' or 'EM'.
# @param check_null_threshold float a threshold on the log scale to
#   compare likelihood between current estimate and zero the null
#
#' @importFrom Matrix diag
update_each_effect_ss_gxe = function (XtX, XtX_inv, Xty, s_init,
                                  estimate_prior_variance = FALSE,
                                  estimate_prior_method = "optim",
                                  check_null_threshold = 0) {
  if (!estimate_prior_variance)
    estimate_prior_method = "none"

  # Repeat for each effect to update.
  s = s_init
  L = nrow(s$alpha)
  if (L > 0) {
    for (l in 1:L) {
      #message(paste0("L = ", l))
      # Remove lth effect from fitted values.
      s$XtXr = s$XtXr - XtX %*% (rep(s$alpha[l,], 2) * s$mu[l,])

      # Compute residuals.
      if(length(Xty)!=length(s$XtXr)) stop("Length of Xty and XtXr do not match")
      XtR = Xty - s$XtXr #length of 2p
      res = single_effect_regression_ss_gxe(as.matrix(XtR),XtX_inv,s$V[[l]],
              s$sigma2,s$pi,estimate_prior_method,check_null_threshold)

      # Update the variational estimate of the posterior mean.
      s$mu[l,]    = as.vector(res$mu) # length = 2p
      s$alpha[l,] = res$alpha
      s$mu2[[l]]   = res$mu2 # a list with matrix of 2p*2p
      s$V[[l]]    = res$V
      s$lbf[l]    = res$lbf_model
      s$lbf_variable[l,] = res$lbf
      s$KL[l]     = -res$lbf_model +
        SER_posterior_e_loglik_ss_gxe(attr(XtX,"dKtK"),XtR,s$sigma2,
                                  rep(s$alpha[l,], 2) * res$mu, rep(s$alpha[l,], 2) * res$mu2)

      s$XtXr = s$XtXr + XtX %*% (rep(s$alpha[l,], 2) * s$mu[l,])
    }
  }
  s$XtXr = unname(as.matrix(s$XtXr))
  return(s)
}
