#' @title Single effect regression
#'
#' @param Xty A 2p-vector.
#'
#' @param XtX_inv a return object of function S_inverse_crossprod
#'
#' @importFrom stats uniroot
#' @importFrom stats optim
#' @importFrom mvtnorm dmvnorm
#'
#' @keywords internal
#'
single_effect_regression_ss_gxe =
  function (Xty, XtX_inv, V = matrix(c(1, 0.5, 0.5, 1), ncol = 2), residual_variance = 1, prior_weights = NULL,
            optimize_V = c("none", "optim", "uniroot", "EM", "simple"),
            check_null_threshold = 0) {
  optimize_V = match.arg(optimize_V)
  #betahat = (1/dXtX) * Xty
  p = length(XtX_inv$dXtX)
  #betahat = cbind(XtX_inv$Kty[1:p],XtX_inv$Kty[(p+1):(2*p)])
  betahat = XtX_inv$Kty
  #shat2 = residual_variance/dXtX
  #shat2_tmp = residual_variance*cbind(XtX_inv$dXtX, XtX_inv$dXtZ, XtX_inv$dXtZ, XtX_inv$dZtZ)
  #shat2_list = lapply(1:p, function(i) {matrix(shat2_tmp[i, ], nrow = 2, byrow = TRUE)})
  shat2_mat = residual_variance*rbind(cbind(Diagonal(x = XtX_inv$dXtX), Diagonal(x = XtX_inv$dXtZ)),
                                  cbind(Diagonal(x = XtX_inv$dXtZ), Diagonal(x = XtX_inv$dZtZ)))

  if (is.null(prior_weights))
    prior_weights = rep(1/p,p)

  if (optimize_V != "EM" && optimize_V != "none")
    V = optimize_prior_variance(optimize_V,betahat,shat2_mat,prior_weights,
                                alpha = NULL,post_mean2 = NULL,V_init = V,
                                check_null_threshold = check_null_threshold)

  # log(po) = log(BF * prior) for each SNP
  #lbf = dnorm(betahat,0,sqrt(V + shat2),log = TRUE) -
  #      dnorm(betahat,0,sqrt(shat2),log = TRUE)
  # Load necessary libraries

  lbf = mapply(function(i) {
            mvtnorm::dmvnorm(betahat[c(i, i+p)], mean = rep(0,2), sigma = as.matrix((V+shat2_mat[c(i, i+p), c(i, i+p)])), log = TRUE) -
    	    mvtnorm::dmvnorm(betahat[c(i, i+p)], mean = rep(0,2), sigma = as.matrix((shat2_mat[c(i, i+p), c(i, i+p)])), log = TRUE)
  	  }, 1:p)
  lpo = lbf + log(prior_weights + sqrt(.Machine$double.eps))

  # Deal with special case of infinite shat2 (e.g., happens if X does
  # not vary).
  infinite_ind = unique(c(is.infinite(residual_variance*XtX_inv$dXtX),
			  is.infinite(residual_variance*XtX_inv$dXtZ),
			  is.infinite(residual_variance*XtX_inv$dZtZ)))
  lbf[c(infinite_ind, infinite_ind+p)] = 0
  lpo[c(infinite_ind, infinite_ind+p)] = 0
  maxlpo = max(lpo)

  # w is proportional to
  #
  #   posterior odds = BF * prior,
  #
  # but subtract max for numerical stability.
  w_weighted = exp(lpo - maxlpo)
  weighted_sum_w = sum(w_weighted)
  alpha = w_weighted / weighted_sum_w

  #post_var = (1/V + dXtX/residual_variance)^(-1) # Posterior variance.
  #post_var = lapply(1:p, function(i) {V %*% solve(shat2_mat[c(i, i+p), c(i, i+p)] + V) %*% shat2_mat[c(i, i+p), c(i, i+p)]})
  #post_mean = (1/residual_variance) * post_var * Xty
  #post_mean_tmp = lapply(1:p, function(i) {(1/residual_variance) * post_var[[i]] %*% Xty[c(i, (i+p))]}) # Xj, Zj
  #post_mean_tmp1 <- unlist(post_mean_tmp)
  #dim(post_mean_tmp1) <- c(2, p)
  #post_mean <- as.vector(t(post_mean_tmp1)) # X, Z, length of 2*p
  #post_mean2 = post_var + post_mean^2 # Second moment.
  #post_mean_2 = lapply(1:p, function(i) {tcrossprod(post_mean_tmp[[i]]) + post_var[[i]]})


  V_spread = V %x% diag(p)
  shat2plusV = shat2_mat + V_spread
  shat2plusV_inv = S_inverse_crossprod(diag(shat2plusV[1:p, 1:p]), diag(shat2plusV[1:p, (p+1):(2*p)]), diag(shat2plusV[(p+1):(2*p), (p+1):(2*p)]), rep(NA, 2*p), inv_only = TRUE)
  shat2plusV_inv1 = rbind(cbind(Diagonal(x = shat2plusV_inv$dXtX), Diagonal(x = shat2plusV_inv$dXtZ)),
                          cbind(Diagonal(x = shat2plusV_inv$dXtZ), Diagonal(x = shat2plusV_inv$dZtZ)))
  post_var = V_spread %*% shat2plusV_inv1 %*% shat2_mat
  post_mean = 1/residual_variance * post_var %*% Xty
  post_mean_2 = rbind(cbind(Diagonal(x = post_mean[1:p]^2), Diagonal(x = post_mean[1:p] * post_mean[(p+1):(2*p)])),
		      cbind(Diagonal(x = post_mean[1:p] * post_mean[(p+1):(2*p)]), Diagonal(x = post_mean[(p+1):(2*p)]^2))) + post_var
  lbf_model = maxlpo + log(weighted_sum_w) # Analogue of loglik in the
                                           # non-summary case.

  if (optimize_V == "EM")
    V = optimize_prior_variance(optimize_V,betahat,shat2,prior_weights,alpha,
          post_mean2,check_null_threshold = check_null_threshold)
  return(list(alpha = alpha,mu = post_mean,mu2 = post_mean_2,lbf = lbf,
              V = V,lbf_model = lbf_model))
}
