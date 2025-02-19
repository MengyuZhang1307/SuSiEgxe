
optimize_prior_variance = function (optimize_V, betahat, shat2, prior_weights,
                                    alpha = NULL, post_mean2 = NULL,
                                    V_init = NULL, check_null_threshold = 0) { # V_init is a matrix
  V = unlist(V_init)
  lV_repara = c(log(V[1]), log(V[4]), V[2]/sqrt(V[1]*V[4]))
  p = length(betahat)/2
  if (optimize_V != "simple") {
    if(optimize_V == "optim") {
#c(log(c(max(c((betahat[1:p]^2-diag(shat2)[1:p]), 1)), max(c((betahat[(p+1):(2*p)]^2-diag(shat2)[(p+1):(2*p)]), 1)))),-0.1)
      lV = optim(par = c(log(c(2, 4)),-0.02/sqrt(2*4)), # log a, log c, and rho
          fn = neg.loglik.logscale,betahat = betahat,shat2 = shat2,
          prior_weights = prior_weights,method = "L-BFGS-B",
          lower = c(-30, -30, -1),
          upper = c(15, 15, 1)
	       )$par
      ## if the estimated one is worse than current one, don't change it.
      if(neg.loglik.logscale(lV, betahat = betahat,shat2 = shat2,prior_weights = prior_weights) >
         neg.loglik.logscale(lV_repara, betahat = betahat,
                             shat2 = shat2, prior_weights = prior_weights)){
        lV = lV_repara
      }
      #V = exp(lV)
      log_a = lV[1]
      log_c = lV[2]
      rho = lV[3]

      a = exp(log_a)
      c = exp(log_c)
      b = rho * sqrt(a * c)
      V = matrix(c(a, b, b, c), ncol = 2)
    } else if (optimize_V == "uniroot")
      stop("Current version does not support uniroot")
      #V = est_V_uniroot(betahat,shat2,prior_weights)
    else if (optimize_V == "EM")
      stop("Current version does not support EM")
      #V = sum(alpha * post_mean2)
    else
     stop("Invalid option for optimize_V method")
  }

  # Set V exactly 0 if that beats the numerical value by
  # check_null_threshold in loglik. check_null_threshold = 0.1 is
  # exp(0.1) = 1.1 on likelihood scale; it means that for parsimony
  # reasons we set estiate of V to zero, if its numerical estimate is
  # only "negligibly" different from zero. We use a likelihood ratio
  # of exp(check_null_threshold) to define "negligible" in this
  # context. This is fairly modest condition compared to, say, a
  # formal LRT with p-value 0.05. But the idea is to be lenient to
  # non-zeros estimates unless they are indeed small enough to be
  # neglible. See more intuition at
  # https://stephens999.github.io/fiveMinuteStats/LR_and_BF.html
  if (loglik(c(-Inf, -Inf, 0),betahat,shat2,prior_weights) +
      check_null_threshold >= loglik(lV,betahat,shat2,prior_weights))
    V = matrix(rep(0, 4), ncol = 2)
  return(V)
}

# In these functions, s2 represents residual_variance, and shat2 is an
# estimate of it.

# The log likelihood function for SER model (based on summary data
# betahat, shat2) as a function of prior variance V.
#
#' @importFrom Matrix colSums
#' @importFrom stats dnorm
#' @importFrom mvtnorm dmvnorm
loglik = function (lV, betahat, shat2, prior_weights) { # V is a vector
   p = length(betahat)/2
   log_a = lV[1]
   log_c = lV[2]
   rho = lV[3]

   a = exp(log_a)
   c = exp(log_c)
   b = rho * sqrt(a * c)
   V = matrix(c(a, b, b, c), ncol = 2)
   #V = as.matrix(rbind(V[1:2], V[2:3]))
  #log(bf) for each SNP
  #lbf = dnorm(betahat,0,sqrt(V + shat2),log = TRUE) -
  #      dnorm(betahat,0,sqrt(shat2),log = TRUE)
  lbf = mapply(function(i) {
            mvtnorm::dmvnorm(betahat[c(i, i+p)], mean = rep(0,2), sigma = as.matrix((V+shat2[c(i, i+p), c(i, i+p)])), log = TRUE) -
    	    mvtnorm::dmvnorm(betahat[c(i, i+p)], mean = rep(0,2), sigma = as.matrix((shat2[c(i, i+p), c(i, i+p)])), log = TRUE)
   	  }, 1:p)
  lpo = lbf + log(prior_weights + sqrt(.Machine$double.eps))

  # Deal with special case of infinite shat2 (e.g., happens if X does
  # not vary).
  infinite_ind = unique(c(is.infinite(diag(shat2[1:p, 1:p])),
			  is.infinite(diag(shat2[1:p, (p+1):(2*p)])),
			  is.infinite(diag(shat2[(p+1):(2*p), (p+1):(2*p)]))))
  lbf[c(infinite_ind, infinite_ind+p)] = 0
  lpo[c(infinite_ind, infinite_ind+p)] = 0

  maxlpo = max(lpo)
  w_weighted = exp(lpo - maxlpo)
  weighted_sum_w = sum(w_weighted)
  return(log(weighted_sum_w) + maxlpo)
}

neg.loglik.logscale = function(lV,betahat,shat2,prior_weights) { #lV is a vector
  return(-loglik(lV,betahat,shat2,prior_weights))
}


#' @importFrom Matrix colSums
#' @importFrom stats dnorm
loglik.grad = function(V, betahat, shat2, prior_weights) {

  # log(bf) for each SNP.
  lbf = dnorm(betahat,0,sqrt(V + shat2),log = TRUE) -
        dnorm(betahat,0,sqrt(shat2),log = TRUE)
  lpo = lbf + log(prior_weights + sqrt(.Machine$double.eps))

  # Deal with special case of infinite shat2 (e.g., happens if X does
  # not vary).
  lbf[is.infinite(shat2)] = 0
  lpo[is.infinite(shat2)] = 0

  maxlpo = max(lpo)
  w_weighted = exp(lpo - maxlpo)
  weighted_sum_w = sum(w_weighted)
  alpha = w_weighted / weighted_sum_w
  return(sum(alpha * lbf.grad(V,shat2,betahat^2/shat2)))
}

# Define loglikelihood and gradient as function of lV:=log(V)
# to improve numerical optimization
negloglik.grad.logscale = function (lV, betahat, shat2, prior_weights)
  -exp(lV) * loglik.grad(exp(lV),betahat,shat2,prior_weights)

# Vector of gradients of logBF_j for each j, with respect to prior
# variance V.
lbf.grad = function (V, shat2, T2) {
  l = 0.5*(1/(V + shat2)) * ((shat2/(V + shat2))*T2 - 1)
  l[is.nan(l)] = 0
  return(l)
}

lbf = function (V, shat2, T2) {
  l = 0.5*log(shat2/(V + shat2)) + 0.5*T2*(V/(V + shat2))
  l[is.nan(l)] = 0
  return(l)
}
