#' @title Initialize a susie object


# Set default susie initialization.
init_setup_gxe = function (n, p, L, scaled_prior_variance, residual_variance,
                       prior_weights, null_weight, varY, standardize) {
  if (!is.numeric(unlist(scaled_prior_variance)) || sum(unlist(scaled_prior_variance)[c(1,4)] < 0) > 0)
    stop("Scaled prior variance should be positive number")
  #if (sum(scaled_prior_variance > 1)>0 && standardize)
  #  stop("Scaled prior variance should be no greater than 1 when ",
  #       "standardize = TRUE")
  if(is.null(residual_variance))
    residual_variance = varY
  if(is.null(prior_weights)){
    prior_weights = rep(1/p,p)
  }else{
    if(all(prior_weights == 0)){
      stop("Prior weight should greater than 0 for at least one variable.")
    }
    prior_weights = prior_weights / sum(prior_weights)
  }
  if(length(prior_weights) != p)
    stop("Prior weights must have length p")
  if (p < L)
    L = p

  if (L > 0) {
    if(!is.list(scaled_prior_variance))
      stop("When L > 0, prior_variance should be a list")
    v = lapply(scaled_prior_variance, function(x) x * varY)
  } else {v = scaled_prior_variance*varY}

  s = list(alpha  = matrix(1/p,nrow = L,ncol = p),
           mu     = matrix(0,nrow = L,ncol = 2*p),
           mu2    = replicate(L, matrix(0, nrow = 2*p, ncol = 2*p), simplify = FALSE), # L*each element is a 2px2p matrix
           Xr     = rep(0,n),
           KL     = rep(as.numeric(NA),L),
           lbf    = rep(as.numeric(NA),L),
           lbf_variable = matrix(as.numeric(NA),L,p),
           sigma2 = residual_variance,
           V      = v,
           pi     = prior_weights)
  if (is.null(null_weight))
    s$null_index = 0
  else
    s$null_index = p
  class(s) = "susie"
  return(s)
}

# Update a susie fit object in order to initialize susie model.
init_finalize_gxe = function (s, X = NULL, Xr = NULL) {
  if(is.matrix(s$V))
    s$V = replicate(nrow(s$alpha), s$V, simplify = FALSE) # a list

  # Check sigma2.
  if (!is.numeric(s$sigma2))
    stop("Input residual variance sigma2 must be numeric")

  # Avoid problems with dimension if input is a 1 x 1 matrix.
  s$sigma2 = as.numeric(s$sigma2)
  if (length(s$sigma2) != 1)
    stop("Input residual variance sigma2 must be a scalar")
  if (s$sigma2 <= 0)
    stop("Residual variance sigma2 must be positive (is your var(Y) zero?)")

  # check prior variance
  if (!is.numeric(unlist(s$V)))
    stop("Input prior variance must be numeric")
  #if (!all(unlist(s$V)[c(1,4)] >= 0))
  #  stop("prior variance must be non-negative")
  if (!all(dim(s$mu) == dim(s$mu2)))
    stop("dimension of mu and mu2 in input object do not match")
  if (nrow(s$mu) != nrow(s$alpha) & ncol(s$mu) != 2*ncol(s$alpha))
    stop("dimension of mu and alpha in input object do not match")
  if (nrow(s$alpha) != length(s$V))
    stop("Input prior variance V must have length of nrow of alpha in ",
         "input object")

  # Update Xr.
  #if (!missing(Xr))
  #  s$Xr = Xr
  #if (!missing(X))
  #  s$Xr = compute_Xb(X,colSums(s$mu * s$alpha))

  # Reset KL and lbf.
  s$KL = rep(as.numeric(NA),nrow(s$alpha))
  s$lbf = rep(as.numeric(NA),nrow(s$alpha))
  class(s) = "susie"
  return(s)
}