#' @title susie_gxe
#'
#' @param XtX A p by p matrix \eqn{X'X}
#'
#' @param XtZ A p by p matrix \eqn{X'Z}
#'
#' @param ZtZ A p by p matrix \eqn{Z'Z}
#'
#' @param Xty A p-vector \eqn{X'y}
#'
#' @param yty A scalar \eqn{y'y} in which y is centered to have mean
#'   zero.
#'
#' @param n The sample size.
#'
#' @param X_colmeans A p-vector of column means of \code{X}. If both
#'   \code{X_colmeans} and \code{y_mean} are provided, the intercept
#'   is estimated; otherwise, the intercept is NA.
#'
#' @param y_mean A scalar containing the mean of \code{y}. If both
#'   \code{X_colmeans} and \code{y_mean} are provided, the intercept
#'   is estimated; otherwise, the intercept is NA.
#'
#' @param maf Minor allele frequency; to be used along with
#'   \code{maf_thresh} to filter input summary statistics.
#'
#' @param maf_thresh Variants having a minor allele frequency smaller
#'   than this threshold are not used.
#'
#' @param r_tol Tolerance level for eigenvalue check of positive
#'   semidefinite matrix of R.
#'
#' @param check_input If \code{check_input = TRUE},
#'   \code{susie_suff_stat} performs additional checks on \code{XtX} and
#'   \code{Xty}. The checks are: (1) check that \code{XtX} is positive
#'   semidefinite; (2) check that \code{Xty} is in the space spanned by
#'   the non-zero eigenvectors of \code{XtX}.
#'
#' @param check_prior If \code{check_prior = TRUE}, it checks if the
#'   estimated prior variance becomes unreasonably large (comparing with
#'   10 * max(abs(chi-square))^2).
#'
#' @param n_purity Passed as argument \code{n_purity} to
#'   \code{\link{susie_get_cs}}.
#'
#' @param ... Additional arguments to provide backward compatibility
#'   with earlier versions of \code{susie_suff_stat}.
#'
#' @export
#'
susie_suff_stat_gxe = function (XtX, XtZ, ZtZ, Xty, yty, n,
                            X_colmeans = NA, y_mean = NA,
                            maf = NULL, maf_thresh = 0, L = 10,
                            scaled_prior_variance = matrix(c(0.5, 0.1, 0.1, 0.5), ncol = 2),
                            residual_variance = NULL,
                            estimate_residual_variance = TRUE,
                            estimate_prior_variance = TRUE,
                            estimate_prior_method = c("optim","EM","simple"),
                            check_null_threshold = 0, prior_tol = 1e-9,
                            r_tol = 1e-08, prior_weights = NULL,
                            null_weight = 0, standardize = TRUE,
                            max_iter = 100, s_init = NULL, coverage = 0.95,
                            min_abs_corr = 0.5, tol = 1e-3,
                            verbose = FALSE, track_fit = FALSE,
                            check_input = FALSE, refine = FALSE,
                            check_prior = FALSE, n_purity = 100, ...) {

  # Check for use of arguments that are now deprecated.
  args <- list(...)
  if (any(is.element(names(args),c("bhat","shat","R","var_y"))))
    stop("susie_suff_stat no longer accepts inputs bhat, shat, R or var_y; ",
         "these inputs are now accepted by susie_rss instead")

  # Process input estimate_prior_method.
  estimate_prior_method = match.arg(estimate_prior_method)

  if (missing(n))
    stop("n must be provided")
  if (n <= 1)
    stop("n must be greater than 1")

  # Check sufficient statistics.
  missing_XtX = c(missing(XtX), missing(XtZ), missing(ZtZ), missing(Xty), missing(yty))

  if (all(missing_XtX))
    stop("Please provide all of XtX, XtZ, ZtZ, Xty, yty, n")

  if (ncol(XtX) > 1000 & !requireNamespace("Rfast",quietly = TRUE))
    warning_message("For large R or large XtX, consider installing the ",
                    "Rfast package for better performance.", style="hint")

  # Check input XtX.
  if (2*ncol(XtX) != length(Xty) || 2*ncol(XtZ) != length(Xty) || 2*ncol(ZtZ) != length(Xty))
    stop(paste0("The dimension of XtX (",nrow(XtX)," by ",ncol(XtX),
                ") or XtZ (",nrow(XtZ)," by ",ncol(XtZ),
                ") or ZtZ (",nrow(ZtZ)," by ",ncol(ZtZ),
                ")does not agree with expected (",length(Xty)/2," by ",
                length(Xty)/2,")"))
  if (!is_symmetric_matrix(XtX)) {
    warning_message("XtX is not symmetric; forcing XtX to be symmetric by ",
            "replacing XtX with (XtX + t(XtX))/2")
    XtX = XtX + t(XtX)
    XtX = XtX/2
  }
  if (!is_symmetric_matrix(XtZ)) {
    warning_message("XtZ is not symmetric; forcing XtZ to be symmetric by ",
            "replacing XtZ with (XtZ + t(XtZ))/2")
    XtZ = XtZ + t(XtZ)
    XtZ = XtZ/2
  }
  if (!is_symmetric_matrix(ZtZ)) {
    warning_message("ZtZ is not symmetric; forcing ZtZ to be symmetric by ",
            "replacing ZtZ with (ZtZ + t(ZtZ))/2")
    ZtZ = ZtZ + t(ZtZ)
    ZtZ = ZtZ/2
  }

  # MAF filter.
  if (!is.null(maf)) {
    if (length(maf) != length(Xty))
      stop(paste("The length of maf does not agree with expected",length(Xty)))
    id = which(maf > maf_thresh)
    XtX = XtX[id,id]
    XtZ = XtZ[id,id]
    ZtZ = ZtZ[id,id]
    Xty = Xty[c(id,(id+nrow(XtX)))]
  }

  if (any(is.infinite(Xty)))
    stop("Input Xty contains infinite values")
  if (!(is.double(XtX) & is.matrix(XtX)) & !inherits(XtX,"CsparseMatrix"))
    stop("Input XtX must be a double-precision matrix, or a sparse matrix")
  if (anyNA(XtX))
    stop("Input XtX matrix contains NAs")
  if (!(is.double(XtZ) & is.matrix(XtZ)) & !inherits(XtZ,"CsparseMatrix"))
    stop("Input XtZ must be a double-precision matrix, or a sparse matrix")
  if (anyNA(XtZ))
    stop("Input XtZ matrix contains NAs")
  if (!(is.double(ZtZ) & is.matrix(ZtZ)) & !inherits(ZtZ,"CsparseMatrix"))
    stop("Input ZtZ must be a double-precision matrix, or a sparse matrix")
  if (anyNA(ZtZ))
    stop("Input ZtZ matrix contains NAs")

  # Replace NAs in Xty with zeros.
  if (anyNA(Xty)) {
    warning_message("NA values in Kty are replaced with 0")
    Xty[is.na(Xty)] = 0
  }

  KtK = rbind(
    cbind(XtX, XtZ),
    cbind(XtZ, ZtZ)
  )

  if (check_input) {

    # Check whether KtK is positive semidefinite.
    semi_pd = check_semi_pd(KtK,r_tol)
    if (!semi_pd$status)
      stop("KtK is not a positive semidefinite matrix")

    # Check whether Xty in space spanned by the non-zero eigenvectors of KtK
    proj = check_projection(semi_pd$matrix,Xty)
    if (!proj$status)
      warning_message("Kty does not lie in the space of the non-zero eigenvectors ",
              "of KtK")
  }

  if (is.numeric(null_weight) && null_weight == 0)
    null_weight = NULL
  #if (!is.null(null_weight)) { #tbd
  #  if (!is.numeric(null_weight))
  #    stop("Null weight must be numeric")
  #  if (null_weight < 0 || null_weight >= 1)
  #    stop("Null weight must be between 0 and 1")
  #  if (is.null(prior_weights))
  #    prior_weights = c(rep(1/ncol(XtX)*(1-null_weight),ncol(XtX)),null_weight)
  #  else
  #    prior_weights = c(prior_weights*(1 - null_weight),null_weight)
  #  XtX = cbind(rbind(XtX,0),0)
  #  Xty = c(Xty,0)
  #}

  p = ncol(XtX)

  if (standardize) {
    dKtK = diag(KtK)
    csd = sqrt(dKtK/(n-1))
    csd[csd == 0] = 1
    KtK = t((1/csd) * KtK) / csd
    Xty = Xty / csd
  } else
    csd = rep(1,length = p*2)

  #I <- diag(p)
  #diag_ind_mat <- kronecker(matrix(1, nrow = 2, ncol = 2), I) # 2x2 block matrix with
  #attr(KtK,"d") = KtK * diag_ind_mat
  attr(KtK,"dXtX") = diag(XtX)
  attr(KtK,"dXtZ") = diag(XtZ)
  attr(KtK,"dZtZ") = diag(ZtZ)
  attr(KtK, "dKtK") = rbind(cbind(Diagonal(x = attr(KtK,"dXtX")), Diagonal(x = attr(KtK,"dXtZ"))),
               cbind(Diagonal(x = attr(KtK,"dXtZ")), Diagonal(x = attr(KtK,"dZtZ"))))
  attr(KtK,"scaled:scale") = csd

  # Check that X_colmeans has length 1 or p.
  if (length(X_colmeans) == 1)
    X_colmeans = rep(X_colmeans,p*2)
  if (length(X_colmeans) != p*2)
    stop("The length of X_colmeans does not agree with number of variables")

  # Initialize susie fit.
  # beta version: residual_variance  = NULL, prior_weights = NULL, null_weights = NULL
  s = init_setup_gxe(0,p,L,scaled_prior_variance,residual_variance,prior_weights,
                 null_weight,yty/(n-1),standardize)
  s$Xr = NULL
  s$XtXr = rep(0,2*p)

  if (!missing(s_init)&& !is.null(s_init)) {
    if (!inherits(s_init,"susie"))
      stop("s_init should be a susie object")
    if (max(s_init$alpha) > 1 || min(s_init$alpha) < 0)
      stop("s_init$alpha has invalid values outside range [0,1]; please ",
           "check your input")

    # First, remove effects with s_init$V = 0
    s_init = susie_prune_single_effects(s_init)
    num_effects = nrow(s_init$alpha)
    if(missing(L)){ # if L is not specified, we set it as in s_init.
      L = num_effects
    }else if(min(p, L) < num_effects){
      warning_message(paste("Specified number of effects L =",min(p, L),
                    "is smaller than the number of effects",num_effects,
                    "in input SuSiE model. The initialized SuSiE model will have",
                    num_effects,"effects."))
      L = num_effects
    }
    # expand s_init if L > num_effects.
    s_init = susie_prune_single_effects(s_init, min(p, L), s$V)
    s = modifyList(s,s_init)
    s = init_finalize(s,X = XtX)
    s$XtXr = s$Xr
    s$Xr = NULL
  } else
    s = init_finalize_gxe(s)

  # Initialize elbo to NA.
  elbo = rep(as.numeric(NA),max_iter + 1)
  elbo[1] = -Inf;
  max_pip = rep(as.numeric(NA),max_iter + 1)
  max_pip[1] = 1;
  tracking = list()
  s$pip_tmp = rep(0, p)

  KtK_inv = S_inverse_crossprod(attr(KtK,"dXtX"), attr(KtK,"dXtZ"), attr(KtK,"dZtZ"), Xty)

  bhat = KtK_inv$Kty
  #shat = cbind(sqrt(s$sigma2*KtK_inv$dXtX), sqrt(s$sigma2*KtK_inv$dZtZ))
  #if(length(bhat) != length(shat))
  #  stop("bhat does not match the lenght of shat before interation in susie_ss_gxe.R")
  #z = bhat/shat
  chisq_tmp = bhat * crossprod(as.matrix(attr(KtK,"dKtK")), as.matrix(bhat))/s$sigma2
  chisq = chisq_tmp[1:p]+chisq_tmp[(p+1):(2*p)]
  chisqm = max(abs(chisq[!is.nan(chisq)]))

  Sys.time()
  for (i in 1:max_iter) {
    print(paste0("iter: ", i))
    if (track_fit)
      tracking[[i]] = susie_slim(s)
    system.time({s = update_each_effect_ss_gxe(KtK,KtK_inv,Xty,s,estimate_prior_variance,
                              estimate_prior_method,check_null_threshold)})
    if(check_prior){
      if(any(unlist(s$V) > 100*(chisqm^2))){
        stop('The estimated prior variance is unreasonably large.
       This is usually caused by mismatch between the summary statistics and the LD matrix.
             Please check the input.')
      }
    }

    if (verbose) # Need KL
      {print(paste0("objective: ",get_objective_ss_gxe(KtK,Xty,s,yty,n)))
       print(paste0("Max pip difference:", max_pip[i+1]))}

    # Compute objective before updating residual variance because part
    # of the objective s$kl has already been computed under the
    # residual variance before the update.
    elbo[i+1] = get_objective_ss_gxe(KtK,Xty,s,yty,n)
    if(is.infinite(elbo[i+1])){
      stop('The objective becomes infinite. Please check the input.')
    }

    pip = susie_get_pip(s,prune_by_cs = FALSE,prior_tol = prior_tol)
    max_pip[i+1] = max(abs(pip-s$pip_tmp))
    if (max_pip[i+1] < tol) {
      s$converged = TRUE
      s$pip_tmp = NULL
      break
    } else {s$pip_tmp = pip}

    #if ((elbo[i+1] - elbo[i]) < tol) {
    #  s$converged = TRUE
    #  break
    #}
    if (estimate_residual_variance) { # TRUE for in-sample R
      est_sigma2 = estimate_residual_variance_ss(KtK,Xty,s,yty,n)
      if (est_sigma2 < 0)
        stop("Estimating residual variance failed: the estimated value ",
             "is negative")
      s$sigma2 = est_sigma2
      if (verbose)
        print(paste0("objective:",get_objective_ss(KtK,Xty,s,yty,n)))
        print(paste0("Max pip difference:", max_pip[i+1]))
    }
  }
  Sys.time()

  elbo = elbo[2:(i+1)] # Remove first (infinite) entry, and trailing NAs.
  max_pip = max_pip[2:(i+1)]
  s$max_pip = max_pip
  s$elbo = elbo
  s$niter = i

  if (is.null(s$converged)) {
    warning(paste("IBSS algorithm did not converge in",max_iter,"iterations!
                  Please check consistency between summary statistics and LD matrix.
                  See https://stephenslab.github.io/susieR/articles/susierss_diagnostic.html"))
    s$converged = FALSE
  }

  s$X_column_scale_factors = attr(KtK,"scaled:scale") # hard code, scale = 1

  # Compute intercept.
  s$intercept = y_mean -
    sum(X_colmeans * (colSums(rep(s$alpha, each = 2) * s$mu)/s$X_column_scale_factors))


  if (track_fit)
    s$trace = tracking

  # SuSiE CS and PIP.
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
    if(any(!(diag(KtK) %in% c(0,1)))){
      s$sets = susie_get_cs(s,coverage = coverage,Xcorr = muffled_cov2cor(KtK),
                            min_abs_corr = min_abs_corr,
                            check_symmetric = FALSE,
                            n_purity = n_purity)
    }else
      s$sets = susie_get_cs(s,coverage = coverage,Xcorr = KtK,
                            min_abs_corr = min_abs_corr,
                            check_symmetric = FALSE,
                            n_purity = n_purity)
    s$pip = susie_get_pip(s,prune_by_cs = FALSE,prior_tol = prior_tol)
  }

  if (!is.null(colnames(KtK))) {
    variable_names = colnames(KtK)
    if (!is.null(null_weight)) {
      variable_names[length(variable_names)] = "null"
      names(s$pip) = variable_names[-p]
    } else
      names(s$pip) = variable_names[1:p]
    colnames(s$alpha) = variable_names[1:p]
    colnames(s$mu)    = variable_names
    #colnames(s$mu2)   = variable_names
    colnames(s$lbf_variable) = variable_names[1:p]
  }

  if (refine) {
    stop("do not support refine = TRUE")
    if (!missing(s_init) && !is.null(s_init))
      warning("The given s_init is not used in refinement")
    if (!is.null(null_weight) && null_weight != 0) {
      ## if null_weight is specified
      ## we remove the extra 0 column
      XtX = XtX[1:(p-1), 1:(p-1)]
      Xty = Xty[1:(p-1)]
      pw_s = s$pi[-s$null_index]/(1 - null_weight)
    } else
      pw_s = s$pi
    conti = TRUE
    while (conti & length(s$sets$cs)>0) {
      m = list()
      for(cs in 1:length(s$sets$cs)) {
        pw_cs = pw_s
        pw_cs[s$sets$cs[[cs]]] = 0
        if(all(pw_cs == 0)){
          break
        }
        s2 = susie_suff_stat_gxe(XtX = XtX, XtZ = XtZ, ZtZ = ZtZ, Xty = Xty, yty = yty, n = n, L = L,
            X_colmeans = X_colmeans, y_mean = y_mean,
            prior_weights = pw_cs, s_init = NULL,
            scaled_prior_variance = scaled_prior_variance,
            residual_variance = residual_variance,
            estimate_residual_variance = estimate_residual_variance,
            estimate_prior_variance = estimate_prior_variance,
            estimate_prior_method = estimate_prior_method,
            check_null_threshold = check_null_threshold, prior_tol = prior_tol,
            r_tol = r_tol, max_iter = max_iter,
            null_weight = null_weight, standardize = standardize,
            coverage = coverage, min_abs_corr = min_abs_corr, tol = tol,
            verbose = FALSE, track_fit = FALSE, check_input = FALSE,
            refine = FALSE)
        sinit2 = s2[c("alpha","mu","mu2")]
        class(sinit2) = "susie"
        s3 = susie_suff_stat(XtX = XtX, XtZ = XtZ, ZtZ = ZtZ, Xty = Xty, yty = yty, n = n, L = L,
            X_colmeans = X_colmeans, y_mean = y_mean,
            prior_weights = pw_s, s_init = sinit2,
            scaled_prior_variance = scaled_prior_variance,
            residual_variance = residual_variance,
            estimate_residual_variance = estimate_residual_variance,
            estimate_prior_variance = estimate_prior_variance,
            estimate_prior_method = estimate_prior_method,
            check_null_threshold = check_null_threshold, prior_tol = prior_tol,
            r_tol = r_tol, max_iter = max_iter, null_weight = null_weight,
            standardize = standardize, coverage = coverage,
            min_abs_corr = min_abs_corr, tol = tol, verbose = FALSE,
            track_fit = FALSE, check_input = FALSE, refine = FALSE)
        m = c(m,list(s3))
      }
      if(length(m) == 0){
        conti = FALSE
      }else{
        elbo = sapply(m, function(x) susie_get_objective(x))
        if ((max(elbo) - susie_get_objective(s)) <= 0)
          conti = FALSE
        else
          s = m[[which.max(elbo)]]
      }
    }
  }

  return(s)
}
