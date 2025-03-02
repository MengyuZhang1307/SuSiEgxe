# @title Check whether A is positive semidefinite
# @param A a symmetric matrix
# @return a list of result:
# \item{matrix}{The matrix with eigen decomposition}
# \item{status}{whether A is positive semidefinite}
# \item{eigenvalues}{eigenvalues of A truncated by r_tol}
check_semi_pd = function (A, tol) {
  attr(A,"eigen") = eigen(A,symmetric = TRUE)
  v = attr(A,"eigen")$values
  v[abs(v) < tol] = 0
  return(list(matrix      = A,
              status      = !any(v < 0),
              eigenvalues = v))
}

# @title Check whether b is in space spanned by the non-zero eigenvectors
#   of A
# @param A a p by p matrix
# @param b a length p vector
# @return a list of result:
# \item{status}{whether b in space spanned by the non-zero
#  eigenvectors of A}
# \item{msg}{msg gives the difference between the projected b and b if
#   status is FALSE}
check_projection = function (A, b) {
  if (is.null(attr(A,"eigen")))
    attr(A,"eigen") = eigen(A,symmetric = TRUE)
  v = attr(A,"eigen")$values
  B = attr(A,"eigen")$vectors[,v > .Machine$double.eps]
  msg = all.equal(as.vector(B %*% crossprod(B,b)),as.vector(b),
                  check.names = FALSE)
  if (!is.character(msg))
    return(list(status = TRUE,msg = NA))
  else
    return(list(status = FALSE,msg = msg))
}

# Check for symmetric matrix.
#' @keywords internal
is_symmetric_matrix = function (x) {
  if (requireNamespace("Rfast",quietly = TRUE))
    return(Rfast::is.symmetric(x))
  else
    return(Matrix::isSymmetric(x))
}

S_inverse_crossprod <- function(A, B, D, r, inv_only = FALSE) {
  # A, B, D are vectors representing the diagonal elements of diagonal matrices A, B, D
  # B = C by assumption (so use B for both B and C)
  # Inverses of A and D (element-wise inversion)
  A_inv <- 1 / A
  D_inv <- 1 / D
  # Compute Schur complement: S = D - B * A_inv * B (element-wise operation)
  S <- D - B * A_inv * B
  # Check if Schur complement is invertible (no zeros in S)
  if (any(S == 0)) {
    stop("Schur complement is singular; cannot compute the inverse.")
  }
  # Compute S inverse (element-wise)
  S_inv <- 1 / S
  # Decompose vector r into two parts corresponding to the block matrix
  r1 <- r[1:length(A)]  # Top part of r (corresponding to A and B)
  r2 <- r[(length(A) + 1):length(r)]  # Bottom part of r (corresponding to B and D)
  # Compute the block inverse using element-wise operations
  M11 <- A_inv + A_inv * B * S_inv * B * A_inv  # Top-left block
  M12 <- -A_inv * B * S_inv
  M22 <- S_inv                                  # Bottom-right block
  # Construct the inverse block matrix M_inv
  # M_inv <- rbind(
  #  cbind(Diagonal(x = M11), Diagonal(x = M12)),
  #  cbind(Diagonal(x = M12), Diagonal(x = M22))
  #)
  crossprod_top = crossprod_bottom = rep(0, length(A))
  if(!inv_only) {
    # Compute each part of the cross product without matrix multiplication
    # Top part: A_inv * r1 + A_inv * B * S_inv * B * A_inv * r1 - A_inv * B * S_inv * r2
    crossprod_top <- M11 * r1 + M12 * r2
    # Bottom part: -S_inv * B * A_inv * r1 + S_inv * r2
    crossprod_bottom <- M12 * r1 + M22 * r2
  }

  return(list(dXtX = M11, dXtZ = M12, dZtZ = M22, Kty= c(crossprod_top, crossprod_bottom)))

}

# Subsample and compute min, mean, median and max abs corr.
#
#' @importFrom stats median
get_purity = function (pos, X, Xcorr, squared = FALSE, n = 100,
                       use_rfast) {
  if (missing(use_rfast))
    use_rfast = requireNamespace("Rfast",quietly = TRUE)
  if (use_rfast) {
    get_upper_tri = Rfast::upper_tri
    get_median    = Rfast::med
  } else {
    get_upper_tri = function (R) R[upper.tri(R)]
    get_median    = stats::median
  }
  if (length(pos) == 1)
    return(c(1,1,1))
  else {

    # Subsample the columns if necessary.
    if (length(pos) > n)
      pos = sample(pos,n)

    if (is.null(Xcorr)) {
      X_sub = X[,pos]
      X_sub = as.matrix(X_sub)
      value = abs(get_upper_tri(muffled_corr(X_sub)))
    } else
      value = abs(get_upper_tri(Xcorr[pos,pos]))
    if (squared)
      value = value^2
    return(c(min(value),
             sum(value)/length(value),
             get_median(value)))
  }
}

#' @title Get credible sets
#'
#' @param X n by p matrix of values of the p variables (covariates) in
#'   n samples. When provided, correlation between variables will be
#'   computed and used to remove CSs whose minimum correlation among
#'   variables is smaller than \code{min_abs_corr}.
#'
#' @param Xcorr p by p matrix of correlations between variables
#'   (covariates). When provided, it will be used to remove CSs whose
#'   minimum correlation among variables is smaller than
#'   \code{min_abs_corr}.
#'
#' @param coverage A number between 0 and 1 specifying desired
#'   coverage of each CS.
#'
#' @param min_abs_corr A "purity" threshold for the CS. Any CS that
#'   contains a pair of variables with correlation less than this
#'   threshold will be filtered out and not reported.
#'
#' @param dedup If \code{dedup = TRUE}, remove duplicate CSs.
#'
#' @param squared If \code{squared = TRUE}, report min, mean and
#' median of squared correlation instead of the absolute correlation.
#'
#' @param check_symmetric If \code{check_symmetric = TRUE}, perform a
#'   check for symmetry of matrix \code{Xcorr} when \code{Xcorr} is
#'   provided (not \code{NULL}).
#'
#' @param n_purity The maximum number of credible set (CS) variables
#'   used in calculating the correlation (\dQuote{purity})
#'   statistics. When the number of variables included in the CS is
#'   greater than this number, the CS variables are randomly subsampled.
#'
#' @param use_rfast Use the Rfast package for the purity calculations.
#'   By default \code{use_rfast = TRUE} if the Rfast package is
#'   installed.
#'
#'
#' @export
#'
susie_get_cs = function (res, X = NULL, Xcorr = NULL, coverage = 0.95,
                         min_abs_corr = 0.5, dedup = TRUE, squared = FALSE,
                         check_symmetric = TRUE, n_purity = 100, use_rfast) {
  if (!is.null(X) && !is.null(Xcorr))
    stop("Only one of X or Xcorr should be specified")
  if (check_symmetric) {
    if (!is.null(Xcorr) && !is_symmetric_matrix(Xcorr)) {
      warning_message("Xcorr is not symmetric; forcing Xcorr to be symmetric",
                  "by replacing Xcorr with (Xcorr + t(Xcorr))/2")
      Xcorr = Xcorr + t(Xcorr)
      Xcorr = Xcorr/2
    }
  }

  null_index = 0
  include_idx = rep(TRUE,nrow(res$alpha))
  if (!is.null(res$null_index)) null_index = res$null_index
  if (is.numeric(res$V)) include_idx = unlist(lapply(res$V, function(x) {x[1,1] > 1e-9 & x[2,2] > 1e-9}))
  # L x P binary matrix.
  status = in_CS(res$alpha,coverage)

  # L-list of CS positions.
  cs = lapply(1:nrow(status),function(i) which(status[i,]!=0))
  claimed_coverage = sapply(1:length(cs),
                            function (i) sum(res$alpha[i,][cs[[i]]]))
  include_idx = include_idx * (lapply(cs,length) > 0)

  # FIXME: see issue 21
  # https://github.com/stephenslab/susieR/issues/21
  if (dedup)
    include_idx = include_idx * (!duplicated(cs))
  include_idx = as.logical(include_idx)
  if (sum(include_idx) == 0)
    return(list(cs = NULL,
                coverage = NULL,
                requested_coverage = coverage))
  cs = cs[include_idx]
  claimed_coverage = claimed_coverage[include_idx]

  # Compute and filter by "purity".
  if (missing(use_rfast))
    use_rfast = requireNamespace("Rfast",quietly = TRUE)
  if (is.null(Xcorr) && is.null(X)) {
    names(cs) = paste0("L",which(include_idx))
    return(list(cs = cs,
                coverage = claimed_coverage,
                requested_coverage = coverage))
  } else {
    purity = NULL
    for (i in 1:length(cs)) {
      if (null_index > 0 && null_index %in% cs[[i]])
        purity = rbind(purity,c(-9,-9,-9))
      else
        purity =
          rbind(purity,
            matrix(get_purity(cs[[i]],X,Xcorr,squared,n_purity,use_rfast),1,3))
    }
    purity = as.data.frame(purity)
    if (squared)
      colnames(purity) = c("min.sq.corr","mean.sq.corr","median.sq.corr")
    else
      colnames(purity) = c("min.abs.corr","mean.abs.corr","median.abs.corr")
    threshold = ifelse(squared,min_abs_corr^2,min_abs_corr)
    is_pure = which(purity[,1] >= threshold)
    if (length(is_pure) > 0) {
      cs        = cs[is_pure]
      purity    = purity[is_pure,]
      row_names = paste0("L",which(include_idx)[is_pure])
      names(cs) = row_names
      rownames(purity) = row_names

      # Re-order CS list and purity rows based on purity.
      ordering = order(purity[,1],decreasing = TRUE)
      return(list(cs       = cs[ordering],
                  purity   = purity[ordering,],
                  cs_index = which(include_idx)[is_pure[ordering]],
                  coverage = claimed_coverage[ordering],
                  requested_coverage=coverage))
    } else
      return(list(cs = NULL,coverage = NULL,requested_coverage = coverage))
  }
}

#' @title Get PIPs
#'
#' @param prune_by_cs Whether or not to ignore single effects not in
#'   a reported CS when calculating PIP.
#'
#' @param prior_tol Filter out effects having estimated prior variance
#'   smaller than this threshold.
#'
#' @export
#'
susie_get_pip = function (res, prune_by_cs = FALSE, prior_tol = 1e-9) {

  if (inherits(res,"susie")) {

    # Drop null weight columns.
    if (!is.null(res$null_index) && res$null_index > 0)
      res$alpha = res$alpha[,-res$null_index,drop=FALSE]

    # Drop the single-effects with estimated prior of zero.
    if (is.numeric(unlist(res$V)))
      include_idx = which(unlist(lapply(res$V, function(x) {x[1,1] > prior_tol & x[2,2] > prior_tol})))
    else
      include_idx = 1:nrow(res$alpha)

    # Only consider variables in reported CS.
    # This is not what we do in the SuSiE paper.
    # So by default prune_by_cs = FALSE means we do not run the
    # following code.
    if (!is.null(res$sets$cs_index) && prune_by_cs)
      include_idx = intersect(include_idx,res$sets$cs_index)
    if (is.null(res$sets$cs_index) && prune_by_cs)
      include_idx = numeric(0)

    # now extract relevant rows from alpha matrix
    if (length(include_idx) > 0)
      res = res$alpha[include_idx,,drop = FALSE]
    else
      res = matrix(0,1,ncol(res$alpha))
  }

  return(as.vector(1 - apply(1 - res,2,prod)))
}

# cov2cor function with specified warning muffled.
#
#' @importFrom stats cov2cor
#' @keywords internal
muffled_cov2cor = function (x)
  withCallingHandlers(cov2cor(x),
    warning = function(w) {
      if (grepl("had 0 or NA entries; non-finite result is doubtful",
                w$message))
          invokeRestart("muffleWarning")
      })

# Find how many variables in the CS.
# x is a probability vector.
#' @keywords internal
n_in_CS_x = function (x, coverage = 0.9)
  sum(cumsum(sort(x,decreasing = TRUE)) < coverage) + 1

# Return binary vector indicating if each point is in CS.
# x is a probability vector.
#' @keywords internal
in_CS_x = function (x, coverage = 0.9) {
  n = n_in_CS_x(x,coverage)
  o = order(x,decreasing = TRUE)
  result = rep(0,length(x))
  result[o[1:n]] = 1
  return(result)
}

# Returns an l-by-p binary matrix indicating which variables are in
# susie credible sets.
#' @keywords internal
in_CS = function (res, coverage = 0.9) {
  if (inherits(res,"susie"))
    res = res$alpha
  return(t(apply(res,1,function(x) in_CS_x(x,coverage))))
}

# @title Utility function to display warning messages as they occur
# @param ... warning message
# @param style either "warning" or "hint"
#'@importFrom crayon combine_styles
warning_message = function(..., style=c("warning", "hint")) {
  style = match.arg(style)
  if (style=="warning" && getOption("warn")>=0) {
    alert <- crayon::combine_styles("bold", "underline", "red")
    message(alert("WARNING:"), " ", ...)
  } else {
    alert <- crayon::combine_styles("bold", "underline", "magenta")
    message(alert("HINT:"), " ", ...)
  }
}
