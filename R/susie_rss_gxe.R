#' @title Sum of Single Effects (SuSiE) Regression using Summary Statistics from Gene by Environment (GxE) Interaction
#'
#' @description \code{susie_rss_gxe} performs variable selection under a
#'   sparse Bayesian multiple linear regression of \eqn{Y} on \eqn{X}
#'   using the effect sizes estimates and SE from standard univariate regression
#'   of \eqn{Y} on each column of \eqn{X}, an estimate, \eqn{R}, of
#'   the correlation matrix for the columns of \eqn{X} (\dQuote{LD matrix}). See
#'   \dQuote{Details} for other ways to call \code{susie_rss}
#'
#' @details With the inputs \code{bhat}, \code{bhat_gxe}, \code{shat}, \code{shat_gxe}, and \code{R},
#'\code{susie_rss_gxe} calls \code{susie_suff_stat_gxe} with constructed \code{XtX}, \code{XtZ}, \code{ZtZ}
#' and \code{Xty = dS_inv$Kty}, and with \code{residual_variance = 1}. The
#' underlying assumption of performing the analysis in this way is
#' that the sample size is large (\emph{i.e.}, infinity), and/or the
#' effects are small. 
#'
#' The \code{estimate_residual_variance} setting is \code{FALSE} by
#' default, which is recommended when the LD matrix is estimated from
#' a reference panel.
#'
#' @param z p-vector of z-scores. Not supported in GxE.
#'
#' @param R p x p correlation matrix.
#'
#' @param n The sample size.
#'
#' @param bhat Alternative summary data giving the estimated effects
#'   (a vector of length p). This, together with \code{shat},\code{shat_gxe},\code{bhat_gxe}, may be
#'   provided instead of \code{z}.
#'   
#' @param bhat_gxe Alternative summary data giving the estimated effects
#'   (a vector of length p). This, together with \code{bhat}, \code{shat},\code{shat_gxe}, \code{covhat}, may be
#'   provided instead of \code{z}.
#'
#' @param shat Alternative summary data giving the standard errors of
#'   the estimated effects (a vector of length p). This, together with
#'   \code{bhat}, \code{bhat_gxe},\code{shat_gxe}, \code{covhat}, may be provided instead of \code{z}.
#'   
#' @param shat_gxe Alternative summary data giving the standard errors of
#'   the estimated effects (a vector of length p). This, together with
#'   \code{bhat}, \code{shat},\code{shat_gxe}, \code{covhat}, may be provided instead of \code{z}.
#'
#' @param covhat Alternative summary data giving the covariance of the estimated effects
#'   (a vector of length p). This, together with
#'   \code{bhat} and \code{shat}, \code{bhat_gxe},\code{shat_gxe}, may be provided instead of \code{z}.
#'
#' @param var_y The sample variance of y, defined as \eqn{y'y/(n-1)}.
#'   When the sample variance is not provided, the coefficients
#'   (returned from \code{coef}) are computed on the
#'   \dQuote{standardized} X, y scale. Not supported in GxE.
#'
#' @param z_ld_weight This parameter is included for backwards
#'   compatibility with previous versions of the function, but it is no
#'   longer recommended to set this to a non-zero value. When
#'   \code{z_ld_weight > 0}, the matrix \code{R} is adjusted to be
#'   \code{cov2cor((1-w)*R + w*tcrossprod(z))}, where \code{w =
#'   z_ld_weight}.
#'
#' @param prior_variance The prior variance(s) for the non-zero
#'   noncentrality parameterss \eqn{\tilde{b}_l}. It is either a 2x2 matrix (prior variance of gene main effect, gene by E and prior covariance of gene and gene by E interaction),
#'   or a list of length L and each element is a 2x2 matrix . When the \code{susie_suff_stat_gxe} option
#'   \code{estimate_prior_variance} is set to \code{TRUE} this simply provides an initial value for the
#'   prior variance.
#'
#' @param estimate_residual_variance The default is FALSE, the
#'   residual variance is fixed to 1 or variance of y. 
#'
#' @param check_prior When \code{check_prior = TRUE}, it checks if the
#'   estimated prior variance becomes unreasonably large (comparing with
#'   100 * max(abs(chi-square))^2).
#'
#' @param \dots Other parameters to be passed to
#' \code{\link{susie_suff_stat_gxe}}.
#'
#' @return A \code{"susie"} object with the following
#'   elements:
#'
#' \item{alpha}{An L by p matrix of posterior inclusion probabilites.}
#'
#' \item{mu}{An L by 2p matrix of posterior means, conditional on
#'   inclusion.}
#'
#' \item{mu2}{A list of L 2p by 2p matrices of posterior second moments,
#'   conditional on inclusion.}
#'
#' \item{lbf}{log-Bayes Factor for each single effect.}
#'
#' \item{lbf_variable}{log-Bayes Factor for each variable and single effect.}
#'
#' \item{V}{Prior variance of the non-zero elements of b.}
#'
#' \item{elbo}{The value of the variational lower bound, or
#'   \dQuote{ELBO} (objective function to be maximized), achieved at
#'   each iteration of the IBSS fitting procedure.}
#'
#' \item{sets}{Credible sets estimated from model fit; see
#'   \code{\link{susie_get_cs}} for details.}
#'
#' \item{pip}{A vector of length p giving the (marginal) posterior
#'   inclusion probabilities for all p covariates.}
#'
#' \item{niter}{Number of IBSS iterations that were performed.}
#'
#' \item{converged}{\code{TRUE} or \code{FALSE} indicating whether
#'   the IBSS converged to a solution within the chosen tolerance
#'   level.}
#'
#' @references
#' G. Wang, A. Sarkar, P. Carbonetto and M. Stephens (2020). A simple
#' new approach to variable selection in regression, with application
#' to genetic fine-mapping. \emph{Journal of the Royal Statistical
#' Society, Series B} \bold{82}, 1273-1300 \doi{10.1101/501114}.
#'
#' Y. Zou, P. Carbonetto, G. Wang, G and M. Stephens
#' (2022). Fine-mapping from summary data with the \dQuote{Sum of
#' Single Effects} model. \emph{PLoS Genetics} \bold{18},
#' e1010299. \doi{10.1371/journal.pgen.1010299}.

#' @export
#'


susie_rss_gxe = function (z, R, n, bhat, bhat_gxe, shat, shat_gxe, covhat = NULL, var_y,
                      z_ld_weight = 0,
                      estimate_residual_variance = FALSE,
                      prior_variance = matrix(c(50, 1, 1, 50), ncol = 2),
                      check_prior = TRUE, ...) {

  if (estimate_residual_variance)
    warning_message("For estimate_residual_variance = TRUE, please check ",
                    "that R is the \"in-sample\" LD matrix; that is, the ",
                    "correlation matrix obtained using the exact same data ",
                    "matrix X that was used for the other summary ",
                    "statistics. Also note, when covariates are included in ",
                    "the univariate regressions that produced the summary ",
                    "statistics, also consider removing these effects from ",
                    "X before computing R.",style = "hint")

  # Check input R.
  if (missing(z))
    {p <- length(bhat)
  } else {
    p <- length(z)}
  if (nrow(R) != p)
      stop(paste0("The dimension of R (",nrow(R)," x ",ncol(R),") does not ",
                  "agree with expected (",p," x ",p,")"))

  # Check input n.
  if (!missing(n))
    if (n <= 1)
      stop("n must be greater than 1")

  # Check inputs z, bhat and shat. Note that bhat is no longer used
  # after this step.
  if (sum(c(missing(z),missing(bhat) || missing(shat))) != 1)
    stop("Please provide either z or (bhat, shat, covhat), but not both")
  if (missing(z)) {
    if (length(shat) == 1)
      shat = rep(shat,length(bhat))
    if (length(bhat) != length(shat))
      stop("The lengths of bhat and shat do not agree")
    if(length(bhat) != length(covhat))
      stop("The lengths of bhat and covhat do not agree")
    if (length(bhat_gxe) != length(shat_gxe))
      stop("The lengths of bhat_gxe and shat_gxe do not agree")
    if (anyNA(bhat) || anyNA(shat) || anyNA(covhat) || anyNA(bhat_gxe) || anyNA(shat_gxe))
      stop("bhat, shat, covhat, bhat_gxe, shat_gxe cannot have missing values")
    if (any(shat <= 0) || any(shat_gxe <= 0))
      stop("shat and shat_gxe cannot have zero or negative elements")
  }
  #if (length(z) < 1 & is.null(covhat))
  #  stop("Input vector z should have at least one element")
  #z[is.na(z)] = 0

  # When n is provided, compute the PVE-adjusted z-scores.
  #if (!missing(n) & is.null(covhat)) {
  #  adj = (n-1)/(z^2 + n - 2)
  #  z   = sqrt(adj) * z
  #}

  # Modify R by z_ld_weight; this modification was designed to ensure
  # the column space of R contained z, but susie_suff_stat does not
  # require this, and is no longer recommended.
  #if (z_ld_weight > 0) {
  #  warning_message("As of version 0.11.0, use of non-zero z_ld_weight is no longer ",
  #          "recommended")
  #  R = muffled_cov2cor((1-z_ld_weight)*R + z_ld_weight*tcrossprod(z))
  #  R = (R + t(R))/2
  #}

  dS_inv = S_inverse_crossprod(shat^2, covhat, shat_gxe^2, c(bhat, bhat_gxe)) #block diag(S_inv), M_inv, and crossprod

  XtX = sqrt(tcrossprod(dS_inv$dXtX)) * R
  XtZ = sqrt(tcrossprod(dS_inv$dXtZ)) * R
  ZtZ = sqrt(tcrossprod(dS_inv$dZtZ)) * R

  

  # Call susie_suff_stat. We call susie_suff_stat in two different
  # ways depending on whether n is provided.
  if (missing(n)) { # do not need N?

    # The sample size (n) is not provided, so use unadjusted z-scores.
    # The choice of n=2, yty=1 is mostly arbitrary except in that it
    # ensures var(y) = yty/(n-1) = 1, and because of this
    # scaled_prior_variance = prior_variance.
    warning_message("Providing the sample size (n), or even a rough estimate of n, ",
            "is highly recommended. Without n, the implicit assumption is ",
            "n is large (Inf) and the effect sizes are small (close to zero).")
    s = susie_suff_stat_gxe(XtX = XtX, XtZ = XtZ, ZtZ = ZtZ, Xty = dS_inv$Kty,n = 2,yty = 1,
                        scaled_prior_variance = prior_variance,
                        estimate_residual_variance = estimate_residual_variance,
                        standardize = FALSE,check_prior = check_prior,...)
  } else { 
    stop("n is not used")
    # The sample size (n) is provided, so use PVE-adjusted z-scores.
    if (!missing(shat) & !missing(var_y)) {

      # var_y, shat (and bhat) are provided, so the effects are on the
      # *original scale*.
      XtXdiag = var_y * adj/(shat^2)
      XtX = t(R * sqrt(XtXdiag)) * sqrt(XtXdiag)
      XtX = (XtX + t(XtX))/2
      Xty = z * sqrt(adj) * var_y / shat
    } else {

      # The effects are on the *standardized* X, y scale.
      XtX = (n-1)*R
      Xty = sqrt(n-1)*z
      var_y = 1
    }
    s = susie_suff_stat(XtX = XtX,Xty = Xty,n = n,yty = (n-1)*var_y,
                        estimate_residual_variance = estimate_residual_variance,
                        check_prior = check_prior,...)
  }
  return(s)
}