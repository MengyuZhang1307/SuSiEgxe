% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/susie_ss_gxe.R
\name{susie_suff_stat_gxe}
\alias{susie_suff_stat_gxe}
\title{susie_gxe}
\usage{
susie_suff_stat_gxe(
  XtX,
  XtZ,
  ZtZ,
  Xty,
  yty,
  n,
  X_colmeans = NA,
  y_mean = NA,
  maf = NULL,
  maf_thresh = 0,
  L = 10,
  scaled_prior_variance = matrix(c(0.5, 0.1, 0.1, 0.5), ncol = 2),
  residual_variance = NULL,
  estimate_residual_variance = TRUE,
  estimate_prior_variance = TRUE,
  estimate_prior_method = c("optim", "EM", "simple"),
  check_null_threshold = 0,
  prior_tol = 1e-09,
  r_tol = 1e-08,
  prior_weights = NULL,
  null_weight = 0,
  standardize = TRUE,
  max_iter = 100,
  s_init = NULL,
  coverage = 0.95,
  min_abs_corr = 0.5,
  tol = 0.001,
  verbose = FALSE,
  track_fit = FALSE,
  check_input = FALSE,
  refine = FALSE,
  check_prior = FALSE,
  n_purity = 100,
  ...
)
}
\arguments{
\item{XtX}{A p by p matrix \eqn{X'X}}

\item{XtZ}{A p by p matrix \eqn{X'Z}}

\item{ZtZ}{A p by p matrix \eqn{Z'Z}}

\item{Xty}{A p-vector \eqn{X'y}}

\item{yty}{A scalar \eqn{y'y} in which y is centered to have mean
zero.}

\item{n}{The sample size.}

\item{X_colmeans}{A p-vector of column means of \code{X}. If both
\code{X_colmeans} and \code{y_mean} are provided, the intercept
is estimated; otherwise, the intercept is NA.}

\item{y_mean}{A scalar containing the mean of \code{y}. If both
\code{X_colmeans} and \code{y_mean} are provided, the intercept
is estimated; otherwise, the intercept is NA.}

\item{maf}{Minor allele frequency; to be used along with
\code{maf_thresh} to filter input summary statistics.}

\item{maf_thresh}{Variants having a minor allele frequency smaller
than this threshold are not used.}

\item{r_tol}{Tolerance level for eigenvalue check of positive
semidefinite matrix of R.}

\item{check_input}{If \code{check_input = TRUE},
\code{susie_suff_stat} performs additional checks on \code{XtX} and
\code{Xty}. The checks are: (1) check that \code{XtX} is positive
semidefinite; (2) check that \code{Xty} is in the space spanned by
the non-zero eigenvectors of \code{XtX}.}

\item{check_prior}{If \code{check_prior = TRUE}, it checks if the
estimated prior variance becomes unreasonably large (comparing with
10 * max(abs(chi-square))^2).}

\item{n_purity}{Passed as argument \code{n_purity} to
\code{\link{susie_get_cs}}.}

\item{...}{Additional arguments to provide backward compatibility
with earlier versions of \code{susie_suff_stat}.}
}
\description{
susie_gxe
}
