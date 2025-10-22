library(Matrix)
library(mvtnorm)
library(parallel)
#devtools::install_github("MengyuZhang1307/SuSiEgxe", ref = "main")
library(SuSiEgxe)

#### This is an example for loci 2:26692756:T_G

#### Example 1: Fine-mapping with SuSiEgxe in cross population meta analysis (CPMA)
# input parameters

sumstatdir <- system.file("extdata", "CMA.DBP.CURSMK.COMBINED.2df.2_26673079_26714801.rsid.txt", package = "SuSiEgxe")
LDdir <- system.file("extdata", "CMA.DBP.CURSMK.COMBINED.2df.2_26673079_26714801.LD_mat.txt", package = "SuSiEgxe")
outrds <- "CMA.DBP.CURSMK.COMBINED.2df.2_26673079_26714801.rds"
L <- 5
beta <- TRUE
print(beta)

joined_data_p = data.table::fread(sumstatdir)
ld_matrix = as.matrix(read.table(LDdir, header = TRUE))
joined_data_p = joined_data_p[which(joined_data_p$refsnp_id %in% colnames(ld_matrix)),]
R <- ld_matrix[match(joined_data_p$refsnp_id, colnames(ld_matrix)), 
                       match(joined_data_p$refsnp_id, colnames(ld_matrix))]
summstats = joined_data_p
if(nrow(joined_data_p) != nrow(R)) stop("Dimention does not match")
#L = 5
#beta = TRUE
print(dim(R))
print(dim(joined_data_p))
Sys.time()
if(!beta) {
  prior_var = replicate(L, matrix(c(10, -2, -2, 10), ncol = 2), simplify = FALSE)
 
  system.time({results = susie_rss_gxe(R = R, bhat = summstats$Effect, 
				     bhat_gxe = summstats$IntEffect, shat = summstats$StdErr, 
				     shat_gxe = summstats$IntStdErr, covhat = summstats$IntCov,
				     z_ld_weight = 0, estimate_residual_variance = FALSE,
                                     prior_variance = prior_var,
	                             check_prior = TRUE, L = L, verbose = T, max_iter = 600, min_abs_corr = 0)})
  saveRDS(results, outrds)
} else {
  prior_var = replicate(L, matrix(c(10, -0.02, -0.02, 10), ncol = 2), simplify = FALSE)
  system.time({results = susie_rss_gxe(R = R, bhat = summstats$Effect, 
				     bhat_gxe = summstats$IntEffect, shat = summstats$StdErr, 
				     shat_gxe = summstats$IntStdErr, covhat = summstats$IntCov, 
				     z_ld_weight = 0, estimate_residual_variance = FALSE,
                                     prior_variance = prior_var, estimate_prior_variance = FALSE,
	                             check_prior = TRUE, L = L, verbose = T, max_iter = 600, min_abs_corr = 0)})
  saveRDS(results, outrds)
}
Sys.time()
#WARNING: Providing the sample size (n), or even a rough estimate of n, is highly recommended. Without n, the implicit assumption is n is large (Inf) and the effect sizes are small (close to zero).
#[1] "iter: 1"
#[1] "objective: 794.048683925696"
#[1] "Max pip difference:0.882959103310611"
#[1] "iter: 2"
#[1] "objective: 834.429873639985"
#[1] "Max pip difference:0"

library(susieR)
summary(results)

## Expected output
#Variables in credible sets:
#
# variable variable_prob cs
#      229     0.8829591  1
#       97     0.7876524  1
#        3     0.4704407  1
#      283     0.2925315  1
#      334     0.1986305  1
#      228     0.1464314  1
#       48     0.1344961  1
#      273     0.1192192  1
#      252     0.1155109  1
#
#Credible sets summary:
#
# cs cs_log10bf cs_avg_r2 cs_min_r2                        variable
#  1   69.17326 0.9017696 0.8049858 3,48,97,228,229,252,273,283,334




#### Example 2: Fine-mapping with SuSiEgxe in Eastern Asian ancestry (EAS)
# input parameters
sumstatdir <- system.file("extdata", "EAS.DBP.CURSMK.COMBINED.2df.2_26681997_26710019.rsid.txt", package = "SuSiEgxe")
LDdir <- system.file("extdata", "EAS.DBP.CURSMK.COMBINED.2df.2_26681997_26710019.LD_mat.txt", package = "SuSiEgxe")
outrds <- "EAS.DBP.CURSMK.COMBINED.2df.2_26681997_26710019.rds"
L <- 5
beta <- TRUE
print(beta)

joined_data_p = data.table::fread(sumstatdir)
ld_matrix = as.matrix(read.table(LDdir, header = TRUE))
joined_data_p = joined_data_p[which(joined_data_p$refsnp_id %in% colnames(ld_matrix)),]
R <- ld_matrix[match(joined_data_p$refsnp_id, colnames(ld_matrix)), 
                       match(joined_data_p$refsnp_id, colnames(ld_matrix))]
summstats = joined_data_p
if(nrow(joined_data_p) != nrow(R)) stop("Dimention does not match")
#L = 5
#beta = TRUE
print(dim(R))
print(dim(joined_data_p))
Sys.time()
if(!beta) {
  prior_var = replicate(L, matrix(c(10, -2, -2, 10), ncol = 2), simplify = FALSE)
 
  system.time({results = susie_rss_gxe(R = R, bhat = summstats$Effect, 
				     bhat_gxe = summstats$IntEffect, shat = summstats$StdErr, 
				     shat_gxe = summstats$IntStdErr, covhat = summstats$IntCov,
				     z_ld_weight = 0, estimate_residual_variance = FALSE,
                                     prior_variance = prior_var,
	                             check_prior = TRUE, L = L, verbose = T, max_iter = 600, min_abs_corr = 0)})
  saveRDS(results, outrds)
} else {
  prior_var = replicate(L, matrix(c(10, -0.02, -0.02, 10), ncol = 2), simplify = FALSE)
  system.time({results = susie_rss_gxe(R = R, bhat = summstats$Effect, 
				     bhat_gxe = summstats$IntEffect, shat = summstats$StdErr, 
				     shat_gxe = summstats$IntStdErr, covhat = summstats$IntCov, 
				     z_ld_weight = 0, estimate_residual_variance = FALSE,
                                     prior_variance = prior_var, estimate_prior_variance = FALSE,
	                             check_prior = TRUE, L = L, verbose = T, max_iter = 600, min_abs_corr = 0)})
  saveRDS(results, outrds)
}
Sys.time()
#WARNING: Providing the sample size (n), or even a rough estimate of n, is highly recommended. Without n, the implicit assumption is n is large (Inf) and the effect sizes are small (close to zero).
#[1] "iter: 1"
#[1] "objective: 253.657624243725"
#[1] "Max pip difference:0.616429413667341"
#[1] "iter: 2"
#[1] "objective: 255.297819339268"
#[1] "Max pip difference:0"
library(susieR)
summary(results)

### Expected output
#Variables in credible sets:
#
# variable variable_prob cs
#        6     0.6164294  1
#       31     0.6049876  1
#       28     0.5527195  1
#        2     0.3049689  1
#       32     0.2856208  1
#       33     0.2688363  1
#        4     0.2191402  1
#       39     0.2107923  1
#       42     0.2051069  1
#       48     0.1924634  1
#       67     0.1700362  1
#       24     0.1396807  1
#       65     0.1058730  1
#
#Credible sets summary:
#
# cs cs_log10bf cs_avg_r2 cs_min_r2                            variable
#  1   22.23549  0.964939 0.9146111 2,4,6,24,28,31,32,33,39,42,48,65,67