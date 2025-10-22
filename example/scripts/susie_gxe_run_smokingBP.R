library(Matrix)
library(mvtnorm)
library(parallel)
suppressMessages(library(optparse))
devtools::install_github("MengyuZhang1307/SuSiEgxe", ref = "main")
library(SuSiEgxe)

#### This is an example for loci 2:26692756:T_G

#### Example 1: Fine-mapping with SuSiEgxe in cross population meta analysis (CPMA)
# input parameters

sumstatdir <- system.file("example/data", "CMA.DBP.CURSMK.COMBINED.2df.2_26673079_26714801.rsid.txt", package = "SuSiEgxe")
LDdir <- system.file("example/data", "CMA.DBP.CURSMK.COMBINED.2df.2_26673079_26714801.LD_mat.txt", package = "SuSiEgxe")
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


#### Example 1: Fine-mapping with SuSiEgxe in Eastern Asian ancestry (EAS)
# input parameters
sumstatdir <- system.file("example/data", "EAS.DBP.CURSMK.COMBINED.2df.2_26681997_26710019.rsid.txt", package = "SuSiEgxe")
LDdir <- system.file("example/data", "EAS.DBP.CURSMK.COMBINED.2df.2_26681997_26710019.LD_mat.txt", package = "SuSiEgxe")
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