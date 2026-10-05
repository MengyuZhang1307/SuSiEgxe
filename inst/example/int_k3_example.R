suppressPackageStartupMessages({
  library(SuSiEgxe)   # SuSiEgxe >= 0.2.0
  library(susieR)
  library(ggplot2)
})

## ---- Data -------------------------------------------------------------------
f <- system.file("example", "int_k3_example.rds", package = "SuSiEgxe") # three causal variants with an interaction effect only in the exposed group
dat <- readRDS(f)
ss  <- dat$sumstats   # summary statistics (meta-analysis of the 4 cohorts)
R   <- dat$R          # in-sample LD: signed correlation, same SNP order as `ss`

print(dat$description)
print(dat$causal)     # the three true causal SNPs (bG = 0: no effect in the unexposed)
print(head(ss[, c("snp", "bhat", "shat", "bhat_gxe", "shat_gxe", "covhat", "bhat_marg", "shat_marg")], 3))

## ---- Checks -----------------------------------------------------------------
cat("Range of off-diagonal LD (signed):", range(R[upper.tri(R)]), "\n")

## ---- SuSiEgxe ---------------------------------------------------------------
L <- 5
prior <- replicate(L, matrix(c(10, -0.02, -0.02, 10), 2, 2), simplify = FALSE)

fit_gxe <- susie_rss_gxe(
  R = R,
  bhat = ss$bhat, bhat_gxe = ss$bhat_gxe,
  shat = ss$shat, shat_gxe = ss$shat_gxe, covhat = ss$covhat,
  estimate_residual_variance = FALSE,
  prior_variance = prior, estimate_prior_variance = FALSE,
  check_prior = TRUE, L = L, max_iter = 600,
  coverage = 0.95, min_abs_corr = 0.5)

cat("\n== SuSiEgxe ==\n")
print(data.frame(converged = fit_gxe$converged, iterations = fit_gxe$niter,
                 credible_sets = length(fit_gxe$sets$cs)))

cs_table <- function(fit, snp = ss$snp, causal = dat$causal$snp) {
  cs <- fit$sets$cs
  if (length(cs) == 0) return("no credible set")
  data.frame(
    CS = paste0("CS", seq_along(cs)), single_effect = names(cs), size = lengths(cs),
    min_abs_r = sapply(cs, function(v) round(min(abs(R[v, v])), 2)),
    SNPs = sapply(cs, function(v) paste(head(snp[v], 6), collapse = ", ")),
    causal_inside = sapply(cs, function(v) {
      h <- intersect(snp[v], causal)
      if (length(h)) paste(h, collapse = ", ") else "none"
    }),
    row.names = NULL)
}
print(cs_table(fit_gxe))

## ---- Marginal SuSiE ---------------------------------------------------------
fit_marg <- susie_rss(
  R = R, n = max(ss$N),
  bhat = ss$bhat_marg, shat = ss$shat_marg,
  estimate_residual_variance = TRUE, check_prior = TRUE,
  L = L, max_iter = 500, coverage = 0.95, min_abs_corr = 0.5)

cat("\n== Marginal SuSiE ==\n")
print(data.frame(converged = fit_marg$converged, iterations = fit_marg$niter,
                 credible_sets = length(fit_marg$sets$cs)))
print(cs_table(fit_marg))

## ---- The three causal SNPs --------------------------------------------------
idx <- dat$causal$position
cat("\n== Causal SNPs ==\n")
print(data.frame(
  snp = dat$causal$snp,
  chisq_2df = round(ss$chisq2df[idx], 1),
  z_marginal = round((ss$bhat_marg / ss$shat_marg)[idx], 2),
  z_interaction = round((ss$bhat_gxe / ss$shat_gxe)[idx], 2),
  SuSiEgxe_CS = sapply(idx, function(j) {
    k <- which(sapply(fit_gxe$sets$cs, function(v) j %in% v))
    if (length(k)) paste0("CS", k[1]) else "-"
  }),
  SuSiEgxe_PIP = round(fit_gxe$pip[idx], 3),
  SuSiEgxe_PIP_in_CS = round(susie_get_pip(fit_gxe, prune_by_cs = TRUE)[idx], 3),
  SuSiE_PIP = round(fit_marg$pip[idx], 3)))

cat("\n== PIP sums ==\n")
print(c(reported = sum(fit_gxe$pip),
        only_effects_with_a_credible_set = sum(susie_get_pip(fit_gxe, prune_by_cs = TRUE))))

## ---- Locus plot -------------------------------------------------------------
label <- c(gxe = "SuSiEgxe (fixed prior; 2 df test)", marg = "SuSiE (marginal test)")
pos <- ss$position
stat <- rbind(
  data.frame(method = label[["gxe"]],  pos = pos,
             logp = -pchisq(ss$chisq2df, 2, lower.tail = FALSE, log.p = TRUE) / log(10)),
  data.frame(method = label[["marg"]], pos = pos,
             logp = -(pnorm(-abs(ss$bhat_marg / ss$shat_marg), log.p = TRUE) + log(2)) / log(10)))
stat$cs <- "none"
for (m in names(label)) {
  cs <- list(gxe = fit_gxe, marg = fit_marg)[[m]]$sets$cs
  for (i in seq_along(cs)) stat$cs[stat$method == label[[m]] & stat$pos %in% cs[[i]]] <- paste0("CS", i)
}
stat$method <- factor(stat$method, levels = unname(label))
cs_levels <- c("none", paste0("CS", 1:6)); stat$cs <- factor(stat$cs, levels = cs_levels)
cols <- c(none = "grey70", CS1 = "#D55E00", CS2 = "#0072B2", CS3 = "#009E73",
          CS4 = "#CC79A7", CS5 = "#E69F00", CS6 = "#56B4E9")
causal <- do.call(rbind, lapply(unname(label), function(l)
  data.frame(method = l, pos = idx, logp = stat$logp[stat$method == l][idx])))
causal$method <- factor(causal$method, levels = unname(label))

p <- ggplot(stat, aes(pos, logp)) +
  geom_vline(xintercept = seq(50, 250, by = 50) + 0.5, colour = "grey85", linewidth = 0.3) +
  geom_hline(yintercept = -log10(5e-8), colour = "red", linetype = "dashed", linewidth = 0.4) +
  geom_point(data = causal, shape = 21, size = 5, stroke = 0.9, colour = "black", fill = NA) +
  geom_point(aes(colour = cs), size = 1.8) +
  scale_colour_manual(values = cols, breaks = cs_levels[-1], drop = TRUE, name = "credible set") +
  facet_wrap(~ method, ncol = 1) +
  labs(x = "SNP position in the locus (grey lines: LD blocks)", y = expression(-log[10](p)),
       title = "Three causal SNPs with an interaction effect only in the exposed group",
       subtitle = "Black circles: causal SNPs. Red line: p = 5e-8.") +
  theme_bw(base_size = 12) + theme(strip.text = element_text(face = "bold"))

out_png <- file.path(getwd(), "int_k3_locus.png")
ggsave(out_png, p, width = 9, height = 6.5, dpi = 150)

