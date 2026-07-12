## Table 3: accuracy comparison across spatial confounding strength u.
## DGP as in the paper (theta_U = 0.4); nsim = 1000 in the paper (slow);
## set NSIM lower for a quick check.
suppressMessages(library(spaci)); library(parallel)
NSIM <- as.integer(Sys.getenv("NSIM", "1000")); n <- 250
dir.create("results", showWarnings = FALSE)
cells <- c(2.0, 1.5, 1.0, 0.5)
res <- lapply(cells, function(u) {
  M <- do.call(rbind, mclapply(1:NSIM, function(r) {
    sim <- simulate_spatial_causal(n = n, delta_u = u, theta_spatial = 0.4,
             gamma_interference = 1.5, tau_exp = 0.1, seed = 1e5*u*10 + r)
    ate <- spatial_ate(sim$Y, sim$Z, sim$X, sim$coords, tau = 0.1,
                       caliper = 0.25, seed = r)
    setNames(ate$ATT, ate$Method)
  }, mc.cores = max(1, detectCores() - 1)))
  data.frame(u = u, method = colnames(M),
             mean = colMeans(M, na.rm = TRUE),
             bias = colMeans(M, na.rm = TRUE) - 2,
             mse  = colMeans((M - 2)^2, na.rm = TRUE))
})
out <- do.call(rbind, res); print(out, row.names = FALSE)
saveRDS(out, "results/table3_accuracy.rds")
