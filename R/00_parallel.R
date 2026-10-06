## Parallel helper shared by the simulation scripts (01-06).
## parallel::mclapply forks, which Windows cannot do, so the scripts use a
## socket cluster on every platform instead. The workers are fresh R sessions:
## they are given the calling session's RNG kind (scripts 02-06 set
## "L'Ecuyer-CMRG", and every set.seed() inside a replicate depends on it), the
## spaci package and the script's global objects. Each replicate sets its own
## seed, so results do not depend on the number of cores.
## Set the environment variable CORES to choose the number of workers.
par_lapply <- function(X, FUN, cores = as.integer(Sys.getenv("CORES", NA))) {
  if (is.na(cores)) cores <- max(1L, parallel::detectCores() - 1L, na.rm = TRUE)
  cores <- min(cores, length(X))
  if (cores < 2) return(lapply(X, FUN))
  cl <- parallel::makeCluster(cores); on.exit(parallel::stopCluster(cl))
  parallel::clusterCall(cl, function(kind, lib) {
    .libPaths(lib); do.call(RNGkind, as.list(kind))
    suppressMessages(library(spaci)); NULL
  }, RNGkind(), .libPaths())
  parallel::clusterExport(cl, ls(.GlobalEnv), envir = .GlobalEnv)
  parallel::parLapply(cl, X, FUN)
}
