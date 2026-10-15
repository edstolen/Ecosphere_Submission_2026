#Build and compile model -----

m <- nimbleModel(
  code       = code,
  constants  = constants,
  data       = data_list,
  inits      = inits(),
  check      = FALSE  # faster build; data & types already aligned
)

Cmodel <- compileNimble(m, showCompilerOutput = T)

conf <- configureMCMC(m, monitors = c(
  "beta_h",
  "beta_h_mean",       
  "beta_h_sd", 
  "happy.alpha0",
  "happy.sigma",
  "psi",
  "log_sigma_int_mean",
  "log_sigma_int_sd",
  "alpha0_base",
  "alpha0_dev"),
  thin      = 1,
  monitors2 = c("u", "z"),
  thin2     = 10,
  enableWAIC = F)

## RW samplers for scalar nodes
for (nm in c("happy.alpha0","happy.sigma", "log_sigma_int_mean")) {
  conf$removeSamplers(nm)
  conf$addSampler(target = nm, type = "RW",
                  control = list(adaptive = TRUE, adaptScaleInterval = 100))
}

## Slice for positive scale
conf$removeSamplers("log_sigma_int_sd")
conf$addSampler(target = "log_sigma_int_sd", type = "slice")

## Slice for psi
for(t in 1:3) {
  conf$removeSamplers(paste0("psi[", t, "]"))
  conf$addSampler(target = paste0("psi[", t, "]"), 
                  type = "slice")
}

# Block sample the means
conf$removeSamplers("beta_h_mean[1]")
conf$removeSamplers("beta_h_mean[2]")
conf$addSampler(
  target = c("beta_h_mean[1]", "beta_h_mean[2]"),
  type = "RW_block",
  control = list(adaptive = TRUE, scale = 0.2)
)

# Slice sampler for SD
conf$removeSamplers("beta_h_sd")
conf$addSampler("beta_h_sd", type = "slice")

# Season-specific deviations (small blocks)
for (t in 1:3) {
  conf$removeSamplers(paste0("beta_h[", t, ",1]"))
  conf$removeSamplers(paste0("beta_h[", t, ",2]"))
  conf$addSampler(
    target = paste0("beta_h[", t, ",1:2]"),
    type = "RW_block",
    control = list(adaptive = TRUE, adaptScaleInterval = 100, scale = 0.15)
  )
}


## Block sampler per season (strongly recommended: parameters are correlated)
for (t in 1:3) {
  conf$removeSamplers(sprintf("alpha0_base[%d]", t))
  conf$removeSamplers(sprintf("alpha0_dev[%d,1]", t))
  conf$removeSamplers(sprintf("alpha0_dev[%d,2]", t))
  conf$addSampler(
    target  = c(sprintf("alpha0_base[%d]", t),
                sprintf("alpha0_dev[%d,1]", t),
                sprintf("alpha0_dev[%d,2]", t)),
    type    = "RW_block",
    control = list(adaptive = TRUE, adaptScaleInterval = 100, scale = 0.25)
  )
}

# Block RW for ACs
for (t in 1:3) {
  for (i in 1:M_aug) {
  nm <- sprintf("u[%d,1:2,%d]", i, t)
  conf$removeSamplers(nm)
  conf$addSampler(target = nm, type = "RW_block",
                  control = list(adaptive = TRUE, adaptScaleInterval = 100, scale = 50))
  }
}

## Build MCMC
Rmcmc <- buildMCMC(conf)
Cmcmc <- compileNimble(Rmcmc, project = m)


## RUN ----
library(posterior)
# reset inits if outliers fail
# Cmodel$setInits(inits())          
# Cmcmc$run(niter, reset = TRUE)

start   <- Sys.time()
Cmcmc$run(300)        # first chunk
end     <- Sys.time()
print(end - start)
# Cmodel$setInits(inits()) # reset inits
# Cmodel$setData(data_list)
# Cmcmc$run(300, reset = TRUE, resetMV = TRUE)  

mvburnin <- as.matrix(Cmcmc$mvSamples)
print(summarize_draws(mvburnin), n=30)
mcmc_trace(mvburnin)

start   <- Sys.time()
Cmcmc$run(1000, reset = FALSE, resetMV = TRUE)        # continue
end     <- Sys.time()
print(end - start)