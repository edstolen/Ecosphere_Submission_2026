## Nimble model

code <- nimbleCode({
  
  beta0 <- 0

  ## Hierarchical habitat effects
  # Population-level means (shared across seasons)
  beta_h_mean[1] ~ dnorm(0, sd = 0.3)
  beta_h_mean[2] ~ dnorm(0, sd = 0.3)
  # Season-specific variation
  beta_h_sd ~ T(dnorm(0, sd = 0.2), 0, )
  # Season-specific effects (deviations from mean)
  for (t in 1:3) {
    beta_h[t,1] ~ dnorm(beta_h_mean[1], sd = beta_h_sd)
    beta_h[t,2] ~ dnorm(beta_h_mean[2], sd = beta_h_sd)
    beta_h[t,3] <- - (beta_h[t,1] + beta_h[t,2])
  }

  ## log_sigma_int: RE at individual level
  log_sigma_int_sd   ~ T(dnorm(0, sd = 0.15), 0, )   # weakly informative - pulled towards 0
  log_sigma_int_mean ~ dnorm(log(15), sd = 0.7) # weakly informative around 15m                  
  for (t in 1:3) {
    for (i in 1:M_aug) {
      log_sigma_int_re[i,t] ~ dnorm(0, sd = log_sigma_int_sd)
      log_sigma_int[i,t]    <- log_sigma_int_mean + log_sigma_int_re[i,t]
    }
  }
  
  ## alpha0.int: time×habitat contrasts
  for (t in 1:3) {
    ## Baseline detection at d=0: weakly informative
    alpha0_base[t] ~ dnorm(0, sd = 1.5)
    ## Habitat-specific contrasts: weakly informative
    alpha0_dev[t,1] ~ dnorm(0, sd = 2)
    alpha0_dev[t,2] ~ dnorm(0, sd = 2)
    alpha0_dev[t,3] <- - (alpha0_dev[t,1] + alpha0_dev[t,2])
    
    for (j in 1:J) {
      alpha0.int[j,t] <- alpha0_base[t] + alpha0_dev[t, trap_hab[j]]
    }
  }
  
  ## Behavior terms
  happy.sigma  ~ dnorm(0, sd = 0.5)
  happy.alpha0 ~ dnorm(0, sd = 0.5)
  
  ## Inclusion indicators
  for(t in 1:3) {
    psi[t] ~ dbeta(2, 2)  # mildly informative around 0.5 (dbeta(1, 1) is flat prior)
    for(i in 1:M_aug) {
      z[i,t] ~ dbern(psi[t])
    }
  }
  
  # Season-indexed normalizer 
  for (t in 1:3) {
    hab_contrib[1:nSuit, t] <- H_s[1:nSuit,1] * beta_h[t,1] +
        H_s[1:nSuit,2] * beta_h[t,2] +
        H_s[1:nSuit,3] * beta_h[t,3]
    lp_s[1:nSuit, t] <- beta0 + hab_contrib[1:nSuit, t]
    w_s[1:nSuit, t]   <- exp(lp_s[1:nSuit, t])
    sumExpEta_t[t]    <- exp(beta0) * sum(w_s[1:nSuit, t])
    logNormConst[t]   <- log(sumExpEta_t[t])
  }
  
  # AC density call (pass gates)
  for (t in 1:3) for (i in 1:M_aug) {
    
    u[i,1:2,t] ~ dAC_inhomog_hab(
      beta_h       = beta_h[t, 1:3],                             # length 3
      beta0        = beta0,                                      # scalar
      Habitat_s    = Habitat_s[1:nSuit],                         # <-- slice
      gridXY_s     = coordsHabitatGridCenter_s[1:nSuit, 1:2],    # <-- slice (nSuit x 2)
      nSuit        = nSuit,
      logNormConst = logNormConst[t],
      # off_s        = off_s[1:nSuit],                             # <-- slice
      cell_half_x  = cell_half_x,
      cell_half_y  = cell_half_y,
      xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax
    )
    
  }
  
  ## Capture model: multinomial over available traps within R_local
  r2 <- R_local * R_local
  
  for(t in 1:3) {
    for(i in 1:M_aug) {
      for(k in 1:n_nights) {
        Ycat[i,k,t] ~ dMultiCapture(
          u = u[i,1:2,t],
          alpha0_int = alpha0.int[1:J,t],
          log_sigma_int = log_sigma_int[i,t],
          happy_alpha0 = happy.alpha0,
          happy_sigma = happy.sigma,
          Cap = Cap.already[i,k,t],
          trapCoords = trapCoords[1:J,1:2],
          trap_avail = trap_avail[1:J,t,k],
          r2 = r2,
          zPresent = z[i,t],
          J = J
        )
      }
    }
  }
})

