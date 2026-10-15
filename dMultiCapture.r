## ---------- Capture density ----------
dMultiCapture <- nimbleFunction(
  run = function(x = integer(0),
                 u = double(1),
                 alpha0_int = double(1),     
                 log_sigma_int = double(0),
                 happy_alpha0 = double(0),
                 happy_sigma = double(0),
                 Cap = double(0),
                 trapCoords = double(2),     # J x 2
                 trap_avail = double(1),     # length J for (t,k)
                 r2 = double(0),             # R_local^2
                 zPresent = integer(0),      # 0/1 scalar
                 J = integer(0),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    
    ## Make sure index is in Bounds
    if((x < 1L) | (x > (J + 1L))) {
      if(log == 1) return(-Inf) else return(0.0)
    }
    
    ## No-capture case (exits function due to return)
    if(zPresent == 0L) {
      if(x == (J + 1L)) {
        if(log == 1) return(0.0) else return(1.0)
      }
      if(log == 1) return(-Inf) else return(0.0)
    }
    
    ## Behavior-adjusted parameters
    log_sigma <- log_sigma_int + happy_sigma * Cap
    sigma <- exp(log_sigma)
    sigma2 <- sigma * sigma
    
    ## Weights
    denom <- 1.0    # +1 for no-capture
    w_x <- 0.0
    
    for(j in 1:J) {
      if(trap_avail[j] > 0.5) {
        dx <- u[1] - trapCoords[j,1]
        dy <- u[2] - trapCoords[j,2]
        d2 <- dx*dx + dy*dy
        if(d2 <= r2) {
          log_alpha0 <- alpha0_int[j] + happy_alpha0 * Cap
          alpha0 <- 1.0 / (1.0 + exp(-log_alpha0))
          wj <- alpha0 * exp(-d2 / (2.0 * sigma2))
          denom <- denom + wj
          if(x == j) w_x <- wj
        }
      }
    }
    
    if(x == (J + 1L)) {  # no-capture
      if(log == 1) return(-log(denom)) else return(1.0 / denom)
    } else {
      if(w_x == 0.0) {
        if(log == 1) return(-Inf) else return(0.0)
      }
      if(log == 1) return(log(w_x) - log(denom)) else return(w_x / denom)
    }
  }
)

registerDistributions(list(
  dSoftmaxCapture = list(
    BUGSdist = "dSoftmaxCapture(u,alpha0_int,log_sigma_int,happy_alpha0,happy_sigma,Cap,trapCoords,trap_avail,r2,zPresent,J)",
    types = c(
      "value = integer(0)",
      "u = double(1)",
      "alpha0_int = double(1)",     
      "log_sigma_int = double(0)",
      "happy_alpha0 = double(0)",
      "happy_sigma = double(0)",
      "Cap = double(0)",
      "trapCoords = double(2)",    
      "trap_avail = double(1)",    
      "r2 = double(0)",             
      "zPresent = integer(0)",     
      "J = integer(0)"
    ),
    discrete = TRUE
  ))
)
