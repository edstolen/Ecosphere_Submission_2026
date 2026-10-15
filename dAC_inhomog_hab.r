## ---------- compute candidate AC cell probability ----------
dAC_inhomog_hab <- nimbleFunction(
  run = function(x = double(1),
                 beta_h = double(1),          
                 beta0 = double(0),
                 Habitat_s = double(1),
                 gridXY_s = double(2),
                 nSuit = integer(0),
                 logNormConst = double(0),
                 cell_half_x = double(0),  
                 cell_half_y = double(0),  
                 xmin = double(0), xmax = double(0),
                 ymin = double(0), ymax = double(0),
                 log = integer(0, default = 0)) {
    
    returnType(double(0))
    
    ## Reject cells outside of global bounds
    if((x[1] < xmin) | (x[1] > xmax) |
       (x[2] < ymin) | (x[2] > ymax)) {
      if(log == 1) return(-Inf)
      else return(0.0)
    }
    
    ## Identify nearest suitable cell
    s_local <- 1L
    dx <- x[1] - gridXY_s[1,1]; dy <- x[2] - gridXY_s[1,2]
    dmin <- dx*dx + dy*dy
    
    if(nSuit > 1L) {
      for(s in 2:nSuit) {
        dx <- x[1] - gridXY_s[s,1]; dy <- x[2] - gridXY_s[s,2]
        dtmp <- dx*dx + dy*dy
        if(dtmp < dmin) {
          dmin <- dtmp
          s_local <- s
        }
      }
    }
    
    ## Make sure coordinates are inside the cell rectangle
    if((abs(x[1] - gridXY_s[s_local, 1]) > cell_half_x) |
       (abs(x[2] - gridXY_s[s_local, 2]) > cell_half_y)) {
      
      if(log == 1) return(-Inf)
      else return(0.0)
    }
    
    ## Get habitat class index
    hab <- Habitat_s[s_local]
    h_idx <- 1L
    if(hab == 2.0) h_idx <- 2L
    if(hab == 3.0) h_idx <- 3L
    
    ## Linear predictor 
    lpIdx <- beta0 + beta_h[h_idx]
    
    ## normalize
    logProb <- lpIdx - logNormConst
    
    if(log == 1) return(logProb)
    else return(exp(logProb))
  }
)

registerDistributions(list(
  dAC_inhomog_hab = list(
    BUGSdist = "dAC_inhomog_hab(beta_h, beta0, Habitat_s, gridXY_s, nSuit, logNormConst, cell_half_x, cell_half_y, xmin, xmax, ymin, ymax)",
    types = c(
      "value = double(1)",
      "beta_h = double(1)",
      "beta0 = double(0)",
      "Habitat_s = double(1)",
      "gridXY_s = double(2)",
      "nSuit = integer(0)",
      "logNormConst = double(0)",
      "cell_half_x = double(0)",     
      "cell_half_y = double(0)",     
      "xmin = double(0)", "xmax = double(0)",
      "ymin = double(0)", "ymax = double(0)"
    ),
    discrete = FALSE
  )
))
