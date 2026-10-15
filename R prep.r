### R-side setup

## Dimensions and time settings
M_aug     <- N.total                     # augmented and observed individuals
n_nights  <- 3                           # capture occasions
J         <- nrow(trapCoords)            # number of traps
S         <- nrow(coordsHabitatGridCenter) # number of cells in mask

## Habitat and grid variables
isSuitable01 <- as.numeric(Cells.Habitat.extended$Suitable)
sSuit        <- which(isSuitable01 > 0.5)
nSuit        <- length(sSuit)
gridXY    <- as.matrix(coordsHabitatGridCenter[,1:2])
gridXY_s  <- as.matrix(coordsHabitatGridCenter[sSuit, 1:2])        
Habitat.integer   <- Cells.Habitat.extended$Habitat + 1 # must be >0
Habitat_int <- as.integer(Habitat.integer)
Habitat_s <- as.numeric(Habitat_int[sSuit]) # use only cells in suitable habitat                        

## Bounds and radius
xmin <- min(coordsHabitatGridCenter[,1]); xmax <- max(coordsHabitatGridCenter[,1])
ymin <- min(coordsHabitatGridCenter[,2]); ymax <- max(coordsHabitatGridCenter[,2])
R_local_m <- 1500L                       # distance beyond which traps are ignored

## Ensure non-scalar arguments to custom dists are 'double'
storage.mode(gridXY_s) <- "double"
storage.mode(Habitat_s)                 <- "double"
storage.mode(trapCoords)                <- "double"
storage.mode(trap.avail)                <- "double"

## Ensure trap habitat codes are integer 1..3
trapHabIdx <- as.integer(trapHabIdx)
if (!all(trapHabIdx %in% c(1L,2L,3L))) {
  stop("trapHabIdx must be in {1,2,3}.")
}

## Habitat matrix
H_s <- cbind(
  as.numeric(Habitat_s == 1),
  as.numeric(Habitat_s == 2),
  as.numeric(Habitat_s == 3)
)
storage.mode(H_s) <- "double"

dx <- diff(sort(unique(gridXY_s[,1])))
dy <- diff(sort(unique(gridXY_s[,2])))
cell_half_x <- if (length(dx) > 0) min(dx)/2 else 0
cell_half_y <- if (length(dy) > 0) min(dy)/2 else 0


## inits function
inits <- (function() {
  
  #Activity centers (snap to suitable cells)
  u <- array(NA_real_, dim = c(M_aug, 2, 3))
  for(t in 1:3) {
    N_t <- nrow(AC.ls[[t]])
    add <- M_aug - N_t
    
    pick  <- sample(sSuit, size = add, replace = TRUE)
    u_aug <- gridXY[pick, , drop = FALSE]
    
    AC_fixed <- AC.ls[[t]]
    for(i in 1:nrow(AC_fixed)) {
      s_local <- 1L
      dx <- AC_fixed[i,1] - gridXY_s[1,1]
      dy <- AC_fixed[i,2] - gridXY_s[1,2]
      dmin <- dx*dx + dy*dy
      if(nSuit > 1L) {
        for(s in 2:nSuit) {
          dx <- AC_fixed[i,1] - gridXY_s[s,1]
          dy <- AC_fixed[i,2] - gridXY_s[s,2]
          dtmp <- dx*dx + dy*dy
          if(dtmp < dmin) {
            dmin <- dtmp
            s_local <- s
          }
        }
      }
      AC_fixed[i,] <- gridXY_s[s_local,]
    }
    u[,1:2,t] <- rbind(AC_fixed, u_aug)
  }
  
  # Other priors 
  happy.alpha0     <- rnorm(1, 0, 1)
  happy.sigma      <- rnorm(1, 0, 1)
  beta_h_mean <- rnorm(2, 0, 0.2)  # Small initial values
  beta_h_sd <- abs(rnorm(1, 0.15, 0.05))
  beta_h <- matrix(NA_real_, nrow = 3, ncol = 3)
  for(t in 1:3) {
    beta_h[t, 1] <- rnorm(1, beta_h_mean[1], beta_h_sd)
    beta_h[t, 2] <- rnorm(1, beta_h_mean[2], beta_h_sd)
    beta_h[t, 3] <- -(beta_h[t, 1] + beta_h[t, 2])
  }
  
  psi <- runif(3, 0.1, 0.9)
  N1 <- nrow(AC.ls[[1]]); N2 <- nrow(AC.ls[[2]]); N3 <- nrow(AC.ls[[3]])
  z <- cbind(
    c(rep(1L, N1), rbinom(M_aug - N1, 1, 0.3)),
    c(rep(1L, N2), rbinom(M_aug - N2, 1, 0.3)),
    c(rep(1L, N3), rbinom(M_aug - N3, 1, 0.3))
  )
 
  # Movement random effects (centered; shared mean/sd across seasons)
  log_sigma_int_mean <- rnorm(1, mean = log(15), sd = 0.10)                 # ≈ 2.708
  log_sigma_int_sd   <- abs(rnorm(1, mean = 0.30,  sd = 0.10)) + 1e-3       # positive
  # Individual×season random effects: small spread at start
  log_sigma_int_re <- matrix(rnorm(M_aug * 3, mean = 0, sd = 0.20),
                             nrow = M_aug, ncol = 3)
  # Deterministic composition: log_sigma_int = mean + re
  log_sigma_int <- log_sigma_int_mean + log_sigma_int_re   # broadcasts over matrix
  
  # assignment for alpha0.int
  alpha0_base <- rnorm(3, 0, 1.5)
  alpha0_dev  <- matrix(NA_real_, 3, 3)
  alpha0_dev[,1] <- rnorm(3, 0, 2)
  alpha0_dev[,2] <- rnorm(3, 0, 2)
  alpha0_dev[,3] <- - (alpha0_dev[,1] + alpha0_dev[,2])
  
  
  # Return 
  list(
    happy.alpha0     = happy.alpha0,
    happy.sigma      = happy.sigma,
    beta_h_mean = beta_h_mean,
    beta_h_sd   = beta_h_sd,
    beta_h      = beta_h,
    
    psi              = psi,
    z                = z,
    u                = u,
    
    log_sigma_int_mean = log_sigma_int_mean,  # scalar shared across seasons
    log_sigma_int_sd   = log_sigma_int_sd,    # scalar shared across seasons
    log_sigma_int_re   = log_sigma_int_re,    # M_aug × 3
    log_sigma_int      = log_sigma_int,       # M_aug × 3 (deterministic OK to pass)
    
    alpha0_base      = alpha0_base,
    alpha0_dev       = alpha0_dev
    
  )
})

## Constants and data
 
constants <- list(
  J          = J,
  M_aug      = M_aug,
  n_nights   = n_nights,
  R_local    = as.integer(R_local_m),
  xmin       = as.numeric(xmin), 
  xmax = as.numeric(xmax),
  ymin = as.numeric(ymin), 
  ymax = as.numeric(ymax),
  nSuit      = as.integer(nSuit),
  trap_hab   = trapHabIdx,
  cell_half_x = as.numeric(cell_half_x),
  cell_half_y = as.numeric(cell_half_y)
)

data_list <- list(
  coordsHabitatGridCenter_s = gridXY_s,   # double(2)
  Habitat_s   = Habitat_s,  # double(1)
  H_s  = H_s,
  trapCoords  = trapCoords,
  trap_avail  = trap.avail,
  Cap.already = Cap.array.collapsed.shifted.augment,
  Ycat        = Ycat.aug
)
