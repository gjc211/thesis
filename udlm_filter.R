######################################################################################
################## Forward filtering #################################################
######################################################################################


filter.dlm = function(y, F, G, V, W, m0, C0) {
  
  ## Expects the following inputs:
  ## y  - the observations, vector of length tt
  ## F  - the forecast matrix, array of dimension (nn,tt)
  ## G  - the evolution matrix, array of dimension (nn,nn,tt)
  ## V  - the observation error matrix, vector of length tt
  ## W  - the evolution covariance matrix, array of dimension (nn,nn,tt)
  ## m0 - the expectation of the state vector at t = 0, vector of length nn
  ## C0 - the covariance matrix of the state vector at t = 0, matrix (nn,nn)
  
  ## Dimensions
  nn = nrow(F)  ## Length of state vector
  tt = ncol(F)  ## Length of observation time series
  
  ## Create arrays to store results
  a  = array(NA, dim = c(nn,   tt))  ## state prediction mean
  R  = array(NA, dim = c(nn,nn,tt))  ## state prediction scale
  m  = array(NA, dim = c(nn,   tt))  ## state posterior mean
  C  = array(NA, dim = c(nn,nn,tt))  ## state posterior scale
  f  = numeric(tt)                   ## forecast mean
  Q  = numeric(tt)                   ## forecast scale
  
  ## Initialisation
  mt = m0
  Ct = C0
  
  ## Loop over time
  for (t in 1:tt) {
    
    ## Prediction step
    at = G[,,t] %*% mt
    Rt = G[,,t] %*% Ct %*% t(G[,,t]) + W[,,t]
    ft = t(F[,t]) %*% at
    Qt = t(F[,t]) %*% Rt %*% F[,t] + V[t]
    
    ## Update step   
    # if observation is missing, posterior mean becomes prediction mean
    # 
    if (is.na(y[t])) {
      mt = at
      Ct = Rt
    } else {
      A  = Rt %*% F[,t] / Qt[1,1]
      e  = y[t] - ft
      mt = at + A %*% e[1,1]
      Ct = Rt - A %*% t(A) * Q[1]
    }
    
    ## Store results
    m[, t] = mt  # state posterior mean
    C[,,t] = Ct  # state posterior scale
    a[, t] = at  # state prediction mean
    R[,,t] = Rt  # state prediciton scale
    f[  t] = ft  # forecast mean 
    Q[  t] = Qt  # forecast scale
    
  } ## t
  
  ## Return results
  list(m = m, C = C, a = a, R = R, f = f, Q = Q,
       y = y, F = F, G = G, V = V, W = W, m0 = m0, C0 = C0)
  
}

