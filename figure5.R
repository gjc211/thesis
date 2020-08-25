#####
load("fig5.RData") # refer to figure8.R
##### 

# forward filtering
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


# DLM inputs with slight over-reading dealt with through F
F <- array(1.03, dim = c(1,    length(y1)))
G <- array(1, dim = c(1, 1, length(y1)))
V <- c(rep(0.5, length = length(y1)))
W <- array(0.009, dim = c(1, 1, length(y1)))
m0 <- c(rep(0.5, length = 1))
C0 <- array(0.5, dim = c(1,   1))

# more over-reading below 0.1 g/m^3
F <- array(1.02, dim = c(1,    length(y1)))
for(i in 1:length(y1)){
  if(exp(2.303*y1_ge[i]) < 0.08){
    F[i] <- 1.06
  }
}


# defining a quadratic function to F to deal with lag
quadratic <- function(a,x){
  return(exp(a*x) - (a*x))
}
quadratic <- quadratic(a=0.5/1500,x=c(-200:1300))
quadratic <- (quadratic-min(quadratic))/(max(quadratic)-min(quadratic)) * (0.95-0.85) + 0.85
F[9000:10500] <- 1/quadratic[1:1501]

# dealing with problems caused by negative values
y2_ge <- y1_ge + abs(min(y1_ge)) + 1
out1_ge <- filter.dlm(y2_ge,F,G,V,W,m0,C0)
out2_ge <- as.numeric(out1_ge$m) - abs(min(y1_ge)) - 1

# plotting
plot(b878$time[8000:14000],exp(2.303*out2_ge[8000:14000]),type="l",log="y",col="darkorange3",ylab="Absolute Humidity",xlab="Minutes from Midnight")
lines(b878$time[8000:14000],exp(2.303*y1_ge[8000:14000]),log="y")
legend("bottomright", legend=c("'True' humidity", "GE"),
       col=c("darkorange3", "black"), lty=1, box.lty=0)

