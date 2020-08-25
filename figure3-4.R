################################################


# read in first flight file
faam_mat <- readMat("c007_dv.mat")
df_mat <- data.frame(Time=faam_mat$Time,Buck=faam_mat$buck,GE=faam_mat$ge,WVSS2F=faam_mat$wvss2f,WVSS2R=faam_mat$wvss2r,amb_press=faam_mat$PS.RVSM,amb_temp=faam_mat$TAT.ND.R)
df_mat_full <- data.frame(time=df_mat$Time,buck=log10(df_mat$Buck),wvss2f=log10(df_mat$WVSS2F),wvss2r=log10(df_mat$WVSS2R),ge=log10(df_mat$GE),press=df_mat$amb_press,temp=df_mat$amb_temp,ind=c(rep(1:10,1478),1:6))

# simple NA imputation

y <- df_mat_full$buck
y_wf <- df_mat_full$wvss2f
y_wr <- df_mat_full$wvss2r
y_ge <- df_mat_full$ge
pr <- df_mat_full$press

for(i in 2:length(y)+1){
  if(is.na(y[i])){
    y[i] <- y[i-1] }
  
  if(is.na(y_wf[i])){
    y_wf[i] <- y_wf[i-1] }
  
  if(is.na(y_wr[i])){
    y_wr[i] <- y_wr[i-1] }
  
  if(is.na(y_ge[i])){
    y_ge[i] <- y_ge[i-1] }
  
}
y <- y[1:14786]
y_wf <- c(y_wf[2],y_wf[2:14786])
y_wr <- c(y_wr[2],y_wr[2:14786])
y_ge <- y_ge[1:14786]
#####

################## figure 3 ####################

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


# DLM inputs

y <- y[1000:2500]
F <- array(1, dim = c(1,    length(y)))
G <- array(1, dim = c(1, 1, length(y)))
V <- c(rep(0.5, length = length(y)))
W <- array(0.009, dim = c(1, 1, length(y)))
m0 <- c(rep(0.5, length = 1))
C0 <- array(0.5, dim = c(1,   1))

# forward filtering output
out <- filter.dlm(y,F,G,V,W,m0,C0)

# plotting
plot(df_mat_full$time[1000:2500],exp(2.303*y),log="y",type="l",ylab="Absolute Humidity",xlab="Minutes from Midnight")
lines(df_mat_full$time[1000:2500],exp(2.303*as.numeric(out$m)),log="y",col="deepskyblue4")
legend("bottomright", legend=c("'True' humidity", "Buck"),
       col=c("deepskyblue4", "black"), lty=1, box.lty=0)



################## figure 4 ####################


# DLM inputs
F <- array(1, dim = c(1,    length(y)))
G <- array(1, dim = c(1, 1, length(y)))
V <- c(rep(0.5, length = length(y)))
W <- array(0.009, dim = c(1, 1, length(y)))
m0 <- c(rep(0.5, length = 1))
C0 <- array(0.5, dim = c(1,   1))

# conditioning on ambient pressure
pr <- pr[1000:2500]

for(i in 1:length(y)){
  if(pr[i] < 465.4) {
    W[,,i] <- 0.00005
  } else{
    W[,,i] <- 0.007
  }
}

# forward filtering output
out <- filter.dlm(y,F,G,V,W,m0,C0)

# plotting
plot(df_mat_full$time[1000:2500],exp(2.303*y),log="y",type="l",ylab="Absolute Humidity",xlab="Minutes from Midnight")
lines(df_mat_full$time[1000:2500],exp(2.303*as.numeric(out$m)),log="y",col="deepskyblue4")
legend("bottomright", legend=c("'True' humidity", "Buck"),
       col=c("deepskyblue4", "black"), lty=1, box.lty=0)





