###################################################################
################### figure 9 ######################################
###################################################################

#####
load("fig8.RData") # refer to figure1.R
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


# DLM inputs
F <- array(1, dim = c(1,    length(y1)))
G <- array(1, dim = c(1, 1, length(y1)))
V <- c(rep(0.5, length = length(y1)))
W <- array(0.05, dim = c(1, 1, length(y1)))
m0 <- c(rep(0.5, length = 1))
C0 <- array(0.5, dim = c(1,   1))

##### 


#### buck (W below) ####

for(i in 1:length(y1)){
  if(exp(y1[i]*2.303) < 0.7){
    V[i] <- 1.5
  }
}

###
for(i in 1:length(y1)){
  if(exp(y1[i]*2.303) > 0.03){
    W[,,i] <- 0.045 
  }
  if(exp(y1[i]*2.303) > 0.7){
    W[,,i] <- 0.005 
  }
}

W[,,c(5600:5900,2500:2800)] <- 0.045
out1 <- filter.dlm(y1,F,G,V,W,m0,C0)



#### flush ####

W <- array(0.05, dim = c(1, 1, length(y1)))
for(i in 1:length(y1)){
  if(exp(y1[i]*2.303) > 0.03){
    W[,,i] <- 0.045 
  }
  if(exp(y1[i]*2.303) > 0.7){
    W[,,i] <- 0.005 
  }
}
W[,,c(5600:5900,2500:2800)] <- 0.045


V <- c(rep(0.5, length = length(y1)))
for(i in 1:length(y1)){
  if(exp(y1[i]*2.303) < 0.7){
    V[i] <- 0.6
  }
  if(exp(y1[i]*2.303) < 0.008){
    V[i] <- 0.7
  }
}

out1_wf <- filter.dlm(y1_wf,F,G,V,W,m0,C0)


#### rosemount #### 

W <- array(0.05, dim = c(1, 1, length(y1)))

V <- c(rep(0.5, length = length(y1)))
for(i in 1:length(y1)){
  if(exp(y1[i]*2.303) < 0.7){
    V[i] <- 0.6
  }
  if(exp(y1[i]*2.303) < 0.03){
    V[i] <- 0.7
  }
}

out1_wr <- filter.dlm(y1_wr,F,G,V,W,m0,C0)



#### ge ####

quadratic <- function(a,x){
  return(exp(a*x) - (a*x) - 1)
}
quadratic <- quadratic(a=0.0009,x=c(-200:900))
quadratic <- (quadratic-min(quadratic))/(max(quadratic)-min(quadratic)) * (0.9995-0.95) + 0.95
G[,,600:1700] <- quadratic[1:1101]

quadratic <- function(a,x){
  return(exp(a*x) - (a*x) - 1)
}
quadratic <- quadratic(a=0.0009,x=c(-200:1300))
quadratic <- (quadratic-min(quadratic))/(max(quadratic)-min(quadratic)) * (0.9995-0.95) + 0.95
G[,,9000:10500] <- quadratic[1:1501]

y2_ge <- y1_ge + abs(min(y1_ge)) + 1
out1_ge <- filter.dlm(y2_ge,F,G,V,W,m0,C0)
out2_ge <- as.numeric(out1_ge$m) - abs(min(y1_ge)) - 1

##### 
# pooling

# weights
#####
w_b <- rep(0.2,length(y1))
w_wf <- rep(0.2,length(y1))
w_wr <- rep(0.25,length(y1))
w_ge <- rep(0.35,length(y1))

for(i in 1:length(y1)){
  if(exp(y1_wf[i]*2.303) < 0.008){
    w_b[i] <- 0.40
    w_wf[i] <- 0.10
    w_wr[i] <- 0.10
    w_ge[i] <- 0.40
  }
}

# weighting decreases gradually with diff from average of others

diff <- vector("numeric", length(y1))
for(i in 1:length(y1)){
  if(y1_wf[i] < 0.03){
    
    diff[i] <- mean(as.numeric(out1$m)[i],out2_ge[i]) - as.numeric(out1_wf$m)[i]
    w_wf[i] <- w_wf[i]-((abs(diff[i])/15))*0.5
    w_wr[i] <- w_wr[i]+((abs(diff[i])/15))*(0.5*(1/6))
    w_b[i] <- w_b[i]+((abs(diff[i])/15))*(0.5*(2.5/6))
    w_ge[i] <- w_ge[i]+((abs(diff[i])/15))*(0.5*(2.5/6))
  }
}

# weight buck 0 between 5600 and 5900 and between 2300 and 2800

w_wf[c(5600:5900,2500:2800)] <- w_wf[c(5600:5900,2500:2800)]+(w_b[c(5600:5900,2500:2800)]/3)
w_wr[c(5600:5900,2500:2800)] <- w_wr[c(5600:5900,2500:2800)]+(w_b[c(5600:5900,2500:2800)]/3)
w_ge[c(5600:5900,2500:2800)] <- w_ge[c(5600:5900,2500:2800)]+(w_b[c(5600:5900,2500:2800)]/3)
w_b[c(5600:5900,2500:2800)] <- 0



# sd posterior

sd_b <- as.numeric(out1$C)
sd_ge <- as.numeric(out1_ge$C)
sd_wf <- as.numeric(out1_wf$C)
sd_wr <- as.numeric(out1_wr$C)

# mean posterior

m_b <- as.numeric(out1$m)
m_ge <- out2_ge
m_wf <- as.numeric(out1_wf$m)
m_wr <- as.numeric(out1_wr$m)


#####

seq <- seq(-10,10,length=5000)
pool <- array(0, dim = c(1,5000,length(y1)))
mean <- vector("numeric", length(y1))
for(i in 1:length(y1)){
  
  pool[,,i] <- 1 * dnorm(seq, m_b[i], sd_b[i])^w_b[i] * 
    dnorm(seq, m_wf[i], sd_wf[i])^w_wf[i] *
    dnorm(seq, m_wr[i], sd_wr[i])^w_wr[i] * 
    dnorm(seq, m_ge[i], sd_ge[i])^w_ge[i]
  
  mean[i] <- seq[pool[,,i] == max(pool[,,i])]
  
}

plot(b878$time,exp(2.303*y1_wf),log="y",type="l",col="grey",xlab="Minutes from Midnight",ylab="Absolute Humidity")
legend("bottomright", legend=c("'True' humidity", "Hygrometers", "95% Credible Intervals"),
       col=c("blue", "grey", "darkblue"), lty=c(1,1,2), box.lty=0, cex=0.8)
lines(b878$time,exp(2.303*y1),log="y",col="grey")
lines(b878$time,exp(2.303*y1_wr),log="y",col="grey")
lines(b878$time,exp(2.303*y1_ge),log="y",col="grey")
lines(b878$time,exp(2.303*mean),log="y",lwd=1,col="blue")
#####


######################################################

# CI
#####
upper <- vector("numeric", length(y1))
lower <- vector("numeric", length(y1))

for (i in 1:length(mean)){
  ci <- data.frame(seq,cumsum=cumsum(pool[,,i]))
  upper[i] <- ci[which.min(abs(0.975*sum(pool[,,i])-
                                 cumsum(pool[,,i]))),]$seq
  lower[i] <- ci[which.min(abs(0.025*sum(pool[,,i])-
                                 cumsum(pool[,,i]))),]$seq
  
}

lines(b878$time,exp(2.303*upper),log="y",col="darkblue",lty=2)
lines(b878$time,exp(2.303*lower),log="y",col="darkblue",lty=2)

#######################################################

save(b878,m_b,m_wf,m_wr,m_ge,sd_b,sd_wf,sd_wr,sd_ge,file="fig2.RData")
save(b878,y1,y1_wf,y1_wr,y1_ge,mean,upper,lower, file="fig5.RData")



#####################################################################
################### figure 9 ########################################
#####################################################################



# weight rosemount 0 between 9900 and 12500

w_b[c(9900:12500)] <- w_b[c(9900:12500)]+(w_wr[c(9900:12500)]/2)
w_ge[c(9900:12500)] <- w_ge[c(9900:12500)]+(w_wr[c(9900:12500)]/2)
w_wr[c(9900:12500)] <- 0

# buck rosemount 0 between 13200 and 14800

w_wr[c(13200:14800)] <- w_wr[c(13200:14800)]+(w_b[c(13200:14800)]/3)
w_wf[c(13200:14800)] <- w_wf[c(13200:14800)]+(w_b[c(13200:14800)]/3)
w_ge[c(13200:14800)] <- w_ge[c(13200:14800)]+(w_b[c(13200:14800)]/3)
w_b[c(13200:14800)] <- 0



#####

seq <- seq(-10,10,length=5000)
pool <- array(0, dim = c(1,5000,length(y1)))
mean <- vector("numeric", length(y1))
for(i in 1:length(y1)){
  
  pool[,,i] <- 1 * dnorm(seq, m_b[i], sd_b[i])^w_b[i] * 
    dnorm(seq, m_wf[i], sd_wf[i])^w_wf[i] *
    dnorm(seq, m_wr[i], sd_wr[i])^w_wr[i] * 
    dnorm(seq, m_ge[i], sd_ge[i])^w_ge[i]
  
  mean[i] <- seq[pool[,,i] == max(pool[,,i])]
  
}

plot(b878$time,exp(2.303*y1_wf),log="y",type="l",col="grey",xlab="Minutes from Midnight",ylab="Absolute Humidity")
legend("bottomright", legend=c("'True' humidity", "Hygrometers", "95% Credible Intervals"),
       col=c("blue", "grey", "darkblue"), lty=c(1,1,2), box.lty=0, cex=0.8)
lines(b878$time,exp(2.303*y1),log="y",col="grey")
lines(b878$time,exp(2.303*y1_wr),log="y",col="grey")
lines(b878$time,exp(2.303*y1_ge),log="y",col="grey")
lines(b878$time,exp(2.303*mean),log="y",lwd=1,col="blue")
#####



######################################################

# CI
#####
upper <- vector("numeric", length(y1))
lower <- vector("numeric", length(y1))

for (i in 1:length(mean)){
  ci <- data.frame(seq,cumsum=cumsum(pool[,,i]))
  upper[i] <- ci[which.min(abs(0.975*sum(pool[,,i])-
                                 cumsum(pool[,,i]))),]$seq
  lower[i] <- ci[which.min(abs(0.025*sum(pool[,,i])-
                                 cumsum(pool[,,i]))),]$seq
  
}

#for(i in 1:length(lower)){
#  if((2.303*exp(lower[i])) < 0.05){
#    lower[i] <- NA
#  }
#}

lines(b878$time,exp(2.303*upper),log="y",col="darkblue",lty=2)
lines(b878$time,exp(2.303*lower),log="y",col="darkblue",lty=2)

#######################################################




