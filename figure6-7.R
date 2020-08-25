################## 
load("fig2.RData") # refer to figure8.R
################## 

####################### equal weighting #########################


seq <- seq(-5,1,length=5000)

# integrate to one 
f <- function(x){ 1 * dnorm(x, m_b[3000], sd_b[3000])^0.25 * dnorm(x, m_wf[3000], sd_wf[3000])^0.25 * 
    dnorm(x, m_wr[3000], sd_wr[3000])^0.25 * dnorm(x, m_ge[3000], sd_ge[3000])^0.25 }
k1 <- (integrate(f, -Inf, Inf))$value

# logarithmic pooling
pool <- (1 * dnorm(seq, m_b[3000], sd_b[3000])^0.25 * dnorm(seq, m_wf[3000], sd_wf[3000])^0.25 * 
           dnorm(seq, m_wr[3000], sd_wr[3000])^0.25 * dnorm(seq, m_ge[3000], sd_ge[3000])^0.25) / k1

# plotting
plot(seq,dnorm(seq, m_wf[3000], sd_wf[3000]),xlim=c(-1.5,1),type="l",xaxt='n',yaxt="n",xlab="",ylab="",col="gray85")
legend(0,2.5, legend=c("Posterior estimates of 'true' humidity", "given each hygrometer's measurements", "'Pooled' posteriors","95% credible intervals"),
       col=c("gray85", "white", "darkred", "darkred"), lty=c(1,1,1,2), box.lty=0)
lines(seq,dnorm(seq, m_b[3000], sd_b[3000]),col="gray85")
lines(seq,dnorm(seq, m_wr[3000], sd_wr[3000]),col="gray85")
lines(seq,dnorm(seq, m_ge[3000], sd_ge[3000]),col="gray85")
lines(seq,pool,col="darkorange3")

# posterior mean and credible intervals
abline(v=seq[pool == max(pool)],col="darkorange3",lty=5)

ci <- data.frame(seq,cumsum=cumsum(pool))
upp <- ci[which.min(abs(0.975*sum(pool)-ci$cumsum)),]$seq
dwn <- ci[which.min(abs(0.025*sum(pool)-ci$cumsum)),]$seq
abline(v=upp,col="darkorange3",lty=2)
abline(v=dwn,col="darkorange3",lty=2)




################# removing buck weighting #######################


f <- function(x){ 1 * dnorm(x, m_b[3000], sd_b[3000])^0 * dnorm(x, m_wf[3000], sd_wf[3000])^(0.25+(0.25/3)) * 
    dnorm(x, m_wr[3000], sd_wr[3000])^(0.25+(0.25/3)) * dnorm(x, m_ge[3000], sd_ge[3000])^(0.25+(0.25/3)) }
k1 <- (integrate(f, -Inf, Inf))$value


pool <- (1 * dnorm(seq, m_b[3000], sd_b[3000])^0 * dnorm(seq, m_wf[3000], sd_wf[3000])^(0.25+(0.25/3)) * 
           dnorm(seq, m_wr[3000], sd_wr[3000])^(0.25+(0.25/3)) * dnorm(seq, m_ge[3000], sd_ge[3000])^(0.25+(0.25/3))) / k1

plot(seq,dnorm(seq, m_wf[3000], sd_wf[3000]),xlim=c(-1.5,1),type="l",xaxt='n',yaxt="n",xlab="",ylab="",col="gray85")
legend(0,2.5, legend=c("Posterior estimates of 'true' humidity", "given each hygrometer's measurements", "'Pooled' posteriors","95% credible intervals"),
       col=c("gray85", "white", "darkorange3", "darkorange3"), lty=c(1,1,1,2), box.lty=0)
lines(seq,dnorm(seq, m_b[3000], sd_b[3000]),col="gray85")
lines(seq,dnorm(seq, m_wr[3000], sd_wr[3000]),col="gray85")
lines(seq,dnorm(seq, m_ge[3000], sd_ge[3000]),col="gray85")
lines(seq,pool,col="darkorange3")
abline(v=seq[pool == max(pool)],col="darkorange3",lty=5)

ci <- data.frame(seq,cumsum=cumsum(pool))
upp <- ci[which.min(abs(0.975*sum(pool)-ci$cumsum)),]$seq
dwn <- ci[which.min(abs(0.025*sum(pool)-ci$cumsum)),]$seq
abline(v=upp,col="darkorange3",lty=2)
abline(v=dwn,col="darkorange3",lty=2)


