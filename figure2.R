################## 
load("fig2.RData") # refer to figure8-9.R
################## 

seq <- seq(-1.5,1,length=5000)

# distributions integrate to one 
f <- function(x){ (1 * ((2.5/10) * dnorm(x, m_ge[3000], sd_ge[3000]) + (2.5/10) * 
                          dnorm(x, (m_wf[3000]-0.3), sd_wf[3000]) + (2.5/10) *
                          dnorm(x, m_b[3000], sd_b[3000])))}
k1 <- (integrate(f, -Inf, Inf))$value

f <- function(x){ 1 * dnorm(x, (m_wf[3000]-0.3), sd_wf[3000])^0.25 * 
    dnorm(x, m_ge[3000], sd_ge[3000])^0.25 *
    dnorm(x, m_b[3000], sd_b[3000])^0.25}
k2 <- (integrate(f, min(seq), max(seq)))$value

# linear opinion pooling
par(mfrow=c(1,1))
pool <- (1 * ((2.5/10) * dnorm(seq, m_ge[3000], sd_ge[3000]) + (2.5/10) * 
                dnorm(seq, (m_wf[3000]-0.3), sd_wf[3000]) +
                (2.5/10) * dnorm(seq, m_b[3000], sd_b[3000]))) / k1

# plotting
plot(seq,dnorm(seq, (m_wf[3000]-0.3), sd_wf[3000]),xlim=c(-1.5,1),yaxt="n",xaxt="n",type="l",xlab="",ylab="",col="grey")
legend(0.2,3, legend=c("Linear pooling", "Logarithmic pooling"),
       col=c("darkorange3","deepskyblue4"), lty=1, box.lty=0)
lines(seq,(dnorm(seq, m_b[3000], sd_b[3000]))/k1,col="grey")
lines(seq,dnorm(seq, m_ge[3000], sd_ge[3000]),col="grey")
lines(seq,pool,col="darkorange4")

# logarithmic pooling

pool <- (1 * dnorm(seq, (m_wf[3000]-0.3), sd_wf[3000])^0.25 * 
           dnorm(seq, m_ge[3000], sd_ge[3000])^0.25 * 
           dnorm(seq, m_b[3000], sd_b[3000])^0.25) / k2 
lines(seq,pool,col="deepskyblue4")


