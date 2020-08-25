library(R.matlab)
B878 <- readMat("b878_dv.mat")
b878 <- data.frame(time=B878$Time,buck=log10(B878$buck),wvss2f=log10(B878$wvss2f),wvss2r=log10(B878$wvss2r),ge=log10(B878$ge),pr=B878$PS.RVSM)

# simple NA imputation
##### 

y1 <- b878$buck
y1_wf <- b878$wvss2f
y1_wr <- b878$wvss2r
y1_ge <- b878$ge
pr1 <- b878$pr
y1_tdb <- exp(2.303*log10(B878$TDEW.CR2))

for(i in 2:length(y1)+1){
  if(is.na(y1[i])){
    y1[i] <- y1[i-1] }
  
  if(is.na(y1_wf[i])){
    y1_wf[i] <- y1_wf[i-1] }
  
  if(is.na(y1_wr[i])){
    y1_wr[i] <- y1_wr[i-1] }
  
  if(is.na(y1_ge[i])){
    y1_ge[i] <- y1_ge[i-1] }
  
  if(is.na(y1_tdb[i])){
    y1_tdb[i] <- y1_tdb[i-1] }
  
}
y1 <- y1[1:18559]
y1_wf <- c(y1_wf[3],y1_wf[3],y1_wf[3:18559])
y1_wr <- c(y1_wr[3],y1_wr[3],y1_wr[3:18559])
y1_ge <- y1_ge[1:18559]
y1_tdb <- y1_tdb[1:18559]


#####
plot(b878$time,exp(2.303*y1_wf),log="y",type="l",col="darkred",xlab="Minutes from Midnight",ylab="Absolute Humidity")
legend("bottomright", legend=c("Buck", "WVSS2F", "WVSS2R","GE"),
       col=c("darkblue", "darkred", "darkorange3","darkslategray4"), lty=1, box.lty=0, cex=0.8)
lines(b878$time,exp(2.303*y1),log="y",col="darkblue")
lines(b878$time,exp(2.303*y1_wr),log="y",col="darkorange3")
lines(b878$time,exp(2.303*y1_ge),log="y",col="darkslategray4")

#############

save(y1,y1_wf,y1_wr,y1_ge, file="fig8.RData")



