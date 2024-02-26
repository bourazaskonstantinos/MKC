# install.packages("ks")
# install.packages("mvtnorm")
# install.packages("plotrix")

##################################################################
#################   Libraries   ##################################
##################################################################


library(ks)
library(mvtnorm)
library(plotrix)
##############################################################################

####################################################################
#################   gravel data   ##################################
####################################################################


ICX=c(5.4,3.2,5.2,3.5,2.9,4.6,4.4,5.0,8.4,4.2,3.8,4.3,3.7,3.8,2.6,
      2.7,7.9,6.6,4.0,2.5,3.8,2.8,2.9,3.3,7.2,7.3,7.0,6.0,7.4,6.8, 
      6.3,6.1,6.6,6.2,6.5,6.0,4.8,4.9,5.8,7.2,5.6,6.9,7.4,8.9,10.9,
      8.2,6.7,5.9,8.7,6.4,8.4,9.6,5.1,5.0,5.0,5.9)
ICY=c(93.6,92.6,91.7,86.9,90.4,92.1,91.5,90.3,85.1,89.7,92.5,91.8,91.7,90.3,94.5,
      94.5,88.7,84.6,90.7,90.2,92.7,91.5,91.8,90.6,87.3,79.0,82.6,83.5,83.6,84.8,
      87.1,87.2,87.3,84.8,87.4,86.8,88.8,89.8,86.9,83.8,89.2,84.5,84.4,84.3,82.2,
      89.8,90.4,90.1,83.6,88.0,84.7,80.6,93.0,91.4,86.2,87.2)


d=cbind(ICX,ICY)


######################################################################
#################   MKC procedure   ##################################
######################################################################

# The process takes about 20-30 minutes to run due to the large number of iterations (It)
# for the estimation the adaptive decision limit hn.


p=2

f0=function(x,y,sH) {
  
  exp(-t(y-x)%*%sH%*%(y-x) / 2) }

f1=function(x,y,sH,rn) {
  
  exp(-t(y-x-rn)%*%sH%*%(y-x-rn) / 2) }

mn=d[1:3,]
wmn=apply(d[1:3,],2,mean)
wCn=cov(d[1:3,])
N=3
Hn=(4/(p+2))^(2/(p+4))*N^(-2/(p+4))*wCn # bandwidth
sHn=solve(Hn)

dn=d[N+1,] - wmn
k=1
SN=c()
Sn=0



It=100000 # iterations for the derivation of the adaptive decision limit hn

Sht=rep(0,It) 
HT=rep(NA,16)

for (n in 4:(dim(d)[1]-1)) {
  
  if (n>=16) {
    
    set.seed(1+n)
    nn=sample(1:N,It, replace = T)
    set.seed(2+2*n)
    nnd=sample(1:N,It, replace = T)
    
    simd=matrix(NA, It, 2)
    set.seed(3+3*n)
    for (ii in 1:N) {
      kk=which(nn==ii)
      simd[kk,] = rmvnorm(length(kk),mn[ii,],Hn) }
    
    dk=matrix(NA, It, 2)
    set.seed(4+4*n)
    for (ii in 1:N) {
      kkd=which(nnd==ii)
      dk[kkd,] = rmvnorm(length(kkd),mn[ii,],Hn) }
    
    ddk=t(apply(dk, 1, function(x) x-wmn))
    
    log_fun=function(jj) {
      
      lf = log ( sum(apply(mn,1,f1,y=simd[jj,],sH=sHn,rn=ddk[jj,])) ) - log ( sum(apply(mn,1,f0,y=simd[jj,],sH=sHn)) )
      
      return(lf)
    }
    
    
    Lht=apply(matrix(1:It,1,It),2,log_fun)    
    Sht=ifelse(Sht+Lht>0,Sht+Lht,0)
    ht=quantile(Sht,(1-1/500))
    HT=c(HT,ht)
    print(ht)    }
  
  
  
  Ln = log ( sum(apply(mn,1,f1,y=d[n+1,],sH=sHn,rn=dn)) ) - log ( sum(apply(mn,1,f0,y=d[n+1,],sH=sHn)) )
  Sn=max(0,Sn+Ln)
  print(d[n+1,])
  print(c(n+1,Ln,Sn))
  SN=c(SN,Sn)
  
  if (Sn==0) {
    
    mn=d[1:n,]
    wmn=apply(d[1:n,],2,mean)
    wCn=cov(d[1:n,])
    
    Hn=(4/(p+2))^(2/(p+4))*(n)^(-2/(p+4))*wCn # bandwidth
    sHn=solve(Hn)
    
    dn=d[n+1,] - wmn
    
    N=n
    k=1
  } else {
    
    k=k+1
    dn=apply(d[(n-k+2):(n+1),], 2, mean) - wmn 
  }
  
  
  
}


#################################################################################
#################   Graphical representation   ##################################
#################################################################################






dev.new( width=6, height=9)


SNf=c(NA,NA,NA,NA,SN)
hn=unname(HT)


par(mfrow=c(2,1),oma=c(0,0.2,0,0.2), mar=c(3.5,3.5,2,0.5))
par(bg="white")

par(fig=c(0,1,1/3,1))

plot(1,1,xlim=c(min(ICX),max(ICX)), ylim=c(min(ICY),max(ICY)), xlab="",ylab="")
abline(h=c(80,85,90,95),v=c(4,6,8,10),col="gray75")
mtext("Medium %",side=2, line = 2.2 , cex = 1.2 )
mtext("Large %",side=1, line =2.2 , cex = 1.2 )
mtext("Observations",side=3, line =0.5 , cex = 1.8 )
points(ICX[1:16],ICY[1:16], pch=1, cex=1.25, col=rgb(46, 119, 191, alpha = 125, max=255))
points(ICX[17:21],ICY[17:21], pch=16, cex=1.25, col="dodgerblue3")
points(ICX[22:29],ICY[22:29], pch=16, cex=1.25, col="firebrick3")
#points(ICX[25],ICY[25], pch=16, cex=1.75, col="firebrick4")
points(ICX[30:56],ICY[30:56], pch=16, cex=1.25, col=rgb(238, 44, 44, alpha = 85, max=255))
#points(ICX[29:56],ICY[29:56], pch=16, cex=1.25, col=rgb(238, 44, 44, alpha = 65, max=255))
points(4.405,90.767, pch=4, col="dodgerblue4", cex=2, lwd=4)
points(5.633,86.078, pch=4, col="firebrick4", cex=2, lwd=4)
#draw.circle(ICX[30], ICY[30], 0.16, border="firebrick2", lwd=2 )
#points(ICX[31:56], ICY[31:56], pch=4, col="firebrick2", cex=1.3, lwd=2)


legend("topright", legend=c("Monitoring data", "IC data", "OOC data until T", "Post-alarm data",
                            expression(paste(bar(y)[{IC}])), expression(paste(bar(y)[{OOC}]))), lwd=c(2,0,0,0,3,3), 
       col=c(rgb(46, 119, 191, alpha = 125, max=255),"dodgerblue3","firebrick3", col=rgb(238, 44, 44, alpha = 85, max=255), "dodgerblue4", "firebrick4"), 
       bty ="n", cex=0.9, lty=0,
       pch=c(1,16,16,16,4,4), 
       y.intersp = 1.1)


#########################################################


par(fig=c(0,1,0,0.32), new = T)


mg=rgb(157,220,93, max=255)
mg2=rgb(0,255,0,25, max=255)
mgd=rgb(104,184,21, max=255)
mbd=rgb(16,81,143, max=255)
mar=rgb(238,48,167, alpha=150, max=255)
mar2=rgb(174,55,129, alpha=255, max=255)
br=rgb(238,44,44, alpha=150, max=255)


plot(100,100,main="",xlim=c(2.5,56),ylim=c(-1.5, 35), axes=F, yaxs = "i", xlab="", ylab="")
axis(1, at=seq(5,55,by=10), labels=seq(5,55,by=10)) 
axis(2, at=seq(0,35,by=5), labels=seq(0,35,by=5))
mtext("n",side=1, line=2.2, cex=1.2)
mtext(expression(paste(S[n])),side=2,  line=2.2, cex=1.2)
mtext("MKC",side=3, line =0.5 , cex = 1.8 )
polygon( c(15.5, 17:55, 56.5, 56.5, 55:17, 15.5), c(hn[16:56],rep(-1.25,41)), col=mg2, border = mg2)
points(c(1:15,15.5,17:55, 56.5), hn,col="firebrick1", lwd=2, lty=1, type="l")


points(1:16,SNf[1:16],type="b", lwd=2, col=rgb(46, 119, 191, alpha = 125, max=255), lty=3, pch=1)
points(17:56,SNf[17:56],type="b", lwd=2, col="dodgerblue3", lty=2, pch=16)
points(30:56,SNf[30:56], lwd=2, col="firebrick3", cex=1.25, pch=4)
points(22,SNf[22], lwd=2, col="firebrick4", cex=1.5, pch=16)
draw.circle(30,SNf[30],4/5,border="firebrick1",lwd=3)

text(28.5,8, "T=30", col = "firebrick1", cex=1.25, adj=1)
text(30.5,1, expression(hat(tau)==22), col = "firebrick4", cex=1.25)

legend("topleft", legend=c(expression(paste(h[n])), "Monitoring data",
                           expression(paste(S[n])), "Alarm", "T", expression(hat(tau))), lwd=c(3,2,2,2,2,2,2), 
       col=c("firebrick2",rgb(46, 119, 191, alpha = 125, max=255),"dodgerblue3","firebrick3", "firebrick1", "firebrick4"), 
       bty ="n", cex=0.9, 
       pt.cex=c(0,1,1,1.25,1.25,1.25), lty=c(1,3,1,0,0,0), pch=c(-1,1,16,4,1,16), 
       y.intersp = 1.1)



#############################################################################
#################   post-alarm inference   ##################################
#############################################################################

### IC mean
apply(d[1:21,],2,mean)

### OOC mean (until the alarm)
apply(d[22:30,],2,mean)

### OOC mean (until the end of the sample)
apply(d[22:56,],2,mean)


### OOC mean excluding the observations 22-24 (until the alarm)
apply(d[25:30,],2,mean)

### OOC mean excluding the observations 22-24 (until the end of the sample)
apply(d[25:56,],2,mean)



 