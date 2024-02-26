# The first dataset came from an aluminum electrolytic capacitor (AEC) manufacturing process.

library(expm)
library(stats)
library(mvtnorm)
library(psych)
library(MVN)

d=read.table("C:/Users/pc/Desktop/submission pattern recognition/water_data.txt", header = T)
d=as.matrix(d)
p=dim(d)[2] # the dimension

#

constnorm=1/(2*pi)^(p/2)

f0=function(x,y,sH, dH) {
  
  exp(-t(y-x)%*%sH%*%(y-x) / 2) / sqrt(dH) }

f1=function(x,y,sH,dH,rn) {
  
  exp(-t(y-x-rn)%*%sH%*%(y-x-rn) / 2) / sqrt(dH) }

pI=function(data, mm0) {
  
  mn=apply(data[1:mm0,], 2, mean)
  CCn=lapply( 1:mm0, function(x) data[x,]%*%t(data[x,]))
  sumn=Reduce("+",CCn)
  Cn=sumn/(mm0) -  mn%*%t(mn)
  Cn=cov(data[1:mm0,])
  return(list(mn, Cn, sumn))
  
}
m0=288 # training dataset

ic0=d[1:m0,]
ie=pI(ic0,m0)
mn=colMeans(ic0)
Cn=cov(ic0)
sumn=ie[[3]]
Gn=(4/(p+2))^(2/(p+4))*m0^(-2/(p+4))*Cn # bandwidth
set.seed(1)
KN=kmeans(ic0,2, nstart = 20, iter.max = 50)
cl=KN$cluster
dj=as.list(rep(0,2))
for (I in 1:2) {
  dj[[I]]=which(cl==I)
}

if (length(dj[[1]])<8) {mv=1} else{ 
  mv=min(as.numeric(as.character(mvn(ic0[dj[[1]],], mvnTest="mardia", covariance = FALSE)$multivariateNormality[1:2,3])))}

if (length(dj[[2]])<8) {mv=c(mv,1)} else{ 
  mv=c(mv,min(as.numeric(as.character(mvn(ic0[dj[[2]],], mvnTest="mardia", covariance = FALSE)$multivariateNormality[1:2,3]))))}

if(min(mv)<0.10) {aw=0} else {aw=1}


while (aw==0) {
  
  lik=which(mv<0.10)
  ind=length(dj)
  
  for (iw in lik) {
    set.seed(1)
    KN0=kmeans(ic0[dj[[iw]],],2, nstart = 20, iter.max = 50)
    cl0=KN0$cluster
    
    dj1=dj[[iw]][which(cl0==1)]
    dj2=dj[[iw]][which(cl0==2)]
    
    dj[[iw]]=dj1
    if (length(dj1)<8) {mv[iw]=1} else{
      mv[iw]=min(as.numeric(as.character(mvn(ic0[dj[[iw]],], mvnTest="mardia", covariance = FALSE)$multivariateNormality[1:2,3])))}
    dj[[length(dj)+1]]=dj2
    if (length(dj2)<8) {mv[length(dj)]=1} else{
      mv[length(dj)]=min(as.numeric(as.character(mvn(ic0[dj[[length(dj)]],], mvnTest="mardia", covariance = FALSE)$multivariateNormality[1:2,3])))}  
    
  }
  
  if (min(mv)>0.10) {aw=1}  
}

wn=sapply(dj, length)

N=length(wn)

wij=outer(wn/sum(wn),wn/sum(wn))
wij[lower.tri(wij, diag = FALSE)]=0

Ssj=as.list(rep(0,N))
kmn=matrix(NA,N,p)
for (I in 1:N) {
  
  if (length(dj[[I]])>1) {
    
    Ssj[[I]]= cov(ic0[dj[[I]],] )   
    kmn[I,]=apply(ic0[dj[[I]],],2,mean) } else { kmn[I,]=ic0[dj[[I]],]}
  
  
  
}


Sgj=lapply(Ssj, function(x) x+Gn)

{Fc = Cn}

RF1=matrix(0,N,N)
for (jj in 1:(N)) {
  for (ii in 1:jj) {
    Sumij=(Sgj[[ii]] + Ssj[[jj]])
    Aij = solve(Sumij)
    Dij = kmn[ii,]-kmn[jj,]
    Bij = Aij%*%(diag(p) - 2*Dij%*%t(Dij)%*%Aij)
    Cij = Aij%*%(diag(p) - Dij%*%t(Dij)%*%Aij)
    
    mij = t(Dij)%*%Aij%*%Dij
    
    rf=wij[ii,jj]*exp(-mij/2) * sqrt(det(Aij)) * constnorm *
      ( 2*tr(Fc%*%Aij%*%Fc%*%Bij) + (tr(Fc%*%Cij)^2) )
    
    RF1[ii,jj]=rf
    
  }}
Rf=2*sum(RF1)-tr(RF1)
Ht=Cn*(m0^(-1) /( sqrt(4*pi)^p * Rf * p ))^(2/(p+4))

Sm0=lapply(Ssj, function(x) x+Ht)


sSm0=lapply(Sm0, solve)
dSm0=lapply(Sm0, det)



Sn=0
SN=c()
k=1
dn=d[m0+1,]-mn

for (n in (m0+1):(dim(d)[1]-1) ) {
  
  Ln = log ( sum(mapply(f1,x=asplit(kmn[1:N,],1),sH=sSm0, dH=dSm0, MoreArgs = list(y=d[n+1,], rn=dn))*wn)/sum(wn) ) - log ( sum(mapply(f0,x=asplit(kmn[1:N,],1),sH=sSm0, dH=dSm0, MoreArgs = list(y=d[n+1,]))*wn)/sum(wn) )
  
  
  if (is.nan(Ln)) { 
    
    if ((d[n+1,]-dn)^T%*%solve(Cn)%*%(d[n+1,]-dn) < (d[n+1,]-mn)^T%*%solve(Cn)%*%(d[n+1,]-mn)) {
      Ln=5} else (Ln=0) }
  
  Sn=max(0,Sn+Ln)
  SN=c(SN,Sn)
  
  if (Sn==0) {
    
    if (k==1) {
      
      tn=d[n,]-mn
      Cn=((n-2)*Cn/(n-1)+tn%*%t(tn)/n ) 
      Gn=(4/(p+2))^(2/(p+4))*n^(-2/(p+4))*Cn # bandwidth
      
      a=unname(which.min(apply((t(t(kmn)-d[n,]))^2,1,sum)))
      kmn[a,]=(kmn[a,]*wn[a]+d[n,])/(wn[a]+1)
      wn[a]=wn[a]+1
      wij=outer(wn/sum(wn),wn/sum(wn))
      wij[lower.tri(wij, diag = FALSE)]=0
      
      dj[[a]]=c(dj[[a]],n)
      Ssj[[a]]=cov(d[dj[[a]], ]) 
      
    } else{
      
      
      mn=(n-k)/n*mn + apply(d[(n-k+1):n,], 2, sum)/n
      CCk=lapply( (n-k+1):n, function(x) d[x,]%*%t(d[x,]))
      sumn=sumn+Reduce("+",CCk)
      Cn=sumn/n -  mn%*%t(mn)
      Gn=(4/(p+2))^(2/(p+4))*n^(-2/(p+4))*Cn # bandwidth
      
      
      for (j in k:1) {
        
        a=unname(which.min(apply((t(t(kmn)-d[(n-j+1),]))^2,1,sum)))
        kmn[a,]=(kmn[a,]*wn[a]+d[(n-j+1),])/(wn[a]+1)
        wn[a]=wn[a]+1
        wij=outer(wn/sum(wn),wn/sum(wn))
        wij[lower.tri(wij, diag = FALSE)]=0
        
        dj[[a]]=c(dj[[a]],n-j+1)
        Ssj[[a]]=cov(d[dj[[a]], ]) 
        
      }
      
      Sgj=lapply(Ssj, function(x) x+Gn)
      
      
      {Fc = Cn}
      
      RF1=matrix(0,N,N)
      for (jj in 1:(N)) {
        for (ii in 1:jj) {
          Sumij=(Sgj[[ii]] + Ssj[[jj]])
          Aij = solve(Sumij)
          Dij = kmn[ii,]-kmn[jj,]
          Bij = Aij%*%(diag(p) - 2*Dij%*%t(Dij)%*%Aij)
          Cij = Aij%*%(diag(p) - Dij%*%t(Dij)%*%Aij)
          
          mij = t(Dij)%*%Aij%*%Dij
          
          rf=wij[ii,jj]*exp(-mij/2) * sqrt(det(Aij)) * constnorm *
            ( 2*tr(Fc%*%Aij%*%Fc%*%Bij) + (tr(Fc%*%Cij)^2) )
          
          RF1[ii,jj]=rf
          
        }}
      Rf=2*sum(RF1)-tr(RF1)
      Ht=Cn*(n^(-1)  /( sqrt(4*pi)^p * Rf * p ))^(2/(p+4))

      Sm0=lapply(Ssj, function(x) x+Ht)
      sSm0=lapply(Sm0, solve)
      dSm0=lapply(Sm0, det)

    }
    
    dn=d[n+1,] - mn
    k=1
    
  } else {
    
    k=k+1
    dn=apply(d[(n-k+2):(n+1),], 2, mean) - mn }     
  
}

SNf=c(NA,SN)
OC_IC_diff_mean=c(-1.36405584, -0.56991579,  0.07425124, -0.70378801)

mg=rgb(157,220,93, max=255)
mg2=rgb(0,255,0,25, max=255)
mgd=rgb(104,184,21, max=255)
mbd=rgb(16,81,143, max=255)
mar=rgb(238,48,167, alpha=150, max=255)
mar2=rgb(174,55,129, alpha=255, max=255)
br=rgb(238,44,44, alpha=150, max=255)

hs=6.75

####################################################################

dev.new(width=8, height=6)

par(fig=c(0,0.5,0.75,1), oma=c(0,0.2,0,0.2), mar=c(3.5,3.5,2,0.5))


plot(d[,1], col=rgb(46, 119, 191, alpha = 125, max=255), type="l",
     xlab="",ylab="",lwd=1,axes=F, yaxs = "i", ylim = c(-2.5,3))
axis(1) ; axis(2)
mtext("Time",side=1, line=2.2, cex=1)
mtext("Values",side=2,  line=2.2, cex=1)
mtext(expression(bold(Y)[1]), side=3,  line=0.5, cex=1.3)
text(140,2.75,"Training data")
text(325,2.75,"Test data")
points(289:349,d[289:349,1], col="dodgerblue2", type="l")
points(350:360,d[350:360,1], col="firebrick2", type="l")
abline(v=288.5, lwd=2)

####################################################################

par(fig=c(0.5,1,0.75,1), oma=c(0,0.2,0,0.2), mar=c(3.5,3.5,2,0.5), new=T)


plot(d[,2], col=rgb(46, 119, 191, alpha = 125, max=255), type="l",
     xlab="",ylab="",lwd=1,axes=F, yaxs = "i", ylim = c(-8.2,6.8))
axis(1) ; axis(2)
mtext("Time",side=1, line=2.2, cex=1)
mtext("Values",side=2,  line=2.2, cex=1)
mtext(expression(bold(Y)[2]), side=3,  line=0.5, cex=1.3)
text(140,5,"Training data")
text(325,5,"Test data")
points(289:349,d[289:349,2], col="dodgerblue2", type="l")
points(350:360,d[350:360,2], col="firebrick2", type="l")
abline(v=288.5, lwd=2)

####################################################################

par(fig=c(0,0.5,0.5,0.75), oma=c(0,0.2,0,0.2), mar=c(3.5,3.5,2,0.5), new=T)


plot(d[,3], col=rgb(46, 119, 191, alpha = 125, max=255), type="l",
     xlab="",ylab="",lwd=1,axes=F, yaxs = "i", ylim = c(-7,5.3))
axis(1) ; axis(2)
mtext("Time",side=1, line=2.2, cex=1)
mtext("Values",side=2,  line=2.2, cex=1)
mtext(expression(bold(Y)[3]), side=3,  line=0.5, cex=1.3)
text(140,4.5,"Training data")
text(325,4.5,"Test data")
points(289:349,d[289:349,3], col="dodgerblue2", type="l")
points(350:360,d[350:360,3], col="firebrick2", type="l")
abline(v=288.5, lwd=2)

####################################################################

par(fig=c(0.5,1,0.5,0.75), oma=c(0,0.2,0,0.2), mar=c(3.5,3.5,2,0.5), new=T)


plot(d[,4], col=rgb(46, 119, 191, alpha = 125, max=255), type="l",
     xlab="",ylab="",lwd=1,axes=F, yaxs = "i", ylim = c(-5.9,3.5))
axis(1) ; axis(2)
mtext("Time",side=1, line=2.2, cex=1)
mtext("Values",side=2,  line=2.2, cex=1)
mtext(expression(bold(Y)[3]), side=3,  line=0.5, cex=1.3)
text(140,3.25,"Training data")
text(325,3.25,"Test data")
points(289:349,d[289:349,4], col="dodgerblue2", type="l")
points(350:360,d[350:360,4], col="firebrick2", type="l")
abline(v=288.5, lwd=2)

####################################################################

par(fig=c(0,1,0.19,0.5),oma=c(0,0.2,0,0.2), mar=c(3.5,3.5,2,0.5), new=T)


plot(27000,27000,main="",xlim=c(1,72),ylim=c(-0.5,26.5), axes=F, yaxs = "i", xlab="", ylab="")
axis(1, at=c(1,12,24,36,48,60,72), labels=c(1,12,24,36,48,60,72)) 
axis(2, at=seq(0,25,by=5), labels=seq(0,25,by=5))
mtext("n",side=1, line=2.1, cex=1.4)
mtext(expression(paste(S[n])),side=2,  line=2.1, cex=1.4)
mtext("MKC",side=3, line =0.3 , cex = 2 )

rect(0.75,-0.5,72.25,hs, col=mg2, border = mg2)
segments(0.75,hs,72.25,hs,col="firebrick1", lwd=2, lty=1)

points(1:72,SNf[1:72],type="b", lwd=2, col="dodgerblue2", lty=1, pch=16, cex=1.4)
points(65:72,SNf[65:72], lwd=3, col="firebrick3", cex=1.8, pch=4)
draw.circle(65,SNf[65],2.25/5,border="firebrick1",lwd=2)
points(62,SNf[62], lwd=2, col="firebrick4", cex=1.8, pch=16)




legend("topleft", legend=c("h", 
                           expression(paste(S[n])), "Alarm", expression(T==65), 
                           expression(hat(tau)==62)), 
       lwd=c(3,2,2,2,2), col=c("firebrick2","dodgerblue2","firebrick3", "firebrick1", "firebrick4"), 
       bty ="n", cex=1.1, 
       pt.cex=c(0,1,1.25,1.25,1.225), lty=c(1,1,0,0,0), pch=c(-1,16,4,1,1), 
       y.intersp = 1)


####################################################################

par(fig=c(0,1,0,0.25),oma=c(0,0.2,0.2,0.2), mar=c(3.5,3.5,4,0.5), new=T)

barplot(rep(NA,4), ylim=c(-1.5,0.5), las=1, 
        names.arg = expression(Y[1], Y[2], Y[3], Y[4]), axes=F, 
        cex.names = 1.2) 
axis(2, at=c(-1.5,-1,-0.5,0,0.5))
segments(0,c(-1.5,-1,-0.5,0,0.5),5.75,c(-1.5,-1,-0.5,0,0.5), lty=2, col="gray75")
barplot(OC_IC_diff_mean, ylim=c(-1,2), las=1, axes=F, 
        col=c("firebrick4", "firebrick1", "green3", "firebrick3"), 
        add = T, cex.names = 1.2) 

mtext(expression(paste("(",hat(mu)[OOC]-hat(mu)[IC],") / diag(",hat(Sigma)[IC]^{1/2},")" )),side=3,  line=1, cex=1.5)


#############################################################################
#################   post-alarm inference   ##################################
#############################################################################

### IC mean
apply(d[1:349,],2,mean)

### OOC mean (until the alarm)
apply(d[350:353,],2,mean)

### OOC mean (until the end of the sample)
apply(d[350:360,],2,mean)



(apply(d[350:353,],2,mean)-apply(d[1:349,],2,mean))/apply(d[1:349,],2,sd)


