library(mvtnorm)
library(MASS)
library(pracma)

#---------Task a Question 1--------

param<-read.csv(choose.files(),header=F)
mu<-param[,1]
sig<-matrix(nrow = 10,ncol=10)
for (i in 2:11) {
  sig[,i-1]<-param[,i]
}
set.seed(123)
p<-10
n<-70
data1<-rmvnorm(n,mu,sig)
mu0<-c(-0.9,0,-0.9,1,-1,1,-2,0,0,0)
xavg<-colMeans(data1)
s<-cov(data1)
Thotteling<-n*t((xavg-mu0))%*%inv(s)%*%(xavg-mu0)
c<-(n-p)/(p*(n-1))
(statisti<-c*Thotteling)
Fcritical<-qf(0.95,p,n-p)
statisti>Fcritical
(pval<-1-pf(statisti,p,n-p))


#---------Task b Question 1--------
set.seed(523511)
power<-0
n<-200
mu<-rep(1,10)
Fcritical<-qf(0.95,p,n-p)
c<-(n-p)/(p*(n-1))
for (i in 1:100) {
  data<-mvrnorm(n,mu,sig)
  xavg<-colMeans(data)
  s<-cov(data)
  Thotteling<-n*t((xavg-mu0))%*%solve(s)%*%(xavg-mu0)
  (statisti<-c*Thotteling)
  if(statisti>Fcritical){
    power<-power+1
  }
}
(power<-power/100)


#---------Task c Question 1--------
n<-70
set.seed(123)
data1<-mvrnorm(n,mu,sig)
B<-matrix(NA,2,4)
B[1,]<-c(1,-1,0,0)
B[2,]<-c(0,0,1,-1)
mu0<-rep(0,2)
xavg<-colMeans(data1)[c(1,2,4,5)]
Sx<-cov(data1)[c(1,2,4,5),c(1,2,4,5)]
Sy<-B%*%Sx%*%t(B)
Thotteling<-n*t(B%*%xavg)%*%solve(Sy)%*%B%*%xavg
c<-(n-2)/(2*(n-1))
(statisti<-c*Thotteling)
Fcritical<-qf(0.95,2,n-2)
(statisti>Fcritical)
(pval<-1-pf(statisti,2,n-2))

#---------Task a Question 2--------
n<-100
p<-3
set.seed(124)
mu<-c(0,1,0)
sig<-matrix(c(16,1,1,1,16,2,1,2,16), nrow = 3, ncol = 3)
data<-mvrnorm(n,mu,sig)
s<-cov(data)
s0<-matrix(c(s[1,1],0,0,0,s[2,2],0,0,0,s[3,3]),3,3)
R<-cor(data)
(statisti<--n*log(det(R)))
chiCritical<-qchisq(0.95,p*(p-1)/2)
(statisti>chiCritical)
(pval<-1-pchisq(statisti,3))


#---------Task b Question 2--------
mu<-c(0,0,0)
sig<-matrix(c(2,0,1,0,2,0,1,0,2),3,3)
power<-0
n<-200
p<-3
chiCritical<-qchisq(0.95,p*(p-1)/2)
for (i in 1:100) {
  data<-mvrnorm(n,mu,sig)
  R<-cor(data)
  statisti<--n*log(det(R))
  if(statisti>chiCritical)
    power<-power+1
}
(power<-power/100)


#---------Task b Question 3--------
n<-44
p<-2
mu<-c(1,1)
sig<-diag(2)
chiCritical<-qchisq(0.95,2)
power<-0
for (i in 1:50000) {
 data<-mvrnorm(n,mu,sig)
 xavg<-colMeans(data)
 statisti<-n*t(xavg)%*%xavg
 if(statisti>chiCritical)
   power<-power+1
}
(power<-power/50000)

#---------Task d Question 3--------
n<-44
p<-2
mu<-c(1,1)
sig<-diag(2)
chiCritical<-qchisq(0.95,2)
power<-0
for (i in 1:50000) {
  data<-mvrnorm(n,mu,sig)
  xavg<-colMeans(data)
  statisti<--2*n*log(1- n*(sum(xavg^2))/sum(data^2))
  if(statisti>chiCritical)
    power<-power+1
}
(power<-power/50000)


