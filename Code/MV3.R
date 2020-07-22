library(mvtnorm)
library(MASS)
library(pracma)


#------Question 1------------
#-------task a---------
update_h <- function(hPrev, R){ #function for updating value in each iteration
  (psi<-1-hPrev)
  a<-R-diag(psi,p,p)
  eigenv<-eigen(a)
  (eigenvalues<-eigenv$values)
  (k<-1)
  (eigenvec<-eigenv$vectors)
  lambda<-matrix(nrow = p, ncol=k)
  for (i in c(1:k)) {
    (i_max<-which (eigenvalues==max(eigenvalues)))
    (lambda[,i]<-sqrt(eigenvalues[i_max])*as.matrix(eigenvec)[,i_max])
    (eigenvalues<-eigenvalues[-i_max])
    (eigenvec<-as.matrix(eigenvec)[,-i_max])
  }
  return (lambda)
}



R<-matrix(c(1,0.34,0.26,0.34,1,0.12,0.26,0.12,1), nrow = 3)
p<-3
n<-200
epsilon<-10^-20
hPrev<- c(0.34,0.34,0.26)
flag<-TRUE
count<-0
while(flag){
  count<-count+1
  (lambda<-update_h(hPrev, R))
  (hCurr<-rowSums(lambda^2))
  if(sum(abs(hPrev-hCurr)< epsilon)==3){ #stop condition
    flag<-FALSE
  }else{
    hPrev<-hCurr
  }
}
(count)
(psi<-diag(1-hCurr,p,p))
(lambda)

#-------task b---------
likelihoodFunc<- function(psi, lambda, r){ #function we want to minimize - we remove the minus
  psi<-diag(psi)
  A<-lambda%*%t(lambda)+psi
  return((log(det(A))+tr(solve(A)%*%R)))
}

R<-matrix(c(1,0.34,0.26,0.34,1,0.12,0.26,0.12,1), nrow = 3)
epsilon<-10^-15
psi0<-rep(1,p)
psi0<-c(0.66,0.66,0.74)
lambda_t<-c(0,1,1)
flag<-TRUE
count<-0
psi_t<-psi0
while (flag) { 
  count<-count+1
  result<-optim(par = c(0,1,2), fn=likelihoodFunc, psi=psi0, r=R)
  (lambda_t<-result$par)
  result<-optim(par = c(0,1,2), fn=likelihoodFunc, lambda=lambda_t, r=R)
  (psi_t<-result$par)
  (x<-as.numeric(t(psi_t-psi0)%*%(psi_t-psi0))) 
  if(x< epsilon | count>500){ #stop condition
    flag<-FALSE
  }else{
    psi0<-psi_t}
  (flag)
}
psi<-diag(psi_t,p)
(psi)
t(lambda_t)


  
#------Question 2--------
distance<-function(x,y){ #distance function
  return(as.numeric(t(x-y)%*%(x-y)))
}

ss<-function(points, avg=c(0,0)){ #function that calculate sum of squares
  sumi<-0
  if(Norm(avg)==0)
    avg<-colMeans(points)
  for (i in 1:nrow(points)) {
    sumi<-sumi+as.numeric(t(points[i,]-avg)%*%(points[i,]-avg))
  }
  return(sumi)
}
set.seed(111)
p1<-c(1.9,0.64) #initial parameters
p2<-c(0.87,-1.2)
p3<-c(0.3,0.09)
p4<-c(-2.08,-0.45)
p5<-c(0.85,0.46)
a<-runif(1,0,1)
b<-runif(1,1,5)
c<-rnorm(1,5,1)
points<-rbind(p1,p2,p3,p4,p5)
avgPoi<-colMeans(points)
centers<-points
cluster<-list(c(1),c(2),c(3),c(4),c(5))
flag<-T
k<-5
(measure<--b*ss(centers,avgPoi)+c*k)
clusterBef<-cluster


while(flag && k>1){ #stop condition
    distMat<-matrix(NA, nrow = k,ncol = k) #distance matrix
    for (i in 1:k) {
      for (j in 1:k) {
        distMat[i,j]<-ifelse(i==j,10000 ,distance(centers[i,],centers[j,]))
      }
    }
    temp<- which(distMat==min(distMat), arr.ind = T)[,1]
    clusterBef<-cluster
    cluster[[temp[1]]]<-rbind(cluster[[temp[1]]],cluster[[temp[2]]])
    cluster<-cluster[-cluster[[temp[2]]]]
    centers<-matrix(NA, nrow=length(cluster), ncol = 2)
    ssw<-0
    for (i in 1:length(cluster)) {
      if(length(cluster[[i]])==1)
        centers[i,]<-points[cluster[[i]],]
      else{
        centers[i,]<-colMeans(points[cluster[[i]],])
        ssw<-ssw+ss(points[cluster[[i]],])
      }
    }
    ssb<-ss(centers,avgPoi)
    print(a*ssw-b*ssb+c*nrow(centers))
    if(measure<a*ssw-b*ssb+c*nrow(centers))
      flag<-F
    else{
      measure<-a*ssw-b*ssb+c*nrow(centers)
      k<-k-1    
    }
}
(clusterBef)


#-----Question 3-------
#task a
data<-read.csv(choose.files(),header=F) #uploading the data

optimise(fBin, interval = c(0.00001^100,10^100), maximum = TRUE, tol = 0.00001)


likelihood<-function(theta){
  -1*(theta*sum(data[,1]))+3*sum(data[,2])*log(theta)+
    sum(log(1+theta*data[38:100,1]+(theta*data[38:100,1])^2/2))
}

fff<-function(theta){
  a<-sum(data[,1])
  k<-sum(data[,2])
  b<-k*log((theta^3)/2)
  sum<-0
  for (i in 1:k){ 
    sum<-sum+log(data[i,1]^2)}
  for (i in 38:100){ 
    sum<-sum+log(1+theta*data[i,1]+(theta*data[i,1])^2/2)}
  
  return(-1*theta*a+b+sum)
  
}
#View(data)
optimise(fff, interval = c(0.00001^100,10^100), maximum = TRUE, tol = 0.00001)

optimise(likelihood, interval = c(0.00001,10^100), maximum = TRUE, tol = 0.00001)


f<-function(theta){
  c<- -1*sum(data[,1])+3*37/theta
  sum<-0
  for (i in 38:100) {
    sum<- sum+ (data[i,1]+data[i,1]^2*theta)/(1+theta*data[i,1]+(theta*data[i,1])^2/2)
  }
  return(sum+c)
}

ff<-function(theta){
-sum(data[,1])+3*37/theta+sum(data[38:100,1]+theta*(data[38:100,1])^2/(1+data[38:100,1]*(theta+theta^2/2)))
}

uniroot.all(f, c(0.00001,100000) )


#-----task d--------
theta<-runif(1,0,10)
flag<-TRUE
epsilon<-10^-10
count<-0
while (flag) {
  count<-count+1
  temp<-3*nrow(data)/
    (sum(data[1:37,1])+
       (1/theta)*(sum(((theta*data[38:100,1])^3+3*((theta*data[38:100,1])^2)+6*theta*data[38:100,1]+6)/
                    (2+2*(theta*data[38:100,1])+(theta*data[38:100,1])^2))))
  if(abs(theta-temp)<epsilon)
    flag<-FALSE
  else
    theta<-temp
}
(theta)


