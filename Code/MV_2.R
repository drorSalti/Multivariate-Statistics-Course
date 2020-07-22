library(mvtnorm)
library(MASS)
library(pracma)


#------Question 1------------
data<-read.csv(choose.files(),header=T) #uploading the data
data<-as.matrix(data[,-1])
group1<-as.matrix(data[1:115,]) #data of group 1
group2<-as.matrix(data[116:230,]) #data of group 2

(mu1<-colMeans(group1))  #estimate expectation of group 1
(mu2<-colMeans(group2)) #estimate expectation of group 2

s1<-matrix(0,10,10)
for (i in 1:115) {
  r<-(group1[i,]-mu1)
  s1<-s1+r%*%t(r)
}
(s1<-s1/115) #estimate covariance matrix for group 1
write.csv(s1, file = "cov_matrix1.csv") #for nice display from excel

s2<-matrix(0,10,10)
for (i in 1:115) {
  r<-(group2[i,]-mu2)
  s2<-s2+r%*%t(r)
}
(s2<-s2/115) #estimate covariance matrix for group 2
write.csv(s2, file = "cov_matrix2.csv") #for nice display from excel


sigma1<-cov(group1)*114/115    #estimate covariance of group 1 - for checking ourselves
sigma2<-cov(group2)*114/115   #estimate covariance of group 2

#-------Question 2--------
n1<-115
n2<-115
p<-10
spooled<-((n1-1)*cov(group1)+(n2-1)*cov(group2))/(n1+n2-2)
(c<-n1*n2/(n1+n2))
(tSquareSt<-t(mu1-mu2)%*%solve(spooled)%*%(mu1-mu2)*c)
(fSt<-tSquareSt*(n1+n2-p-1)/((n1+n2-2)*p))
(Fcritical<-qf(0.95,p,n1+n2-p-2))
fSt>Fcritical

#----Question 3------
#-----Task 1-----
mu1<-rep(0,10)
mu2<-rep(1,10)
s1<-matrix(0,10,10)
s2<-matrix(0,10,10)
for (i in 1:115) {
  r<-(group1[i,]-mu1)
  s1<-s1+r%*%t(r)
}
for (i in 1:115) {
  r<-(group2[i,]-mu2)
  s2<-s2+r%*%t(r)
}
s<-(s1+s2)/(n1+n2-2)
write.csv(s, file = "cov_matrix3.csv") #for nice display from excel

#----Task c-------
classify<-function(x, mu1, mu2, S){
  side1<-t(mu1-mu2)%*%solve(S)%*%x
  side2<-0.5*t(mu1+mu2)%*%solve(S)%*%(mu1-mu2)
  return(side1>side2) #function that clasify to group 1
}
ans<-rep(NA,230)
for (i in 1:230) {
  ans[i]<-classify(data[i,],mu1,mu2,s)
}
((sum(ans[1:115])+115-sum(ans[116:230]))/230) #rate of right classifications
length(ans[116:230])


#-------Question 4------
#stan_data<-scale(group1)
sigma1<-cov(group1)*114/115
eigenvalues<-eigen(sigma1)
totEigenValue<-sum(eigenvalues$values)
sum(eigenvalues$values[1:4])/totEigenValue
#values<-0
#p<-0


PCA<-prcomp(group1) #finding eigenvalues and eigenvectors of covariance matrix
eigenValues<-as.vector(PCA$sdev^2) #eigenvalues
totEigenValue<-sum(eigenValues) #total sum of eigenvalues
p<-0
val<-0
eigenVectors<-rep(NA,10)
while(val/totEigenValue<0.8) { #adding components until we define 80%
  p<-p+1
  val<-val+eigenValues[p]
  eigenVectors[p]<-list(PCA$rotation[,p])}
eigenVectors1<-matrix(NA,10,5)
for (i in 1:5) {  eigenVectors1[,i]<-eigenVectors[[i]]} #matrix T 
newdataset<-t(eigenVectors1)%*%t(group1) #new data set


#----Question 5-----
sig5<-matrix(c(1,0.8,0.8,1), 2,2)
mu1<-c(0,0)
mu2<-c(1,0)
mu3<-c(0,1)
a1<-solve(sig5)%*%(mu1-mu2)
b1<-(0.5*t(mu1-mu2)%*%solve(sig5)%*%(mu1+mu2))/a1[2]
a2<-solve(sig5)%*%(mu1-mu3)
b2<-(0.5*t(mu1-mu3)%*%solve(sig5)%*%(mu1+mu3))/a2[2]
a3<-solve(sig5)%*%(mu2-mu3)
b3<-(0.5*t(mu2-mu3)%*%solve(sig5)%*%(mu2+mu3))/a3[2]
plot(0:5,0:5, type = "n")
abline(b1,-a1[1]/a1[2],col='blue')
abline(b2,-a2[1]/a2[2],col='green')
abline(b3,-a3[1]/a3[2],col='red')


#------Question 6-------
sig5<-matrix(c(1,0.8,0.8,1), 2,2)
mu1<-c(0,0)
mu2<-c(1,0)
mu3<-c(0,1)

classifyFunc <- function(x1,x2){ #Function that classify
  if(x1 < (2.2222*x2+1.3889)/2.7778 & x1 >(2.7778*x2-1.3889)/2.2222)
    return(1)
  if(x1>=(2.2222*x2+1.3889)/2.7778 & x1> x2)
    return(2)
  if(x1<= x2 &  x1 <= (2.7778*x2-1.3889)/2.2222 )
    return(3)
}
set.seed(11)
data1<- rmvnorm(30,mu1,sig5)
data2<- rmvnorm(30,mu2,sig5)
data3<- rmvnorm(30,mu3,sig5)
sum<- 0
for (i in 1:30) { 
  if(classifyFunc(data1[i,1],data1[i,2])==1) 
     sum<- sum+1
}
for (i in 1:30) { 
  if(classifyFunc(data2[i,1],data2[i,2])==2) 
     sum<- sum+1
}
for (i in 1:30) { 
  if(classifyFunc(data3[i,1],data3[i,2])==3) 
    sum<- sum+1
}
(sum/90)
#------Question 7-------
mu1<-colMeans(data1)
mu2<-colMeans(data2)
mu3<-colMeans(data3)
sig<-(cov(data1)+cov(data2)+cov(data3))/3
a1<-solve(sig)%*%(mu1-mu2)
b1<-(0.5*t(mu1-mu2)%*%solve(sig)%*%(mu1+mu2))/a1[2]
a2<-solve(sig)%*%(mu1-mu3)
b2<-(0.5*t(mu1-mu3)%*%solve(sig)%*%(mu1+mu3))/a2[2]
a3<-solve(sig)%*%(mu2-mu3)
b3<-(0.5*t(mu2-mu3)%*%solve(sig)%*%(mu2+mu3))/a3[2]
plot(-1:5,-1:5, type = "n")
abline(b1,-a1[1]/a1[2],col='blue')
abline(b2,-a2[1]/a2[2],col='green')
abline(b3,-a3[1]/a3[2],col='red')

#-----Question 8-------
set.seed(123)
data8<-matrix(NA,100,3)
lambda1<-1
lambda2<-3
for (i in 1:100) {
  if(runif(1)<0.5)
    data8[i,]<-c(1,rexp(1,lambda1),rexp(1, lambda2))
  else
    data8[i,]<-c(2,rexp(1,lambda2),rexp(1, lambda1))
}
ans<-rep(NA,100)
for (i in 1:100) {
  if(data8[i,2]>data8[i,3])
    ans[i]<-1
  else
    ans[i]<-2
}
sum(ans==data8[,1])

set.seed(123)
data8<-matrix(NA,100,3)
lambda1<-1
lambda2<-3
for (i in 1:100) {
  if(runif(1)<0.25)
    data8[i,]<-c(1,rexp(1,lambda1),rexp(1, lambda2))
  else
    data8[i,]<-c(2,rexp(1,lambda2),rexp(1, lambda1))
}
ans<-rep(NA,100)
for (i in 1:100) {
  if((data8[i,2]-data8[i,3])>log(3,base = exp(1))/(lambda2-lambda1))
    ans[i]<-1
  else
    ans[i]<-2
}
sum(ans==data8[,1])


