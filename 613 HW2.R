#Exercise1 Data creation
x1<-as.vector(runif(10000,1,3))# Generate a vector of a uniform distribution
x2<-as.vector(rgamma(10000,3,2))# Generate a vector of a gamma distribution
x3<-as.vector(rbinom(10000,1,0.3))#Generate a vector of a binomial distribution
eps<-as.vector(rnorm(10000,2,1))# Generate a vector of a normal distribution
Y<-as.vector(0.5+1.2*x1-0.9*x2+0.1*x3+eps)# Generate Y variable
ydum<-as.vector(rep(0,10000)) # Generate a vector of 10000 zeros
for (i in 1:10000) {
  if( Y[i]> mean(Y)){
    ydum[i]=1
  }
}# Generate ydum.


#Exercise2 OLS
#2.1
cor.test(Y,x1)
# The correlation between Y and x1 is 0.4967823. 
#The difference between the correlation and 1.2 is because Y are also affected by x2 and x3.
#2.3
x0<-as.vector(rep(1,10000))
X<-matrix(c(x0,x1,x2,x3),10000,4)
beta<-as.vector(solve(t(X)%*%X)%*%t(X)%*%Y)
# The coefficients on this regression
#2.4
#Caculate the standard errors
#Using the formula
beta<-matrix(c(beta),4,1)
res<-Y-X%*%beta
s2<-t(res)%*%res/9995
#Using bootstrap with 49 and 499 replications respectively
B<-49
betasample<-as.data.frame(matrix(numeric(0),ncol = 4))
for (i in 1:B) {
  data <- cbind(Y,X)
  sampleX<- data[sample(10000,10000,replace = T),]
  beta<-matrix(c(solve(t(sampleX[,2:5])%*%sampleX[,2:5])%*%t(sampleX[,2:5])%*%sampleX[,1]),1,4)
  betasample<-rbind(betasample,beta)
}
seboot1<-data.frame(apply(betasample, 2,sd))

B<-499
betasample<-as.data.frame(matrix(numeric(0),ncol = 4))
for (i in 1:B) {
  data <- cbind(Y,X)
  sampleX<- data[sample(10000,10000,replace = T),]
  beta<-matrix(c(solve(t(sampleX[,2:5])%*%sampleX[,2:5])%*%t(sampleX[,2:5])%*%sampleX[,1]),1,4)
  betasample<-rbind(betasample,beta)
}
seboot2<-apply(betasample, 2,sd)
#Exercise 3
#3.1
ydum<-matrix(ydum,10000,1)
beta_func<-function(beta_probit){
  return(sum(ydum*log(pnorm(X%*%beta_probit)))+sum((1-ydum)*log(1-pnorm(X%*%beta_probit))))
}# Generate the likelihood function
# test:a = matrix(c(0.5, 0.5, 0.5, 0.3), 4, 1)
#  beta_func(a)
#3.2
i<-matrix(c(0.5, 0.5, 0.5, 0.3), 4, 1)
a<-diag(0.00001,4,4)
i_new<-matrix(c(0.49999,0.49999,0.49999,0.29999),4,1)#Generate the start point of beta
while (beta_func(i)<beta_func(i_new)) {
  i<-i_new
  i_combo<-matrix(c(rep(i,4)),4,4)
  i_combo_new<-i_combo+a
  b0<-(beta_func(i_combo_new[,1])-beta_func(i_combo[,1]))/0.00001
  b1<-(beta_func(i_combo_new[,2])-beta_func(i_combo[,2]))/0.00001
  b2<-(beta_func(i_combo_new[,3])-beta_func(i_combo[,3]))/0.00001
  b3<-(beta_func(i_combo_new[,4])-beta_func(i_combo[,4]))/0.00001
  dk<-matrix(c(b0,b1,b2,b3),4,1)
  i_new<-i+0.00001*dk
}
print(i)
# The i when the loop stops is the maximum likelihood
#3.3
#The difference between the result from 3.2 and beta is in the coeffience of x0=0.5
# There is no difference among other coefficients

#Exercise 4
probit<-glm(ydum~0+X, family = binomial(link = "probit"))
logit<-glm(ydum~0+X, family = binomial(link = "logit"))
linear<-lm(ydum~0+X)
beta_probit<-as.vector(probit$coefficients)
beta_logit<-as.vector(logit$coefficients)
beta_linear<-as.vector(linear$coefficients)
beta_compare<-data.frame(beta_probit,beta_logit,beta_linear)
summary(probit)
summary(logit)
summary(linear)
# There are literally no difference of coefficients between logit and probit models.

#Exercise 5
#5.1
Xbeta<-mean(X%*%beta_probit)
fideriv<-(pnorm(Xbeta+0.00001)-pnorm(Xbeta))/0.00001
margief_probit<-data.frame(fideriv*beta_probit)

# The outcome of marginal effects of probit models
Xbeta2<-mean(X%*%beta_logit)
logit_func<-exp(Xbeta2)/(1+exp(Xbeta2))
logit_func1<-exp(Xbeta2+0.00001)/(1+exp(Xbeta2+0.00001))
fideriv2<-(logit_func1-logit_func)/0.00001
margief_logit<-data.frame(fideriv2*beta_logit)
# The outcome of marginal effects of logit models
#5.2
#Using delta method to compute the standard deviations.
V_probit <- vcov(probit)
beta_probit_new<-data.frame(beta_probit,beta_probit,beta_probit,beta_probit)
beta_probit_new<-beta_probit_new+diag(0.00001,4)
M1 <- function(beta_probit) mean(dnorm(X%*%beta_probit)*beta_probit[1])
g11<-(M1(beta_probit_new[,1])-M1(beta_probit))/0.00001
g12<-(M1(beta_probit_new[,2])-M1(beta_probit))/0.00001
g13<-(M1(beta_probit_new[,3])-M1(beta_probit))/0.00001
g14<-(M1(beta_probit_new[,4])-M1(beta_probit))/0.00001
g1<-matrix(c(g11,g12,g13,g14),1,4)
# The first row of Jacobian matrix
M2 <- function(beta_probit) mean(dnorm(X%*%beta_probit)*beta_probit[2])
g21<-(M2(beta_probit_new[,1])-M2(beta_probit))/0.00001
g22<-(M2(beta_probit_new[,2])-M2(beta_probit))/0.00001
g23<-(M2(beta_probit_new[,3])-M2(beta_probit))/0.00001
g24<-(M2(beta_probit_new[,4])-M2(beta_probit))/0.00001
g2<-matrix(c(g21,g22,g23,g24),1,4)
M3 <- function(beta_probit) mean(dnorm(X%*%beta_probit)*beta_probit[3])
g31<-(M3(beta_probit_new[,1])-M3(beta_probit))/0.00001
g32<-(M3(beta_probit_new[,2])-M3(beta_probit))/0.00001
g33<-(M3(beta_probit_new[,3])-M3(beta_probit))/0.00001
g34<-(M3(beta_probit_new[,4])-M3(beta_probit))/0.00001
g3<-matrix(c(g31,g32,g33,g34),1,4)
M4 <- function(beta_probit) mean(dnorm(X%*%beta_probit)*beta_probit[4])
g41<-(M4(beta_probit_new[,1])-M4(beta_probit))/0.00001
g42<-(M4(beta_probit_new[,2])-M4(beta_probit))/0.00001
g43<-(M4(beta_probit_new[,3])-M4(beta_probit))/0.00001
g44<-(M4(beta_probit_new[,4])-M4(beta_probit))/0.00001
g4<-matrix(c(g41,g42,g43,g44),1,4)

J_probit<-rbind(g1,g2,g3,g4)
delta_probit<-t(J_probit)%*%V_probit%*%J_probit# The result of probit model
sd_probit<-diag(delta_probit)
V_logit <- vcov(logit)
beta_logit_new<-data.frame(beta_logit,beta_logit,beta_logit,beta_logit)
beta_logit_new<-beta_logit_new+diag(0.00001,4)
M1 <- function(beta_logit) mean(dlogis(X%*%beta_logit)*beta_logit[1])
g11<-(M1(beta_logit_new[,1])-M1(beta_logit))/0.00001
g12<-(M1(beta_logit_new[,2])-M1(beta_logit))/0.00001
g13<-(M1(beta_logit_new[,3])-M1(beta_logit))/0.00001
g14<-(M1(beta_logit_new[,4])-M1(beta_logit))/0.00001
g1<-matrix(c(g11,g12,g13,g14),1,4)
# The first row of Jacobian matrix
M2 <- function(beta_logit) mean(dlogis(X%*%beta_logit)*beta_logit[2])
g21<-(M2(beta_logit_new[,1])-M2(beta_logit))/0.00001
g22<-(M2(beta_logit_new[,2])-M2(beta_logit))/0.00001
g23<-(M2(beta_logit_new[,3])-M2(beta_logit))/0.00001
g24<-(M2(beta_logit_new[,4])-M2(beta_logit))/0.00001
g2<-matrix(c(g21,g22,g23,g24),1,4)
M3 <- function(beta_logit) mean(dlogis(X%*%beta_logit)*beta_logit[3])
g31<-(M3(beta_logit_new[,1])-M3(beta_logit))/0.00001
g32<-(M3(beta_logit_new[,2])-M3(beta_logit))/0.00001
g33<-(M3(beta_logit_new[,3])-M3(beta_logit))/0.00001
g34<-(M3(beta_logit_new[,4])-M3(beta_logit))/0.00001
g3<-matrix(c(g31,g32,g33,g34),1,4)
M4 <- function(beta_logit) mean(dlogis(X%*%beta_logit)*beta_logit[4])
g41<-(M4(beta_logit_new[,1])-M4(beta_logit))/0.00001
g42<-(M4(beta_logit_new[,2])-M4(beta_logit))/0.00001
g43<-(M4(beta_logit_new[,3])-M4(beta_logit))/0.00001
g44<-(M4(beta_logit_new[,4])-M4(beta_logit))/0.00001
g4<-matrix(c(g41,g42,g43,g44),1,4)

J_logit<-rbind(g1,g2,g3,g4)
delta_logit<-t(J_logit)%*%V_logit%*%J_logit
sd_logit<-diag(delta_logit)


#Using bootstrap to compute the standard errors
B<-499
mesample<-as.data.frame(matrix(numeric(0),ncol = 4))
for (i in 1:B) {
  data <- cbind(ydum,X)
  sampleX<- data[sample(10000,10000,replace = T),]
  probitboot<-glm(sampleX[,1]~0+sampleX[,2:5], family = binomial(link = "probit"))
  coeboot<-probitboot$coefficients
  meanboot<-mean(sampleX[,2:5]%*%coeboot)
  fiderivboot<-(pnorm(meanboot+0.00001)-pnorm(meanboot))/0.00001
  margief<-matrix(fideriv*coeboot,1,4)
  mesample<-rbind(mesample, margief)
}
mead<-apply(mesample,2,sd)
# The standard deviations of probit model
B<-499
mesample<-as.data.frame(matrix(numeric(0),ncol = 4))
for (i in 1:B) {
  data <- cbind(ydum,X)
  sampleX<- data[sample(10000,10000,replace = T),]
  probitboot<-glm(sampleX[,1]~0+sampleX[,2:5], family = binomial(link = "logit"))
  coeboot<-probitboot$coefficients
  meanboot<-mean(sampleX[,2:5]%*%coeboot)
  fiderivboot<-(pnorm(meanboot+0.00001)-pnorm(meanboot))/0.00001
  margief<-matrix(fideriv*coeboot,1,4)
  mesample<-rbind(mesample, margief)
}
mead2<-apply(mesample,2,sd)
# The standard deviations of logit model



