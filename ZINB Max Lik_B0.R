#Jan 30, 2018 ZINB Model Simulation
#Mixed Model

library(pscl)
library(MASS)

#1) Generate data
n=1000
X1 <- runif(n, min=0, max=1)
B0 <- -0.5
B1 <- 2
size=3 #dispersion parameter

Y=ifelse(rbinom(n=n, size = 1, prob = (exp(B0+B1*X1)/(1+(exp(B0+B1*X1)))))==0, 0, rnbinom(n=n, size=size, mu=exp(B0+B1*X1)))
d <- data.frame(x=X1,y=Y)
hist(d$y)


d.no0 <- d
d.no0$y <- ifelse(d.no0$y>0, d.no0$y, NA)
d.0 <- d
d.0$y <- ifelse(d.0$y!=0, NA, 0)


#2) prob of getting 0 counts, 
zinb<- function(d, theta){
  B0 <- theta[1]
  B1 <- theta[2]
  size <- theta[3]
  Z0 <- theta[4]
  Z1 <- theta[5]
  
  mu=exp(B0+B1*d$x) #count process
  false.0 <- (exp(Z0+Z1*d$x)/(1+(exp(Z0+Z1*d$x))))#zero mass
  count.0 <- (size/(mu+size))^size
  logC <- log((1-false.0)*dnbinom(d$y, size=size, mu=mu, log=F))
  log0 <- log(false.0+(1-false.0)*count.0)
  z <- d$y==0
  c <- d$y>0
  logLik <- sum(log0[z])+sum(logC[c])
  return(-logLik)
}

op1 <- optim(c(3,1,2,0,0), zinb, d=d, method="BFGS")
op1


#GLM for ZINB
summary(zeroinfl(y~x|x, link="logit", dist="negbin", data = d)) #theta=dispersion parameter
