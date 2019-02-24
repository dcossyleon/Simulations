#Jan 30, 2018 ZINB Model Simulation
#Mixed Model

library(pscl)
library(MASS)

#1) Generate data
n=10000
X1 <- runif(n, min=0, max=1)
B0 <- 4
B1 <- 2
size=3 #DISPERSION PARAMETER
C0 <- 1.6
C1 <- 5

Y=ifelse(rbinom(n=n, size = 1, prob = (exp(C0+C1*X1)/(1+(exp(C0+C1*X1)))))== 0, 0, rnbinom(n=n, size=size, mu=exp(B0+B1*X1)))
d <- data.frame(x=X1,y=Y)
hist(d$y)

#GLM for ZINB
summary(zeroinfl(y~x|x, link = "logit", dist="negbin", data = d)) #theta= dispers.

#glm binom
d.bin <- d
d.bin$y <- ifelse(d.bin$y>0, 1, 0)
summary(glm(y~x, family = binomial, data=d.bin))
glm.bin<- glm(y~x, family = binomial, data=d.bin)
logLik(glm.bin)
#glm neg binom
summary(glm.nb(y~x, data=d))


#2) prob of getting 0 counts, 
zero.mass<- function(y, theta){
  B0 <- theta[1]
  B1 <- theta[2]
  size <- theta[3]
  C0 <- theta[4]
  C1 <- theta[5]
  mu=exp(B0+B1*X1)
  false.0 <- (exp(C0+C1*X1)/(1+(exp(C0+C1*X1))))
  count.0 <- (size/(mu+size))^size
  nlog.L <- -(sum(log(false.0+(1-false.0)*count.0)))
  return(nlog.L)
}

#3) then prob of getting non-zero counts
count.process<- function(y, theta){
  B0 <- theta[1]
  B1 <- theta[2]
  size <- theta[3]
  C0 <- theta[4]
  C1 <- theta[5]
  mu=exp(B0+B1*X1)
  false.0 <- (exp(C0+C1*X1)/(1+(exp(C0+C1*X1))))
  nlog.L <- -sum(log((1-false.0)*dnbinom(d$y, size=size, mu=exp(B0+B1*X1), log=F)))
  return(nlog.L)
}

op1 <- optim(c(.5,.5,.5,.5,.5), zero.mass, y=d, method="BFGS")
op2 <- optim(c(4,2,3,1.6,5), count.process, y=d, method="BFGS")

