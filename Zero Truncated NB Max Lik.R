#Jan 25, 2018 Zero Truncated NB Model Simulation

library(aster)
library(VGAM)
stats::resid #is this right?


#1) Generate data
  n=1000
  X1 <- runif(n, min=0, max=1)
  size=3 #DISPERSION PARAMETER
  k=0 #truncation limit
  B0 <- 4
  B1 <- 2
  Y=rktnb(n=n, size=size, k=k,  mu=exp(B0+B1*X1), xpred=1)  #k is value for truncation
  d <- data.frame(x=X1,y=Y)
  hist(d$y)

summary(vglm(y~x, family=posnegbinomial, data=d)) #how do i find size?, *Intercept2 is junk?

#2)
  opt.input<- function(y, theta){
    B0 <- theta[1]
    B1 <- theta[2]
    size <- theta[3]
    mu=exp(B0+B1*X1)
    p0 <-(size/(mu+size)^size)
    nlog.L <- -(sum(dnbinom(d$y, size=size, mu=exp(B0+B1*X1),log=TRUE)-(log(1-p0))))
    return(nlog.L)
  }
        
optim(c(1,1,1), opt.input, y=d, method="CG")
