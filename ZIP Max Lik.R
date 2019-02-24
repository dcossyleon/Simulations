#ZIP Max Lik

#1) Generate data
n=1000
X1 <- runif(n, min=0, max=1)
B0 <- -0.5
B1 <- 2
Y=ifelse(rbinom(n=n, size = 1, prob = (exp(B0+B1*X1)/(1+(exp(B0+B1*X1)))))==0, 0, rpois(n=n, lambda=exp(B0+B1*X1)))
d <- data.frame(x=X1,y=Y)
hist(d$y)


zip <- function(d,theta){
  zB0 <- theta[1] #zero mass
  zB1 <- theta[2]
  cB0 <- theta[3] #count process
  cB1 <- theta[4] 
  lambda=exp(cB0+cB1*d$x)
  mu=exp(zB0+zB1*d$x)/(1+(exp(zB0+zB1*d$x)))
  log0 <- log(mu+(1-mu)*exp(-lambda)) #zero mass
  logP  <- log((1-mu)*dpois(d$y, lambda=lambda, log=FALSE)) #count process
  z <- d$y==0
  c <- d$y>0
  logLik <- sum(log0[z])+sum(logP[c]) #compare to my vector of logical values, and only give me back the ones that are TRUE. banking on position being the same.
  return(-logLik)
}


op1 <- optim(c(0,0,0,0), zip, d=d, method="BFGS")
op1
#if optim won't run, the estimates may just be way too big for R to handle
#try changing to smaller parameters
#has to be run in the same optim, because each function is not mutually exclusive of the others

#GLM for ZIP
summary(zeroinfl(y~x|x, link="logit", dist="poisson", data = d))


