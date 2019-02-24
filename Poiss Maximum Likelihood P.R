#Jan 18, 2018

#Poisson Density Function
pois <- function(y, lambda) {
  ((lambda^(y))*exp(1)^(-lambda))/factorial(y)
} #What the prob of getting a "y"i exactly, if our lambda =__?

data <- rpois(n=100, lambda=5)

#Poisson Likelihood function
  poisson.lik<-function(mu,y){
    n<-length(y)
    logl <-sum(y)*log(mu)-n*mu
    return(-logl)
  }
  
  logl <- sum(y*log(mu)-mu-log(factorial(y))) #can also use this as logl
  optim(0.1,poisson.lik, y=data, method="BFGS") #Optim $value is NEGATIVE MAX LOG LIKELIHOOD


#Use this, Hasse's much faster way
  opt.input<- function(lambda, y){
    nlog.L <- -sum(dpois(y, lambda, log=TRUE))
    return(nlog.L)
  }
  
  optim(0.1,opt.input, y=data, method="BFGS")

#GLM
  summary(glm(data~1, family="poisson"))

#How to use:
  #- given a data set, it's going to run through every possible value of lambda, and see which lambda value is most probable to fit the data
  #- for all possible lambdas (x axis), which has the highest likelihood (y axis) to have produced the data we observe? That's the MLE.
  #- Wrong way: Given the mean of our data, what's the probability of getting value x? 
#^^I know this is wrong way to think about it, but would these give you the same answer? Are they equivalent?



#Finding MLE without Optim: Take L function, have it go from 0 to infinity(?) and report the max value...
#then find the corresponding x value for that max y value
data <- rpois(n=100, lambda=5)

pois <- function(y, lambda) {
  ((lambda^(y))*exp(1)^(-lambda))/factorial(y)
}

zp <- function(y,)

L <- do.call(rbind, lapply(0.1:100, function(x) cbind(x,prod(pois(data, x))))) #returns all Ls
ML <- max(L[,2]) #returns Maximum Likelihood value
MLE <- subset(L[,1],L[,2]==ML) #returns Max Lik Estimate for lambda
MLE


#Use this
opt.input<- function(lambda, y){
   nlog.L <- -sum(dpois(y, lambda, log=TRUE))
  return(nlog.L)
}
optim(0.1,opt.input, y=data, method="BFGS")

#GLM
summary(glm(data~1, family="poisson"))
