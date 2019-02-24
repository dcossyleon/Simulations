#Jan 19, 2018
#Maximum likelihood normal distrib

#Normal Distribution
data <- rnorm(n=100, mean=50, sd=30)


opt.input <- function(data, theta){
  b0 <- theta[1]
  b1 <- theta[2]
  std <- theta[3]    
  n.logl <- -sum(dnorm(data$y, mean=(b0+b1*data$x), sd=std, log=T))
  return(n.logl)
}

n=10000
x <- runif(n, min=0, max=100)
B0 <- 0
B1 <- 5
STD <- 2
y=B0+B1*x+rnorm(n,0,STD)
  
summary(lm(y~x)) #it works!
d <- data.frame(x=x,y=y)

optim(c(1,1,1), opt.input, data=d, method="BFGS")







opt.input<- function(y, theta){
  mean <- theta[1]
  sd <- theta[2]
  nlog.L <- -sum(dnorm(y, mean, sd, log=TRUE))
  return(nlog.L)
}

optim(c(44,2.5), opt.input, y=data, method="BFGS")


#__

normal.lik1 <- function(theta,y){
  mu <- theta[1]
  sigma2 <- theta[2]
  n <- length(y)
  log1 <- -.5*n*log(2*pi)-0.5*n*log(sigma2) - (1/(2*sigma2))*sum((y-mu)**2)
  return(-log1)
}  #sigma2= variance

optim(c(30,2), normal.lik1, y=data, method="BFGS")


normal.lik2 <- function(theta,y){
  mu <- theta[1]
  sigma <- theta[2]
  n <- length(y)
  z <- (y-mu)/sigma
  logl <- -n*log(sigma)- sum(log(dnorm(z)))
  return(-logl)
}

optim(c(20,2.5), normal.lik2, y=data, method="L-BFGS-B")

#__



