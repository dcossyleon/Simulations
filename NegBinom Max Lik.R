#Jan 22, 2018
#Maximum likelihood Neg. Binomial distrib

library(MASS)
#Neg Binomial Distribution

#1) Generage data
n=1000
X1 <- runif(n, min=0, max=1) #try normal with mu 0
k=2 #DISPERSION PARAMETER
B0 <- 0
B1 <- 2
Y=rnbinom(n=n, size=k, mu=exp(B0+B1*X1))
d <- data.frame(x=X1,y=Y)
hist(d$y)

summary(glm.nb(y~x, data=d))



opt.input<- function(y, theta){
  B0 <- theta[1]
  B1 <- theta[2]
  size <- theta[3]
  nlog.L <- -sum(log(dnbinom(d$y, size=size, mu=exp(B0+B1*X1), log=FALSE)))
  return(nlog.L)
}

ob <- optim(c(.5, .5, .5), opt.input, y=d, method="BFGS")
ob1 <- optim(c(.5, .5, .5), opt.input, y=d, method="BFGS")


#Why it doesn't like some of the Xs; -->R doesn't like too big of numbers for NBinom or Poisson
  
#why am i getting errors for dnbinom
  #try to get a df with the warning messages

#find a wordy version of optim's iterations output
  
  