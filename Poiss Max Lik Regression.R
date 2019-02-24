#Jan19,2018
#Poiss Maximum Likelihood Regression


#1) Generate data
  n=3000
  X1 <- runif(n, min=0, max=1)
  X2 <- runif(n, min=0, max=1)
  B0 <- 3
  B1 <- 2
  B2 <- .05
  Y <- rpois(n=n, lambda=exp(B0+B1*X1+B2*X2))
  d <- data.frame(x1=X1, x2=X2, y=Y)

s2 <- (glm(Y~X1+X2, family="poisson")) 
s1<- (glm(Y~X1, family="poisson")) 

#2) Misspecify the model
  opt.input <- function(data, theta){
    b0 <- theta[1]
    b1 <- theta[2]
    n.logl <- -sum(dpois(data$y, lambda=exp(b0+b1*data$x1), log=T))
    return(n.logl)
  }

  P1 <- optim(c(1,1), opt.input, data=d, method="BFGS")
  num <-   P1$value


#3) Correctly specify model: Poisson regression#2  MLE
  opt.input2 <- function(data, theta){
    b0 <- theta[1]
    b1 <- theta[2]
    b2 <- theta[3]
    n.logl <- -sum(dpois(data$y, lambda=exp(b0+b1*data$x1+b2*data$x2), log=T))
    return(n.logl)
  }
  
  P2<- optim(c(1,1,1), opt.input2, data=d, method="BFGS")
  num2 <- P2$value

#Likelihood Ratio Test
  #LRT=2(Log Like of full model - Log like of reduced/null model)
  #test statistic is approximately a chi-squared distribution with degrees of freedom equal to df(alt)-df(null)
  LRT <- 2*(-(num2)-(-num))
  
