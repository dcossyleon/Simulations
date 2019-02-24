#Jan 22, 2018
#Maximum likelihood Binomial distrib

#1) Generate binomial data
n=1000
size=1
X1 <- runif(n, min=0, max=1)
B0 <- 1
B1 <- 2
Y <- rbinom(n=n, size=size, prob=(exp(B0+B1*X1)/(1+(exp(B0+B1*X1))))) #not good
Y2 <-rbinom(n=n, size=size, prob=(1/(1+exp(-(B0+B1*X1)))))
Y3 <-rbinom(n=n, size=size, prob=(1/(1+exp(-(Z)))))

Z <- B0+B1*X1
PR <- 1/(1+(exp(-Z)))
Y1 <- rbinom(n, 1, PR)

d <- data.frame(x1=X1, y=Y)
hist(d$y)

#GLM
summary(glm(Y~X1, family = binomial))
summary(glm(Y2~X1, family = binomial))
summary(glm(Y3~X1, family = binomial))
summary(glm(Y1~X1, family = binomial))



#log likelihood
opt.input.binom <- function(data, theta){
  b0 <- theta[1]
  b1 <- theta[2]
  nlog.L <- -sum(dbinom(data$y, size=size, prob=(exp(b0+b1*data$x)/(exp(b0+b1*data$x)+1)), log=T))
  return(nlog.L)
}

optim(c(0,1),opt.input.binom, data=d, method="BFGS")




