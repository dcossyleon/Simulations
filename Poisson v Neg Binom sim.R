#1/15/2018
#Poisson vs Neg Binomial Sim

library(MASS)

random <- function(x) {
  x <- rpois(n=100, 2)
  y <- rpois(n=100, 2.3)
  df<- data.frame(event=c(x,y), group=c(rep("one", 100), rep("two", 100)))
  p1 <- summary(glm(event~group, family=poisson, data=df))$coefficient[2,4]
  p2 <- summary(glm.nb(event~group, data=df))$coefficient[2,4]
  p <- cbind(p1,p2)
  return(p)
}
random()

sim <- do.call(rbind, lapply(1:500, function(x) random()))

table(sim[,1]<0.05)

