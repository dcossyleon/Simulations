#Bootstrapping: read something formal

n=1000
size=1
X1 <- runif(n, min=0, max=1)
B0 <- 1
B1 <- .1
Y <- rbinom(n=n, size=size, prob=(exp(B0+B1*X1)/(1+(exp(B0+B1*X1))))) 
df <- data.frame(x=X1,y=Y)

summary(glm(y~x, data=df, family=binomial)) #Generate data, and then run a glm


boot <- function(x){
  new <- df[sample(nrow(df), replace=T),]  #sample from your data with replacement
  int<- summary(glm(y~x, family = binomial, data=new))$coefficient[1]
  b1<- summary(glm(y~x, family = binomial, data=new))$coefficient[2]
  d <- data.frame(int,b1)
  return(d)
}


summary(glm(y~x, data=df, family=binomial))


practice <- do.call(rbind, lapply(1:10, function(x) boot()))

mean(practice$b1) #should be equivalent to the glm output's estimate
hist(practice$b1) #significant estimates should not overlap zero

#95CI

quantile(x=practice$b1,probs=c(.025,.975)) #includes 0, so point estimate is not significant
