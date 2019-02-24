#Poisson sim 1/15/2018
#Count data simulation

x <- rpois(n=10000, 2)
table(x)
hist(x)
mean(x)
var(x)

y <- rpois(n=10000, 3)
table(y)
hist(y)
mean(y)

df<- data.frame(event=c(x,y), group=c(rep("one", 10000), rep("two", 10000)))


summary(glm(event~group, family=poisson, data=df))
