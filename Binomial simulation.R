#binomial simulation 

x <- rbinom(n=10000, size=1, prob=.3)
table(x)

y <- rbinom(n=10000, size=1, prob=.6)
table(y)

z <- rbinom(n=10000, size=1, prob=.85)


df<- data.frame(event=c(x,y), group=c(rep("one", 10000), rep("two", 10000)))
summary(glm(event~group, family=binomial, data=df))


df1<- data.frame(event=c(x,y,z), group=c(rep("one", 10000), rep("two", 10000), rep("three", 10000)))
summary(glm(event~group, family=binomial, data=df1))


df2<- data.frame(event=c(x,y,z), group=c(rep(1, 10000), rep(2, 10000), rep(3, 10000)))
M2 <- glm(event~group, family=binomial, data=df2)
summary(glm(event~group, family=binomial, data=df2))


#Plot the sigmoid
MyData <- data.frame(group = seq(from = 0, to = 3,by = 1))
Pred <- predict(M2, newdata = MyData, type = "response")
plot(y = df2$event, x= df2$group)
lines(MyData$group, Pred)
