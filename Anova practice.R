#practice


y <- rnorm(n=100, mean=c(6,2,4,6,1), sd=1)
group <- rep(LETTERS[1:5], times=20)
sex <- rbinom(n = 100, size = 1, prob = 0.5)
var <- rnorm(n=100, mean = 0, sd=1)


df <- data.frame(y, group, sex, var)


fit <- lm(y~group+sex+var, data=df)
summary(fit)


#Redundant var, if Y then 1, otherwise 0
df$dummyA <- ifelse(df$group=="A", 1, 0)


fit1 <- lm(y~group + sex +var +dummyA -1, data=df) 
summary(fit1)

#ANOVA
# 
# y1 <- rnorm(n=16, mean=c(6,2,4,6), sd=1)
# y <- round(y)
y <- c(6,2,4,7,5,1,5,4,6,0,4,6,7,3,4,5)
group <- rep(LETTERS[1:4], times=4)
num<- rep(1:4, times=4)

df <- data.frame(y, group, num)




#Reg
fit <- lm(y~group, data=df)
summary(fit)


g <- aggregate(y~group, FUN=sd, data=df)
g$div4 <- (g$y)/sqrt(16)
g

#No int
fit.0 <- lm(y~group+0, data=df)
summary(fit.0)


fit.aov<- aov(y~group, data=df)
summary(fit.aov)


aggregate(y~group, FUN=mean, data=df)
mean(df$y)


sum(((6-mean(df$y))^2)*4,((1.5-mean(df$y))^2)*4,((4.25-mean(df$y))^2)*4,((5.5-mean(df$y))^2)*4)

#R2
R2=1-[1-()]




