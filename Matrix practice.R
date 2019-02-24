#Learning matrices and regression

x <- rnorm(100, mean=0, sd=1)
y= 1+ 5*x + rnorm(100, mean=0, sd=1)
ones <- rep(1, times=100)

y<- as.matrix(y)
d <- as.matrix(data.frame(ones, x))

b <- (solve(t(d)%*%(d)))%*%(t(d)%*%y)

lm(y~x)
summary(lm(y~x))

#fitted values
d%*%b


#Multiply d by itself
d%*%t(d) #gives me the squared diag, product off diags
t(d)%*%d #gives me: N, sum of x, and sum of squared x's

t(d)%*%y #gives me sum of y and the sum of x*y



#Sum of y
  s <- t(y)%*%ones

#mult by 1/n
  mu <- s*(1/length(y))
  
#repeat the mean value n times
I <- as.matrix(diag(100))
ones.square <- ones%*%t(ones)



#matrix practice
fam5 <- matrix(c(1, 2, 3,
                 2, 3, 2,
                 1, 1, 1), nrow=3)

dz5 <- matrix(c(1, 2,
                2, 3,
                1, 1), nrow=3)
              
fam5%*%t(fam5)
fam5







 
 
vec <- matrix(c(1, 2, 2, 3), nrow=4)

vec%*%t(vec)

t(y)%*%I # gives me back y
I%*%y # also gives me back y


t(y)%*%ones.square #gives me back the Sum Square in every cell. not ideal

#Vector of Sum of all y
ones.square%*%y


#to get mean, divide by 1/n
ones.square%*%y*(1/(length(y)))

#Deviation score
y-(ones.square%*%y*(1/(length(y))))

#sum and square the devs
t(y-(ones.square%*%y*(1/(length(y)))))%*%(y-(ones.square%*%y*(1/(length(y)))))

#Var
(t(y-(ones.square%*%y*(1/(length(y)))))%*%(y-(ones.square%*%y*(1/(length(y))))))/(length(y)-1)
var(y) #check w R

#Sd
sqrt((t(y-(ones.square%*%y*(1/(length(y)))))%*%(y-(ones.square%*%y*(1/(length(y))))))/(length(y)-1))

sd(y) #dbl check with R




#Covariance of x and y

#Deviance scores of x
#x-[(square matrix of ones)*(x)]
dev.x <- x-(ones.square%*%x*(1/(length(x))))

#Deviance scores of y
dev.y <- y-(ones.square%*%y*(1/(length(y))))

#Multiplying these two matrices will give me the sum of the products (the order does not matter)
t(dev.x)%*%dev.y
t(dev.y)%*%dev.x #same (the order does not matter)

#divide by N-1
(t(dev.y)%*%dev.x)/(length(y)-1) #ANSWER!

cov(x,y) #confirm with r's answer



#Correlation, r_xy= covariance/(s_x*s_y)

((t(dev.y)%*%dev.x)/(length(y)-1))/
  ((sqrt((t(y-(ones.square%*%y*(1/(length(y)))))%*%(y-(ones.square%*%y*(1/(length(y))))))/(length(y)-1)))*
  (sqrt((t(x-(ones.square%*%x*(1/(length(x)))))%*%(x-(ones.square%*%x*(1/(length(x))))))/(length(y)-1))))



combo <- cbind(dev.x,dev.y)
(t(combo)%*%combo)*(1/(length(y)-1))

#check with 
var(combo)
