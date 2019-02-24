library(MASS) #load all packages before starting. All these packages need to be installed (using install.packages('MASS') for example) before loaded
library(reshape2)
library(regress)
library(gap)

mzcor <- 0.6 #Set the twin correlations. These values will apear in the 'mvrnorm' function. 
dzcor <- 0.5

mz <- as.data.frame(mvrnorm(100, mu=c(0, 0), Sigma=matrix(c(1,mzcor,mzcor,1), nrow=2), empirical=T)) #This is the code to simulate the data. 'mu' will set the means, 'Sigma' is the variance-covariance matrix. The diagonal of this matrix is the variance (in this case 1 meaning that the variance is the same as the standard deviation), and when the diagonal is 1 the off diagonal values are correlations (in this case the twin correlations). The 'emprirical=T' tells the function that no error should be added so that the mean==0, SD==1, and twin correlation==mzcor. 
dz <- as.data.frame(mvrnorm(100, mu=c(0, 0), Sigma=matrix(c(1,dzcor,dzcor,1), nrow=2), empirical=T))

mz$fam <- seq(1:dim(mz)[1]) #adds a variable to keep track of families
mz <- melt(mz, id='fam') # The 'melt' function will will convert from wide to long format, meaning each twin will be listed on a seperate row. This long format is necessary for running the mixed model below. 
mz$zyg=1 #Adds a column indicating zygosity. 

dz$fam <- c((dim(dz)[1]+1):(dim(dz)[1]*2)) #Same as above but foro dz twins. 
dz <- melt(dz, id='fam')
dz$zyg=2

tw <- rbind(mz,dz) #This combines mz and dz data into one data set. 




tw <- tw[order(tw$fam),] #Orders the data by family id to make the different matricies look 'nice'.
family <- as.data.frame(tw$fam) #Saves a vector of family ids to be used with the apply function below. 
tw$id <- c(1:dim(tw)[1]) #Creates a new column with unique values for each twin. This info needs to be added to the matricies so the model knows what phenotype value belongs to what individual. 

a <- apply(family, 1, function(x) ifelse(tw$fam==x & tw$zyg==1, 1, ifelse(tw$fam==x & tw$zyg==2, 0.5, 0))) #creates the additive genetic relationship matrix. 
diag(a) <- 1 #The apply function will add 0.5 to the diagonal for dz and is corrected by this line. 
rownames(a) <-  tw$id #Adds the unique ids rows and columns. 
colnames(a) <-  tw$id

c <- apply(family, c(1), function(x) ifelse(tw$fam==x, 1, 0)) #Creates the C matrix.
rownames(c) <-  tw$id
colnames(c) <-  tw$id

d <- ifelse(a==0.5, 0.25,a) #Changes A to D by changing 0.5 values to 0.25 and keeping the rest


fite <- regress(value~1, data=tw) # The different models 
fitae <- regress(value~1, ~a, data=tw)
fitace <- regress(value~1, ~a+c, data=tw)
fitade <- regress(value~1, ~a+d, data=tw)
data.frame('a2'=(fitae$sigma[1]/sum(fitae$sigma)), 'e2'= (fitae$sigma[2]/sum(fitae$sigma))) # Calculates the amount of variance explained by the different relationship matricies. Check these if results look weird because of negative values. 
data.frame('a2'=(fitace$sigma[1]/sum(fitace$sigma)), 'c2'= (fitace$sigma[2]/sum(fitace$sigma)),'e2'= (fitace$sigma[3]/sum(fitace$sigma)))
data.frame('a2'=(fitade$sigma[1]/sum(fitade$sigma)), 'd2'= (fitade$sigma[2]/sum(fitade$sigma)),'e2'= (fitade$sigma[3]/sum(fitade$sigma)))

1-pchisq(2*(fitace$llik-fitae$llik),1) #Model comparison. This is what is called a log likelihood test and it could be good for you to read up on this a little bit. Basically you contrast different nested models to see the significance of a specific variance component, in this case C. 


