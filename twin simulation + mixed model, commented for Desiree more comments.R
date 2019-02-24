library(MASS) #load all packages before starting. All these packages need to be installed (using install.packages('MASS') for example) before loaded
library(reshape2)
library(regress)
library(gap)

mzcor <- 0.6 #Set the twin correlations. These values will apear in the 'mvrnorm' function. #can play around with these. 
dzcor <- 0.5

#Data simulation; this part of the script will eventually be replaced with our real sample data.
mz <- as.data.frame(mvrnorm(100, mu=c(0, 0), Sigma=matrix(c(1,mzcor,mzcor,1), nrow=2), empirical=T)) 
					#This is the code to simulate the data. "100" is just to tell it how many observations  n to include. 'mu' will set the means, 'Sigma' is the variance-covariance matrix. The diagonal of this matrix is the variance (in this case 1 meaning that the variance is the same as the standard deviation), and when the diagonal is 1 the off diagonal values are correlations (in this case the twin correlations). The 'emprirical=T' tells the function that no error should be added so that the mean==0, SD==1, and twin correlation==mzcor. 
					#the mu=(0,0) is just arbitrary to say pull these simulation values from a population that overall has mu of 0 for variable 1 and 0 for variable 2. 
					#Sigma is building the actual matrix...it's a 4-cell matrix, and the (1,mzcor,mzcor,1) is just saying the values to put in the matrix from left to right, top to bottom. We could have used other numbers here for the standard dev. but then we'd have to standardize it before using the other numbers as correlations.
					#empirical can have two arguments (T or F....true or false (Default is F). with True as the argument, this means that it should actually use these values, and not just sample from a population that has these values on average.

dz <- as.data.frame(mvrnorm(100, mu=c(0, 0), Sigma=matrix(c(1,dzcor,dzcor,1), nrow=2), empirical=T))

mz$fam <- seq(1:dim(mz)[1]) #adds a variable to keep track of families
mz <- melt(mz, id='fam') # The 'melt' function will will convert from wide to long format, meaning each twin will be listed on a seperate row. This long format is necessary for running the mixed model below. 
mz$zyg=1 #Adds a column indicating zygosity. 
#this is just so that we can identify later, what the appropriate expected correlations are for each kin relationships. (i.e. we can later tell R that, if the zygosity=1, this means R should use the correlation for the mz. )
#I will have to add an identifier for siblings, half siblings, and cousins, not-related,....and mother (?)

dz$fam <- c((dim(dz)[1]+1):(dim(dz)[1]*2)) #Same as above but for dz twins. 
dz <- melt(dz, id='fam')
dz$zyg=2

tw <- rbind(mz,dz) #This combines mz and dz data into one data set. 
#"value", when you type "tw"fit, is the phenotype score


#Part 2: Create the matricies

tw <- tw[order(tw$fam),] #Orders the data by family id to make the different matricies look 'nice'.
family <- as.data.frame(tw$fam) #Saves a vector of family ids to be used with the apply function below. 
tw$id <- c(1:dim(tw)[1]) #Creates a new column with unique values for each twin. This info needs to be added to the matricies so the model knows what phenotype value belongs to what individual. 

a <- apply(family, 1, function(x) ifelse(tw$fam==x & tw$zyg==1, 1, ifelse(tw$fam==x & tw$zyg==2, 0.5, 0))) #creates the additive genetic relationship matrix. 
#says: when tw$zyg==1 (when the twins are monozygotic (when zygosity ID is 1)), they will be 100% correlated (1). If else, if the relationships between the first individual and the second individual is zygosity=2, then the value should be 0.5. And if there is no zygosity (from different family), then make 0. 

diag(a) <- 1 #The apply function will add 0.5 to the diagonal for dz and is corrected by this line. 
rownames(a) <-  tw$id #Adds the unique ids rows and columns. 
colnames(a) <-  tw$id

c <- apply(family, c(1), function(x) ifelse(tw$fam==x, 1, 0)) #Creates the C matrix. #shared environment.
rownames(c) <-  tw$id
colnames(c) <-  tw$id

d <- ifelse(a==0.5, 0.25,a) #Changes A to D by changing 0.5 values to 0.25 and keeping the rest


fite <- regress(value~1, data=tw) # The different models. #the E model is not taking into account any of the matrix stuff.  because it's the non-shared environment, or in this case, the error.
fitae <- regress(value~1, ~a, data=tw)
fitace <- regress(value~1, ~a+c, data=tw)
fitade <- regress(value~1, ~a+d, data=tw)


data.frame('a2'=(fitae$sigma[1]/sum(fitae$sigma)), 'e2'= (fitae$sigma[2]/sum(fitae$sigma))) # Calculates the amount of variance explained by the different relationship matricies. Check these if results look weird because of negative values. 
data.frame('a2'=(fitace$sigma[1]/sum(fitace$sigma)), 'c2'= (fitace$sigma[2]/sum(fitace$sigma)),'e2'= (fitace$sigma[3]/sum(fitace$sigma)))
data.frame('a2'=(fitade$sigma[1]/sum(fitade$sigma)), 'd2'= (fitade$sigma[2]/sum(fitade$sigma)),'e2'= (fitade$sigma[3]/sum(fitade$sigma)))


#The model comparing whether there is a significant difference between the fit of the model with additive genetics and dominance genetics, for example, vs the 
1-pchisq(2*(fitace$llik-fitae$llik),1) #Model comparison. This is what is called a log likelihood test and it could be good for you to read up on this a little bit. Basically you contrast different nested models to see the significance of a specific variance component, in this case C. 


