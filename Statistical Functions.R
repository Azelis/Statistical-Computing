#Justinas Azelis
#MATH6166
#Statistical Computing

##### Random walk Metropolis-Hastings algorithm #####
### 1-A part###----------------------------------

set.seed(100)

#Creating function to trasfer value to logarithmic value. Will be used
#for rwmet function

#rwmet - Random walk Metropolis-Hastings algorithm
  
gValueLog <- function(x,feta){
  g_value_log <- feta * x - exp(feta*x)
  g_value_log
}

##Writing the main function for A part with arguments M-sample size and
#phi - quantity phi controlling the variance of the proposal distribution
rwmet <-function(M, phi){
  
  
  feta <-6#Student ID
  
  #Creating variables
  x <- c()
  y <- c()
  alfa <- c()
  miu <- c()
  
  #C part
  acceptance = 0
  
  ##Simulation part
  x[1] <- 1 #As starting value0
  
  #x0 is not included in the chain of values, so x0 will be excluded
  #M+1, because first value will be excluded and I need to have M sample size
  for(i in 1:(M+1)){
    
    #(a), normal distribution with mean x[i-1 (in this case I used as i, because
    # all the values in my list are moved + 1 in order to exclude X0 from the
    #chain of values)] and variance exp(phi)
    y <- rnorm(1,mean = x[i], sd = sqrt(exp(phi)))
    
    #(b), Calculating acceptence probability
    #Equation g_x from the hint. gValueLog is a function, to calculate
    #logarithmic values instead of log(), because log() result is not so
    #accurate as gValueLog()
    g_x <- exp(gValueLog(y,feta) - gValueLog(x[i],feta))
    
    #Getting acceptence probability
    alfa <- min(1,g_x)
    
    #(c)
    #Getting 1 uniform distribution in the interval 0-1
    miu <- runif(1,0,1)
    
    #Comparing the values if I need to accept the proposal or to reject it based
    #on alfa from (b) part and uniform distribution
    if (alfa >= miu){
      acceptance <- acceptance + 1 #For C part to calculate total number of
      #acceptance
      x[i+1] <- y
    }else{
      x[i+1] <- x[i]
    }
  }
  
  #Result the values of 1 - vector of samples and 2-acceptance rate
  list(x[2:(M+1)],acceptance/M)
}
#Result of rwmet is a list of [[1]] - B, D parts, vector of samples;
#[[2]] - C part acceptance rate


### 1-B part###----------------------------------
#Giving values based on the question
N <- 100#function to generate
M <- 10000#sample times
phi <- -4 #phi
n_means <- c()

#Calculating expectation N times rwmet[[1]]
for (j in 1:N){  
  n_means[j] <- mean(rwmet(M,phi)[[1]])
}

#Calculating the variance of expectation of samples
var_B <- var(n_means)
var_B#0.0001178915


### 1-C part###----------------------------------

#sample of size
M <- 10000

#Calcualting equally-spaced sequence of 100 times from -5 to 10
phi_acc <- seq(from = -5, to = 10, length.out = 100)

#Calculating acceptance rate rwmet[[2]]
acc_result <- 0
for (i in 1:100){
  acc_result[i] <- rwmet(M,phi_acc[i])[[2]]
}

#Creating a plot to show acceptance and line under certain condition
#the optimal acceptance rate for a random walk Metropolis-Hastings algorithm
plot(x = phi_acc, y = acc_result, main="Acceptance rate vs phi values", 
     ylab = "Acceptance rate", xlab = "phi values" )
abline(h = 0.234, col ="blue")
#Legend will be added after when optimal phi will be found, because
#optimal phi will have to be specified in the legend as well

#Selecting the most optimal phi value (matching which circle(phi value) is the
#closest one to line value = 0.234))
key_optimal <- which(abs(acc_result-0.234)==min(abs(acc_result-0.234)))
phi_optimal <- phi_acc[key_optimal]#33 number of phi: -0.1515152
phi_optimal#-0.1515152
abline(v = phi_optimal, col ="green", lty = 2)

#Creating the legend based on the values in plot
legend(6, 0.88, lty = c(1, 2), cex=0.8,
       legend = c("Acceptance rate", "Optimal phi"),
       col=c("blue", "green"))


### 1-D part###----------------------------------
N <- 100 
M <- 10000
n_means1 <- c()

#Calculating expectation N times rwmet[[1]]
for (j in 1:N){
  n_means1[j] <- mean(rwmet(M,phi_optimal)[[1]])
}
var_D <- var(n_means1)
var_D#2.923068e-05
var_D - var_B#-8.866078e-05; optimal phi has lower variance
#It means that a found optimal phi gives values with lower variablity from the mean


##### Question2 #####
### 2-A part###----------------------------------

#Writing a function to return the maximum likelihood estimates 
#of beta for a logistic regression model. X-matrix nxp;
#Y- nx1 vector; itmax-maximum number of iterations of the 
#IWLS algorithm; eps - the tolerance
iwls <- function(X,Y, itmax = 25,eps = 10^(-9)){
  
  #length of matrix rows
  n <- length(X[,1])
  
  #calculating beta value
  beta <- as.vector(solve(t(X) %*% X) %*% t(X) %*% Y)
  
  #creating variables to support the function
  cond <- 1
  miu <- as.vector(rep(0,n))
  miu_up <- as.vector(rep(0,n))
  z_vector <- as.vector(rep(0,n))
  w <- rep(0,n)
  j <- 1
  
  while(cond >= eps){
    
    for (i in 1:n){
      
      #selecting every new row vector from matrix
      x <- as.vector(X[i,])
      
      #calculating miu values
      miu[i] <- exp(x %*% beta)/(1+exp(x %*% beta))
      miu_up[i] <- 1/((miu[i])*(1-miu[i]))
      
      #calculating z vector nx1
      z_vector[i] <- (Y[i]-miu[i]) * miu_up[i]
      
      #calculating w nxn diagonal matrix
      w[i] <- 1/(miu[i] * (1-miu[i]) * (miu_up[i] ^ 2) )    
    }
    
    w_matrix <- diag(w)
    
    #formula to meet the convergence
    beta_hat <- beta + (solve(t(X) %*% w_matrix %*% X)) %*% t(X) %*% w_matrix %*% z_vector
    
    #Checking the condition to compare with tolerance level
    cond <- max(abs(beta - beta_hat))
    
    #giving a new beta for a new intertion
    beta <- beta_hat
    
    #calculating how many iterations are in already
    j <- j + 1
    
    #condition if interations reached maximum 25 value (25 rounds)
    if (j >= itmax){stop("itmax reached maximum value")}
  }
  beta_hat
}


### 2-B part###----------------------------------
#importing data into R
diab <- read.table("diab.txt", header = TRUE)

#creating a matrix for analysis
df_B <- as.matrix(diab, nrow=nrow(diab), ncol = 9)

#first column have to be covered by 1 as intercept
df_B[,1] <- 1
colnames(df_B)[1] <- 'intercept'

#moving values to matrix, except 'out' column
for (i in 2:9){
  df_B[,i] <- diab[,i-1]
  colnames(df_B)[i] <- colnames(diab)[i-1]
}

#creating y vector
df_Y = diab[,9]

#using iwls to get Coefficients for variables
iwls(df_B, df_Y)
#intercept -8.4046963669
#preg       0.1231822984
#gluc       0.0351637146
#bp        -0.0132955469
#skin       0.0006189644
#insulin   -0.0011916990
#bmi        0.0897009700
#ped        0.9451797406
#age        0.0148690047

#Creating glm model based on binomial family
glm(out ~., data = diab, family = "binomial")

#comparing the differences between glm ir iwls
glm(out ~., data = diab, family = "binomial")$coefficients -
iwls(df_B, df_Y)[,1]
#Intercept 9.716672e-13
#preg -1.236511e-14
#gluc -2.775558e-15 
#bp 1.139713e-15
#skin 1.734723e-18
#insulin 1.511378e-16
#bmi -1.504352e-14
#ped -1.418865e-13
#age -1.765948e-15

#There are minor differences based on results, but they are not significant