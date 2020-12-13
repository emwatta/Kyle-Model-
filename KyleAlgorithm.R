cat("\014")
remove()
rm(list = ls(all = TRUE)) 
rm(list = ls(all.names = TRUE)) 
rm(list = ls())
graphics.off()


library(dplyr)
library(rootSolve)
library(ggplot2)
library(tidyverse)

######## Kyle Model Simulation ##############
#############################################
myImplySigma0_funct <- function(targetSigma0, inputSigmaN,
                                p_start, steps, sigma_u_sq, delta, alpha, returnvalue){
  
  Sigma.0 <- targetSigma0
  Sigma.N <- inputSigmaN
  
  p <- numeric(steps+1) ### price vector
  delta.t <- 1/steps
  v <- rnorm(n =1, mean = p_start, sd = sqrt(Sigma.0)) 
  delta_p <- numeric(steps)
  delta.u <- rnorm(n = steps, mean = 0, sd = sqrt(sigma_u_sq*delta.t)) 
  delta_x <- numeric(steps)
  Sigma <- numeric(steps+1)
  
  ## Nth period:
  #### Lambda-solver function #
  lambdafunction <- function(lambda) lambda -(((1-2*alpha*lambda)*Sigma.N)/(2*delta.t*sigma_u_sq*lambda*(1-alpha*lambda)))
  solver <- uniroot.all(lambdafunction, c(0, 100))
  lambda_N <- solver[1]
  
  #### Beta, sigma function ##
  beta_N <- ((1-2*alpha*lambda_N)/(delta.t*(2*lambda_N*(1-alpha*lambda_N))))
  beta_N
  
  Sigma_N1 <- Sigma.N/((1-delta.t*beta_N*lambda_N))
  Sigma_N1
  
  ###################
  ## Put the values into the vectors on the n-th place ##
  
  lambda <- numeric(steps+1)
  lambda[length(lambda)] <- lambda_N
  lambda
  
  beta <- numeric(steps+1)
  beta[length(beta)] <- beta_N
  beta
  
  Sigma <- numeric(steps+1)
  Sigma[length(Sigma)] <- Sigma.N
  Sigma
  
  Sigma[steps+1] <- Sigma.N
  p <- numeric(steps+1)
  alpha <- numeric(steps+1)

  
  # Put values into the vector ##
  lambda
  lambda_i <- nth(lambda, steps+1)
  delta <- numeric(steps+1)
  alpha_i <- nth(alpha, steps+1)
  delta_i <- nth(delta, steps+1)
  
  ####### N-1 period: ##
  # alpha,Sigma:
  for (i in steps:1){
    
    alpha[i] <- 1/(4*lambda[i+1]*(1-alpha[i+1]*lambda[i+1]))
    delta[i] <- delta[i+1] + alpha[i+1]*(lambda[i+1])^2*sigma_u_sq*delta.t
    Sigma[i] <- Sigma[i+1]/((1-delta.t*beta[i+1]*lambda[i+1]))
    
    ## lambda:
    alpha_i <- alpha[i] 
    Sigma_i <- Sigma[i] 
    
    
    lambda_func <- function(lambda) lambda -(((1-2*alpha_i*lambda)*Sigma_i)/(2*delta.t*sigma_u_sq*lambda*(1-alpha_i*lambda)))
    solver <- uniroot.all(lambda_func, c(0, 100))
    lambda_N <- solver[1]
    
    #put the lambda_n value into lambda vector
    lambda[i] <- lambda_N
    
    
    # beta, Sigma
    beta[i] <- (1-2*alpha[i]*lambda[i])/
      (2*lambda[i]*(1-alpha[i]*lambda[i])*delta.t)
    
    Sigma[i] <- Sigma[i+1]/((1-delta.t*beta[i+1]*lambda[i+1]))
  }
  
  ## After the iteration, we calculate delta_x, delta_p, p: ##
  ################ profits and expected profits
  p_start
  p[1] <- p_start
  p
  exp_profits <- numeric(steps)
  profits_periode <- numeric(steps)
  PROFITS <- numeric(steps)
  EXPECTED_PROFITS <- numeric(steps)

  for (j in 1:steps){
    delta_x[j] <- delta.t * (v-p[j])*beta[j+1]
    delta_p[j] <- (delta.u[j]+ delta_x[j])*lambda[j+1]
    p[j+1] <- p[j]+delta_p[j]
    exp_profits[j] <- alpha[j]*(v-p[j])^2+delta[j]
    profits_periode[j] <- (v-p[j+1])*delta_x[j]
  }
  
  #new_profits = cumsum(profits_periode)
  for (m in 1:steps){
    PROFITS[m] <- sum(profits_periode[m:1])
    EXPECTED_PROFITS[m] <- sum(exp_profits[m:1])
  }
  
  
  ## calculate y = u + x (delta.u + delta_x)
  y <- numeric(steps)
  
  for (i in 1:steps){
    y[i]=delta.u[i]+delta_x[i]
  }
  
  mySigma0 <- Sigma[1]
  mySigma0
  
  if(returnvalue == "err"){  
    return(mySigma0 - targetSigma0)
  }
  if(returnvalue == "list"){  
    return(list(alpha = alpha, 
                profits_periode = profits_periode,
                exp_profits = exp_profits,
                Sigma = Sigma,
                beta = beta, 
                delta_p= delta_p,
                delta_x = delta_x,
                delta=delta,
                p = p,
                v = v,
                delta.t = delta.t,
                delta.u = delta.u,
                lambda=lambda,
                PROFITS = PROFITS,
                y = y,
                EXPECTED_PROFITS = EXPECTED_PROFITS,
                steps=steps
    ))
  }
}


#### KYLE-function #### 
my_Kyle_funct <- function(Sigma.0, p_start, steps, sigma_u_sq, delta, alpha){
  
  my_implied_SigmaN_lst <- uniroot(f = myImplySigma0_funct,
                                   interval = c(0.00000000001, 10),
                                   targetSigma0 = Sigma.0,
                                   p_start = p_start,
                                   steps = steps,
                                   sigma_u_sq = sigma_u_sq,
                                   delta = delta,
                                   alpha = alpha,
                                   returnvalue = "err")
  
  my_Kyle_results_lst <- myImplySigma0_funct(targetSigma0 = Sigma.0, 
                                             inputSigmaN = my_implied_SigmaN_lst$root,
                                             p_start = p_start, 
                                             steps = steps, 
                                             sigma_u_sq = sigma_u_sq, 
                                             delta = delta, 
                                             alpha = alpha,
                                             returnvalue = "list")
  
  return(my_Kyle_results_lst)
}

#### Inputs of your variables (change inputs)
Sigma.0 = 0.5 
p_start = 8 
steps = 50
sigma_u_sq = 0.3  
delta = 0
alpha = 0

#####################################################################
mykyleresults <- my_Kyle_funct(Sigma.0 = Sigma.0,
                                p_start = p_start,
                                steps = steps,
                                sigma_u_sq = sigma_u_sq,
                                delta = delta,
                                alpha = alpha)

mykyleresults

### outputs/vectors of variables: 
delta <- mykyleresults$delta
lambda <- mykyleresults$lambda
alpha <- mykyleresults$alpha
beta <- mykyleresults$beta
Sigma <- mykyleresults$Sigma
p <- mykyleresults$p
v <- mykyleresults$v
delta.t <- mykyleresults$delta.t
steps <- mykyleresults$steps  
delta_x <- mykyleresults$delta_x
delta_p <- mykyleresults$delta_p
delta.u <- mykyleresults$delta.u

############################################################
names(mykyleresults)

timedata <- data.frame(1:(steps+1))
names(timedata)[1] <- "time"
time <- timedata

data_1 <- data.frame(delta_p, delta_x, delta.u)
data_2 <- rep(NA, time=steps)
data_2 <- rbind( data_1 , data_2)

data_3 <- data.frame(Sigma, alpha, beta, delta, lambda, delta, p, time)
data <- data.frame(data_2, data_3)
data

#### PLOTS #########
############################################################

### alpha ###
plot_alpha <- ggplot(data = data, aes(x = time, y = alpha))+
  geom_point(shape=16) + labs(x = "time", y = "alpha")
plot_alpha

### delta ###
plot_delta <- ggplot(data = data, aes(x = time, y = delta)) +
  geom_point(shape=16) + labs(x = "time", y = "delta")
plot_delta

### beta ###
plot_beta <- ggplot(data = data, aes(x = time, y = beta)) +
  geom_point(shape=16) + labs(x = "time", y = "beta")
plot_beta

### lambda ###
plot_lambda <- ggplot(data = data, aes(x = time, y = lambda)) +
  geom_point(shape=16) + labs(x = "time", y = "lambda")
plot_lambda

### Sigma ###
plot_Sigma <- ggplot(data = data, aes(x = time, y = Sigma)) +
  geom_point(shape=16) + labs(x = "time",y = "Sigma")
plot_Sigma

##### Plot - Outputs ####
plot_alpha
plot_delta
plot_Sigma
plot_beta
plot_lambda

