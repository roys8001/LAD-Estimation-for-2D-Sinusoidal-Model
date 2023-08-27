library(writexl)
########## Two dimensional sinusoidal model ##################

## data generation
generate_data <- function(par, t, s, distr, sd = 0.1, df = 1){
  data <- matrix(0, nrow = t, ncol = s)
  
  if(distr == "normal"){  
    for(i in 1:t){
      for(j in 1:s){
        data[i,j] <- par[1] * cos(par[3] * i + par[4] * j) + par[2] * sin(par[3] * i + par[4] * j) + rnorm(1, sd = sd)
      }
    }
  }else if(distr == "t"){
    for(i in 1:t){
      for(j in 1:s){
        data[i,j] <- par[1] * cos(par[3] * i + par[4] * j) + par[2] * sin(par[3] * i + par[4] * j) + rt(1, df = df)
      }
    }
  }else if(distr == "slash"){
    for(i in 1:t){
      for(j in 1:s){
        data[i,j] <- par[1] * cos(par[3] * i + par[4] * j) + par[2] * sin(par[3] * i + par[4] * j) + (rnorm(1) / runif(1))
      }
    }
  }

  ## returning the generated data
  return(data)
}

## phi function (written a more general function to incorporate LAD as a special case of RQ estimation)
phi.fn <- function(beta = 0.5, lambda){
  if(lambda >= 0)
    return(beta * lambda)
  else
    return((beta - 1) * lambda)
}

## objective function for LAD
obj.fn.lad <- function(par, beta = 0.5, y, t, s){
  x <- 0
  for (i in 1:t){
    for (j in 1:s){
      x <- x + phi.fn(beta, y[i, j] - par[1] * cos(par[3] * i + par[4] * j) - par[2] * sin(par[3] * i + par[4] * j))
    }
  }
  ## returning the objective function 
  return(x / (t * s))
  
}

## objective function for LSE
obj.fn.ls <- function(par, y, t, s){
  x <- 0
  for (i in 1:t) {
    for (j in 1:s){
      x <- x + (y[i, j] - par[1] * cos(par[3] * i + par[4] * j) - par[2] * sin(par[3] * i + par[4] * j))^2
    }
  }
  ## returning the LAD objective function 
  return(x)
  
}

## LAD Estimation 
LAD_est_2d <- function(original, initial, t, s, rep, distr, sd = 0.1, df = 1){
  ## declaring the variables to store the LAD and LS estimates 
  est.lad <- matrix(0, nrow = rep, ncol = 4)
  est.lse <- matrix(0, nrow = rep, ncol = 4)
  for(k in 1:rep){
    Y <- generate_data(par = original, t = t, s = s, distr, sd, df)
    LAD <- optim(initial, obj.fn.lad , y = Y, t = t, s = s, method = "Nelder-Mead")
    LS <- optim(initial, obj.fn.ls , y = Y, t = t, s = s, method = "Nelder-Mead")
    est.lad[k, ] <- LAD$par
    est.lse[k, ] <- LS$par
  }
  ## returning the LAD and LS estimates
  return(list("lad" = est.lad, 'lse' = est.lse))
}

## Functions to calculate asymptotic variance
asym_var <- function(beta = 0.5, distr, sd = 0.1, df = 1, t, s, a.0, b.0){
  ## precision matrix \Sigma^{-1}
  sigma <- matrix(c(1/2, 0, b.0/4, b.0/4, 0, 1/2, -a.0/4, -a.0/4,
                    b.0/4, -a.0/4, (a.0^2+b.0^2)/6, (a.0^2+b.0^2)/8, b.0/4,
                    -a.0/4, (a.0^2+b.0^2)/8, (a.0^2+b.0^2)/6), 
                  nrow = 4, ncol = 4, byrow = TRUE)
  precision_mat <- solve(sigma)
  ## Extracting the diagonal elements of the precision matrix
  diagonal <- diag(precision_mat)
  
  ## Calculation of g_0
  if(distr == "normal"){
    g_0 <- 1 / sqrt(2 * pi * sd)
  }else if(distr == "t"){
    g_0 <- 1 / pi
  }else if(distr == "slash"){
    g_0 <- 1 / (2 * sqrt(2 * pi))
  }
  
  temp <- (beta * (1 - beta) / g_0^2) * diagonal
  v1 <- temp[1] / (t * s)
  v2 <- temp[2] / (t * s)
  v3 <- temp[3] / (t^3 * s)
  v4 <- temp[4] / (s^3 * t)
  
  ## returning the asymptotic variance
  return(c(v1,v2,v3,v4))
}

##################################################################
######################## Experiment ############################## 
##################################################################

###############  Normal(0, 0.1^2) Distribution ###################

############# Combination for (t, s) =  (25, 25) #################

## Calculation of asymtotic variance
asyvar_25_normal <- asym_var(beta = 0.5, distr = "normal", sd = 0.1, df = 1, t = 25, s = 25, a.0 = 2.4, b.0 = 1.4)

## Comaparison of LAD and LS Estimates
fit1_normal <- LAD_est_2d(original = c(2.4, 1.4, 0.4, 0.6), initial = c(2.4, 1.4, 0.4, 0.6), t = 25, s = 25, rep = 1000, distr = "normal")

## Average estimate for LAD
LAD_ae_25_normal <- apply(fit1_normal$lad, 2, mean)
## Mean Squared Error for LAD
LAD_mse_25_normal <- apply(fit1_normal$lad, 2, var)
## Average estimate for LS
LSE_ae_25_normal <- apply(fit1_normal$lse, 2, mean)
## Mean Squared Error for LS
LSE_mse_25_normal <- apply(fit1_normal$lse, 2, var)


############# Combination for (t, s) =  (50, 50) #################

## Calculation of asymtotic variance
asyvar_50_normal <- asym_var(beta = 0.5, distr = "normal", sd = 0.1, df = 1, t = 50, s = 50, a.0 = 2.4, b.0 = 1.4)

## Comaparison of LAD and LS Estimates
fit2_normal <- LAD_est_2d(original = c(2.4, 1.4, 0.4, 0.6), initial = c(2.4, 1.4, 0.4, 0.6), t = 50, s = 50, rep = 1000, distr = "normal")

## Average estimate for LAD
LAD_ae_50_normal <- apply(fit2_normal$lad, 2, mean)
## Mean Squared Error for LAD
LAD_mse_50_normal <- apply(fit2_normal$lad, 2, var)
## Average estimate for LS
LSE_ae_50_normal <- apply(fit2_normal$lse, 2, mean)
## Mean Squared Error for LS
LSE_mse_50_normal <- apply(fit2_normal$lse, 2, var)

############# Combination for (t, s) =  (75, 75) #################

## Calculation of asymtotic variance
asyvar_75_normal <- asym_var(beta = 0.5, distr = "normal", sd = 0.1, df = 1, t = 75, s = 75, a.0 = 2.4, b.0 = 1.4)

## Comaparison of LAD and LS Estimates
fit3_normal <- LAD_est_2d(original = c(2.4, 1.4, 0.4, 0.6), initial = c(2.4, 1.4, 0.4, 0.6), t = 75, s = 75, rep = 1000, distr = "normal")

## Average estimate for LAD
LAD_ae_75_normal <- apply(fit3_normal$lad, 2, mean)
## Mean Squared Error for LAD
LAD_mse_75_normal <- apply(fit3_normal$lad, 2, var)
## Average estimate for LS
LSE_ae_75_normal <- apply(fit3_normal$lse, 2, mean)
## Mean Squared Error for LS
LSE_mse_75_normal <- apply(fit3_normal$lse, 2, var)

############# Combination for (t, s) =  (150, 150) #################

## Calculation of asymtotic variance
asyvar_150_normal <- asym_var(beta = 0.5, distr = "normal", sd = 0.1, df = 1, t = 150, s = 150, a.0 = 2.4, b.0 = 1.4)

## Comaparison of LAD and LS Estimates
fit4_normal <- LAD_est_2d(original = c(2.4, 1.4, 0.4, 0.6), initial = c(2.4, 1.4, 0.4, 0.6), t = 150, s = 150, rep = 1000, distr = "normal")

## Average estimate for LAD
LAD_ae_150_normal <- apply(fit4_normal$lad, 2, mean)
## Mean Squared Error for LAD
LAD_mse_150_normal <- apply(fit4_normal$lad, 2, var)
## Average estimate for LS
LSE_ae_150_normal <- apply(fit4_normal$lse, 2, mean)
## Mean Squared Error for LS
LSE_mse_150_normal <- apply(fit4_normal$lse, 2, var)

############# Combination for (t, s) =  (300, 300) #################

## Calculation of asymtotic variance
asyvar_300_normal <- asym_var(beta = 0.5, distr = "normal", sd = 0.1, df = 1, t = 300, s = 300, a.0 = 2.4, b.0 = 1.4)

## Comaparison of LAD and LS Estimates
fit5_normal <- LAD_est_2d(original = c(2.4, 1.4, 0.4, 0.6), initial = c(2.4, 1.4, 0.4, 0.6), t = 300, s = 300, rep = 1000, distr = "normal")

## Average estimate for LAD
LAD_ae_300_normal <- apply(fit5_normal$lad, 2, mean)
## Mean Squared Error for LAD
LAD_mse_300_normal <- apply(fit5_normal$lad, 2, var)
## Average estimate for LS
LSE_ae_300_normal <- apply(fit5_normal$lse, 2, mean)
## Mean Squared Error for LS
LSE_mse_300_normal <- apply(fit5_normal$lse, 2, var)



####################### t(1) Distribution ########################

############# Combination for (t, s) =  (25, 25) #################

## Calculation of asymtotic variance
asyvar_25_t <- asym_var(beta = 0.5, distr = "t", sd = 0.1, df = 1, t = 25, s = 25, a.0 = 2.4, b.0 = 1.4)

## Comaparison of LAD and LS Estimates
fit1_t <- LAD_est_2d(original = c(2.4, 1.4, 0.4, 0.6), initial = c(2.4, 1.4, 0.4, 0.6), t = 25, s = 25, rep = 1000, distr = "t")

## Average estimate for LAD
LAD_ae_25_t <- apply(fit1_t$lad, 2, mean)
## Mean Squared Error for LAD
LAD_mse_25_t <- apply(fit1_t$lad, 2, var)
## Average estimate for LS
LSE_ae_25_t <- apply(fit1_t$lse, 2, mean)
## Mean Squared Error for LS
LSE_mse_25_t <- apply(fit1_t$lse, 2, var)


############# Combination for (t, s) =  (50, 50) #################

## Calculation of asymtotic variance
asyvar_50_t <- asym_var(beta = 0.5, distr = "t", sd = 0.1, df = 1, t = 50, s = 50, a.0 = 2.4, b.0 = 1.4)

## Comaparison of LAD and LS Estimates
fit2_t <- LAD_est_2d(original = c(2.4, 1.4, 0.4, 0.6), initial = c(2.4, 1.4, 0.4, 0.6), t = 50, s = 50, rep = 1000, distr = "t")

## Average estimate for LAD
LAD_ae_50_t <- apply(fit2_t$lad, 2, mean)
## Mean Squared Error for LAD
LAD_mse_50_t <- apply(fit2_t$lad, 2, var)
## Average estimate for LS
LSE_ae_50_t <- apply(fit2_t$lse, 2, mean)
## Mean Squared Error for LS
LSE_mse_50_t <- apply(fit2_t$lse, 2, var)

############# Combination for (t, s) =  (75, 75) #################

## Calculation of asymtotic variance
asyvar_75_t <- asym_var(beta = 0.5, distr = "t", sd = 0.1, df = 1, t = 75, s = 75, a.0 = 2.4, b.0 = 1.4)

## Comaparison of LAD and LS Estimates
fit3_t <- LAD_est_2d(original = c(2.4, 1.4, 0.4, 0.6), initial = c(2.4, 1.4, 0.4, 0.6), t = 75, s = 75, rep = 1000, distr = "t")

## Average estimate for LAD
LAD_ae_75_t <- apply(fit3_t$lad, 2, mean)
## Mean Squared Error for LAD
LAD_mse_75_t <- apply(fit3_t$lad, 2, var)
## Average estimate for LS
LSE_ae_75_t <- apply(fit3_t$lse, 2, mean)
## Mean Squared Error for LS
LSE_mse_75_t <- apply(fit3_t$lse, 2, var)

############# Combination for (t, s) =  (150, 150) #################

## Calculation of asymtotic variance
asyvar_150_t <- asym_var(beta = 0.5, distr = "t", sd = 0.1, df = 1, t = 150, s = 150, a.0 = 2.4, b.0 = 1.4)

## Comaparison of LAD and LS Estimates
fit4_t <- LAD_est_2d(original = c(2.4, 1.4, 0.4, 0.6), initial = c(2.4, 1.4, 0.4, 0.6), t = 150, s = 150, rep = 1000, distr = "t")

## Average estimate for LAD
LAD_ae_150_t <- apply(fit4_t$lad, 2, mean)
## Mean Squared Error for LAD
LAD_mse_150_t <- apply(fit4_t$lad, 2, var)
## Average estimate for LS
LSE_ae_150_t <- apply(fit4_t$lse, 2, mean)
## Mean Squared Error for LS
LSE_mse_150_t <- apply(fit4_t$lse, 2, var)

############# Combination for (t, s) =  (300, 300) #################

## Calculation of asymtotic variance
asyvar_300_t <- asym_var(beta = 0.5, distr = "t", sd = 0.1, df = 1, t = 300, s = 300, a.0 = 2.4, b.0 = 1.4)

## Comaparison of LAD and LS Estimates
fit5_t <- LAD_est_2d(original = c(2.4, 1.4, 0.4, 0.6), initial = c(2.4, 1.4, 0.4, 0.6), t = 300, s = 300, rep = 1000, distr = "t")

## Average estimate for LAD
LAD_ae_300_t <- apply(fit5_t$lad, 2, mean)
## Mean Squared Error for LAD
LAD_mse_300_t <- apply(fit5_t$lad, 2, var)
## Average estimate for LS
LSE_ae_300_t <- apply(fit5_t$lse, 2, mean)
## Mean Squared Error for LS
LSE_mse_300_t <- apply(fit5_t$lse, 2, var)



###################### slash Distribution ########################

############# Combination for (t, s) =  (25, 25) #################

## Calculation of asymtotic variance
asyvar_25_slash <- asym_var(beta = 0.5, distr = "slash", sd = 0.1, df = 1, t = 25, s = 25, a.0 = 2.4, b.0 = 1.4)

## Comaparison of LAD and LS Estimates
fit1_slash <- LAD_est_2d(original = c(2.4, 1.4, 0.4, 0.6), initial = c(2.4, 1.4, 0.4, 0.6), t = 25, s = 25, rep = 1000, distr = "slash")

## Average estimate for LAD
LAD_ae_25_slash <- apply(fit1_slash$lad, 2, mean)
## Mean Squared Error for LAD
LAD_mse_25_slash <- apply(fit1_slash$lad, 2, var)
## Average estimate for LS
LSE_ae_25_slash <- apply(fit1_slash$lse, 2, mean)
## Mean Squared Error for LS
LSE_mse_25_slash <- apply(fit1_slash$lse, 2, var)


############# Combination for (t, s) =  (50, 50) #################

## Calculation of asymtotic variance
asyvar_50_slash <- asym_var(beta = 0.5, distr = "slash", sd = 0.1, df = 1, t = 50, s = 50, a.0 = 2.4, b.0 = 1.4)

## Comaparison of LAD and LS Estimates
fit2_slash <- LAD_est_2d(original = c(2.4, 1.4, 0.4, 0.6), initial = c(2.4, 1.4, 0.4, 0.6), t = 50, s = 50, rep = 1000, distr = "slash")

## Average estimate for LAD
LAD_ae_50_slash <- apply(fit2_slash$lad, 2, mean)
## Mean Squared Error for LAD
LAD_mse_50_slash <- apply(fit2_slash$lad, 2, var)
## Average estimate for LS
LSE_ae_50_slash <- apply(fit2_slash$lse, 2, mean)
## Mean Squared Error for LS
LSE_mse_50_slash <- apply(fit2_slash$lse, 2, var)

############# Combination for (t, s) =  (75, 75) #################

## Calculation of asymtotic variance
asyvar_75_slash <- asym_var(beta = 0.5, distr = "slash", sd = 0.1, df = 1, t = 75, s = 75, a.0 = 2.4, b.0 = 1.4)

## Comaparison of LAD and LS Estimates
fit3_slash <- LAD_est_2d(original = c(2.4, 1.4, 0.4, 0.6), initial = c(2.4, 1.4, 0.4, 0.6), t = 75, s = 75, rep = 1000, distr = "slash")

## Average estimate for LAD
LAD_ae_75_slash <- apply(fit3_slash$lad, 2, mean)
## Mean Squared Error for LAD
LAD_mse_75_slash <- apply(fit3_slash$lad, 2, var)
## Average estimate for LS
LSE_ae_75_slash <- apply(fit3_slash$lse, 2, mean)
## Mean Squared Error for LS
LSE_mse_75_slash <- apply(fit3_slash$lse, 2, var)

############# Combination for (t, s) =  (150, 150) #################

## Calculation of asymtotic variance
asyvar_150_slash <- asym_var(beta = 0.5, distr = "slash", sd = 0.1, df = 1, t = 150, s = 150, a.0 = 2.4, b.0 = 1.4)

## Comaparison of LAD and LS Estimates
fit4_slash <- LAD_est_2d(original = c(2.4, 1.4, 0.4, 0.6), initial = c(2.4, 1.4, 0.4, 0.6), t = 150, s = 150, rep = 1000, distr = "slash")

## Average estimate for LAD
LAD_ae_150_slash <- apply(fit4_slash$lad, 2, mean)
## Mean Squared Error for LAD
LAD_mse_150_slash <- apply(fit4_slash$lad, 2, var)
## Average estimate for LS
LSE_ae_150_slash <- apply(fit4_slash$lse, 2, mean)
## Mean Squared Error for LS
LSE_mse_150_slash <- apply(fit4_slash$lse, 2, var)

############# Combination for (t, s) =  (300, 300) #################

## Calculation of asymtotic variance
asyvar_300_slash <- asym_var(beta = 0.5, distr = "slash", sd = 0.1, df = 1, t = 300, s = 300, a.0 = 2.4, b.0 = 1.4)

## Comaparison of LAD and LS Estimates
fit5_slash <- LAD_est_2d(original = c(2.4, 1.4, 0.4, 0.6), initial = c(2.4, 1.4, 0.4, 0.6), t = 300, s = 300, rep = 1000, distr = "slash")

## Average estimate for LAD
LAD_ae_300_slash <- apply(fit5_slash$lad, 2, mean)
## Mean Squared Error for LAD
LAD_mse_300_slash <- apply(fit5_slash$lad, 2, var)
## Average estimate for LS
LSE_ae_300_slash <- apply(fit5_slash$lse, 2, mean)
## Mean Squared Error for LS
LSE_mse_300_slash <- apply(fit5_slash$lse, 2, var)


####################################################################################
################ Summarizing the simulation results in a data frame ################
####################################################################################
summary <- data.frame(Distribution = c(rep("normal(0, 0.1^2)", 25), rep("t(1)", 25), rep("slash", 25)), T = rep(c(25, 50, 75, 150, 300), each = 5, times = 3), S = rep(c(25, 50, 75, 150, 300), each = 5, times = 3), ReportedVals = rep(c("LAD Avg Est", "LAD MSE", "Asymtotic Variance", "LS Avg Est", "LS MSE"), 15),  
                    A = c(LAD_ae_25_normal[1], LAD_mse_25_normal[1], asyvar_25_normal[1], LSE_ae_25_normal[1], LSE_mse_25_normal[1], LAD_ae_50_normal[1], LAD_mse_50_normal[1], asyvar_50_normal[1], LSE_ae_50_normal[1], LSE_mse_50_normal[1], LAD_ae_75_normal[1], LAD_mse_75_normal[1], asyvar_75_normal[1], LSE_ae_75_normal[1], LSE_mse_75_normal[1], LAD_ae_150_normal[1], LAD_mse_150_normal[1], asyvar_150_normal[1], LSE_ae_150_normal[1], LSE_mse_150_normal[1], LAD_ae_300_normal[1], LAD_mse_300_normal[1], asyvar_300_normal[1], LSE_ae_300_normal[1], LSE_mse_300_normal[1], LAD_ae_25_t[1], LAD_mse_25_t[1], asyvar_25_t[1], LSE_ae_25_t[1], LSE_mse_25_t[1], LAD_ae_50_t[1], LAD_mse_50_t[1], asyvar_50_t[1], LSE_ae_50_t[1], LSE_mse_50_t[1], LAD_ae_75_t[1], LAD_mse_75_t[1], asyvar_75_t[1], LSE_ae_75_t[1], LSE_mse_75_t[1], LAD_ae_150_t[1], LAD_mse_150_t[1], asyvar_150_t[1], LSE_ae_150_t[1], LSE_mse_150_t[1], LAD_ae_300_t[1], LAD_mse_300_t[1], asyvar_300_t[1], LSE_ae_300_t[1], LSE_mse_300_t[1], LAD_ae_25_slash[1], LAD_mse_25_slash[1], asyvar_25_slash[1], LSE_ae_25_slash[1], LSE_mse_25_slash[1], LAD_ae_50_slash[1], LAD_mse_50_slash[1], asyvar_50_slash[1], LSE_ae_50_slash[1], LSE_mse_50_slash[1], LAD_ae_75_slash[1], LAD_mse_75_slash[1], asyvar_75_slash[1], LSE_ae_75_slash[1], LSE_mse_75_slash[1], LAD_ae_150_slash[1], LAD_mse_150_slash[1], asyvar_150_slash[1], LSE_ae_150_slash[1], LSE_mse_150_slash[1], LAD_ae_300_slash[1], LAD_mse_300_slash[1], asyvar_300_slash[1], LSE_ae_300_slash[1], LSE_mse_300_slash[1]), 
                    B = c(LAD_ae_25_normal[2], LAD_mse_25_normal[2], asyvar_25_normal[2], LSE_ae_25_normal[2], LSE_mse_25_normal[2], LAD_ae_50_normal[2], LAD_mse_50_normal[2], asyvar_50_normal[2], LSE_ae_50_normal[2], LSE_mse_50_normal[2], LAD_ae_75_normal[2], LAD_mse_75_normal[2], asyvar_75_normal[2], LSE_ae_75_normal[2], LSE_mse_75_normal[2], LAD_ae_150_normal[2], LAD_mse_150_normal[2], asyvar_150_normal[2], LSE_ae_150_normal[2], LSE_mse_150_normal[2], LAD_ae_300_normal[2], LAD_mse_300_normal[2], asyvar_300_normal[2], LSE_ae_300_normal[2], LSE_mse_300_normal[2], LAD_ae_25_t[2], LAD_mse_25_t[2], asyvar_25_t[2], LSE_ae_25_t[2], LSE_mse_25_t[2], LAD_ae_50_t[2], LAD_mse_50_t[2], asyvar_50_t[2], LSE_ae_50_t[2], LSE_mse_50_t[2], LAD_ae_75_t[2], LAD_mse_75_t[2], asyvar_75_t[2], LSE_ae_75_t[2], LSE_mse_75_t[2], LAD_ae_150_t[2], LAD_mse_150_t[2], asyvar_150_t[2], LSE_ae_150_t[2],LAD_ae_300_t[2], LAD_mse_300_t[2], asyvar_300_t[2], LSE_ae_300_t[2], LSE_mse_300_t[2], LSE_mse_150_t[2], LAD_ae_25_slash[2], LAD_mse_25_slash[2], asyvar_25_slash[2], LSE_ae_25_slash[2], LSE_mse_25_slash[2], LAD_ae_50_slash[2], LAD_mse_50_slash[2], asyvar_50_slash[2], LSE_ae_50_slash[2], LSE_mse_50_slash[2], LAD_ae_75_slash[2], LAD_mse_75_slash[2], asyvar_75_slash[2], LSE_ae_75_slash[2], LSE_mse_75_slash[2], LAD_ae_150_slash[2], LAD_mse_150_slash[2], asyvar_150_slash[2], LSE_ae_150_slash[2], LSE_mse_150_slash[2], LAD_ae_300_slash[2], LAD_mse_300_slash[2], asyvar_300_slash[2], LSE_ae_300_slash[2], LSE_mse_300_slash[2]), 
                    lambda = c(LAD_ae_25_normal[3], LAD_mse_25_normal[3], asyvar_25_normal[3], LSE_ae_25_normal[3], LSE_mse_25_normal[3], LAD_ae_50_normal[3], LAD_mse_50_normal[3], asyvar_50_normal[3], LSE_ae_50_normal[3], LSE_mse_50_normal[3], LAD_ae_75_normal[3], LAD_mse_75_normal[3], asyvar_75_normal[3], LSE_ae_75_normal[3], LSE_mse_75_normal[3], LAD_ae_150_normal[3], LAD_mse_150_normal[3], asyvar_150_normal[3], LSE_ae_150_normal[3], LSE_mse_150_normal[3], LAD_ae_300_normal[3], LAD_mse_300_normal[3], asyvar_300_normal[3], LSE_ae_300_normal[3], LSE_mse_300_normal[3], LAD_ae_25_t[3], LAD_mse_25_t[3], asyvar_25_t[3], LSE_ae_25_t[3], LSE_mse_25_t[3], LAD_ae_50_t[3], LAD_mse_50_t[3], asyvar_50_t[3], LSE_ae_50_t[3], LSE_mse_50_t[3], LAD_ae_75_t[3], LAD_mse_75_t[3], asyvar_75_t[3], LSE_ae_75_t[3], LSE_mse_75_t[3], LAD_ae_150_t[3], LAD_mse_150_t[3], asyvar_150_t[3], LSE_ae_150_t[3], LSE_mse_150_t[3], LAD_ae_300_t[3], LAD_mse_300_t[3], asyvar_300_t[3], LSE_ae_300_t[3], LSE_mse_300_t[3], LAD_ae_25_slash[3], LAD_mse_25_slash[3], asyvar_25_slash[3], LSE_ae_25_slash[3], LSE_mse_25_slash[3], LAD_ae_50_slash[3], LAD_mse_50_slash[3], asyvar_50_slash[3], LSE_ae_50_slash[3], LSE_mse_50_slash[3], LAD_ae_75_slash[3], LAD_mse_75_slash[3], asyvar_75_slash[3], LSE_ae_75_slash[3], LSE_mse_75_slash[3], LAD_ae_150_slash[3], LAD_mse_150_slash[3], asyvar_150_slash[3], LSE_ae_150_slash[3], LSE_mse_150_slash[3], LAD_ae_300_slash[3], LAD_mse_300_slash[3], asyvar_300_slash[3], LSE_ae_300_slash[3], LSE_mse_300_slash[3]), 
                    mu = c(LAD_ae_25_normal[4], LAD_mse_25_normal[4], asyvar_25_normal[4], LSE_ae_25_normal[4], LSE_mse_25_normal[4], LAD_ae_50_normal[4], LAD_mse_50_normal[4], asyvar_50_normal[4], LSE_ae_50_normal[4], LSE_mse_50_normal[4], LAD_ae_75_normal[4], LAD_mse_75_normal[4], asyvar_75_normal[4], LSE_ae_75_normal[4], LSE_mse_75_normal[4], LAD_ae_150_normal[4], LAD_mse_150_normal[4], asyvar_150_normal[4], LSE_ae_150_normal[4], LSE_mse_150_normal[4], LAD_ae_300_normal[4], LAD_mse_300_normal[4], asyvar_300_normal[4], LSE_ae_300_normal[4], LSE_mse_300_normal[4], LAD_ae_25_t[4], LAD_mse_25_t[4], asyvar_25_t[4], LSE_ae_25_t[4], LSE_mse_25_t[4], LAD_ae_50_t[4], LAD_mse_50_t[4], asyvar_50_t[4], LSE_ae_50_t[4], LSE_mse_50_t[4], LAD_ae_75_t[4], LAD_mse_75_t[4], asyvar_75_t[4], LSE_ae_75_t[4], LSE_mse_75_t[4], LAD_ae_150_t[4], LAD_mse_150_t[4], asyvar_150_t[4], LSE_ae_150_t[4], LSE_mse_150_t[4], LAD_ae_300_t[4], LAD_mse_300_t[4], asyvar_300_t[4], LSE_ae_300_t[4], LSE_mse_300_t[4], LAD_ae_25_slash[4], LAD_mse_25_slash[4], asyvar_25_slash[4], LSE_ae_25_slash[4], LSE_mse_25_slash[4], LAD_ae_50_slash[4], LAD_mse_50_slash[4], asyvar_50_slash[4], LSE_ae_50_slash[4], LSE_mse_50_slash[4], LAD_ae_75_slash[4], LAD_mse_75_slash[4], asyvar_75_slash[4], LSE_ae_75_slash[4], LSE_mse_75_slash[4], LAD_ae_150_slash[4], LAD_mse_150_slash[4], asyvar_150_slash[4], LSE_ae_150_slash[4], LSE_mse_150_slash[4], LAD_ae_300_slash[4], LAD_mse_300_slash[4], asyvar_300_slash[4], LSE_ae_300_slash[4], LSE_mse_300_slash[4]))

View(summary)
save.image(file = "~/Desktop/Research/2DSinusoidalModel/SimulationResults.RData")
#write_xlsx(table,"~/Desktop/Research/2DSinusoidalModel/SimulationSummary.xlsx")
load(file = "~/Desktop/Research/2DSinusoidalModel/SimulationResults.RData")
