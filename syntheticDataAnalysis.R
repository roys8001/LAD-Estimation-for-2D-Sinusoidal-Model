library(ggplot2)
library(ggpubr)

## data generation
generate_data = function(par, t, s, distr, sd = 0.1, df = 1){
  data = matrix(0, nrow = t, ncol = s)
  
  if(distr == "normal"){  
    for(i in 1:t){
      for(j in 1:s){
        data[i,j] = par[1] * cos(par[3] * i + par[4] * j) + par[2] * sin(par[3] * i + par[4] * j) + rnorm(1, sd = sd)
      }
    }
  }else if(distr == "t"){
    for(i in 1:t){
      for(j in 1:s){
        data[i,j] = par[1] * cos(par[3] * i + par[4] * j) + par[2] * sin(par[3] * i + par[4] * j) + rt(1, df = df)
      }
    }
  }else if(distr == "slash"){
    for(i in 1:t){
      for(j in 1:s){
        data[i,j] = par[1] * cos(par[3] * i + par[4] * j) + par[2] * sin(par[3] * i + par[4] * j) + (rnorm(1) / runif(1))
      }
    }
  }
  
  ## returning the generated data
  return(data)
}

## phi function (written a more general function to incorporate LAD as a special case of RQ estimation)
phi.fn = function(beta = 0.5, lambda){
  if(lambda >= 0)
    return(beta * lambda)
  else
    return((beta - 1) * lambda)
}

## objective function for LAD
obj.fn.lad = function(par, beta = 0.5, y, t, s){
  x = 0
  for (i in 1:t){
    for (j in 1:s){
      x = x + phi.fn(beta, y[i, j] - par[1] * cos(par[3] * i + par[4] * j) - par[2] * sin(par[3] * i + par[4] * j))
    }
  }
  ## returning the objective function 
  return(x / (t * s))
}

## objective function for LSE
obj.fn.ls = function(par, y, t, s){
  x = 0
  for (i in 1:t) {
    for (j in 1:s){
      x = x + (y[i, j] - par[1] * cos(par[3] * i + par[4] * j) - par[2] * sin(par[3] * i + par[4] * j))^2
    }
  }
  ## returning the LAD objective function 
  return(x)
}

## LAD Estimation 
LAD_est_2d = function(original, initial, t, s, distr, sd = 0.1, df = 1){
  ## declaring the variables to store the LAD and LS estimates 
  est.lad = rep(0, length(original))
  est.lse = rep(0, length(original))
  Y = generate_data(par = original, t = t, s = s, distr, sd, df)
  LAD = optim(initial, obj.fn.lad , y = Y, t = t, s = s, method = "Nelder-Mead")
  LS = optim(initial, obj.fn.ls , y = Y, t = t, s = s, method = "Nelder-Mead")
  est.lad= LAD$par
  est.lse= LS$par
  
  ## returning the LAD and LS estimates
  return(list("data" = Y, "lad" = est.lad, 'lse' = est.lse))
}

## Comaparison of LAD and LS Estimates
fit_t = LAD_est_2d(original = c(2.4, 1.4, 0.4, 0.6), initial = c(2.4, 1.4, 0.4, 0.6), t = 100, s = 100, distr = "slash")

## Average estimate for LAD
LADest = fit_t$lad
## Average estimate for LS
LSest = fit_t$lse

###############################################################################################
################################# Synthetic data analysis #####################################
###############################################################################################
t = 100; s = 100

## original data matrix
originalData = fit_t$data

## predicted data from LAD estimation
pred.LAD = matrix(0, nrow = t, ncol = s)
for(i in 1:t){
  for(j in 1:s){
    pred.LAD[i,j] = LADest[1] * cos(LADest[3] * i + LADest[4] * j) + LADest[2] * sin(LADest[3] * i + LADest[4] * j)
  }
}

## predicted data from LAD estimation
pred.LSE = matrix(0, nrow = t, ncol = s)
for(i in 1:t){
  for(j in 1:s){
    pred.LSE[i,j] = LSest[1] * cos(LSest[3] * i + LSest[4] * j) + LSest[2] * sin(LSest[3] * i + LSest[4] * j)
  }
}
## data without noise
data_withoutnoise = matrix(0, t, s)
for(i in 1:t){
  for(j in 1:s){
    data_withoutnoise[i, j] = 2.4 * cos(0.4 * i + 0.6 * j) + 1.4 * sin(0.4 * i + 0.6 * j) 
  }
}

df.data.withoutnoise = reshape2::melt(data_withoutnoise, varnames = c("T","S"),value.name = "value")
df.data.original = reshape2::melt(originalData, varnames = c("T","S"),value.name = "value")
df.data.LAD = reshape2::melt(pred.LAD, varnames = c("T","S"),value.name = "value")
df.data.LSE = reshape2::melt(pred.LSE, varnames = c("T","S"),value.name = "value")

#pdf(file= "plots.pdf", width = 10, height = 10)

### Plotting gray-scale image of original data
ggplot(df.data.original, aes_string(x = "T", y = "S", fill = "value")) +
  geom_raster() +
  scale_x_continuous(name = "T", breaks = seq(0, t, by = 10)) +
  scale_y_continuous(name = "S", breaks = seq(0, s, by = 10)) +
  scale_fill_continuous(high = "white", low = "black", guide = "none")+
  theme_bw(base_size = 14) #+ labs(title = "Original data")

### Plotting gray-scale image of predicted data from LAD estimation
ggplot(df.data.LAD, aes_string(x = "T", y = "S", fill = "value"))+
  geom_raster() +
  scale_x_continuous(name = "T", breaks = seq(0, t, by = 10)) +
  scale_y_continuous(name = "S", breaks = seq(0, s, by = 10)) +
  scale_fill_continuous(high = "white", low = "black", guide = "none")+
  theme_bw(base_size = 14) #+ labs(title = "Predicted data from LAD estimation")

### Plotting gray-scale image of predicted data from LSE estimation
ggplot(df.data.LSE, aes_string(x = "T", y = "S", fill = "value"))+
  geom_raster() +
  scale_x_continuous(name = "T", breaks = seq(0, t, by = 10)) +
  scale_y_continuous(name = "S", breaks = seq(0, s, by = 10)) +
  scale_fill_continuous(high = "white", low = "black", guide = "none")+
  theme_bw(base_size = 14) #+ labs(title = "Predicted data from LS estimation")

### Plotting gray-scale image of noiseless original data 
ggplot(df.data.withoutnoise, aes_string(x = "T", y = "S", fill = "value")) +
  geom_raster() +
  scale_x_continuous(name = "T", breaks = seq(0, t, by = 10)) +
  scale_y_continuous(name = "S", breaks = seq(0, s, by = 10)) +
  scale_fill_continuous(high = "white", low = "black", guide = "none")+
  theme_bw(base_size = 14) #+ labs(title = "Noiseless original data")

# #dev.off()
# figure = ggarrange(p, # First row with line plot
#   # Second row with box and dot plots
#   ggarrange(q, r, ncol = 2, labels = c("B", "C")), 
#   nrow = 2, 
#   labels = "A"       # Label of the line plot
# ) 
# 
# figure



