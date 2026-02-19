#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# 4A - Structural equation modelling
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


## Function
# Moving window SEM function
# the function has to be implemented according to: 
# - the environmental gradient X
# - the number of parameters (P) in the sem model (i.e. arrows)
# - the name of the parameters (P1-3) in the sem model (name of the regression coefficient, variance or covariance specified by an arrow in the sem model)
movingsem <- 	function(data, R, sem.model){
  result <- data.frame(Start=c(1:(nrow(data)-R+1)), av.swe =0, rsq.Enz=0, n0 = 0, com = 0, lat.sd = 0, lat.max = 0, lat.min = 0, ab = 0)
  for(i in 1:(nrow(data)-R+1)) { 
    mat <- data[i:(i+R-1),]
    sem.model = lavaan::sem(full_model,data=vegan::decostand(mat[,-c(16, 17)], 'standardize'))
    sum.sem = lavaan::summary(sem.model, standardize = T, rsq=T, fit.measures=T)
    esti.coef = sum.sem$pe
    esti.coef = as.data.frame(esti.coef)
    latent = as.data.frame(lavPredict(sem.model, type = "lv"))
    
    #Extract the average value of ASV richness
    result$av.X[i] <- mean(mat$n0, na.rm = TRUE)
    
    #Extract variables of interest from SEM
    result$rsq.Enz[i] = esti.coef[47, 5]
    result$tmmn[i] = esti.coef[7, 10]
    result$ph[i] = esti.coef[8, 10]
    result$n0[i] = esti.coef[11, 10]
    result$ab[i] = esti.coef[13, 10]
    result$com[i] = esti.coef[14, 10]
    result$n0enz[i] = esti.coef[15, 10]
    result$lat.sd[i] = sd(latent$enz)
    result$lat.max[i] = max(latent$enz)
    result$lat50[i] = quantile(latent$enz, 0.50)
    result$lat75[i] = quantile(latent$enz, 0.75)
    result$lat95[i] = quantile(latent$enz, 0.95)
    
    #Extract SEM fit indices
    fit <- as.data.frame(sum.sem$fit)
    result$chi[i] = fit[3,]
    result$pvalue[i] = fit[5,]
    result$rmsea[i] = fit[17,]
    result$srmr[i] = fit[25,]
    result$aic[i] = fit[13,]
    result$cfi[i] = fit[9,]
  }
  result
}

## Import data
datsem <- read.csv("SEM_data_cleaned.csv", row.names=1)
datsem$bac_biomass <- datsem$Total_Bacteria_AB_g*0.383*5.60e-7

#order datsem by decreasing ASV richness (for the smoothed appraoch)
datsem <- datsem %>%  arrange(desc(n0))

###
#### Fit SEM model ----
###
## sem model full
full_model <- '
#latent variable definitions
enz =~ BG + NAG + AP + LAP
com =~ dim1 + dim3

# regression
n0 ~ tmmn_lt + ph
TotalAB  ~ com + n0
com ~ n0 + tmmn_lt
enz ~ TotalAB  + com + n0

#residual correlations
BG ~~ NAG
BG ~~ AP
AP ~~ LAP
NAG ~~ AP
BG ~~ dim1
AP ~~ dim3
BG ~~ dim3
NAG ~~ dim1
LAP ~~ dim3
dim1 ~~ dim3
'

## Fit model on full gradient
fit <- lavaan::sem(full_model,data=vegan::decostand(datsem[,-c(16,17)], 'standardize'))
lavaan::summary(fit, standardize = T, rsq=T, fit.measures=T)

###
#### Individual SEMs ----
###

### BG
bg_model <- '
#latent variable definitions
com =~ dim1 + dim3

# regression
n0 ~ tmmn_lt + ph
TotalAB ~ com + n0
com ~ n0 + tmmn_lt
BG ~ TotalAB + com + n0

#residual correlations
BG ~~ dim1
BG ~~ dim3
dim1 ~~ dim3
'

## Fit model on full gradient
fitbg <- lavaan::sem(bg_model,data=vegan::decostand(datsem[,-c(16, 17)], 'standardize'))
lavaan::summary(fitbg, standardize = T, rsq=T, fit.measures=T)

### AP
ap_model <- '
#latent variable definitions
com =~ dim1 + dim3

# regression
n0 ~ tmmn_lt + ph
TotalAB ~ com + n0
com ~ n0 + tmmn_lt
AP ~ TotalAB + com + n0

#residual correlations
AP ~~ dim3
'

## Fit model on full gradient
fitap <- lavaan::sem(ap_model,data=vegan::decostand(datsem[,-c(16, 17)], 'standardize'))
lavaan::summary(fitap, standardize = T, rsq=T, fit.measures=T)


### LAP
lap_model <- '
#latent variable definitions
com =~ dim1 + dim3

# regression
n0 ~ tmmn_lt + ph
TotalAB ~ com + n0
com ~ n0 + tmmn_lt
LAP ~ TotalAB + com + n0

#residual correlations
LAP ~~ dim3
'

## Fit model on full gradient
fitlap <- lavaan::sem(lap_model,data=vegan::decostand(datsem[,-c(16, 17)], 'standardize'))
lavaan::summary(fitlap, standardize = T, rsq=T, fit.measures=T)

### NAG
nag_model <- '
#latent variable definitions
com =~ dim1 + dim3

# regression
n0 ~ tmmn_lt + ph
TotalAB ~ com + n0
com ~ n0 + tmmn_lt
NAG ~ TotalAB + com + n0

#residual correlations
NAG ~~ dim1
dim1 ~~ dim3
'

## Fit model on full gradient
fitnag <- lavaan::sem(nag_model,data=vegan::decostand(datsem[,-c(16, 17)], 'standardize'))
lavaan::summary(fitnag, standardize = T, rsq=T, fit.measures=T)

###
#### Smoothed SEM modelling ----
###

# Run the moving window SEM function with a window of 40 consecutive samples along the environmental gradient
mov.sem = movingsem(datsem, 85, full_model)

ggplot(mov.sem,aes(x=av.X,y=lat75)) + 
  geom_point(size = 3) + geom_line() + 
  ylab("Enzyme multifunctionality (75th quartile Latent Variable)") + xlab("Averaged ASV richness") + theme_classic() + scale_x_reverse()


###
#### Sensitivity ananlysis | Randomized moving SEM model - 1,000 iterations ----
###

#### loop it
run  <- 1:1000
output <- lapply(X = run, FUN = function(x) {   
  #random order N0
  rand <- sample(nrow(datsem))
  datsem.r <- datsem[rand, ]
  
  #run moving window
  mwsem = movingsem(datsem.r, 85, full_model)
  
  #add iteration
  mwsem <- as.data.frame(mwsem)
  mwsem$run <- run[[x]]
  # Save output  
  my_list <- list(movsem = mwsem)   
  return(my_list)  
})

head(output)

plot(output[[100]]$movsem$lat.max~output[[100]]$movsem$av.X, type="o", las=1, pch=16, cex=2, ylab = "75th quartile Latent Variable Enzymes", xlab = "averaged n0")
test <- summary(lm(output[[100]]$movsem$lat.max~output[[100]]$movsem$av.X))
flat_results <- do.call(rbind, output)
df_list <- as.list(flat_results[, "movsem"])
all_randomSEM <- do.call(rbind, df_list)

### --- plot the results of randomization ---

#Run lm() on every item in the list and extract the slope
random_ses <- sapply(output, function(iteration) {
  model <- summary(lm(lat.max ~ av.X, data = iteration$movsem))
  rho <- sqrt(model$adj.r.squared)*sign(model$coefficients[2])
  return(res(rho, n = 32)$fisher.z) # This grabs just the richness coefficient
})

boxplot(random_ses)
hist(random_ses)

#Get observed slope from non-randomized SEM data
obs_model <- summary(lm(mov.sem$lat.max ~ mov.sem$av.X))
rho_obs <- sqrt(obs_model$adj.r.squared)*sign(obs_model$coefficients[2])
obs_ses <- res(rho_obs, n = 32)$fisher.z


#Create a data frame for plotting
plot_data <- data.frame(Group = "Randomized (Null)", SES = random_ses)

#Boxplot
ggplot(plot_data, aes(x = Group, y = SES)) +
  geom_jitter(color = "gray40", width = 0.05, alpha = 0.2) +
  geom_boxplot(fill = "white", width = 0.1, outlier.alpha = 0.3) +
  # Add the observed SES as a large red point
  geom_point(aes(x = 1, y = obs_ses), color = "red", size = 5) +
  # Label the observed point
  annotate("text", x = 1.2, y = obs_ses, label = "Observed", color = "red", hjust = 0) +
  labs(y = "Standardized Effect Size (SES)",
       x = "") +
  theme_classic()