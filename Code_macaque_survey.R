#############################################################################
#############################################################################
##### ATTENTION: 'Cap.y.km2' are capture data and have been redacted    #####
##### from the csv-files, as the authors do not own this data.          #####
##### However, all models and simulations can run without 'Cap.y.km2'   #####
##### and produce similar results to those presented in the manuscript  #####
##### All code below should run fine after removing 'Cap.y.km2'         #####
##### from the relevant sections                                        #####
#############################################################################
#############################################################################

## Set working directory first

# Clear plots
if(!is.null(dev.list())) dev.off()
# Clear console                                                           
cat("\014") 
# Clean workspace
rm(list=ls())
# Set working directory
setwd("C:/path")
#load "HighstatLibV10.R.txt" from Zuur et al. (2009) to obtain function corvif, pairs, panel.smooth2 and panel.cor
#(only necessary for data exploration, not for running models or simulations)
source("C:/path/HighstatLibV10.R.txt")

#load packages
library(ggplot2)
library(MASS)
library(dplyr)
library(DHARMa)
library(ggpubr)
library(AICcmodavg)
library(data.table)
library(unmarked)
library(nmixgof)
library(parallel)
library(openxlsx)

#######################################################
#######################################################
#######################################################
##########                                 ############
##########        exploration              ############
##########                                 ############
#######################################################
#######################################################
#######################################################
## attach data 
siteCovs <- read.csv("siteCovs_submission.csv",na.strings=c("", "NA"),header=T,sep=",") 
obsCovs <- read.csv("obsCovs_submission.csv",na.strings=c("", "NA"),header=T,sep=",")
attach(siteCovs)
z <- cbind(Max.abundance, Dist_agriculture, Dist_river, Cap.y.km2,Dist_chass)
#now quick visualisation of relationships between variables
pairs(z, lower.panel = panel.smooth2,upper.panel = panel.cor, diag.panel = panel.hist)
#check multicollinearity
corvif(z) #highest vif values are elev_mean and Scrub_frac because they are highly correlated
#remove scrub_frac
z <- cbind(Max.abundance, Dist_agriculture, Dist_river, Dist_roads, Elev_mean, Dist_chass, Cap.y.km2)
corvif(z) #perfect
#visualize relationships
ggplot(siteCovs, aes(x = Dist_agriculture, y = Max.abundance)) +
  geom_point(color = "grey") + geom_smooth(method='lm', color = "black", se = FALSE) + 
  labs(y = "Individuals per 25 ha (n)", x = "Distance from agriculture (m)") +
  theme_classic()
ggplot(siteCovs, aes(x = Dist_chass, y = Max.abundance)) +
  geom_point(color = "grey") + geom_smooth(method='lm', color = "black", se = FALSE) + 
  labs(y = "Individuals per 25 ha (n)", x = "Distance to hunting area (m)") +
  theme_classic()
ggplot(siteCovs, aes(x = Dist_river, y = Max.abundance)) +                                 #don't use length_rivers
  geom_point(color = "grey") + geom_smooth(method='lm', color = "black", se = FALSE) + 
  labs(y = "Individuals per 25 ha (n)", x = "Total river length within 117 ha (m)") +
  theme_classic()
ggplot(siteCovs, aes(x = Cap.y.km2, y = Max.abundance)) +
  geom_point(color = "grey") + geom_smooth(method='lm', color = "black", se = FALSE) + 
  labs(y = "Individuals per 25 ha (n)", x = "Mean captured macaques/year/km2 (n)") +
  theme_classic()

#######################################################
#######################################################
#######################################################
##########                                 ############
##########      N-mixture modelling        ############
##########                                 ############
#######################################################
#######################################################
#######################################################

##############################
### global model fitting #####
##############################
#first standardize all variables (substract mean and divide by SD)
siteCovs <- siteCovs %>% mutate(across(c(Dist_agriculture, Dist_river, Cap.y.km2, Dist_chass), scale))
#siteCovs has one row per site and one column per covariate
#y has one row per site and two columns containing the counts of each visit
#obsCovs should have n = site*visit*observer rows

y <- siteCovs[,c(7:8)]
obs_Covs <- obsCovs[obsCovs$Observer %in% c("Observer2"),]
umfp <- unmarkedFramePCount(y, siteCovs=site_Covs, obsCovs=obs_Covs)
#create model
fmp <- pcount(~ NDVI_mean #with NDVI_mean instead of Temperature you get more sensible estimates
              ~ Dist_agriculture + Dist_chass + Cap.y.km2 + Dist_river, data = umfp, mixture=c("ZIP"))
#evaluate model fit using goodness of fit residual plots from packages nmixgof
residfit(fmp, "marginal") #looks good
residfit(fmp, "site-sum") #looks good
residfit(fmp, "observation") #looks good
#qqplots from package nmixgof
residqq(fmp, "site-sum") #evidence for overdispersion in Poisson and ZIP models
residqq(fmp, "observation") #evidence for overdispersion in Poisson and ZIP models
#check overdispersion with test 
obs.boot <- Nmix.gof.test(fmp, nsim = 1000, ncores = 12)
obs.boot #overdispersion -> estimate of c-hat = 17.9 for Poisson, 5.1 for ZIP and 1.35 for NB 
#(values > 4 indicate potential lack of fit)
#null model test
fmpnull <- pcount(~ NDVI_mean  
                  ~ 1, data = umfp, mixture=c("ZIP"))
LRT(fmp, fmpnull) #significant
#general inference
summary(fmp) #reasonable estimates
confint(fmp, type = 'state') #get predictor confidence intervals for abundance model
#use function at end of script to adjust estimates and SEs by c-hat
c_hat <- 5.1 #for ZIP model
extract_adjust_pcount(fmp, type = "state", c_hat = c_hat) #run function at end of script
extract_adjust_pcount(fmp, type = "det",   c_hat = c_hat)
#get average predictions of both abundance and detection for each plot 
mean(predict(fmp, type="state")$Predicted, na.rm = TRUE) 
mean(predDet <- predict(fmp, type="det")$Predicted, na.rm = TRUE) 

##############################
###  model optimization  #####
##############################
#do a model selection procedure with pcount to find the most plausible candidate models
#always include Dist_agriculture and Capture
fmp1 <- pcount(~ 1 ~ Dist_agriculture + Cap.y.km2 + Dist_chass + Dist_river, data = umfp, mixture = "ZIP")
fmp2 <- pcount(~ 1 ~ Dist_agriculture + Cap.y.km2 + Dist_chass,data = umfp, mixture = "ZIP")
fmp3 <- pcount(~ 1 ~ Dist_agriculture + Cap.y.km2,data = umfp, mixture = "ZIP")
#make model list and compare AIC
models <- fitList(fmp1, fmp2, fmp3)
modSel(models) 
##use function at end of script to calculate QAICc, so scroll down and run function first
models_list <- list(fmp1, fmp2, fmp3)
c_hat_value <- 5.1 #5.1 for ZIP and 1.35 for the NB model
number_of_sites <- nrow(umfp@siteCovs)
compute_QAICc(models_list, c_hat_value, number_of_sites)


##############################
### extrapolation to grid ####
##############################
#load the grid into R and predict for the whole island
grid_data <- read.csv("Grid_submission.csv")
#make sure to scale the predictors here too
grid_data <- grid_data %>% mutate(across(c(Dist_agriculture, Dist_river, Cap.y.km2, Dist_chass), scale))
grid_data$pred_abundance <- predict(fmp1, type = "state", newdata = grid_data)$Predicted
#divide by the area of each cell to get final prediction
grid_data$pred_adjusted <- grid_data$pred_abundance*(grid_data$Area/250000)
grid_data$pred_adjusted <- grid_data$pred_adjusted*(25/(35*(1-0.2)))
sum(grid_data$pred_adjusted, na.rm= TRUE) #for quick result

#use a parametric bootstrap to get CIs for the population estimate, using function parboot
#first build function for extrapolating to full habitat extent grid generated with MaxEnt model
total_abundance_fnalt <- local({
  gd <- grid_data  # this captures grid_data from your global env
  function(fitted_model) {
    lambda_pred <- predict(fitted_model, type = "state", newdata = gd)
    pred_abundance <- lambda_pred$Predicted
    adjusted <- pred_abundance * (gd$Area / 250000)
    adjusted <- adjusted*(25/(35*(1-0.2)))
    sum(adjusted, na.rm = TRUE)
  }
})

#Now parametric bootstrap using parboot to compute confidence intervals for total abundance
set.seed(123)
pb <- parboot(fmp1, statistic = total_abundance_fnalt, nsim = 1000, ncores = 12, parallel = TRUE) 
#1000 simulations takes around 4.5 min with P and 5.5 min with ZIP distribution
#Inspect bootstrap results
hist(pb@t.star, breaks = 30, main = "Bootstrap Population Estimates", xlab = "Total Abundance") #distribution of estimates
mean(pb@t.star) #average estimate
quantile(pb@t.star, probs = c(0.025, 0.975)) #confidence intervals

#also a function for extrapolating to smaller habitat extent (450 km2) --> does not need to access the grid_data
total_abundance_fn2 <- function(fitted_model) {
  pred_abundance <- mean(predict(fitted_model, type="state")$Predicted, na.rm = TRUE)
  pred_abundance <- pred_abundance*(25/(35*(1-0.2)))
  pred_abundance * (45000/25)
}

#to get the mean predicted plot-level abundance 
mean_lambda <- function(fit) {
  mean(predict(fit, type = "state")$Predicted)
}

#to get the mean predicted abundance for each cell for export to QGIS
make_grid_stat <- function(grid) {
  function(fit) {
    pred <- predict(fit, type = "state",  newdata = grid)$Predicted
    pred_adj <- pred * (grid$Area / 250000)
    pred_adj <- pred_adj * (25 / (35 * (1 - 0.2)))
    return(pred_adj)
  }
}
grid_stat <- make_grid_stat(grid_data)
set.seed(123)
pb <- parboot(fmp2, statistic = grid_stat, nsim = 1000, ncores = 12, parallel = TRUE)
grid_data$pred_adjusted <- colMeans(pb@t.star, na.rm = TRUE)
#export to excel so you can visualize in QGIS
library(openxlsx)
write.xlsx(grid_data, 'grid_predictions.xlsx') 
#when plotting these in QGIS, do not forget to multiply by 4 to get n/km2



#########################################
### sample size sensitivity analysis ####
#########################################
source("NMix_SS_analysis.R")

results_all <- NMix_SS_analysis(
  sample_sizes = c(20,30,40,50,60),
  reps_per_n = 500, #Reduce for faster testing. Running with 500 reps_per_n and 100 nsim_parboot takes ~15 hours 
  ncores = 12, #Number of cores for parallel processing, reduce or increase depending on pc
  nsim_parboot = 100,
  umfp = umfp, 
  grid_data = grid_data
)

summary_stats <- aggregate(. ~ SampleSize, data = results_all[, c("SampleSize", "Mean", "LowerCI", "UpperCI")], FUN = mean)
SS_NMix <- ggplot(summary_stats, aes(x = SampleSize, y = Mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 1) +
  labs(x = "Sample size (n)", y = "Estimated population size (N)") +
  theme_classic() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12))



#########################################
#######     power analysis    ###########
#########################################
source("NMix_power_analysis.R")

power_results <- NMix_power_analysis(
  fmp = fmp,
  site_Covs = site_Covs,
  obs_Covs = obs_Covs,
  n_bootstrap = 1000,    #takes about ~30 hours
  ncores = 12
)

# Smooth and plot the results. 
#make labels for the declines
label_positions <- tibble(
  DeclineRate = c(2, 5, 10, 20, 50),
  label = c("2%", "5%", "10%", "20%", "50%"),
  x = c(8, 7.5, 6.1, 4.5, 2.9),   # manually chosen x positions
  y = c(0.08, 0.20, 0.41, 0.55, 0.62)  # manually chosen y positions
)

power1 <- ggplot(power_results, aes(x = Years, y = Power, group = as.factor(DeclineRate))) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey70") +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5),
              se = FALSE, size = 0.8, color = "black") +
  geom_label(data = label_positions, aes(x = x, y = y, label = label), size = 4,
             fill = "white", label.size = 0, label.padding = unit(0.15, "lines"),
             inherit.aes = FALSE) +
  annotate("text", x = 1.5, y = 1.1, label = "Annual decline (%)", hjust = 0, size = 4) +
  scale_x_continuous(breaks = 2:10, limits = c(2,10)) +
  scale_y_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1.1)) +
  labs(x = "Years monitored (n)", y = "Estimated power") + theme_classic() +
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 13),legend.position = "none")


#######################################################
#######################################################
#######################################################
########                                     ##########
######## Modified Horvitz-Thompson estimator ##########
########                                     ##########
#######################################################
#######################################################
#######################################################


###########################################
###########################################
#####    systematic sampling design   #####
###########################################
###########################################
#use auxiliary data from Reinegger et al. (2025) to adjust HT estimator for imperfect detection
#estimated p for selected flight parameters = 0.69 (response scale) with SE = 0.24 (link scale)
#quick and dirty population estimate
C_total <- sum(siteCovs$V1.observer.2,na.rm = TRUE)
area_fraction <- (62 * 0.28) / 608 #in km2, adjusting for ESA
detection_prob <- 0.68703196 #assuming error is similar across all sites
N_total <- C_total / (area_fraction * detection_prob) 

#Now estimate with non-parametric bootstrap and error propagation 
logit <- function(p) log(p / (1 - p))
#estimated detection and SE from auxiliary data 
p_hat <- 0.68703196
se_logit <- 0.2419754
#Convert detection probability to logit scale
logit_p <- logit(p_hat)
#simulate detection on logit scale, then convert to response scale
set.seed(123)
n_bootstrap <- 1000
area_fraction <- (62 * 0.28) / 608 #in km2, adjusting for ESA
bootstrap_estimates <- numeric(n_bootstrap)

for (i in 1:n_bootstrap) {
  # Resample plots
  bootstrap_data <- siteCovs[sample(1:nrow(siteCovs), replace = TRUE), ]
  C_total_bootstrap <- sum(bootstrap_data$V1.observer.2, na.rm = TRUE)
  # Simulate p from normal on logit scale, then inverse-logit back to response scale
  p_sample_logit <- rnorm(1, mean = logit_p, sd = se_logit)
  p_sample <- plogis(p_sample_logit)
  # Horvitz-Thompson estimate adjusted for detection and area
  bootstrap_estimates[i] <- C_total_bootstrap/ (area_fraction * p_sample)
}

#mean boostrap estimate with SD
mean(bootstrap_estimates) #mean
sd(bootstrap_estimates)/sqrt(1000) #SD 
# Calculate confidence intervals (95% CI) from the bootstrap samples
quantile(bootstrap_estimates, 0.025) #lower confidence limit
quantile(bootstrap_estimates, 0.975) #upper confidence limit


#########################################
### sample size sensitivity analysis ####
#########################################
source("ModifiedHT_SS_analysis.R")
bootstrap_results <- ModifiedHT_SS_analysis(siteCovs, 
                                            sample_sizes = c(10, 20, 30, 40, 50, 60), 
                                            n_bootstrap = 1000,
                                            logit_p = logit_p, 
                                            se_logit = se_logit)
#plot means with 95% CIs per sample size
SS_HT <- ggplot(bootstrap_results, aes(x = factor(SampleSize), y = Estimate)) +
  stat_summary(fun = mean, geom = "point", color = "black", size = 3) +
  stat_summary(fun.data = function(x) {
    data.frame(y = mean(x), ymin = quantile(x, 0.025), ymax = quantile(x, 0.975))
  }, geom = "errorbar", width = 0.2, color = "black") +
  labs(x = "Sample size (n)", y = "Estimated population size (N)") +
  theme_classic() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12))


#########################################
#######      power analysis     #########
#########################################
#Get estimated mean and dispersion for the raw counts
nb_fit <- glm.nb(Max.V1 ~ 1, data = y)
#validate GLM
simulationOutput <- simulateResiduals(fittedModel = nb_fit, n = 1000)
plotSimulatedResiduals(simulationOutput = simulationOutput) #perfect
#now save the parameters
mean_hat <- exp(coef(nb_fit))         # Estimated mean, same as observed 14.7
dispersion <- nb_fit$theta            # Estimated dispersion parameter 0.6

source("ModifiedHT_power_analysis.R")

#Run simulation
power_res <- ModifiedHT_power_analysis(n_bootstrap = 1000, 
                                       mean_hat = mean_hat, 
                                       dispersion = dispersion, 
                                       logit_p = logit_p, 
                                       se_logit = se_logit)

#make labels for the declines
label_positions2 <- tibble(
  DeclineRate = c(2, 5, 10, 20, 50),
  label = c("2%", "5%", "10%", "20%", "50%"),
  x = c(8.2, 6.7, 4.9, 2.8, 1.5),   # manually chosen x positions
  y = c(0.17, 0.46, 0.7, 0.7, 0.91)  # manually chosen y positions
)
#smooth and plot results
power2 <- ggplot(power_res, aes(x = Years, y = Power, group = as.factor(DeclineRate))) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey70") +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5),
              se = FALSE, size = 0.8, color = "black") +
  geom_label(data = label_positions2, aes(x = x, y = y, label = label), size = 4,
             fill = "white", label.size = 0, label.padding = unit(0.15, "lines"),
             inherit.aes = FALSE) +
  annotate("text", x = 1.5, y = 1.1, label = "Annual decline (%)", hjust = 0, size = 4) +
  scale_x_continuous(breaks = 2:10, limits = c(2,10)) +
  scale_y_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1.1)) +
  labs(x = "Years monitored (n)", y = "Estimated power") + theme_classic() +
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 13),legend.position = "none")

##############################################
### plot all ss and power analyses together ##
##############################################
void1 <- ggplot() +
  theme_void()
#first only sensitivity or power analyses
ggarrange(void1, void1, SS_NMix, SS_HT, ncol = 2, nrow = 2, labels = c("A", "B", "", ""), 
          label.x = 0.05, heights = c(0.08, 1))
#now everything together
ggarrange(void1, void1, SS_NMix, SS_HT, void1, void1, power1, power2,
          ncol = 2, nrow = 4, labels = c("A", "B", "", "", "C", "D", "", ""), 
          label.x = 0.05, heights = c(0.08, 1, 0.08, 1))









#######################################################
#######################################################
## function for computing QAICc for a list of models ##
#######################################################
#######################################################
compute_QAICc <- function(models, c_hat, n_sites) {
  # models: list of fitted unmarked models
  # c_hat: single numeric value of ĉ applied to all models
  # n_sites: integer, number of sites in dataset
  
  results <- data.frame(
    Model = character(length(models)),
    K = numeric(length(models)),
    logLik = numeric(length(models)),
    c_hat = numeric(length(models)),
    QAICc = numeric(length(models)),
    stringsAsFactors = FALSE
  )
  
  for(i in seq_along(models)) {
    mod <- models[[i]]
    
    K <- length(coef(mod))
    logLik_val <- as.numeric(logLik(mod))
    
    QAICc <- (-2 * logLik_val) / c_hat + 2 * K + (2 * K * (K + 1)) / (n_sites - K - 1)
    
    results[i, ] <- list(
      Model = paste0("Model_", i),
      K = K,
      logLik = logLik_val,
      c_hat = c_hat,
      QAICc = QAICc
    )
  }
  
  results <- results[order(results$QAICc), ]
  rownames(results) <- NULL
  return(results)
}

#######################################################
#######################################################
##     function for adjust SEs and CIs by c-hat      ##
#######################################################
#######################################################
extract_adjust_pcount <- function(model, type = c("state", "det"), c_hat = 1, conf = 0.95) {
  type <- match.arg(type)
  # Z-value for CI
  z <- qnorm((1 + conf) / 2)
  # Extract estimates and SEs
  est <- coef(model, type = type)
  se  <- SE(model, type = type)
  # Adjust SE for overdispersion
  adj_se <- se * sqrt(c_hat)
  # Compute adjusted CIs
  adj_low  <- est - z * adj_se
  adj_high <- est + z * adj_se
  # Return as data frame
  data.frame(
    term     = names(est),
    estimate = est,
    orig_SE  = se,
    adj_SE   = adj_se,
    adj_CI_low  = adj_low,
    adj_CI_high = adj_high,
    row.names = NULL
  )
}

