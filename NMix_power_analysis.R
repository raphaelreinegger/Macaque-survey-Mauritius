# Arguments:
# - fmp: fitted pcount model with Poisson distribution from main script
# - alpha_hat: overdispersion parameter estimated by pcount model fitted with NB distribution
# - site_Covs: data frame with site covariates from main script
# - obs_Covs: data frame with observation covariates from main script
# - n_visits: matrix of sites with either single or two visits. In our case 19 sites (31%) have two visits.
# - decline_rates: annual decline rates (e.g. 0.02 and 0.05)
# - n_years_max: max number of years to simulate
# - ncores: number of parallel cores to be used
# - alpha: significance level for power test

library(unmarked)
library(doParallel)
library(foreach)

NMix_power_analysis <- function(fmp,
                                alpha_hat = exp(-0.179),  # overdispersion parameter from pcount model with NB mixture
                                site_Covs,
                                obs_Covs,
                                n_visits = rep(c(2, 1), times = c(19, 43)),  # 31% of sites revisited
                                decline_rates = c(0.02, 0.05, 0.1, 0.2, 0.5),
                                n_years_max = 10,
                                n_bootstrap = 100,
                                ncores = 12,
                                alpha = 0.05) {
  set.seed(123)
  
  # Set up parallel backend
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  #Extract site-level predictions for abundance and detection probability 
  #from top-ranked Poisson N-mixture model with fitted in the main script (fmp)
  lambda_hat <- predict(fmp, type = "state")$Predicted
  lambda_se  <- predict(fmp, type = "state")$SE
  det_hat    <- predict(fmp, type = "det")$Predicted
  det_se     <- predict(fmp, type = "det")$SE
  n_sites <- length(lambda_hat)
  
  #create data frame to store simulation results
  power_results <- data.frame()
  
  #make sure outer loops stay in order
  for (decline_rate in decline_rates) {
    for (n_years in 1:n_years_max) {
      #start parallel simulation loop
      pvals <- foreach(b = 1:n_bootstrap, .combine = c,
                       .packages = c("unmarked")) %dopar% {
                         #Sample site-specific true abundances and detection probabilities
                         log_lambda <- log(lambda_hat)
                         log_lambda_sim <- rnorm(n_sites, mean = log_lambda, sd = lambda_se / lambda_hat)
                         lambda_site <- exp(log_lambda_sim)
                         logit_det <- qlogis(det_hat)
                         logit_det_sim <- rnorm(n_sites, mean = logit_det, sd = det_se)
                         det_site <- plogis(logit_det_sim)
                         det_site <- pmax(pmin(det_site, 0.999), 0.001) #place limits to avoid silly probabilities
                         #Create data frame to store mean abundances estimated by new N-mixture models
                         mean_abundances <- numeric(n_years + 1)
                         for (year in 0:n_years) {
                           #Expected abundance decline for a specific year within specific scenario
                           lambda_year <- lambda_site * (1 - decline_rate)^year
                           #Simulate the true N per site for new scenario and year
                           #introducing overdispersion similar to actual data
                           true_N <- rnbinom(n_sites, mu = lambda_year, size = 1/alpha_hat)
                           #Simulate observed counts for each visit
                           y_mat <- matrix(NA, nrow = n_sites, ncol = 2)
                           for (j in 1:n_sites) {
                             #count for 1st visit
                             y_mat[j, 1] <- rbinom(1, size = true_N[j], prob = det_site[j])
                             #if site also has 2nd visit, keep same true_N and det_site as first visit
                             if (n_visits[j] == 2) {
                               y_mat[j, 2] <- rbinom(1, size = true_N[j], prob = det_site[j])
                             }
                           }
                           #Create observation covariates matrix 
                           obs_covs_year <- lapply(obs_Covs, function(x) {
                             x_mat <- matrix(x, nrow = n_sites, ncol = 2, byrow = FALSE)
                             as.vector(t(x_mat))
                           })
                           obs_covs_df <- as.data.frame(obs_covs_year)
                           #Build unmarkedFrame and fit new pcount model
                           umf_sim <- unmarkedFramePCount(y = y_mat,
                                                          siteCovs = site_Covs,
                                                          obsCovs = obs_covs_df)
                           fmp_sim <- pcount(~ NDVI_mean ~ Dist_agriculture + Cap.y.km2 + Dist_chass,
                                             data = umf_sim, mixture = "P")
                           #Extract mean estimated abundance across sites
                           lambda_pred <- predict(fmp_sim, type = "state")$Predicted
                           mean_abundances[year + 1] <- mean(lambda_pred)
                         }
                         #Fit LM to the log-transformed mean plot-level abundance with year as predictor
                         df_lm <- data.frame(
                           Year = 0:n_years,
                           LogMeanAbundance = log(mean_abundances + 1)
                         )
                         lm_fit <- lm(LogMeanAbundance ~ Year, data = df_lm)
                         summary(lm_fit)$coefficients[2, 4]
                       }
      #Estimate power (proportion of significantly negative 'Year' effects estimated by LM
      power <- mean(pvals < alpha, na.rm = TRUE)
      power_results <- rbind(power_results, data.frame(
        DeclineRate = decline_rate,
        Years = n_years,
        Power = power
      ))
      cat("Decline rate:", decline_rate, "- Years:", n_years,
          "Power:", round(power, 3), "\n")
    }
  }
  # Stop cluster
  stopCluster(cl)
  return(power_results)
}
