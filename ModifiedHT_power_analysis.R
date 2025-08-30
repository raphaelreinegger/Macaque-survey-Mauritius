# Arguments:
# - mean_hat: estimated mean abundance from GLM in main script
# - dispersion: dispersion parameter from GLM (NB) in main script
# - logit_p: estimated detection probability on logit scale (in main script)
# - se_logit: SE of logit_p (in main script)
# - n_plots: number of survey plots
# - n_years_max: max number of years to simulate
# - n_bootstrap: bootstrap reps per scenario
# - decline_rates: decline rates per year
# - total_area: total study area (km²)
# - plot_area: area of a single plot (km²)
# - alpha: significance level for power test

ModifiedHT_power_analysis <- function(mean_hat,
                           dispersion,
                           logit_p,
                           se_logit,
                           n_plots = 62,
                           n_years_max = 10,
                           n_bootstrap = 1000,
                           decline_rates = c(0.02, 0.05, 0.10, 0.2, 0.5),
                           total_area = 700,
                           plot_area = 0.25,
                           alpha = 0.05) {
  set.seed(123)
  
  #area fraction per plot
  area_fraction <- (n_plots * plot_area) / total_area
  
  #data.frame for results
  power_results <- data.frame()
  
  #start simulation
  for (decline_rate in decline_rates) {
    for (n_years in 1:n_years_max) {
      pvals <- numeric(n_bootstrap)
      for (i in 1:n_bootstrap) {
        #simulate detection probability based on previous study (Reinegger et al. 2025)
        p_sample <- rnorm(1, logit_p, se_logit)
        detection_prob <- plogis(p_sample)
        #create data frame to store HT estimators for each year within specific scenario
        df <- data.frame()
        for (year in 0:n_years) {
          #expected abundance decline for specific year within scenario
          mean_year <- mean_hat * (1 - decline_rate)^year
          #draw new counts from negative binomial distribution parameterized by mean_year and dispersion
          counts <- rnbinom(n_plots, mu = mean_year, size = dispersion)
          #compute modified horvitz-thompson estimator
          N_hat <- sum(counts) / (area_fraction * detection_prob)
          #store estimate
          df_year <- data.frame(
            HT = N_hat,
            Year = year
          )
          df <- rbind(df, df_year)
        }
        #fit linear model on log-transformed HT estimator to estimate effect of Year
        df$logHT <- log(pmax(df$HT, 1))  # small floor to avoid -Inf
        lm_model <- lm(logHT ~ Year, data = df)
        pvals[i] <- summary(lm_model)$coefficients[2, 4]
      }
      #estimate power (proportion of simulations where year has significantly negative effect)
      power <- mean(pvals < alpha)
      power_results <- rbind(power_results, data.frame(
        Years = n_years,
        DeclineRate = decline_rate * 100,
        Power = power
      ))
    }
  }
  return(power_results)
}
