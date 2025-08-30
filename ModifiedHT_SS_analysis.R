# Arguments:
# - siteCovs: data frame containing sites with their respective counts
# - sample_sizes: different sample size for resampling 
# - n_bootstrap: number of bootstrap iterations
# - logit_p: estimated p from auxiliary data (on logit scale) from main script
# - se_logit: SE of the same estimated p from main script
# - plot_size: size of a single plot (km²)
# - study_area: total study area (km²)
# - plot_size: size of a single plot (km²)

# Function to bootstrap population estimates across sample sizes
ModifiedHT_SS_analysis <- function(siteCovs, 
                                   sample_sizes, 
                                   n_bootstrap = 1000, 
                                   logit_p, 
                                   se_logit, 
                                   plot_size = 0.25, 
                                   study_area = 700) {
  set.seed(123)
  
  #create list to store results
  results <- vector("list", length(sample_sizes) * n_bootstrap)
  idx <- 1
  for (n in sample_sizes) {
    for (b in 1:n_bootstrap) {
      #Resample sites
      sampled_data <- siteCovs[sample(1:nrow(siteCovs), size = n, replace = TRUE), ]
      C_total <- sum(sampled_data$V1.observer.2, na.rm = TRUE)
      #Sample detection probability from logit-normal distribution
      p_sample <- rnorm(1, mean = logit_p, sd = se_logit)
      prob_detection <- plogis(p_sample)
      #Population estimate
      N_total <- C_total / (prob_detection * (n * plot_size / study_area))
      results[[idx]] <- data.frame(SampleSize = n, Estimate = N_total)
      idx <- idx + 1
    }
  }
  return(do.call(rbind, results))
}

