# Arguments:
# - sample_sizes: different sample size for resampling 
# - reps_per_n: number of times sites are resampled
# - ncores: number of cores to be used for parallel processing
# - nsim_parboot: number of parametric bootstrap iterations for estimating population size at grid level
# - umfp: fitted umfp object in the original pcount model 'fmp'
# - grid_data: 25-ha resolution grid covering all macaque habitat in Mauritius

require(unmarked)
require(foreach)
require(doParallel)

NMix_SS_analysis    <- function(sample_sizes = c(10,20,30,40,50,60),
                                reps_per_n = 500,
                                ncores = 12,
                                nsim_parboot = 1000,
                                umfp,
                                grid_data) {
  
  #Prepare full dataset based on originally fitted umfp object in pcount model 'fmp'
  full_y <- getY(umfp)
  site_covs_full <- siteCovs(umfp)
  obs_covs_full <- obsCovs(umfp)
  num_occasions <- ncol(full_y)
  #Identify replicated vs single-visit sites
  num_visits_per_site <- rowSums(!is.na(full_y))
  replicated_sites <- which(num_visits_per_site > 1)
  single_visit_sites <- which(num_visits_per_site == 1)
  #Allow parallel processing
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl), add = TRUE)  #ensures cluster closes when function exits
  #Create data frame to store results
  results <- data.frame()
  
  set.seed(123)
  #Loop over the different sample sizes
  for (n in sample_sizes) {
    cat("Sample size:", n, "\n")
    n_rep <- round(n * 0.31)  #makes sure 31% of reps are always two-visit sites, reflecting the sampling design
    n_single <- n - n_rep
    results_n <- foreach(rep = 1:reps_per_n, .combine = rbind,
                         .packages = c("unmarked")) %dopar% {
                           #Check if enough sites are available
                           if (n_rep > length(replicated_sites) || n_single > length(single_visit_sites)) return(NULL)
                           #Draw sites randomly without replacement: 31% has two be 2-visit sites, rest 1-visit
                           selected_rep <- sample(replicated_sites, n_rep, replace = FALSE)
                           selected_single <- sample(single_visit_sites, n_single, replace = FALSE)
                           selected <- c(selected_rep, selected_single)
                           #Subset data, making sure the right covariate values are matched with the right y value
                           y_sub <- full_y[selected, , drop = FALSE]
                           site_covs_sub <- site_covs_full[selected, , drop = FALSE]
                           obs_covs_sub <- lapply(obs_covs_full, function(x) {
                             x_mat <- matrix(x, nrow = nrow(full_y), ncol = num_occasions, byrow = FALSE)
                             x_sub <- x_mat[selected, , drop = FALSE]
                             as.vector(t(x_sub))
                           })
                           obs_covs_df <- as.data.frame(obs_covs_sub)
                           #Create new unmarked frame with the resampled sites
                           umf_sub <- try(unmarkedFramePCount(
                             y = y_sub,
                             siteCovs = site_covs_sub,
                             obsCovs = obs_covs_df
                           ), silent = TRUE)
                           if (inherits(umf_sub, "try-error")) return(NULL)
                           #Fit pcount model
                           fmp_sub <- try(pcount(~ NDVI_mean ~ Dist_agriculture + Dist_chass + Dist_roads +
                                                   Cap.y.km2 + Elev_mean + Dist_river,
                                                 data = umf_sub, mixture = "P"), silent = TRUE)
                           if (inherits(fmp_sub, "try-error")) return(NULL)
                           #Grid extrapolation function
                           total_abundance_fn <- function(fitted_model) {
                             lambda_pred <- predict(fitted_model, type = "state", newdata = grid_data)
                             pred_abundance <- lambda_pred$Predicted
                             adjusted <- pred_abundance * (grid_data$Area / 250000)
                             sum(adjusted, na.rm = TRUE)
                           }
                           #Use parametric bootstrap to extrapolate abundance to grid
                           pb <- try(parboot(fmp_sub,
                                             statistic = total_abundance_fn,
                                             nsim = nsim_parboot,
                                             parallel = TRUE, ncores = 1),
                                     silent = TRUE)
                           if (inherits(pb, "try-error")) return(NULL)
                           #Summarize results
                           mean_est <- mean(pb@t.star, na.rm = TRUE)
                           lower_ci <- quantile(pb@t.star, 0.025, na.rm = TRUE)
                           upper_ci <- quantile(pb@t.star, 0.975, na.rm = TRUE)
                           data.frame(SampleSize = n, Replicate = rep,
                                      Mean = mean_est, LowerCI = lower_ci, UpperCI = upper_ci)
                         }
    results <- rbind(results, results_n)
  }
  return(results)
}