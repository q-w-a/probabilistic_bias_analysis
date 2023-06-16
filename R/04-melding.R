
get_melded <- function(alpha_mean = 0.9,
                       alpha_sd = 0.04,
                       alpha_bounds = NA,
                       beta_mean = .15,
                       beta_sd =.09,
                       beta_bounds = NA,
                       s_untested_mean = .025,
                       s_untested_sd = .0225,
                       s_untested_bounds = NA,
                       p_s0_pos_mean = .4,
                       p_s0_pos_sd = .1225,
                       p_s0_pos_bounds = NA,
                       pre_nsamp = 1e4,
                       post_nsamp = 1e3,
                       include_corrected = TRUE,
                       bde = FALSE) {
  
  
  given_args <- as.list(environment())
  cat("Arguments to get_melded:\n")
  print(given_args)
  
  
  theta <- tibble(alpha = sample_gamma_density(pre_nsamp,
                                               mean = alpha_mean,
                                               sd = alpha_sd,
                                               bounds = alpha_bounds),
                  beta= sample_beta_density(pre_nsamp,
                                            mean = beta_mean,
                                            sd = beta_sd,
                                            bounds = beta_bounds),
                  P_S_untested = sample_beta_density(pre_nsamp,
                                                     mean = s_untested_mean,
                                                     sd = s_untested_sd,
                                                     bounds = s_untested_bounds)) %>%
    mutate(phi_induced = est_P_A_testpos(P_S_untested = P_S_untested,
                                         alpha = alpha,
                                         beta=beta))
  
  # theta contains values sampled from alpha, beta, P_S_untested, and M(theta) = phi_induced
  # induced phi
  phi <- theta$phi_induced
  
  # approximate induced distribution via a density approximation
  if(bde) {
    message("Beta Kernel Density Estimation")
    phi_induced_density <- kdensity::kdensity(x = phi,  bw="SJ", adjust = 1, start = "beta", kernel = "beta", support = c(0,1))
    phi_sampled_density <- phi_induced_density(phi) 
  }
  else {
    phi_induced_density <- density(x = phi, n = pre_nsamp, adjust = 2, kernel = "gaussian")
    indexes <- findInterval(phi, phi_induced_density$x)
    phi_sampled_density <- phi_induced_density$y[indexes]
  }
  
  dp_s0_pos <- function(x) {
    
    beta_density(x,
                 mean=p_s0_pos_mean,
                 sd = p_s0_pos_sd,
                 bounds=p_s0_pos_bounds)
  }
  
  
  
  weights <- (phi_sampled_density/ dp_s0_pos(phi))^(.5)
  
  
  post_samp_ind <-sample.int(n=pre_nsamp,
                             size=post_nsamp,
                             prob=1/weights,
                             replace=TRUE)
  
  
  pi_samp <- cbind(theta[post_samp_ind,],
                   P_A_testpos =  phi[post_samp_ind]) %>%
    select(-phi_induced)
  
  
  return(list(post_melding = pi_samp, pre_melding = theta))
  
  
}



get_corrected_counts <- function(county_df, 
                                 melded_df, 
                                 post_nsamp,
                                 num_reps_counts = 1e3) {
  
  # melded_df <- melded_df %>%
  #   rename(Z_S = alpha,
  #          Z_A = beta)
  
  corrected <- pmap_df(county_df, ~ {
    process_priors_per_county(
      priors = melded_df,
      county_df = list(...),
      nsamp = post_nsamp) %>%
      generate_corrected_sample(., num_reps = num_reps_counts) %>%
      summarize_corrected_sample() })
  
  corrected %>%
    left_join(dates)
}