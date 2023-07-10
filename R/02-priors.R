
get_priors <- function() { 
  
  # list(
  #   alpha_mean = .95,
  #   alpha_sd = 0.08,
  #   alpha_bounds = NA,
  #   # alpha_bounds = c(.8,1),
  #   beta_mean = .15,
  #   beta_sd =.09,
  #   beta_bounds = NA,
  #   #  beta_bounds = c(0.002, 0.4),
  #   s_untested_mean = .03,
  #   s_untested_sd = .0225,
  #   #  s_untested_bounds = c(0.0018, Inf),
  #   s_untested_bounds = NA,
  #   p_s0_pos_mean = .4,
  #   p_s0_pos_sd = .1225,
  #   p_s0_pos_bounds = NA,
  #   #  p_s0_pos_bounds = c(.25, .7),
  #   pre_nsamp = 1e6,
  #   post_nsamp = 1e5)
  
  list(
    alpha_mean = .95,
    alpha_sd = 0.08,
    alpha_bounds = NA,
    # alpha_bounds = c(.8,1),
    beta_mean = .15,
    beta_sd =.09,
    beta_bounds = NA,
    #  beta_bounds = c(0.002, 0.4),
    s_untested_mean = .03,
    s_untested_sd = .0225,
    #  s_untested_bounds = c(0.0018, Inf),
    s_untested_bounds = NA,
    p_s0_pos_mean = .55,
    p_s0_pos_sd = .12,
    p_s0_pos_bounds = NA,
    #  p_s0_pos_bounds = c(.25, .7),
    pre_nsamp = 1e6,
    post_nsamp = 1e5)
  
  
}

