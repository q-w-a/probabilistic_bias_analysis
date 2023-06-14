



###############################################################
# BETA PARAMETERS FROM DESIRED MEAN AND VARIANCE
###############################################################
get_beta_params <- function(mu, sd){
  var = sd ^ 2
  a = ((1 - mu) / var - 1 / mu) * mu ^ 2
  b = a * (1 / mu - 1)
  return (list(a = a, b = b))
}


# vectorized version to get each parameter separately
get_shape1 <- function(mu, sd) {
  var = sd^2
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(alpha)
}

get_shape2 <- function(mu, sd) {
  var = sd^2
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(beta)
}




###############################################################
# BETA DENSITY WITH DESIRED MEAN AND VARIANCE
###############################################################
beta_density <- function(x, mean, sd, bounds=NA) {
  shape_params <-  get_beta_params(
    mu = mean,
    sd = sd)
  
  if(!length(bounds) == 1){
    # message("here")
    dtrunc(x,
           spec = "beta",
           a = bounds[1],
           b = bounds[2],
           shape1 = shape_params$a,
           shape2 = shape_params$b) %>%
      return()
  }else{
    dbeta(x,
          shape1 = shape_params$a,
          shape2 = shape_params$b)  %>%
      return()
  }
}




###############################################################
# SAMPLE FROM BETA DENSITY WITH DESIRED MEAN AND VARIANCE
###############################################################

sample_beta_density <- function(n, mean, sd, bounds = NA) {
  
  shape_params <-  get_beta_params(
    mu = mean,
    sd = sd)
  
  rbeta(n,
        shape1 = shape_params$a,
        shape2 = shape_params$b)
  
  if(!length(bounds) == 1){
    # message("here")
    rtrunc(n,
           spec = "beta",
           a = bounds[1],
           b = bounds[2],
           shape1 = shape_params$a,
           shape2 = shape_params$b) %>%
      return()
  }else{
    rbeta(n,
          shape1 = shape_params$a,
          shape2 = shape_params$b)  %>%
      return()
  }
}




###############################################################
# GAMMA PARAMETERS FROM DESIRED MEAN AND VARIANCE
###############################################################
get_gamma_params <- function(mu, sd) {
  var = (mu/sd)^2
  shape = (mu/sd)^2
  scale = sd^2/mu
  return(params = list(shape = shape,
                       scale = scale))
}


###############################################################
# GAMMA DENSITY WITH DESIRED MEAN AND VARIANCE
###############################################################
gamma_density <- function(x, mean, sd, bounds=NA) {
  
  shape_params <-  get_gamma_params(
    mu = mean,
    sd = sd)
  
  if(!length(bounds) == 1){
    #message("here")
    dtrunc(x,
           spec = "gamma",
           a = bounds[1],
           b = bounds[2],
           shape = shape_params$shape,
           scale = shape_params$scale) %>%
      return()
  }else{
    dgamma(x,
           shape = shape_params$shape,
           scale = shape_params$scale) %>%
      return()
  }
}


sample_gamma_density <- function(n, mean, sd, bounds = NA) {
  
  shape_params <-  get_gamma_params(
    mu = mean,
    sd = sd)
  
  if(!length(bounds) == 1){
    #message("here")
    rtrunc(n,
           spec = "gamma",
           a = bounds[1],
           b = bounds[2],
           shape = shape_params$shape,
           scale = shape_params$scale) %>%
      return()
  }else{
    rgamma(n,
           shape = shape_params$shape,
           scale = shape_params$scale) %>%
      return()
  }
}





###############################################################
# INDUCED PRIOR ON ASYMPTOMATIC RATE  P(S_0|test+,untested)
###############################################################

# input sampled values of theta and compute M(\theta)
est_P_A_testpos = function(P_S_untested, alpha, beta){
  (beta * (1 - P_S_untested)) / (( beta * (1 - P_S_untested)) + (alpha * P_S_untested))
}



##########################################
# County Correction Functions -----------
##########################################


#' Add sensitivity and specificity
process_priors_per_county <-  function(priors, county_df, nsamp){
  dist_Se <- truncdist::rtrunc(n = nsamp,spec = "beta",a = 0.65,b = 1,
                               shape1 = get_beta_params(mu = 0.8,
                                                        sd = (0.4)^2)$a,
                               shape2 = get_beta_params(mu = 0.8,
                                                        sd = (0.4)^2)$b)
  dist_Sp <- truncdist::rtrunc(n = nsamp,spec = "beta",a = 0.9998,b = 1,
                               shape1 = get_beta_params(mu = 0.99995,
                                                        sd = (0.01)^2)$a,
                               shape2 = get_beta_params(mu = 0.99995,
                                                        sd = (0.01)^2)$b)
  
  priors_out <- priors %>%
    mutate(
      # calculate P(+|S_1, untested) and P(+|S_0, untested)
      P_testpos_S = priors$alpha  * county_df$posrate,
      P_testpos_A = priors$beta  * county_df$posrate,
      empirical_testpos = county_df$posrate,
      population = county_df$population,
      total = county_df$total,
      positive = county_df$positive,
      negative = county_df$negative) %>%
    mutate(Se = dist_Se,
           Sp = dist_Sp,
           biweek = county_df$biweek,
           fips = county_df$fips) %>%
    # compute with constrained priors
    mutate(P_A_testpos = est_P_A_testpos(
      P_S_untested = priors$P_S_untested,
      alpha = priors$alpha,
      beta = priors$beta))
  return(priors_out)
}


#' correct for test inaccuracy
calc_A_star <- function(N, N_tested,
                        N_pos_obs,
                        P_testpos_est,
                        P_S_untested,
                        P_A_testpos,
                        alpha,
                        beta,
                        Se,
                        Sp){
  
  N_untested = N - N_tested
  
  # number symptomatic and asymptomatic among tested
  Npos_tested_S = N_pos_obs * (1 - P_A_testpos)
  Npos_tested_A = N_pos_obs - Npos_tested_S
  
  #  prob testpos among untested
  P_testpos_S = P_testpos_est * alpha
  P_testpos_A = P_testpos_est * beta
  
  # estimate number of positives among untested
  Npos_untested_S = P_S_untested * N_untested * P_testpos_S
  Npos_untested_A = (1 - P_S_untested) * N_untested * P_testpos_A
  
  A_star = Npos_tested_S   + Npos_tested_A +
    Npos_untested_S + Npos_untested_A
  
  # correct for imperfect sensitivity and specificity
  A = (A_star - ((1 - Sp) * N)) / (Se + Sp - 1)
  
  return(max(A,0))
  
}

#' compute corrected counts for randomly sampled combination of parameters
correct_bias <- function(priors_by_county_df){
  # N, N_tested, N_pos_obs, P_testpos_est,
  
  # sample index to draw from distribution
  sample_ind = sample(1:nrow(priors_by_county_df),
                      size = 1,
                      replace=TRUE)
  
  # randomly sample from each distribution
  sampled_priors = priors_by_county_df[sample_ind,]
  
  # corrected case count
  Astar = calc_A_star(N = sampled_priors$population,
                      N_tested =sampled_priors$total,
                      N_pos_obs = sampled_priors$positive,
                      P_testpos_est = sampled_priors$empirical_testpos,
                      P_S_untested = sampled_priors$P_S_untested,
                      P_A_testpos = sampled_priors$P_A_testpos,
                      alpha = sampled_priors$alpha,
                      beta = sampled_priors$beta,
                      Se = sampled_priors$Se,
                      Sp = sampled_priors$Sp
  )
  
  out = bind_cols(sampled_priors, exp_cases=Astar)
  
  return(out)
}


#' generate distribution of corrected counts
generate_corrected_sample = function(priors_by_county_df,
                                     num_reps){
  
  # need to set seed here to ensure that the same random draws are
  # used for a given time period / location with same priors
  set.seed(123)
  
  # perform probabilistic bias analysis
  result = replicate(num_reps, 
                     correct_bias(priors_by_county_df ),
                     simplify=FALSE) %>%
    bind_rows()
  
  return(result)
  
}


#' summarize distribution of corrected counts
summarize_corrected_sample <- function(priors_by_county_df_exp_cases) {
  
  summarized <-  tibble(
    exp_cases_median = median(priors_by_county_df_exp_cases$exp_cases),
    exp_cases_lb = quantile(priors_by_county_df_exp_cases$exp_cases,
                            prob = 0.025) %>% unlist(),
    exp_cases_ub = quantile(priors_by_county_df_exp_cases$exp_cases,
                            prob = 0.975) %>% unlist(),
    biweek = unique(priors_by_county_df_exp_cases$biweek),
    fips = unique(priors_by_county_df_exp_cases$fips),
    empirical_testpos = unique(priors_by_county_df_exp_cases$empirical_testpos),
    population = unique(priors_by_county_df_exp_cases$population),
    negative = unique(priors_by_county_df_exp_cases$negative),
    positive = unique(priors_by_county_df_exp_cases$positive),
    total = unique(priors_by_county_df_exp_cases$total))
  
  return(summarized)
  
}

