
#---------- Version 1 -------------
get_county_v1 <- function(data, params, testing = FALSE) {

    county_testing <- data
    
    county_testing <- county_testing %>%
      select(-date) %>% 
      distinct()
    
    
    # only use a few rows if testing
    county_testing <- if(testing) county_testing %>% 
      filter(biweek >=6) %>% slice_sample(n=5) else county_testing
    
    melded <- do.call(get_melded, params)
    
    corrected <- pmap_df(county_testing, 
                         function(posrate,
                                  population,
                                  total, 
                                  positive,
                                  biweek,
                                  fips,
                                  ...) {
                           process_priors_per_county(
                             priors = melded$post_melding,
                             county_df = tibble(posrate, 
                                                population,
                                                total, 
                                                positive,
                                                biweek,
                                                fips),
                             nsamp = params$post_nsamp) %>%
                             generate_corrected_sample(., num_reps = 1e3) %>%
                             summarize_corrected_sample() })
    
    
    corrected <- corrected %>%
      mutate(version="v1")
    
    return(corrected)
    
    
}


#---------------- Version 2 --------------

get_ctis_biweekly <- function(ctis) {
  
  dates <- readRDS(here("data/data_raw/date_to_biweek.RDS"))
  
  ctis_biweekly <- ctis %>%
    left_join(dates) %>%
    group_by(biweek, state) %>%
    slice_max(n=1, order_by=date) %>%
    select(-date) %>%
    mutate(state =toupper(state)) %>%
    select(state,
           biweek, 
           s_untested_smoothed, 
           beta_estimate_smoothed,
           beta_estimate_spline_smoothed,
           keep) 
  
  return(ctis_biweekly)
  
}



get_corrected_county <- function(data, params, ctis,
                                 vary, spline = FALSE,
                                 testing = FALSE) {
  
  dates <- readRDS(here("data/data_raw/date_to_biweek.RDS"))
  
  county_testing <- data
  
  county_testing <- county_testing %>%
    select(-date) %>% 
    distinct()
  
  
  # only use a few rows if testing
  county_testing <- if(testing) county_testing %>% 
    filter(biweek >=6) %>% slice_sample(n=5) else county_testing
  
  ctis_biweekly <- get_ctis_biweekly(ctis)
  
  county_testing <- county_testing %>% 
    #  mutate(fips = tolower(fips)) %>%
    inner_join(ctis_biweekly, by =c('state'='state', 
                                    'biweek'='biweek')) 
  
  # glimpse(state_testing)
  
  if(vary=="beta" & spline == FALSE) {
    message("v2")
    corrected <- get_v2_county_corrected(county_testing, params) %>%
      mutate(version="v2") }
  
  if(vary=="s_untested" & spline == FALSE) {
    message("v3")
    corrected <- get_v3_county_corrected(county_testing, params) %>%
      mutate(version="v3") }
  
  if(vary=="s_untested_and_beta" & spline == FALSE) {
    message("v4")
    corrected <- get_v4_county_corrected(county_testing, params) %>%
      mutate(version="v4") }
  
  ###########################
  # spline-smoothed beta
  ###########################
  
  if(vary=="beta" & spline == TRUE) {
    message("v5")
    corrected <- get_v5_county_corrected(county_testing, params) %>%
      mutate(version="v5") }
  
  if(vary=="s_untested_and_beta" & spline == TRUE) {
    message("v6")
    corrected <- get_v6_county_corrected(county_testing, params) %>%
      mutate(version="v6") }
  

  return(corrected)
  
}



#---------- Version 2 -------------
# center prior for beta at survey estimates 
get_v2_county_corrected <- function(county_testing, params) {
  
  message("Running version 2")
  
  county_testing <- county_testing %>% 
    # only have CTIS data starting at week 6
    # filter out the beginning dates where beta_estimate_smoothed is NA
    filter(!is.na(beta_estimate_smoothed))
  
  
  corrected <- county_testing %>% 
    arrange(biweek) %>%
    # there will be more than one observation per county since
    # beta estimates are at the state level
    distinct() %>%
    pmap_df( function(fips, 
                      positive, 
                      total, 
                      biweek, 
                      posrate,
                      population,
                      beta_estimate_smoothed,
                      s_untested_smoothed,
                      ...) {
      
      state_data <- tibble(fips, positive, 
                           total,  biweek,
                           posrate, population,
                           beta_estimate_smoothed,
                           s_untested_smoothed)
      # message(paste0("before: ",prior_params$beta_mean))
      params$beta_mean <- state_data$beta_estimate_smoothed
      # message(paste0("after: ",prior_params$beta_mean))
      res <- do.call(get_melded, params)
      constrained <- res$post_melding
      
      # glimpse(constrained)
      
      process_priors_per_county(
        priors = constrained, 
        county_df = state_data,
        nsamp = params$post_nsamp) %>%
        generate_corrected_sample(., num_reps = 1e3) %>%
        summarize_corrected_sample()
    })
  
  return(corrected)
  
}


#---------- Version 3 -------------
# center prior for P(S_1|untested) at survey estimates 
get_v3_county_corrected <- function(county_testing, params) {
  
  message("Running version 3")
  
  county_testing <- county_testing %>% 
    # only have CTIS data starting at week 6
    # filter out the beginning dates where beta_estimate_smoothed is NA
    filter(!is.na(beta_estimate_smoothed))
  
  
  corrected <- county_testing %>% 
    arrange(biweek) %>%
    # there will be more than one observation per county since
    # beta estimates are at the state level
    distinct() %>%
    pmap_df( function(fips, 
                      positive, 
                      total, 
                      biweek, 
                      posrate,
                      population,
                      beta_estimate_smoothed,
                      s_untested_smoothed,
                      ...) {
      
      state_data <- tibble(fips, positive, 
                           total,  biweek,
                           posrate, population,
                           beta_estimate_smoothed,
                           s_untested_smoothed)
      # message(paste0("before: ",prior_params$beta_mean))
      params$s_untested_mean <- state_data$s_untested_smoothed
      # message(paste0("after: ",prior_params$beta_mean))
      res <- do.call(get_melded, params)
      constrained <- res$post_melding
      
      # glimpse(constrained)
      
      process_priors_per_county(
        priors = constrained, 
        county_df = state_data,
        nsamp = params$post_nsamp) %>%
        generate_corrected_sample(., num_reps = 1e3) %>%
        summarize_corrected_sample()
    })
  
  return(corrected)
  
}






#---------- Version 4 -------------
# center prior for BOTH beta and P(S_1|untested) at survey estimates 
get_v4_county_corrected <- function(county_testing, params) {
  
  message("Running version 4")
  
  county_testing <- county_testing %>% 
    # only have CTIS data starting at week 6
    # filter out the beginning dates where beta_estimate_smoothed is NA
    filter(!is.na(beta_estimate_smoothed))
  
  
  corrected <- county_testing %>% 
    arrange(biweek) %>%
    # there will be more than one observation per county since
    # beta estimates are at the state level
    distinct() %>%
    pmap_df( function(fips, 
                      positive, 
                      total, 
                      biweek, 
                      posrate,
                      population,
                      beta_estimate_smoothed,
                      s_untested_smoothed,
                      ...) {
      
      state_data <- tibble(fips, positive, 
                           total,  biweek,
                           posrate, population,
                           beta_estimate_smoothed,
                           s_untested_smoothed)
      # message(paste0("before: ",prior_params$beta_mean))
      params$s_untested_mean <- state_data$s_untested_smoothed
      params$beta_mean <- state_data$beta_estimate_smoothed
      
      # message(paste0("after: ",prior_params$beta_mean))
      res <- do.call(get_melded, params)
      constrained <- res$post_melding
      
      # glimpse(constrained)
      
      process_priors_per_county(
        priors = constrained, 
        county_df = state_data,
        nsamp = params$post_nsamp) %>%
        generate_corrected_sample(., num_reps = 1e3) %>%
        summarize_corrected_sample()
    })
  
  return(corrected)
  
}



#---------- Version 5 -------------
# center prior for BOTH beta and P(S_1|untested) at survey estimates 
get_v5_county_corrected <- function(county_testing, params) {
  
  message("Running version 5")
  
  county_testing <- county_testing %>% 
    # only have CTIS data starting at week 6
    # filter out the beginning dates where beta_estimate_smoothed is NA
    filter(!is.na(beta_estimate_smoothed))
  
  
  corrected <- county_testing %>% 
    arrange(biweek) %>%
    # there will be more than one observation per county since
    # beta estimates are at the state level
    distinct() %>%
    pmap_df( function(fips, 
                      positive, 
                      total, 
                      biweek, 
                      posrate,
                      population,
                      beta_estimate_smoothed,
                      beta_estimate_spline_smoothed,
                      s_untested_smoothed,
                      ...) {
      
      state_data <- tibble(fips, positive, 
                           total,  biweek,
                           posrate, population,
                           beta_estimate_smoothed,
                           beta_estimate_spline_smoothed,
                           s_untested_smoothed)
      # message(paste0("before: ",prior_params$beta_mean))
      # params$s_untested_mean <- state_data$s_untested_smoothed
      params$beta_mean <- state_data$beta_estimate_spline_smoothed
      
      message(paste0("v5 beta mean: ", params$beta_mean ))
      
      # message(paste0("after: ",prior_params$beta_mean))
      res <- do.call(get_melded, params)
      constrained <- res$post_melding
      
      # glimpse(constrained)
      
      process_priors_per_county(
        priors = constrained, 
        county_df = state_data,
        nsamp = params$post_nsamp) %>%
        generate_corrected_sample(., num_reps = 1e3) %>%
        summarize_corrected_sample()
    })
  
  return(corrected)
  
}




#---------- Version 6 -------------
# center prior for BOTH beta (spline smoothed) and loess-smoothed P(S_1|untested) at survey estimates 
get_v6_county_corrected <- function(county_testing, params) {
  
  message("Running version 6")
  
  county_testing <- county_testing %>% 
    # only have CTIS data starting at week 6
    # filter out the beginning dates where beta_estimate_smoothed is NA
    filter(!is.na(beta_estimate_smoothed))
  
  
  corrected <- county_testing %>% 
    arrange(biweek) %>%
    # there will be more than one observation per county since
    # beta estimates are at the state level
    distinct() %>%
    pmap_df( function(fips, 
                      positive, 
                      total, 
                      biweek, 
                      posrate,
                      population,
                      beta_estimate_smoothed,
                      beta_estimate_spline_smoothed,
                      s_untested_smoothed,
                      ...) {
      
      state_data <- tibble(fips, positive, 
                           total,  biweek,
                           posrate, population,
                           beta_estimate_smoothed,
                           beta_estimate_spline_smoothed,
                           s_untested_smoothed)
      # message(paste0("before: ",prior_params$beta_mean))
      params$s_untested_mean <- state_data$s_untested_smoothed
      params$beta_mean <- state_data$beta_estimate_spline_smoothed
      
      message(paste0("v6 beta mean: ", params$beta_mean ))
      message(paste0("v6 s_untested mean: ", params$s_untested_mean ))
      
      
      # message(paste0("after: ",prior_params$beta_mean))
      res <- do.call(get_melded, params)
      constrained <- res$post_melding
      
      # glimpse(constrained)
      
      process_priors_per_county(
        priors = constrained, 
        county_df = state_data,
        nsamp = params$post_nsamp) %>%
        generate_corrected_sample(., num_reps = 1e3) %>%
        summarize_corrected_sample()
    })
  
  return(corrected)
  
}


#----------------- Version 7 ------------------------
#---------- Version 1 -------------
get_county_v7 <- function(data, testing = FALSE) {
  
  params <-  list(
    alpha_mean = .95,
    alpha_sd = 0.08,
    alpha_bounds = NA,
    # alpha_bounds = c(.8,1),
    beta_mean = .18,
    beta_sd =.09,
    beta_bounds = NA,
    #  beta_bounds = c(0.002, 0.4),
    s_untested_mean = .016,
    s_untested_sd = .0225,
    #  s_untested_bounds = c(0.0018, Inf),
    s_untested_bounds = NA,
    p_s0_pos_mean = .4,
    p_s0_pos_sd = .1225,
    p_s0_pos_bounds = NA,
    #  p_s0_pos_bounds = c(.25, .7),
    pre_nsamp = 1e6,
    post_nsamp = 1e5)
  
  
  county_testing <- data
  
  county_testing <- county_testing %>%
    select(-date) %>% 
    distinct()
  
  
  # only use a few rows if testing
  county_testing <- if(testing) county_testing %>% 
    filter(biweek >=6) %>% slice_sample(n=5) else county_testing
  
  melded <- do.call(get_melded, params)
  
  corrected <- pmap_df(county_testing, 
                       function(posrate,
                                population,
                                total, 
                                positive,
                                biweek,
                                fips,
                                ...) {
                         process_priors_per_county(
                           priors = melded$post_melding,
                           county_df = tibble(posrate, 
                                              population,
                                              total, 
                                              positive,
                                              biweek,
                                              fips),
                           nsamp = params$post_nsamp) %>%
                           generate_corrected_sample(., num_reps = 1e3) %>%
                           summarize_corrected_sample() })
  
  
  corrected <- corrected %>%
    mutate(version="v7")
  
  return(corrected)
  
  
}



