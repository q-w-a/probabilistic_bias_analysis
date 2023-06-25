


#---------- Version 1 -------------
get_v1_state <- function(data, params, testing = FALSE) {
  
  state_testing <- data
  
  state_testing <- state_testing %>%
    rename(fips = state) %>%
    select(-date) %>% 
    distinct()
  
  
  # only use a few rows if testing
  state_testing <- if(testing) state_testing %>% 
    filter(biweek >=6) %>% slice_sample(n=5) else state_testing
  
  melded <- do.call(get_melded, params)
  
  corrected <- pmap_df(state_testing, 
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


#------------- Using CTIS Data ----------------
get_smoothed_ctis <- function(ctis_raw,
                              smooth_beta_span = .33,
                              smooth_s_untested_span =.2) {
  
  dates <- readRDS(here("data/data_raw/date_to_biweek.RDS"))
  
  symp <- ctis_raw %>% 
    select(signal, date, value, state = geo_value) %>% 
    pivot_wider(names_from = signal,
                values_from = c(value)) %>%
    mutate(beta_est = smoothed_wscreening_tested_positive_14d/smoothed_wtested_positive_14d)
  
  symp <- symp %>%
    group_by(state) %>%
    mutate(m = sum(!is.na(beta_est)),
           n = n()) %>%
    mutate(prop_with_data = m/n) %>%
    # keep states with less than 40% of observations missing
    mutate(keep = ifelse(prop_with_data > .6, TRUE, FALSE)) %>%
    ungroup() 
  
  # impute missing values
  
  
  # make sure each state has all dates
  symp_all <- symp %>% 
    select(date,state,beta_est) %>% 
    pivot_wider(names_from = date, 
                values_from =beta_est) %>% 
    pivot_longer(cols =2:ncol(.), 
                 values_to = "beta_est",
                 names_to = "date") %>%
    mutate(date=as_date(date))
  
  symp <- symp_all %>%
    left_join(symp)
  
  symp_imputed <- symp %>%
    group_by(state) %>%
    arrange(date) %>%
    mutate(imputed_beta = ifelse(keep,
                                 imputeTS::na_ma(beta_est,
                                          k = 20, 
                                          weighting = "simple"),
                                 beta_est),
           imputed_s_untested  = imputeTS::na_ma(smoothed_wcli,
                                                 k = 20, 
                                                 weighting = "simple")) 
  
  
  symp_smooth <- symp_imputed %>%
    # select(-beta_est) %>%
    # rename(beta_est = imputed_beta) %>%
    group_by(state) %>%
    arrange(date) %>%
    mutate(index = row_number()) %>%
    group_split() %>%
    map_df(~ {
      # message(paste0(unique(.x$state), unique(.x$keep)))
      smoothed_s_untested <- loess(imputed_s_untested~index, 
                                   data = .x,
                                   span = smooth_s_untested_span)
      
      # only add smoothed beta if there are enough beta observations for the state (keep = TRUE)
      if (unique(.x$keep)) {
        smoothed_beta <- loess(imputed_beta~index, 
                               data = .x, 
                               span = smooth_beta_span)
        # 257 corresponds to 2021-12-01
        spline_smoothed_beta <- lm(imputed_beta ~ splines::ns(index,knots = c(118, 257)), 
                                   data = .x )
        
        # include both spline smoothing and loess smoothing for beta
        .x %>%
          mutate(beta_estimate_smoothed = predict(smoothed_beta, newdata = index),
                 beta_estimate_spline_smoothed = predict(spline_smoothed_beta, newx= index),
                 s_untested_smoothed = predict(smoothed_s_untested, newdata = index)) 
      }
      
      else {
        .x %>%
          mutate(s_untested_smoothed = predict(smoothed_s_untested,newdata = index)) 
  
      }
     
    } ) 
  
  symp_smooth %>% 
    select(-c(m,n,prop_with_data, index)) %>%
    return()
  
}



get_corrected_state <- function(data, params,
                                ctis, vary, 
                                spline=FALSE,
                                testing = FALSE) {
  
  dates <- readRDS(here("data/data_raw/date_to_biweek.RDS"))
  
  state_testing <- data
  
  state_testing <- state_testing %>%
    rename(fips = state) %>%
    select(-date) %>% 
    distinct()
  
  
  # only use a few rows if testing
  state_testing <- if(testing) state_testing %>% 
    filter(biweek >=6) %>% slice_sample(n=5) else state_testing
  
  
  ctis_biweekly <- ctis %>%
    left_join(dates) %>%
    group_by(biweek, state) %>%
    slice_max(n=1, order_by=date) %>%
    select(-date) %>%
    mutate(state =toupper(state)) %>%
    select(fips = state,
           biweek, 
           s_untested_smoothed, 
           beta_estimate_smoothed,
           beta_estimate_spline_smoothed,
           keep)  %>%
    ungroup()
  
  glimpse(ctis_biweekly)
    
  state_testing <- state_testing %>% 
  #  mutate(fips = tolower(fips)) %>%
    inner_join(ctis_biweekly, by =c('fips'='fips', 'biweek'='biweek')) 
  
  # glimpse(state_testing)
  
  
  ##########################
  # loess smoothing
  ##########################
  if(vary=="beta" & spline == FALSE) {
    message("v2")
    corrected <- get_v2_corrected(state_testing, params) %>%
    mutate(version="v2") }
  
  if(vary=="s_untested" & spline == FALSE) {
    message("v3")
    corrected <- get_v3_corrected(state_testing, params) %>%
      mutate(version="v3") }
  
  if(vary=="s_untested_and_beta" & spline == FALSE) {
    message("v4")
    corrected <- get_v4_corrected(state_testing, params) %>%
      mutate(version="v4") }
  
  ##########################
  # spline smoothing
  ##########################
  if(vary=="beta" & spline == TRUE) {
    message("v5")
    corrected <- get_v5_corrected(state_testing, params) %>%
      mutate(version="v5") }
  
  if(vary=="s_untested_and_beta" & spline == TRUE) {
    message("v6")
    corrected <- get_v6_corrected(state_testing, params) %>%
      mutate(version="v6") }
  
  

  
  return(corrected)
  
}



#---------- Version 2 -------------
# center prior for beta at survey estimates 
get_v2_corrected <- function(state_testing, params) {
  
  message("Running version 2")
  
  state_testing <- state_testing %>% 
    # only have CTIS data starting at week 6
    # filter out the beginning dates where beta_estimate_smoothed is NA
    filter(!is.na(beta_estimate_smoothed))
  
  
  corrected <- state_testing %>% 
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



#---------- Version 3 ---------------------------------
# center prior for P(S_1|untested) at survey estimates 
get_v3_corrected <- function(state_testing, params) {
  
  message("Running version 3")

  
  state_testing <- state_testing %>% 
    # only have CTIS data starting at week 6
    # filter out the beginning dates where beta_estimate_smoothed is NA
    filter(!is.na(beta_estimate_smoothed))
  
  corrected <- state_testing %>% 
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




#---------- Version 4 ---------------------------------
# center prior for BOTH P(S_1|untested) and beta at survey estimates 
get_v4_corrected <- function(state_testing, params) {
  
  message("Running version 4")
  
  state_testing <- state_testing %>% 
    # only have CTIS data starting at week 6
    # filter out the beginning dates where beta_estimate_smoothed is NA
    filter(!is.na(beta_estimate_smoothed))
  
  
  corrected <- state_testing %>% 
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





#---------- Version 5 ---------------------------------
# center prior for beta at survey estimate (spline smoothed) 
# at survey estimate (loess smoothed) 
get_v5_corrected <- function(state_testing, params) {
  
  glimpse(state_testing)
  nrow(state_testing)
  
  message("Running version 5")

  state_testing <- state_testing %>% 
    # only have CTIS data starting at week 6
    # filter out the beginning dates where beta_estimate_smoothed is NA
    filter(!is.na(beta_estimate_smoothed))
  
  
  corrected <- state_testing %>% 
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
                           s_untested_smoothed,
                           beta_estimate_spline_smoothed)
      # message(paste0("before: ",prior_params$beta_mean))
      params$beta_mean <- state_data$beta_estimate_spline_smoothed
      
      message(paste0("v5 beta mean: ", params$beta_mean ))
      
      
      # params$s_untested_mean <- state_data$s_untested_smoothed
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



#---------- Version 6 ---------------------------------
# center prior for beta at survey estimate (spline smoothed) and P(S_1|untested) 
# at survey estimate (loess smoothed) 
get_v6_corrected <- function(state_testing, params) {
  
  message("Running version 6")
  
  state_testing <- state_testing %>% 
    # only have CTIS data starting at week 6
    # filter out the beginning dates where beta_estimate_smoothed is NA
    filter(!is.na(beta_estimate_smoothed))
  
  
  corrected <- state_testing %>% 
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
                      beta_estimate_spline_smoothed,
                      ...) {
      
      state_data <- tibble(fips, positive, 
                           total,  biweek,
                           posrate, population,
                           beta_estimate_smoothed,
                           beta_estimate_spline_smoothed,
                           s_untested_smoothed)
      # message(paste0("before: ",prior_params$beta_mean))
      params$beta_mean <- state_data$beta_estimate_spline_smoothed
      params$s_untested_mean <- state_data$s_untested_smoothed
      # message(paste0("after: ",prior_params$beta_mean))
      
      
      message(paste0("v6 beta mean: ", params$beta_mean ))
      message(paste0("v6 s_untested mean: ", params$s_untested_mean ))
      
      
      
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



#---------------- Version 7 --------------

get_v7_state <- function(data, testing = FALSE) {
  
  params <-  list(
    alpha_mean = .95,
    alpha_sd = 0.08,
    alpha_bounds = NA,
    # alpha_bounds = c(.8,1),
    beta_mean = .216,
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
  
  state_testing <- data
  
  state_testing <- state_testing %>%
    rename(fips = state) %>%
    select(-date) %>% 
    distinct()
  
  
  # only use a few rows if testing
  state_testing <- if(testing) state_testing %>% 
    filter(biweek >=6) %>% slice_sample(n=5) else state_testing
  
  melded <- do.call(get_melded, params)
  
  corrected <- pmap_df(state_testing, 
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


