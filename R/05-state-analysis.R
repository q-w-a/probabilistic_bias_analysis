


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

objective <- function(shapes, x, prob, ...) {
  fit <- pbeta(x, shapes[1], shapes[2])
  # return sum of squares
  return(sum((fit-prob)^2))
}



get_shape_smoothed <- function(ctis_smoothed, option, quiet=FALSE) {
  
  if (option=="beta") {
    ctis_quantiles <- ctis_smoothed %>%
      summarize(q1 = quantile(beta_estimate_spline_smoothed,.025, na.rm=TRUE),
                median = median(beta_estimate_spline_smoothed, na.rm=TRUE),
                q2 = quantile(beta_estimate_spline_smoothed,.975, na.rm=TRUE),
                sd=sd(beta_estimate_spline_smoothed,na.rm=TRUE))
    

  }
  if (option=="s_untested") {
    ctis_quantiles <- ctis_smoothed %>%
      summarize(q1 = quantile(s_untested_smoothed,.025, na.rm=TRUE),
                median = median(s_untested_smoothed, na.rm=TRUE),
                q2 = quantile(s_untested_smoothed,.975, na.rm=TRUE),
                sd=sd(s_untested_smoothed, na.rm=TRUE))
    
  }
  
  start <- get_beta_params(ctis_quantiles$q1, ctis_quantiles$sd) %>%
    unlist()
  
  optim_results <- nlm(objective, start, 
             x=c(ctis_quantiles$q1, ctis_quantiles$median, ctis_quantiles$q2),
             prob=c(.025,.5,.975), lower=0, upper=1,
             typsize=c(1,1), fscale=1e-14, gradtol=1e-14)
  
  params <- (optim_results$estimate) 
  
  if(!quiet) {
    message(paste0("Input quantiles (0.025, .5, .975):\n",
                   round(ctis_quantiles$q1,4), ", ",
                   round(ctis_quantiles$median,4), ", ",
                   round(ctis_quantiles$q2,4)))
    sim <- rbeta(1e4, params[1], params[2])
    output_quantiles <- list(q1=quantile(sim,.025),
                          median=median(sim),
                          q2=quantile(sim, .975))
    message(paste0("Output quantiles (0.025, .5, .975):\n",
                   round(output_quantiles$q1,4), ", ",
                   round(output_quantiles$median,4), ", ",
                   round(output_quantiles$q2,4)))
    
    
  }
  
  return(params)
  
  
}





get_shape <- function(ctis_smoothed, option, quiet=FALSE) {
  
  dates <- readRDS(here("data/data_raw/date_to_biweek.RDS"))
  
  if (option=="beta") {
    ctis_quantiles <- ctis_smoothed %>% 
      left_join(dates) %>%
      group_by(biweek,state) %>%
      summarize(beta_est=mean(beta_est,na.rm=TRUE)) %>%
      ungroup() %>%
      summarize(q1 = quantile(beta_est,.025, na.rm=TRUE),
                median = median(beta_est, na.rm=TRUE),
                q2 = quantile(beta_est,.975, na.rm=TRUE),
                sd=sd(beta_est,na.rm=TRUE))
    
  }
  if (option=="s_untested") {
    ctis_quantiles <- ctis_smoothed %>% 
      left_join(dates) %>%
      group_by(biweek,state) %>%
      summarize(smoothed_wcli=mean(smoothed_wcli,na.rm=TRUE))  %>%
      ungroup() %>%
      summarize(q1 = quantile(smoothed_wcli,.025, na.rm=TRUE),
                median = median(smoothed_wcli, na.rm=TRUE),
                q2 = quantile(smoothed_wcli,.975, na.rm=TRUE),
                sd=sd(smoothed_wcli, na.rm=TRUE))
    
  }
  
  start <- get_beta_params(ctis_quantiles$median, ctis_quantiles$sd) %>%
    unlist()
  
  optim_results <- nlm(objective, start, 
                       x=c(ctis_quantiles$q1, ctis_quantiles$median, ctis_quantiles$q2),
                       prob=c(.025,.5,.975), lower=0, upper=1,
                       typsize=c(1,1), fscale=1e-14, gradtol=1e-14)
  
  params <- (optim_results$estimate) 
  
  if(!quiet) {
    message(paste0("Input quantiles (0.025, .5, .975):\n",
                   round(ctis_quantiles$q1,4), ", ",
                   round(ctis_quantiles$median,4), ", ",
                   round(ctis_quantiles$q2,4)))
    sim <- rbeta(1e4, params[1], params[2])
    output_quantiles <- list(q1=quantile(sim,.025),
                             median=median(sim),
                             q2=quantile(sim, .975))
    message(paste0("Output quantiles (0.025, .5, .975):\n",
                   round(output_quantiles$q1,4), ", ",
                   round(output_quantiles$median,4), ", ",
                   round(output_quantiles$q2,4)))
    
    
  }
  
  return(params)
  
  
}


#---------------- Version 7 --------------

get_v7_state <- function(data, beta_shape, s_untested_shape, testing = FALSE) {
  
  params <-  list(
    alpha_mean = .95,
    alpha_sd = 0.08,
    alpha_bounds = NA,
    # alpha_bounds = c(.8,1),
    beta_shape1 = beta_shape[1],
    beta_shape2 = beta_shape[2],
    beta_bounds = NA,
    #  beta_bounds = c(0.002, 0.4),
    s_untested_shape1 = s_untested_shape[1],
    s_untested_shape2 = s_untested_shape[2],
    #  s_untested_bounds = c(0.0018, Inf),
    s_untested_bounds = NA,
    p_s0_pos_mean = .4,
    p_s0_pos_sd = .1225,
    p_s0_pos_bounds = NA,
    #  p_s0_pos_bounds = c(.25, .7),
    pre_nsamp = 1e6,
    post_nsamp = 1e5,
    direct_params=TRUE)
  
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


