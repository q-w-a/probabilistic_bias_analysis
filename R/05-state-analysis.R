


#---------- Version 1 -------------
get_v1_state <- function(data, params, testing = TRUE) {
  
  state_testing <- data
  
  state_testing <- state_testing %>%
    rename(fips = state) %>%
    select(-date) %>% 
    distinct()
  
  
  # only use a few rows if testing
  state_testing <- if(testing) state_testing %>% 
    filter(biweek >=6) %>% slice_sample(n=15) else state_testing
  
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
      smoothed_s_untested <- loess(imputed_s_untested~index, 
                                   data = .x,
                                   span = smooth_s_untested_span)
      
      # only add smoothed beta if there are enough beta observations for the state (keep = TRUE)
      if (unique(.x$keep)) {
        smoothed_beta <- loess(imputed_beta~index, 
                               data = .x, 
                               span = smooth_beta_span)
        .x %>%
          mutate(beta_estimate_smoothed = predict(smoothed_beta),
                 s_untested_smoothed = predict(smoothed_s_untested)) 
      }
      
      else {
        .x %>%
          mutate(s_untested_smoothed = predict(smoothed_s_untested)) 
  
      }
     
    } ) 
  
  return(symp_smooth)
  
}



#---------- Version 2 -------------

get_v2_state <- function(data, params, ctis, testing = TRUE) {
  
  state_testing <- data
  
  state_testing <- state_testing %>%
    rename(fips = state) %>%
    select(-date) %>% 
    distinct()
  
  
  # only use a few rows if testing
  state_testing <- if(testing) state_testing %>% 
    filter(biweek >=6) %>% slice_sample(n=15) else state_testing
  
  # only keep states that have enough observations of beta
  ctis <- ctis %>% 
    filter(keep == TRUE) 
  
  ctis_biweekly <- ctis %>%
    left_join(dates) %>%
    group_by(biweek, state) %>%
    slice_max(n=1, order_by=date) %>%
    select(-date) %>%
    mutate(state =toupper(state)) %>%
    select(fips = state,
           biweek, 
           smoothed_s_untested, 
           smoothed_beta) 
    
    
    
    
  state_testing <- state_testing %>% 
  #  mutate(fips = tolower(fips)) %>%
    inner_join(ctis_biweekly, by =c('fips'='fips', 'biweek'='biweek')) %>%
    # only have CTIS data starting at week 6
    # filter out the beginning dates where beta_estimate_smoothed is NA
    filter(!is.na(smoothed_beta))
    
  
  corrected <- state_testing %>% 
    arrange(biweek) %>%
    # there will be more than one observation per county since
    # beta estimates are at the state level
    distinct() %>%
    pmap_df(~{
      
      state_data <- list(...) %>%  as.data.frame()
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
