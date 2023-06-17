
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
           keep) 
  
  return(ctis_biweekly)
  
}



get_corrected_county <- function(data, params, ctis,
                                 vary, testing = FALSE) {
  
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
  
  if(vary=="beta") {
    message("v2")
    corrected <- get_v2_corrected(county_testing, params) %>%
      mutate(version="v2") }
  
  if(vary=="s_untested") {
    message("v3")
    corrected <- get_v3_corrected(county_testing, params) %>%
      mutate(version="v3") }
  
  if(vary=="s_untested_and_beta") {
    message("v4")
    corrected <- get_v4_corrected(county_testing, params) %>%
      mutate(version="v4") }
  
  
  return(corrected)
  
}



#---------- Version 2 -------------
# center prior for beta at survey estimates 
get_v2_corrected <- function(county_testing, params) {
  
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
# center prior for beta at survey estimates 
get_v3_corrected <- function(county_testing, params) {
  
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
# center prior for beta at survey estimates 
get_v4_corrected <- function(county_testing, params) {
  
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


