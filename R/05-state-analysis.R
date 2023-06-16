


#---------- Version 1 -------------
get_v1_state <- function(data, params, testing = FALSE) {
  
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


get_smoothed_ctis <- function(ctis,
                              smooth_beta_span = .33,
                              smooth_s_untested_span =.2) {
  
  dates <- readRDS(here("data/data_raw/date_to_biweek.RDS"))
  
  symp <- ctis %>% 
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
  symp <- symp %>%
    group_by(state) %>%
    arrange(date) %>%
    mutate(imputed_beta = ifelse(
      keep,
      imputeTS::na_ma(beta_est, k = 10, weighting = "simple"),
      beta_est),
      imputed_s_untested  =ifelse(
        keep,
        imputeTS::na_ma(smoothed_wcli, k = 10, weighting = "simple"),
        smoothed_wcli ))
  
  
  # loess smoothing
  
  
  symp_smooth <- symp %>%
    group_by(state) %>%
    group_split() %>%
    map_df(~{
      if(unique(.x$keep)) {
        .x <- .x %>% arrange(date) %>% 
          mutate(index=row_number())
        l <- loess(imputed_beta ~ index,
                   data=.x, 
                   span=smooth_beta_span)
        s <- loess(imputed_s_untested ~ index, data=.x, span=.2)
        
        .x %>% mutate(smoothed_beta =predict(l),
                      smoothed_s_untested = predict(s))
      }
      
      else { 
        .x <-.x %>% arrange(date) %>% mutate(index=row_number())
        s <- loess(imputed_s_untested ~ index, data=.x, span=smooth_s_untested_span)
        
        .x %>% mutate(
          smoothed_s_untested = predict(s))
      }
    })
  
  return(symp_smooth)
  
}



#---------- Version 2 -------------

