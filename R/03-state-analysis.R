

get_v1_state <- function(data, params) {
  
  state_testing <- data
  
  
  state_testing <- state_testing %>%
    rename(fips = state) %>%
    select(-date) %>% 
    distinct()
  
  
  # only use a few rows if testing
  state_testing <- if(testing) state_testing %>% filter(biweek >=6) %>% slice_sample(n=15) else state_testing
  
  results_path <- here("analysis/results/adj_biweekly_state/")
  
  
  melded <- do.call(get_melded, prior_params)
  
  
  
}