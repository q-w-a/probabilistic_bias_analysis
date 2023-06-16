

get_v1_state <- function(data, params, testing = TRUE) {
  
  state_testing <- data
  
  state_testing <- state_testing %>%
    rename(fips = state) %>%
    select(-date) %>% 
    distinct()
  
  
  # only use a few rows if testing
  state_testing <- if(testing) state_testing %>% 
    filter(biweek >=6) %>% slice_sample(n=15) else state_testing
  
  results_path <- here("results/adj_biweekly_state/")
  
  
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
