
library(here)
library(tidyverse)
library(httr)

#------------------ State Level Cleaning -----------------

data_path <- here("data_clean")

#------Census Population Estimates -----

url_2019 <- "https://www2.census.gov/programs-surveys/popest/datasets/2010-2019/counties/totals/co-est2019-alldata.csv"

population_2019 <- read_csv(url_2019) %>%
  mutate(fips_code = paste0(STATE, COUNTY)) %>%
  select(fips_code, population =POPESTIMATE2019)




#------- Covidestim (State Level) -----

# Due to changes in the model to handle the Omicron wave, 
# we have to obtain dates before 2021-12-02 from a different
# source than the latest runs endpoint of the API.


get_covidestim_biweekly <- function() {
  
  # set end date
  end_date <- ymd("2022-02-25")
  
  
  covidestim_api_link <- "https://api2.covidestim.org/latest_runs?geo_type=eq.state&select=*,timeseries(*)"
  covidestim <- GET(covidestim_api_link)
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # access dates BEFORE 2021-12-02
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  last_weeks <-jsonlite::fromJSON(
    httr::content(covidestim,
                  as = "text",
                  encoding = "UTF-8"),
    simplifyVector = TRUE,
    flatten = TRUE)  %>%
    unnest(cols=timeseries, 
           names_repair = "unique") %>%
    mutate( date = ymd(date),
            week = week( date)) %>%
    filter(date <= end_date & year(date) > 2020)
  
  last_weeks <- last_weeks %>%
    select(date, 
           infections.lo = infections_p2_5,
           infections.hi =infections_p97_5,
           infections, created_at,
           state = geo_name) %>%
    mutate(week = week(date),
           created = substr(created_at, 1, 10),
           created = ymd(created))
  
  # remove duplicates from multiple model runs by taking most recent
  last_weeks <- last_weeks %>%
    group_by(week,date) %>%
    slice_max(order_by = created, n=1) %>%
    select(-c(created,created_at)) %>%
    mutate(year = year(date)) %>%
    mutate(week = case_when(
      year == 2022 ~ week + 52,
      year == 2021 ~ week
    )) %>%
    ungroup()
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # access dates AFTER 2021-12-02
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  covidestim_link <- "https://covidestim.s3.us-east-2.amazonaws.com/latest/state/estimates.csv"
  
  covidestim <- read_csv(covidestim_link)
  
  # join data from each source to include dates before and after 2021-12-02
  covidestim <- covidestim %>%
    select(date, contains("infections"), state) %>%
    filter(date <= end_date & year(date) > 2020) %>%
    mutate(week = week(date), year = year(date)) %>%
    mutate(week = case_when(
      year == 2022 ~ week + 52,
      year == 2021 ~ week
    )) %>%
    bind_rows(last_weeks) %>%
    group_by(week, state, year) %>%
    summarize(across(contains("infections"), sum), 
              date = min(date)) %>%
    ungroup()
  
  #~~~~~~~~~~~~~~~~
  # add biweek
  #~~~~~~~~~~~~~~~~
  num_weeks <- covidestim %>% 
    pull(week) %>% 
    unique() %>% 
    length()
  
  num_biweeks <- num_weeks/2
  
  biweek <- tibble(biweek = c(rep(1:num_biweeks, 2))) %>%
    arrange(biweek)
  
  biweek_to_week <- covidestim %>%
    select(week) %>%
    distinct() %>%
    arrange(week) %>%
    cbind(biweek =biweek)
  
  # calculate biweekly sums
  covidestim_biweekly <- covidestim %>%
    left_join(biweek_to_week) %>%
    group_by(biweek, state)  %>%
    mutate(across(
      contains("infections"), 
      sum)) %>%
    ungroup()
  
  # add state codes
  statecodes <- read_csv(
    here("data/demographic/statecodes.csv")) %>%
    rename(state_name = state,
           state = code) %>%
    select(-abbrev)
  
  covidestim_biweekly <- covidestim_biweekly %>%
    rename(state_name = state) %>%
    left_join(statecodes)
  
  return(covidestim_biweekly)
  
}

# covidestim <- get_covidestim_biweekly()
# saveRDS(here("data/data_clean/covidestim_biweekly_all_states.RDS"))






