
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


get_covidestim_state_biweekly <- function(end_date = "2022-02-25") {
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # access dates AFTER 2021-12-02
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  covidestim_api_link <- "https://api2.covidestim.org/latest_runs?geo_type=eq.state&select=*,timeseries(*)"
  covidestim <- GET(covidestim_api_link)
  
  
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
  # access dates BEFORE 2021-12-02
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










#------- Covidestim (County Level) -----

get_covidestim_county_biweekly <- function(end_date = "2022-02-25") {
  legacy_link <- "https://covidestim.s3.us-east-2.amazonaws.com/latest/estimates.csv"
  
  covidestim_county <- read_csv(legacy_link) %>%
    select(fips, date, infections)
  
  
  covidestim_county <- covidestim_county %>%
    filter(date <= end_date &  year(date) > 2020) %>%
    mutate(week = week(date),
           year = year(date)) %>%
    mutate(week = case_when(
      year == 2022 ~ week + 52,
      year == 2021 ~ week
    ))
  
  
  num_weeks <- covidestim_county %>% 
    pull(week) %>% 
    unique() %>% 
    length()
  
  num_biweeks <- num_weeks/2
  
  biweek <- tibble(biweek = c(rep(1:num_biweeks, 2))) %>%
    arrange(biweek)
  
  biweek_to_week <- covidestim_county %>%
    select(week) %>%
    distinct() %>%
    arrange(week) %>%
    cbind(biweek =biweek)
  
  
  covidestim_biweekly_all_counties <- covidestim_county %>%
    left_join(biweek_to_week) %>%
    group_by(biweek, fips)  %>%
    mutate(across(contains("infections"), sum)) %>%
    ungroup() 
}





#--------- Archive Data -----------

archive_data <- function() {
  
  # save county estimates
  legacy_link <- "https://covidestim.s3.us-east-2.amazonaws.com/latest/estimates.csv"
  
  estimates_legacy <- read_csv(legacy_link) 
  
  estimates_legacy <- estimates_legacy %>%
    select(fips, date, infections) %>%
    filter(year(date) >= 2020 & 
             substr(fips,1,2) %in% c("25", "26")) 
  
  saveRDS(
    estimates_legacy,
    here("data/data_raw/covidestim_county_estimates.RDS"))
  
  
  # save state estimates
  
  # access dates AFTER 2021-12-02
  covidestim_api_link <- "https://api2.covidestim.org/latest_runs?geo_type=eq.state&select=*,timeseries(*)"
  covidestim <- GET(covidestim_api_link)
  
  saveRDS(
    covidestim,
    here("data/data_raw/covidestim_state_2021_to_2022.RDS"))
  
  
  # access dates BEFORE 2021-12-02
  covidestim_link <- "https://covidestim.s3.us-east-2.amazonaws.com/latest/state/estimates.csv"
  covidestim_allstates <- read_csv(covidestim_link) %>%
    filter(year(date) == 2021)
  
  saveRDS(
    covidestim_allstates,
    here("data/data_raw/covidestim_state_2021.RDS"))
  
}

# archive_data()


#----------- State-Level Testing -------

get_state_testing <- function() {

    
  dates <- readRDS(here("data/data_raw/date_to_biweek.RDS"))
  
  
  # PAGE: https://www2.census.gov/programs-surveys/popest/datasets/2010-2019/state/detail/
  state_population_link <- "https://www2.census.gov/programs-surveys/popest/datasets/2010-2019/state/detail/SCPRC-EST2019-18+POP-RES.csv"
  state_pop <- read_csv(state_population_link)
  
  
  
  statecodes <- read_csv(paste0(
    here::here(), 
    "/data/demographic/statecodes.csv"))
  
  state_pop <- state_pop %>%
    left_join(statecodes, by = c("NAME" = "state")) %>%
    select(population = POPESTIMATE2019,
           state = code) %>%
    filter(!is.na(state))
  
  
  
  # split into two queries to ensure we obtain all data with the API limit
  
  dat <- httr::GET(URLencode(
    paste0("https://healthdata.gov/resource/j8mb-icvb.json?",
           "$where=date between '2020-12-30T12:00:00' and '2021-10-15T12:00:00'&$limit=50000")))
  
  cdc1 <-jsonlite::fromJSON(
    httr::content(dat,
                  as = "text", 
                  encoding = "UTF-8")) %>%
    as_tibble()
  
  
  dat2 <- httr::GET(URLencode(
    paste0("https://healthdata.gov/resource/j8mb-icvb.json?",
           "$where=date between '2021-10-15T14:00:00' and '2022-02-25T14:00:00'&$limit=50000")))
  
  
  cdc2 <-jsonlite::fromJSON(
    httr::content(dat2,
                  as = "text", 
                  encoding = "UTF-8")) %>%
    as_tibble()
  
  
  cdc <- bind_rows(cdc1, cdc2) %>%
    mutate(date = ymd(substr(date,1,10))) %>%
    mutate(across(c(new_results_reported), as.numeric)) %>%
    filter(!state %in% c("MP", "AS", "GU", "PR", "VI", "MH"))
  
  # overall_outcome is the outcome of the test (Inconclusive, Negative, or Positive)
  # new_results_reported is the number with the given outcome
  cdc_pos <- cdc %>%
    select(-c(fema_region, total_results_reported)) %>%
    pivot_wider(names_from = c("overall_outcome"),
                values_from = c("new_results_reported")) %>%
    mutate(total = Inconclusive + Negative + Positive) %>%
    rename_with(tolower) %>%
    select(state, positive, total, date)
  
  cdc_pos %>%
    left_join(dates) %>%
    filter(!is.na(biweek)) %>%
    group_by(biweek, state) %>%
    mutate(positive = sum(positive, na.rm=TRUE),
           total = sum(total, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(posrate = positive/total) %>%
    left_join(state_pop) 


} 