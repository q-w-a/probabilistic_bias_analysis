
library(here)
library(tidyverse)
library(httr)

#------------------ State Level Cleaning -----------------

# data_path <- here("data/data_clean")

#------Census Population Estimates -----


get_county_pop <- function() {
  url_2019 <- "https://www2.census.gov/programs-surveys/popest/datasets/2010-2019/counties/totals/co-est2019-alldata.csv"
  
  population_2019 <- read_csv(url_2019) %>%
    mutate(fips_code = paste0(STATE, COUNTY)) %>%
    select(fips_code, population =POPESTIMATE2019)
  
  return(population_2019)
}




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
    summarize(across(
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
    ungroup() %>%
    select(-week)
  
  return(covidestim_biweekly_all_counties)
  
}





#--------- Archive Covidestim Data -----------

archive_covidestim_data <- function() {
  
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
    select(state, positive, total, negative, date)
  
  cdc_pos <- cdc_pos %>%
    left_join(dates) %>%
    filter(!is.na(biweek)) %>%
    group_by(biweek, state) %>%
    mutate(positive = sum(positive, na.rm=TRUE),
           total = sum(total, na.rm = TRUE),
           negative=sum(negative)) %>%
    ungroup() %>%
    mutate(posrate = positive/total,
           negative = total) %>%
    left_join(state_pop) 
  
  
  return(cdc_pos)


} 


#--------- Archive County Data -----------


archive_county_data <- function() {
  
  # Michigan
  link <- paste0("https://www.michigan.gov/coronavirus/-/media/Project/Websites/coronavirus/",
                 "Michigan-Data/09-27-2022/Datasets/Diagnostic-Tests-by-Result-and-County-2022-09-27.xlsx?",
                 "rev=7ba61151dcff4b038e33ea9ead95137c&hash=6C1AA97E93A30A3C1F10171233FDF31A")
  httr::GET(link, httr::write_disk(tf <- tempfile(fileext = ".xlsx")))
  df <- readxl::read_excel(tf)
  saveRDS(df, file =here("data/data_raw/mi_county_original.RDS"))
  
  
  # Massachusetts
  mass_link <- "https://www.mass.gov/doc/covid-19-raw-data-april-12-2022/download"
  httr::GET(mass_link, httr::write_disk(tf <- tempfile(fileext = ".xlsx")))
  
  df <- readxl::read_excel(tf, sheet = "County_Weekly")
  df <- df %>% 
    select(contains("date"), 
           total = `Total Tests (Last 14 days)`,
           positive = `Total Positive Tests (Last 14 days)`,
           county_name =County) %>%
    mutate( negative = total - positive,
            state = "MA") 
  saveRDS(df, file =here("data/data_raw/ma_county_original.RDS"))
  
}


# archive_county_data()


#--------- Michigan County Data -----------


get_michigan_county <- function(end_date = "2022-02-25") {
  
  counts_raw <- readRDS(here(
    "data/data_raw/mi_county_original.RDS"))
  county_fips <- read_tsv(here("data/demographic/county_fips.tsv")) %>%
    dplyr::rename_with(.cols =everything(), tolower)
  population_2019 <- get_county_pop()
  
  # basic reformatting; join fips code and population
  counts_raw <- counts_raw %>%
    mutate(week = week(MessageDate),
           date = ymd(MessageDate)) %>%
    dplyr::rename_with(.cols =everything(), tolower) %>%
    mutate(state = "MI") %>%
    left_join(county_fips,
              by = c("county"="name",
                     "state" = "state")) %>%
    select(-messagedate) %>%
    left_join(population_2019, by = c("fips" = "fips_code"))
  
  county <- counts_raw %>% 
    mutate(
      week = week(date),
      year = year(date)) %>%
    filter(year > 2020 & date <= end_date) %>%
    mutate( week = case_when(
      year == 2021 ~ week,
      year == 2022 ~ week + 52)) %>%
    select(fips, county, positive, 
           date, total, week, 
           negative, population) %>%
    mutate(state = "MI")
  
  # add biweek
  num_weeks <- unique(county$week) %>%
    length()
  
  num_biweeks <- num_weeks/2
  
  biweek <- tibble(biweek = c(rep(1:num_biweeks, 2))) %>%
    arrange(biweek)
  
  biweek_to_week <- county %>%
    select(week) %>%
    distinct() %>%
    arrange(week) %>%
    cbind(biweek =biweek)
  
  county_biweekly <- county %>%
    left_join(biweek_to_week) %>%
    group_by(biweek, fips)  %>%
    mutate(across(c(total,positive, negative),
                  sum)) %>%
    ungroup() %>%
    mutate(posrate = positive/ total) %>%
    select(-week) %>%
    filter(!is.na(fips))
  
  
  
  county_biweekly <- county_biweekly %>%
    select(-date) %>%
    distinct()
  
  return(county_biweekly) 
  
}

#---------- Massachusetts County Data -----------


get_mass_county <- function(end_date = "2022-02-25") {
  
  county <- readRDS(here("data/data_raw/ma_county_original.RDS"))
  
  county_fips <- read_tsv(here("data/demographic/county_fips.tsv")) %>%
    dplyr::rename_with(.cols =everything(), tolower)
  
  population_2019 <- get_county_pop()
  
  # Dukes and Nantucket are combined
  together <- county_fips %>% 
    filter(name == "Dukes" | name == "Nantucket") %>%
    pull(fips)
  
  county_fips <- county_fips %>%
    filter(!(fips %in% together)) %>%
    bind_rows(tibble(name = "Dukes and Nantucket", 
                     fips = paste0(together, 
                                   collapse=","),
                     state = "MA"))
  
  
  
  county <- county %>%
    filter(!(is.na(total) | is.na(positive))) %>%
    filter(!county_name %in% c("All of Massachusetts", 
                                        "Unknown County")) %>%
    rename(
           start_period_date = `Start Date`,
           end_period_date = `End Date`) %>%
    mutate(year = year(start_period_date),
           week = week(start_period_date)) %>%
    filter(year> 2020 & start_period_date <= end_date) %>%
    mutate( week = case_when(
      year == 2021 ~ week,
      year == 2022 ~ week + 52)) %>%
    mutate(across(contains("date"),ymd)) %>%
    select(-c(`Report Date`))
  
  county <- county %>%
    mutate(county_name = gsub(" County", "", county_name),
           county_name = gsub(" Counties", "", county_name),
           county_name = trimws(county_name)) %>%
    left_join(county_fips, 
              by = c("state"="state", "county_name"="name"))
  
  
  
  # note that dates are overlapping, so remove the duplicates
  startdate_to_biweek <- county %>%
    filter(!week %% 2 ==0) %>%
    select(-week) %>%
    select(start_period_date) %>%
    distinct() %>%
    arrange(start_period_date) %>%
    mutate(biweek = row_number())
  
  
  
  
  # sum populations for grouped counties
  population_2019 <- population_2019 %>%
    mutate(fips_code = ifelse(fips_code %in% together,
                              paste0(together, 
                                     collapse=","), 
                              fips_code)) %>%
    group_by(fips_code) %>%
    summarize(population = sum(population))
  
  
  # note that summing by biweek is not needed because structure is already in 
  # 2 week interval format
  county  <- county %>%
    filter(!week %% 2 ==0) %>%
    left_join(startdate_to_biweek) %>%
    select(-c(week, end_period_date)) %>%
    left_join(population_2019, by = c("fips" = "fips_code")) %>%
    rename(date = start_period_date) %>%
    mutate(posrate = positive/total) 
  
  
  return(county)
  

  
}



#---------- COVID-19 Trends and Impact Survey Data --------------

save_survey_data <- function() {
  
  screening_data_link <- paste0(
    "https://api.covidcast.cmu.edu/epidata/covidcast/?data_source=fb-survey",
    "&signal=smoothed_wscreening_tested_positive_14d,smoothed_wtested_positive_14d,smoothed_wcli",
    "&geo_type=state&time_type=day&time_values=20210320-20221212&geo_value=*&api_key=TEMP-API-KEY-EXPIRES-2023-06-28")
  
  fb_screening <- httr::GET(screening_data_link)
  
  fb_symptoms <-jsonlite::fromJSON(
    httr::content(fb_screening,
                  as = "text", 
                  encoding = "UTF-8"))$epidata %>%
    mutate(date = lubridate::ymd(time_value),
           week = lubridate::week(date),
           state = geo_value,
           value = value/100,
           stderr = stderr/100) %>% 
    filter(date <= lubridate::ymd("2022-03-01")) %>%
    as_tibble()
  
  
  fb_symptoms %>%
    saveRDS(here(
      "data/data_raw/ctis_all_states.RDS"))
  
}

# save_survey_data()



#----------- Wastewater Data ----------------

get_biweekly_wastewater <- function() {
  
  biobot_link <- "https://raw.githubusercontent.com/biobotanalytics/covid19-wastewater-data/master/wastewater_by_county.csv"
  
  dates <- readRDS(here("data/data_raw/date_to_biweek.RDS"))
  
  w_data <- read_csv(biobot_link)%>% 
    filter(sampling_week >= ymd("2021-03-01") &  
             sampling_week <= ymd("2022-03-01")) %>%
    mutate(fips = as.character(fipscode)) %>%
    select(-fipscode)  %>%
    left_join(dates, by = c("sampling_week" = "date")) %>%
    group_by(fips, biweek, state,name) %>%
    summarize(mean_conc = mean(effective_concentration_rolling_average, na.rm= TRUE)) %>%
    ungroup()
  
  
  
  # fill in missing values with rolling mean
  w_data <- w_data %>% 
    group_by(fips) %>%
    # keep track of min biweek make sure to not fill in values
    # before first data point when using rolled mean
    mutate(min_biweek = min(biweek)) %>%
    ungroup() %>%
    pivot_wider(names_from = biweek,
                values_from = mean_conc) %>% 
    pivot_longer(cols = -c(fips,state, name, min_biweek), 
                 names_to = "biweek", 
                 values_to = "mean_conc") %>%
    mutate(biweek = as.numeric(biweek)) %>%
    group_by(fips) %>%
    arrange(biweek) %>%
    mutate(rolled_mean = RcppRoll::roll_mean(mean_conc,
                                             n = 4,
                                             na.rm = TRUE,
                                             fill = NA),
           mean_conc = ifelse(is.na(mean_conc) & biweek >= min_biweek, 
                              rolled_mean, mean_conc)) %>%
    ungroup()
  
  return(w_data)
  
}


#------------- Variant Data -----------------

get_variant_proportions <- function() {
  
  variant <- "https://data.cdc.gov/resource/jr58-6ysp.json?$limit=50000&$where=week_ending between '2021-02-18T00:00:00.000' and '2022-03-01T00:00:00.000'&usa_or_hhsregion=USA"
  variant <- httr::GET(URLencode(variant))
  
  
  # only go to the 6th because these are already weekly counts
  # get the last couple weeks of 2021
  variant <-jsonlite::fromJSON(
    httr::content(variant,
                  as = "text",
                  encoding = "UTF-8"),
    simplifyVector = TRUE,
    flatten = TRUE)  %>%
    as_tibble() 
  
  saveRDS(variant, 'data/data_raw/variant_prop.RDS')
  
  
  # https://www.cdc.gov/coronavirus/2019-ncov/variants/variant-classifications.html
  variant <- variant %>% 
    mutate(week_end_date = substr(week_ending, 1,10),
           week_end_date=ymd(week_end_date)) %>%
    mutate(variant_category = case_when(
      grepl("B[.]1[.]1[.]529|BA[.]1|BA[.]2|BA[.]4|BA[.]5", variant) ~ "Omicron",
      variant %in% c("B.1.621", "B.1.621.1") ~ "Mu",
      variant == "B.1.1.7" ~ "Alpha",
      grepl("B[.]1[.]351", variant) ~ "Beta",
      grepl("P.1", variant) ~"Gamma",
      grepl("B[.]1[.]617[.]2", variant) ~ "Delta",
      variant %in% c("B.1.427" ,"B.1.429") ~ "Epsilon",
      variant == "B.1.525" ~ "Eta",
      variant == "Other" ~ "Other"
    )) %>%
    mutate(share = as.numeric(share),
           creation_date = ymd(substr(creation_date,1,10))) %>%
    filter(modeltype == "weighted" & time_interval=="weekly") %>%
    group_by(variant, week_end_date, variant_category, time_interval) %>% 
    slice_max(n=1, order_by=creation_date) %>%
    group_by(variant_category, week_end_date, time_interval) %>%
    summarize(share =sum(share, na.rm=TRUE))
  
  
  variant <- variant %>%
    filter(variant_category != "Other") %>%
    select(week_end_date, variant_category, share) %>%
    # week beginning rather than week end
    mutate(week = week_end_date - days(7)) %>%
    filter(!is.na(variant_category)) 
  

  
  return(variant)
  
}



