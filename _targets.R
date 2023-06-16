# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed. # nolint

# Set target options:
tar_option_set(
  # packages needed for targets
  packages = c("tidyverse", "lubridate", 
               "here", "httr", "imputeTS"), 
  format = "rds" # default storage format
  # Set other options as needed.
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multicore")

# Run the R scripts in the R/ folder with your custom functions:
#tar_source("R/priors.R")
tar_source("R/01-data_cleaning.R")
tar_source("R/02-priors.R")
tar_source("R/03-base-functions.R")
tar_source("R/04-melding.R")
tar_source("R/05-state-analysis.R")




# source("other_functions.R") # Source other scripts as needed. # nolint

# Replace the target list below with your own:
list(
  #---------------- data -------------------
  # state-level covidestim
  tar_target(
    name = covidestim_biweekly_state,
    command = get_covidestim_state_biweekly()),
  # county-level covidestim
  tar_target(
    name = covidestim_biweekly_county,
    command = get_covidestim_county_biweekly()),
  # state-level testing data
  tar_target(
    name = tests_biweekly_state,
    command = get_state_testing()),
  # county-level testing data
  tar_target(
    name = mi_biweekly_county,
    command = get_michigan_county()),
  tar_target(
    name = ma_biweekly_county,
    command = get_mass_county()),
  
  #---------------- analysis -------------------
  # set priors
  tar_target(
    name = prior_params,
    command= get_priors()
  ),
  # tar_target(
  #   name = state_v1,
  #   command = get_v1_state(
  #       data = tests_biweekly_state,
  #       params = prior_params)
  # ),
  # get smoothed CTIS survey data
  tar_target(
    name = ctis_smoothed,
    command = get_smoothed_ctis(
      ctis = readRDS(here('data/data_raw/ctis_all_states.RDS')),
      smooth_beta_span = .33,
      smooth_s_untested_span =.2)
  )

)
