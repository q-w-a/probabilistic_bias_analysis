# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
library(crew)
# library(tarchetypes) # Load other packages as needed. # nolint

# Set target options:
tar_option_set(
  # packages needed for targets
  packages = c("tidyverse", "lubridate", 
               "here", "httr", "imputeTS"), 
  controller=crew_controller_local(workers=4),
  format = "rds" # default storage format
  # Set other options as needed.
)

# tar_make_clustermq() configuration (okay to leave alone):
# options(clustermq.scheduler = "multicore")

# Run the R scripts in the R/ folder with your custom functions:
#tar_source("R/priors.R")
tar_source("R/01-data_cleaning.R")
tar_source("R/02-priors.R")
tar_source("R/03-base-functions.R")
tar_source("R/04-melding.R")
tar_source("R/05-state-analysis.R")
tar_source("R/06-county-analysis.R")




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
  tar_target(
    name = variant,
    command = get_variant_proportions()),
  tar_target(
    name = state_deaths,
    command = get_state_deaths()),
  
  #---------------- state-level analysis -------------------
  # set priors
  tar_target(
    name = prior_params,
    command= get_priors()
  ),
  # version 1
  tar_target(
    name = state_v1,
    command = get_v1_state(
        data = tests_biweekly_state,
        params = prior_params)
  ),
  # get smoothed CTIS survey data
  tar_target(
    name = ctis_smoothed,
    command = get_smoothed_ctis(
      ctis_raw = readRDS(here('data/data_raw/ctis_all_states.RDS')),
      smooth_beta_span = .33,
      smooth_s_untested_span =.2)
  ),
  # version 2
  tar_target(
    name = state_v2,
    command = get_corrected_state(
      data = tests_biweekly_state,
      params = prior_params,
      ctis=ctis_smoothed,
      vary = "beta"
    )
  ),
  # version 3
  tar_target(
    name = state_v3,
    command = get_corrected_state(
      data = tests_biweekly_state,
      params = prior_params,
      ctis=ctis_smoothed,
      vary = "s_untested"
    )
  ),
  # version 4
  tar_target(
    name = state_v4,
    command = get_corrected_state(
      data = tests_biweekly_state,
      params = prior_params,
      ctis=ctis_smoothed,
      vary = "s_untested_and_beta"
    )
  ),
  # version 5
  tar_target(
    name = state_v5,
    command = get_corrected_state(
      data = tests_biweekly_state,
      params = prior_params,
      ctis=ctis_smoothed,
      spline=TRUE,
      vary = "beta"
    )
  ),
  # version 6
  tar_target(
    name = state_v6,
    command = get_corrected_state(
      data = tests_biweekly_state,
      params = prior_params,
      ctis=ctis_smoothed,
      spline=TRUE,
      vary = "s_untested_and_beta"
    )
  ),
  # version 7
  tar_target(
    name = state_v7,
    command = get_v7_state(
      data = tests_biweekly_state,
      beta_shape = get_shape(ctis_smoothed,
                             probs=c(.1,.5,.9),
                             option="beta"),
      s_untested_shape = get_shape(ctis_smoothed,
                                   option="s_untested")
    )
  ),
  
  
  #---------------- county-level analysis -------------------
  
  #------ MA ---------
  # version 1
  tar_target(
    name = ma_v1,
    command = get_county_v1(
      data = ma_biweekly_county,
      params = prior_params
    )
  ),
  # version 2
  tar_target(
    name = ma_v2,
    command = get_corrected_county(
      data = ma_biweekly_county,
      params = prior_params,
      ctis=ctis_smoothed,
      vary = "beta"
    )
  ),
  # version 3
  tar_target(
    name = ma_v3,
    command = get_corrected_county(
      data = ma_biweekly_county,
      params = prior_params,
      ctis=ctis_smoothed,
      vary = "s_untested"
    )
  ),
  # version 4
  tar_target(
    name = ma_v4,
    command = get_corrected_county(
      data = ma_biweekly_county,
      params = prior_params,
      ctis=ctis_smoothed,
      vary = "s_untested_and_beta"
    )
  ),
  
  # version 5
  tar_target(
    name = ma_v5,
    command = get_corrected_county(
      data = ma_biweekly_county,
      params = prior_params,
      ctis=ctis_smoothed,
      spline=TRUE,
      vary = "beta"
    )
  ),
  
  # version 6
  tar_target(
    name = ma_v6,
    command = get_corrected_county(
      data = ma_biweekly_county,
      params = prior_params,
      ctis=ctis_smoothed,
      spline=TRUE,
      vary = "s_untested_and_beta"
    )
  ),
  
  # version 7
  tar_target(
    name = ma_v7,
    command = get_county_v7(
      data = ma_biweekly_county,
      beta_shape = get_shape(ctis_smoothed,
                             probs=c(.1,.5,.9),
                             option="beta"),
      s_untested_shape = get_shape(ctis_smoothed,
                                   option="s_untested")
    )
  ),
  
  tar_target(
    name = waste,
    command = get_biweekly_wastewater())
  
  
)

# tar_read(ma_v2)
# tar_read(ma_v3)
# tar_read(ma_v4)


# for testing
# ctis <- tar_read(ctis_smoothed)
# params <-  list(
# alpha_mean = .95,
# alpha_sd = 0.08,
# alpha_bounds = NA,
# # alpha_bounds = c(.8,1),
# beta_mean = .15,
# beta_sd =.09,
# beta_bounds = NA,
# #  beta_bounds = c(0.002, 0.4),
# s_untested_mean = .03,
# s_untested_sd = .0225,
# #  s_untested_bounds = c(0.0018, Inf),
# s_untested_bounds = NA,
# p_s0_pos_mean = .4,
# p_s0_pos_sd = .1225,
# p_s0_pos_bounds = NA,
# #  p_s0_pos_bounds = c(.25, .7),
# pre_nsamp = 1e6,
# post_nsamp = 1e5)
# data <- tar_read(tests_biweekly_state)

