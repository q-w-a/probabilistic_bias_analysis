# Probabalistic Bias Analysis to Approximate COVID-19 Infection Counts

Applying probabilistic bias analysis to approximate the true number of COVID-19 infections, at the state level in all states and at the county level in Massachusetts.


## Target Descriptions

### Data

#### Covidestim

See the documentation [here](covidestim.org) for documentation on Covidestim. Targets for Covidestim data include:

* `covidestim_biweekly_state`: Covidestim infection counts summed for each 2-week interval at the state level
  - `infections`: median
  - `infections.hi`: 97.5th percentile
  - `infections.lo`: 2.5th percentile
* `covidestim_biweekly_county`: Covidestim infection counts summed for each 2-week interval at the county level, only Massachusetts and Michigan
  - `fips`: county FIPS code
  - `infections`: median

#### State-level Testing Data

* `tests_biweekly_state`: PCR test results at the state-level from healthdata.gov endpoint [here](https://healthdata.gov/dataset/COVID-19-Diagnostic-Laboratory-Testing-PCR-Testing/j8mb-icvb)
    * Columns: 
      - `state`
      - `positive` number of positive tests
      - `total` total number of tests
      - `date` first date of 2-week interval
      - `biweek`
      - `posrate` biweekly positivity rate, `positive`/`total`
      - `population` census population size for state
      
#### County-level Testing Data

* `mi_biweekly_county` and `ma_biweekly_county`
  - Columns: 
      -  `fips`: county FIPS code
      - `county`: county name 
      - `positive`: number of positive tests
      - `negative`: number of negative tests
      - `total`: number of total tests 
      - `state`: 2-letter state abbreviation
      - `biweek`: 2-week interval
      - `posrate`: biweekly positivity rate

### Results

#### State Results 

* `state_v1`: priors do not vary by state or date
* `state_v2`: prior for $\beta$ is centered at ratio of screening test positivity to overall test positivity from the COVID-19 Trends and Impact Survey
* `state_v3`: prior for $\Pr(S_1|\text{untested})$ is centered at ratio of screening test positivity to overall test positivity from the COVID-19 Trends and Impact Survey
* `state_v4`: prior for $\Pr(S_1|\text{untested})$ is centered at ratio of screening test positivity to overall test positivity from the COVID-19 Trends and Impact Survey

#### County Results 

##### Massachusetts

State-level [COVID-19 Trends and Impact Survey](https://delphi.cmu.edu/covid19/ctis/) data used to allow priors to vary by date.

* `ma_v1`: priors do not vary by date
* `ma_v2`: prior for $\beta$ is centered at ratio of screening test positivity to overall test positivity from the COVID-19 Trends and Impact Survey
* `ma_v3`: prior for $\Pr(S_1|\text{untested})$ is centered at ratio of screening test positivity to overall test positivity from the COVID-19 Trends and Impact Survey
* `ma_v4`: prior for $\Pr(S_1|\text{untested})$ is centered at ratio of screening test positivity to overall test positivity from the COVID-19 Trends and Impact Survey

# Repository Structure

All analyses can be reproduced by running `tar_make()`, which runs the file containing the full pipeline, `_targets.R`. 

Directories:

- `_targets`: where targets from  `_targets.R` are saved; can be accessed with `tar_read` more quickly rather than reading the RDS files directly from `_targets/objects`
- `R`: scripts forming the core of the analyses, which are run in the pipeline in  `_targets.R`
- `app`: Shiny app for exploring implications of changing priors
