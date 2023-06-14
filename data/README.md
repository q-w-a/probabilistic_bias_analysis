
# File Descriptions

### Demographic 
* `population_2019.RDS` contains the population for each fips code
* `statecodes.csv` contains state names mapped to codes from [the census](https://www2.census.gov/programs-surveys/popest/datasets/2010-2019/counties/totals/)
* `county_fips.tsv` contains fips code mapped to county name and state 
* `state_pop.RDS` contains state population estimates from [the census](https://www2.census.gov/programs-surveys/popest/datasets/2010-2019/state/detail/)

### Covidestim
* `covidestim_original_all.RDS` is an RDS file of the covidestim csv from the link `https://covidestim.s3.us-east-2.amazonaws.com/latest/state/estimates.csv`, unmodified; includes dates through 2021 earlier than 2021-12-02
* `covidestim_last_weeks_all_states.RDS` contains estimates after 2021-12-02 from `https://api2.covidestim.org/latest_runs?geo_type=eq.state&select=*,timeseries(*)`
* `covidestim_biweekly_all_states.RDS` is the cases summed for each 2-week interval, including data from `https://covidestim.s3.us-east-2.amazonaws.com/latest/state/estimates.csv` to obtain dates before 2021-12-02 (latest date available is 2021-11-30) and data from `https://api2.covidestim.org/latest_runs?geo_type=eq.state&select=*,timeseries(*)` to obtain dates after  2021-12-02
* `covidestim_county_estimates.RDS` contains estimates of infections by date and fips from `https://covidestim.s3.us-east-2.amazonaws.com/latest/estimates.csv`, no modifications other than selecting variables of interest
* `covidestim_biweekly_all_counties.RDS` contains biweekly aggregates for all counties, data again from link `https://covidestim.s3.us-east-2.amazonaws.com/latest/estimates.csv`