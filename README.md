# Probabalistic Bias Analysis to Approximate COVID-19 Infection Counts

Applying probablistic bias analysis to approximate the true number of COVID-19 infections, at the state level in all states and at the county level in Massachusetts.


## Target Descriptions

### Data

#### Covidestim

See the documentation [here](covidestim.org) for documentation on Covidestim. Targets for Covidestim data include:

* `covidestim_biweekly_state`: Covidestim infection counts summed for each 2-week interval at the state level
* `covidestim_biweekly_county`: Covidestim infection counts summed for each 2-week interval at the county level

#### State-level Testing Data

* `tests_biweekly_state`: PCR test results at the state-level from healthdata.gov endpoint [here](https://healthdata.gov/dataset/COVID-19-Diagnostic-Laboratory-Testing-PCR-Testing/j8mb-icvb)
* Columns: `state`, `positive` (number of positive tests),   `total` (total number of tests), `date` (first date of 2-week interval), `biweek`, `posrate` (biweekly positivity rate, `positive`/`total`), `population` (census population)



## Repository Structure

- `data_raw`
- `data_clean`
- `scripts`
- `figures`
- `output`
- `results`

