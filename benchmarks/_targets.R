
#### Targets Master Script ####

### Libraries ###

# Organization
library(targets)
library(tarchetypes)

# Data Wrangling
library(tidyverse)

# Modeling
library(tidymodels)
library(tidycausality)

# Source Functions
source("benchmarks/functions/data_prep_f.R")

## Pipeline ##
list(

  # Get dataset for Benchmarking and Analysis
  tar_target(
    name = ad_smart_data,
    command = read.csv("benchmarks/data/AdSmartABdata - AdSmartABdata.csv")
  ),
  tar_target(
    name = online_ad_data,
    command = read_csv("benchmarks/data/online_ad_AB.csv")
  ),
  # Clean data
  tar_target(
    name = clean_smart_data,
    command = clean_ad_smart(ad_smart_data)
  )
)












