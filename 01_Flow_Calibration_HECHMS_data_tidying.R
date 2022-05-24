# Script to import storm flow data and performance metrics and assemble into comomon datasets
# by Annie Holt

# load libraries
library(tidyverse) #for data clearning
library(purrr) #for iterations

#### LISTING FILES #### ------------------------------------------------------------------------------------------------------#

# set working directory
# data is on Teams folder, synced locally
# update path depending on user
setwd("C:/Users/anneh/SCCWRP/OPC Sediment Flux to Coast - SD Creek Modeling/Model/Final")
getwd()

# create lists of file names in each directory desired
# just want text files that have names including storm or performance metrics
# decided to just paste the paths and list of file names together
filenames_flow <- paste(getwd(), 'ModeledFlow', list.files('./ModeledFlow', pattern = '^storm.*txt'), sep = "/")
filenames_metrics <- paste(getwd(), 'PerformanceMetrics', list.files('./PerformanceMetrics', pattern = '^PerformanceMetrics_storm.*txt'), sep = "/")
filenames_metrics_weighted <- paste(getwd(), 'PerformanceMetrics', 'WeightedMetrics', list.files('./PerformanceMetrics/WeightedMetrics', pattern = '^storm.*txt'), sep = "/")


#### CREATING DATA IMPORT FUNCTIONS #### ------------------------------------------------------------------------------------------------------#

# function to read in csv file and add an ID column with the number in the file name
# for example, 'storm_1.txt' will get an ID of '1'
dat_fun <- function(filename){
  
  df <- read_csv(filename) %>% 
    # create ID column based on number in filename
    mutate(StormID = as.numeric(gsub(".*?([0-9]+).*", "\\1", filename)))
  
}

# had to create a second function to deal with the weirdly formatted weighted metrics dataset 
# reads in data into one column, then separates the data correctly
dat_fun_2 <- function(filename){
  
  df <- read_csv(filename,  col_types = cols_only(metric.Weighted = col_character())) %>% 
    # create ID column based on number in filename
    mutate(StormID = as.numeric(gsub(".*?([0-9]+).*", "\\1", filename))) %>% 
    # deal with poor formatting
    rename(baddata = metric.Weighted) %>% 
    # get everything before comma for first column, everything after comma for second column of data
    mutate(metric.num = sub(',.*$','', baddata), metric.Weighted = sub('.*,\\s*', '', baddata)) %>% 
    mutate(metric.Weighted = case_when(metric.Weighted == "NA" ~ as.numeric(NA), 
                                       TRUE ~ as.numeric(metric.Weighted))) %>% 
    select(-baddata)
  
}


#### FINAL RUN #### ------------------------------------------------------------------------------------------------------#

# read in text files and append together into one dataframe using 'map_df()'
# use the data functions created above
ModeledFlow_all_storms <- purrr::map_df(filenames_flow, ~dat_fun(.x))
PerformanceMetrics_all_storms <- purrr::map_df(filenames_metrics,  ~dat_fun(.x))
WeightedMetrics_all_storms <- purrr::map_df(filenames_metrics_weighted, ~dat_fun_2(.x))


####  EXPORT #### ------------------------------------------------------------------------------------------------------#

# save the files in the repo and in the original locations

# # first original directories
# write_csv(ModeledFlow_all_storms, "./ModeledFlow/ModeledFlow_all_storms.csv")
# write_csv(PerformanceMetrics_all_storms, "./PerformanceMetrics/PerformanceMetrics_all_storms.csv")
# write_csv(WeightedMetrics_all_storms, "./PerformanceMetrics/WeightedMetrics/WeightedMetrics_all_storms.csv")
# 
# # changing directory so can save on local GitHub repo
# setwd("C:/Users/anneh/Documents/Repositories/SDCreekModeling/output_01")
# write_csv(ModeledFlow_all_storms, "./ModeledFlow_all_storms.csv")
# write_csv(PerformanceMetrics_all_storms, "./PerformanceMetrics_all_storms.csv")
# write_csv(WeightedMetrics_all_storms, "./WeightedMetrics_all_storms.csv")


