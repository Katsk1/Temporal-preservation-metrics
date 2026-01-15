# Load necessary libraries
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(stringr)
library(gtable)
library(grid)
library(patchwork)
library(scales)
library(reshape2)
library(data.table)
library(RcppHungarian)

source("src/helpers.R")
source("src/analyze_and_plot.R")

# Bring only the data
tmp <- new.env(parent = .GlobalEnv)
res <- sys.source("src/preproc_data_halo.R", envir = tmp) 
wanted <- c("original_data", "synthetic_data") 
list2env(mget(wanted, envir = tmp, inherits = FALSE), envir = .GlobalEnv)
rm(tmp, res)

### Running the metrics ###

start_time <- Sys.time()

# Documentation in the respective .R file

analyze_and_plot_data(
  original = original_data,
  synthetic = synthetic_data,
  id_col = "Subject ID",
  time_col = "Time",
  max_time = 48,
  x_interval = 12,
  bandwidth = 6,
  quantiles = c(0.05, 0.25, 0.5, 0.75, 0.95),
  save_path = "final_results_48h/",  # Set to NULL to print to screen,
  mean_quantile = TRUE,
  cat_prop = TRUE,
  rank_stab = TRUE,
  meas_sim = TRUE,
  meas_n = 2000,
  meas_iter = 100
)

end_time <- Sys.time()

print(end_time - start_time)
