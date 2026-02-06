# Run this script to perform a temporal preservation analysis on the demo data

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

# Import demo data and split it into two 
# (the other part is treated as synthetic data)
demo <- read_csv("data/demo.csv")

sampled_ids <- sample(unique(demo$`Subject ID`), 
                      size = length(unique(demo$`Subject ID`))/2)
demo1 <- demo[demo$`Subject ID` %in% sampled_ids, ]
demo2 <- demo[!(demo$`Subject ID` %in% sampled_ids), ]

# Running the metrics

start_time <- Sys.time()

# In the current implementation, the id column is expected to be "Subject ID"
# and the analysis works only for id, time and time varying variables
analyze_and_plot_data(
  original = demo1[, c("Subject ID", "Time", "Response")],
  synthetic = demo2[, c("Subject ID", "Time", "Response")],
  id_col = "Subject ID",
  time_col = "Time",
  max_time = 5,
  x_interval = 1,
  bandwidth = 2,
  quantiles = c(0.05, 0.25, 0.5, 0.75, 0.95),
  save_path = NULL,  # Set to NULL to print to screen,
  mean_quantile = TRUE,
  cat_prop = TRUE,
  rank_stab = TRUE,
  meas_sim = TRUE,
  meas_n = 150,
  meas_iter = 100
)

end_time <- Sys.time()

print(end_time - start_time)
