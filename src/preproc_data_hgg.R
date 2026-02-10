### This file preprocessess the original and synthetic datasets ###
### Load necessary libraries ###

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

########################### Functions ######################################

### Mask out the columns using the mask columns which end with _m

mask_by_suffix <- function(df, suffix = "_m", drop_mask = FALSE, ignore_case = TRUE) {
  stopifnot(is.data.frame(df))
  
  # Case-insensitive
  nm <- names(df)
  nm_cmp <- if (ignore_case) tolower(nm) else nm
  suff_cmp <- if (ignore_case) tolower(suffix) else suffix
  
  # Finds mask columns and their base names
  is_mask <- endsWith(nm_cmp, suff_cmp)
  mask_cols <- nm[is_mask]
  if (length(mask_cols) == 0) return(df)
  
  strip_regex <- paste0("(", gsub("([\\W])", "\\\\\\1", suffix, perl = TRUE), ")$") 
  if (ignore_case) strip_regex <- paste0("(?i)", strip_regex)
  
  base_names <- sub(strip_regex, "", mask_cols, perl = TRUE)
  
  for (i in seq_along(mask_cols)) {
    mask_col <- mask_cols[i]
    base_col <- base_names[i]
    
    # Skip if the paired base column isn't present
    if (!base_col %in% names(df)) next
    
    m <- suppressWarnings(as.numeric(df[[mask_col]]))
    zero_idx <- !is.na(m) & m == 0
    
    if (any(zero_idx)) {
      df[[base_col]][zero_idx] <- NA
    }
  }
  
  if (drop_mask) {
    df <- df[ , setdiff(names(df), mask_cols), drop = FALSE]
  }
  
  df
}
######################################### Read data and pre-proccesing ####################################################


synthetic <- read_csv("data/train_synth_mimic/Hgym_synth.csv")
original <- read_csv("data/train_synth_mimic/Hgym_real.csv")
synthetic <- synthetic[ , names(original)]
synthetic <- synthetic %>%
  mutate(across(c(3, 4,5, 7, 9,10,13), ~ round(.x, 0)))

synthetic <- synthetic %>%
  mutate(across(c(14, 15), ~ round(.x, 2)))

# Time for each subject from 0-47 to 1-48
original$Time <- original$Time + 1

# Make sure time is sorted per subject in synthetic data
synthetic <- synthetic[order(synthetic$"Subject ID", synthetic$Time), ]
synthetic$Time <- ave(synthetic$Time, synthetic$"Subject ID", FUN = seq_along)


# Keep as character 
original$`Subject ID` <- as.character(original$`Subject ID`)
u_ids <- unique(original$`Subject ID`)           
original$`Subject ID` <- match(original$`Subject ID`, u_ids)
n_patients <- length(u_ids)


## Fluid boluses were binned into one of the four categories of none, low, medium, and high reﬂecting the numeric ranges of 0, [250, 500), [500, 1000), and ≥ 1000, in mL units

original$fluid_boluses <- with(original,
                               ifelse(fluid_boluses == 0, "none",
                                      ifelse(fluid_boluses >= 250 & fluid_boluses < 500,  "low",
                                             ifelse(fluid_boluses >= 500 & fluid_boluses < 1000, "medium",
                                                    ifelse(fluid_boluses >= 1000, "high", NA))))
)


synthetic$fluid_boluses <- with(synthetic,
                                ifelse(fluid_boluses == 0, "none",
                                       ifelse(fluid_boluses >= 250 & fluid_boluses < 500,  "low",
                                              ifelse(fluid_boluses >= 500 & fluid_boluses < 1000, "medium",
                                                     ifelse(fluid_boluses >= 1000, "high", NA))))
)

## vasopressor variable was categorized according to the ranges of 0, (0, 8.4), [8.4, 20.28), and ≥ 20.28,


Q1 <- 8.4
Q3 <- 20.28

original$vasopressors <- with(original,
                              ifelse(vasopressors == 0, "none",
                                     ifelse(vasopressors > 0 & vasopressors < Q1, "low",
                                            ifelse(vasopressors >= Q1 & vasopressors < Q3, "medium",
                                                   ifelse(vasopressors >= Q3, "high", NA))))
)


synthetic$vasopressors <- with(synthetic,
                               ifelse(vasopressors == 0, "none",
                                      ifelse(vasopressors > 0 & vasopressors < Q1, "low",
                                             ifelse(vasopressors >= Q1 & vasopressors < Q3, "medium",
                                                    ifelse(vasopressors >= Q3, "high", NA))))
)



#FiO2 was the fraction of inspired oxygen, and binned the variable into 10 classes according to the ranges of [0, 0.1), [0.1, 0.2), ... , [0.8, 0.9), and [0.9, 1.0].

x <- pmax(pmin(as.numeric(original$FiO2), 1), 0)

brks <- seq(0, 1, by = 0.1)
labs <- paste0("Range-", sprintf("%.1f", head(brks, -1)))  # range-0.0 ... range-0.9

# 1.0 so it falls into [0.9,1.0]
x2 <- pmin(x, 0.9999999)

original$FiO2 <- cut(
  x2, breaks = brks, labels = labs,
  right = FALSE, include.lowest = TRUE
)
# Just to check
round(100 * prop.table(table(original$FiO2)), 2)

x <- pmax(pmin(as.numeric(synthetic$FiO2), 1), 0)
brks <- seq(0, 1, by = 0.1)
labs <- paste0("Range-", sprintf("%.1f", head(brks, -1)))  
x2 <- pmin(x, 0.9999999)
synthetic$FiO2 <- cut(
  x2, breaks = brks, labels = labs,
  right = FALSE, include.lowest = TRUE
)
# Just to check
round(100 * prop.table(table(synthetic$FiO2)), 2)

# GCS has 15 values [1,..,15] therefore 15 classes

x <- as.integer(original$GCS_total)              
x[ x < 1 | x > 15 ] <- NA_integer_         

original$GCS_total <- factor(
  paste0("R-", x),                         
  levels = paste0("R-", 1:15),            
  ordered = TRUE
)


x <- as.integer(synthetic$GCS_total)              
x[ x < 1 | x > 15 ] <- NA_integer_         

synthetic$GCS_total <- factor(
  paste0("R-", x),                        
  levels = paste0("R-", 1:15),  
  ordered = TRUE
)

####################################################################################################

synthetic <- synthetic[ , names(original)]
synthetic_masked <- mask_by_suffix(synthetic, drop_mask = TRUE)
original_masked <- mask_by_suffix(original, drop_mask = TRUE)

export <- list(
  original_masked  = original_masked,
  synthetic_masked = synthetic_masked
)

export  # <-- last evaluated expression



