### This file preprocessess the original and synthetic datasets ###

# Import datasets
synthetic <- read_csv("data/train_synth_mimic/haloDataset.csv")
original <- read_csv("data/train_synth_mimic/traindataset.csv")

# Select wheter the analyses are done for all longitudinal variables 
# (full <- TRUE) or the ones presented in the manuscript (full <- FALSE)
full <- FALSE

# Filter data

clean_class_labels <- function(df) {
  df %>%
    # Remove leading numbers + spaces
    mutate(across(where(is.character), ~ str_replace(., "^\\d+\\s+", ""))) %>%      
    # Trim leading/trailing whitespace
    mutate(across(where(is.character), ~ str_trim(.))) %>%
    # Convert all to lowercase
    mutate(across(where(is.character), ~ str_to_lower(.)))                         
}

# Optional: fix common typos with a named vector dictionary
fix_common_typos <- function(vec) {
  typo_dict <- c(
    "abnorm extensn" = "abnormal extension",
    "abnorm flexion" = "abnormal flexion",
    "inapprop words" = "inappropriate words",
    "incomp sounds" = "incomprehensible sounds",
    "no response-ett" = "no response"
  )
  
  # Vectorized replacement with fallback to original value
  sapply(vec, function(x) {
    # If x matches a typo, replace it, else keep x as is
    if (!is.na(x) && x %in% names(typo_dict)) {
      typo_dict[[x]]
    } else {
      x
    }
  }, USE.NAMES = FALSE)
}

# The original data includes non-longitudinal data, i.e., only one observation per patient
# removing these
original <- original %>%
  filter(`Time Since Last Visit` >= 0 & `Time Since Last Visit` < max(`Time Since Last Visit`, na.rm = TRUE))

# Reformatting the Time variable as cumulative time in hours
original <- original %>%
  group_by(`Subject ID`) %>%
  filter(n() > 1) %>% 
  arrange(`Subject ID`) %>%
  mutate(Time = round(cumsum(`Time Since Last Visit`), 0)) %>%
  ungroup() %>%
  clean_class_labels() %>%
  mutate(across(where(is.character), fix_common_typos))

synthetic <- synthetic %>%
  group_by(`Subject ID`) %>%
  filter(n() > 1) %>% 
  arrange(`Subject ID`) %>%
  mutate(Time = round(cumsum(`Time Since Last Visit`), 0)) %>%
  ungroup() %>%
  clean_class_labels() %>%
  mutate(across(where(is.character), fix_common_typos))


# Compute last time for each subject
subject_end_times <- original %>%
  group_by(`Subject ID`) %>%
  summarise(last_time = max(Time), .groups = "drop")

subject_end_times <- subset(subject_end_times, last_time > 0)

proportions(table(subject_end_times$last_time <= 48))

subject_end_times <- synthetic %>%
  group_by(`Subject ID`) %>%
  summarise(last_time = max(Time), .groups = "drop")

subject_end_times <- subset(subject_end_times, last_time > 0)

proportions(table(subject_end_times$last_time <= 48))

# Subsetting data to include only observations at most 48 hours since admission
original_sub <- subset(original, Time <= 48)
synthetic_sub <- subset(synthetic, Time <= 48)

# Drop rows where all columns except Subject ID and Time are NA
original_sub <- original_sub %>%
  filter(rowSums(!is.na(original_sub[, 3:ncol(original_sub)])) > 0)

synthetic_sub <- synthetic_sub %>%
  filter(rowSums(!is.na(synthetic_sub[, 3:ncol(synthetic_sub)])) > 0)

# Get all variables except ID and time
vars <- setdiff(names(original_sub), c("Subject ID", "Time"))

# Get variable names that contain "Glas"
glas_vars <- vars[grepl("Glas", vars)]

# Convert those variables to factors in the dataframe
# NOTE: all categorical variables with multiple classes (> 2) are assumed to be set to
# factors or characters
original_sub[glas_vars] <- lapply(original_sub[glas_vars], as.factor)
synthetic_sub[glas_vars] <- lapply(synthetic_sub[glas_vars], as.factor)

# Determine min number of subjects across both datasets
n_orig <- n_distinct(original_sub$`Subject ID`)

# Match the number of patients: sample without replacement
set.seed(123)
synthetic_sub <- synthetic_sub %>%
  semi_join(
    synthetic_sub %>% distinct(`Subject ID`) %>% slice_sample(n = n_orig),
    by = "Subject ID"
  )

# Returning data depending on the settings selected in the beginning
columns_full <- colnames(original_sub[,c(43,50,16,22:27,30,35,40,44,45,48,49)])
columns_manuscript <- columns_full[c(1,2,4,5,12,13,15)]

if (full == TRUE) {
  original_data  <- original_sub[, columns_full]
  synthetic_data <- synthetic_sub[, columns_full]
} else {
  original_data  <- original_sub[, columns_manuscript]
  synthetic_data <- synthetic_sub[, columns_manuscript]
}

export <- list(
  original_data  = original_data,
  synthetic_data = synthetic_data
)

export  # <-- last evaluated expression
  

