### Helper: Gaussian Kernel ###
gaussian_kernel <- function(t, time_points, bandwidth) {
  exp(-((t - time_points)^2) / (2 * bandwidth^2)) / bandwidth
}

### Helper: Weighted Quantiles ###
weighted_quantile <- function(values, weights, probs) {
  valid <- !is.na(values) & !is.na(weights)
  values <- values[valid]
  weights <- weights[valid]
  
  if (length(values) < 2 || sum(weights) == 0) {
    return(rep(NA_real_, length(probs)))
  }
  
  ord <- order(values)
  values <- values[ord]
  weights <- weights[ord]
  cum_weights <- cumsum(weights)
  
  approx(cum_weights, values, xout = probs)$y
  
  # cum_values <- aggregate(weights, by = list(category = values), FUN = sum)
  # approx(cumsum(cum_values$x), cum_values$category, xout = quant_levels)$y
}

### Helper: Kernel-smoothed Mean, Quantiles, and Variance ###
compute_smoothed_summary <- function(data, var_name, time_col, quantiles, bandwidth, time_grid, data_type) {
  map_dfr(time_grid, function(t) {
    values <- data$value
    times <- data[[time_col]]
    
    weights <- gaussian_kernel(t, times, bandwidth)
    weights <- weights / sum(weights)
    
    smoothed_mean <- sum(values * weights, na.rm = TRUE)
    q_vals <- weighted_quantile(values, weights, quantiles)
    smoothed_variance <- sum(weights * (values - smoothed_mean)^2, na.rm = TRUE)
    
    tibble(
      time = t,
      variable = var_name,
      data_type = data_type,
      mean = smoothed_mean,
      variance = smoothed_variance,
      !!!set_names(as.list(q_vals), paste0("q", quantiles * 100))
    )
  })
}

### Helper: Check if Variable is Categorical ###
is_categorical <- function(vec) {
  is.character(vec) || is.factor(vec) || (is.numeric(vec) && length(unique(na.omit(vec))) <= 2)
}

### Helper: Kernel-smoothed Class Proportions ###
compute_smoothed_proportions <- function(data, time_col, bandwidth, time_grid, data_type) {
  map_dfr(time_grid, function(t) {
    times <- data[[time_col]]
    classes <- data$value
    
    # Kernel weights and normalization
    weights <- gaussian_kernel(t, times, bandwidth)
    weights <- weights / sum(weights)
    
    # Combine into tibble
    weighted_data <- tibble(
      class = classes,
      weight = weights
    )
    
    # Compute smoothed proportions
    weighted_data %>%
      group_by(class) %>%
      summarise(
        smoothed_proportion = sum(weight, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        time = t,
        data_type = data_type
      ) %>%
      select(time, data_type, class, smoothed_proportion)
  })
}

### Helper: subject specific resicuals ###
compute_subject_smoothed_residuals <- function(data, id_col, time_col, bandwidth, time_grid) {
  data %>%
    group_by(!!sym(id_col)) %>%
    group_map(~ {
      subject_id <- .y[[1]]
      times <- .x[[time_col]]
      values <- .x$value
      
      map_dfr(time_grid, function(t) {
        weights <- gaussian_kernel(t, times, bandwidth)
        weights <- weights / sum(weights)
        
        mu_ij <- sum(weights * values, na.rm = TRUE)
        residuals <- values - mu_ij
        weighted_squared_residuals <- sum(weights * residuals^2, na.rm = TRUE)
        
        tibble(
          time = t,
          subject_id = subject_id,
          weighted_sq_residual = weighted_squared_residuals
        )
      })
    }) %>%
    bind_rows() %>%
    group_by(time) %>%
    summarise(
      smoothed_subject_residual_variance = mean(weighted_sq_residual, na.rm = TRUE),
      .groups = "drop"
    )
}

compute_smoothed_transitions <- function(data, id_col, time_col, value_col, bandwidth, time_grid) {
  # Step 1: Standardize column names using strings
  data <- data %>%
    rename(id = !!id_col, time = !!time_col, value = !!value_col)
  
  # Step 2: Create transitions (including self-transitions)
  transitions <- data %>%
    arrange(id, time) %>%
    group_by(id) %>%
    mutate(
      from = lag(value),
      to = value,
      t_from = lag(time),
      t_to = time
    ) %>%
    filter(!is.na(from) & !is.na(to)) %>%
    ungroup()
  
  # Step 3: Compute smoothed transition probabilities over time
  map_dfr(time_grid, function(t) {
    transitions_at_t <- transitions %>%
      mutate(
        dist = pmax(abs(t - t_from), abs(t - t_to)),
        weight = gaussian_kernel(t, t - dist, bandwidth)
      ) %>%
      filter(!is.na(weight) & weight > 0)
    
    # Ensure subject was observed in both from & to states
    subject_states <- data %>%
      distinct(id, state = value)
    
    valid_transitions <- transitions_at_t %>%
      inner_join(subject_states, by = c("id", "from" = "state")) %>%
      inner_join(subject_states, by = c("id", "to" = "state"))
    
    # Normalize weights within 'from' groups
    valid_transitions <- valid_transitions %>%
      group_by(from) %>%
      mutate(weight = weight / sum(weight, na.rm = TRUE)) %>%
      ungroup()
    
    # Aggregate to get smoothed probabilities
    valid_transitions %>%
      group_by(from, to) %>%
      summarise(prob = sum(weight, na.rm = TRUE), .groups = "drop") %>%
      mutate(time = t)
  })
}

### Helper: comuptes the residuals ###
compute_residuals <- function(data, var_name, id_col, time_col, bandwidth, time_grid, data_type) {
  message("  [*] Computing residuals for ", data_type, " data")
  
  map_dfr(time_grid, function(t) {
    all_times <- data[[time_col]]
    all_values <- data$value
    
    weights <- gaussian_kernel(t, all_times, bandwidth)
    weights <- weights / sum(weights)
    
    smoothed_mean <- sum(all_values * weights, na.rm = TRUE)
    
    exact_data <- data[data[[time_col]] == t, ]
    residuals <- exact_data$value - smoothed_mean
    
    weighted_variance <- sum(weights * (all_values - smoothed_mean)^2, na.rm = TRUE)
    
    tibble(
      time = t,
      variable = var_name,
      data_type = data_type,
      id = exact_data[[id_col]],
      value = exact_data$value,
      mean = smoothed_mean,
      residual = residuals,
      variance = weighted_variance
    )
  })
}

### Helper: Data Processing ###
process_single_dataset <- function(data, var_name, data_type, id_col, time_col) {
  data_long <- data %>%
    select(all_of(id_col), all_of(time_col), all_of(var_name)) %>%
    rename(value = all_of(var_name)) %>%
    mutate(variable = var_name, data_type = data_type) %>%
    na.omit()
  
  return(data_long)
}


### Helper: Check if Variable is Categorical ###
is_categorical <- function(vec) {
  is.character(vec) || is.factor(vec) || (is.numeric(vec) && length(unique(na.omit(vec))) <= 2)
}


### Helper: computes the variogram ###
compute_variogram_fast <- function(data, u_grid, bandwidth, data_type, variable) {
  message("  [*] Computing variogram for ", data_type)
  
  dt <- as.data.table(data)
  setkey(dt, id, time)
  
  # Preallocate list
  pairwise_list <- vector("list", length(unique(dt$id)))
  subject_ids <- unique(dt$id)
  
  for (i in seq_along(subject_ids)) {
    subj <- subject_ids[i]
    subject_data <- dt[id == subj]
    if (nrow(subject_data) < 2) next
    
    comb <- combn(nrow(subject_data), 2)
    delta_t <- abs(subject_data$time[comb[1,]] - subject_data$time[comb[2,]])
    delta_r2 <- (subject_data$residual[comb[1,]] - subject_data$residual[comb[2,]])^2
    pairwise_list[[i]] <- data.table(delta_t, delta_r2)
    
    if (i %% 500 == 0) message("    Processed ", i, "/", length(subject_ids), " subjects")
  }
  
  # Combine all
  all_pairs <- rbindlist(pairwise_list, use.names = TRUE, fill = TRUE)
  
  gamma_vals <- sapply(u_grid, function(u) {
    weights <- gaussian_kernel(u, all_pairs$delta_t, bandwidth)
    weights <- weights / sum(weights)
    0.5 * sum(weights * all_pairs$delta_r2, na.rm = TRUE)
  })
  
  tibble(
    time = u_grid,
    value = gamma_vals,
    data_type = data_type,
    statistic = "Variogram",
    variable = variable
  )
}


### Helper: get legend ###
get_legend <- function(myggplot) {
  gtable::gtable_filter(ggplotGrob(myggplot), "guide-box")
}

### Helper: rank-order stability ###
compute_smoothed_ranks <- function(data, var_name, time_col, quantiles, bandwidth, time_grid, data_type) {
  bind_rows(lapply(time_grid, function(t) {
    values <- data$value
    times <- data[[time_col]]
    
    # Kernel weights centered at current time grid point
    weights <- gaussian_kernel(t, times, bandwidth)
    weights <- weights / sum(weights)
    
    # Compute smoothed quantiles
    q_vals <- weighted_quantile(values, weights, quantiles)
    breaks <- sort(c(-Inf, q_vals, Inf))
    
    # For each observation, check if it's observed at this time point
    data_at_t <- data %>% filter(!!sym(time_col) == t)
    
    if (nrow(data_at_t) == 0) return(NULL)
    
    # Assign quantile ranks
    data_at_t <- data_at_t %>%
      mutate(
        smoothed_rank = findInterval(value, vec = breaks, rightmost.closed = TRUE),
        time = t,
        variable = var_name,
        data_type = data_type
      )
    
    return(data_at_t)
  }))
}

### Helper: Rank-order variability ###
bounded_rank_variability <- function(ranks, max_change = 5) {
  mean(diff(ranks)) / max_change
}


### Helper: creates measurement matrix ###
create_measurement_matrix_fast <- function(df, time_grid) {
  colnames(df) <- c("id", "time", "value")

  # Save all unique IDs in the full df (not filtered)
  all_ids <- unique(df$id)
  
  # Filter non-NA values only for matrix creation
  df_non_na <- df[!is.na(df$value), ]
  
  # Use factor with levels including all subjects to keep all rows
  df_non_na$id <- factor(df_non_na$id, levels = all_ids)
  
  # Force time column to match the fixed grid
  df_non_na$time <- factor(df_non_na$time, levels = time_grid)
  
  # Create table
  M <- with(df_non_na, table(id, time))
  M <- (M > 0) * 1L  # binarize
  
  # Ensure all time_grid columns are present
  M <- M[, as.character(time_grid), drop = FALSE]
  
  return(M)
}


### Helper: calculates the measurement similarity score ###
measurement_similarity <- function(M_orig, M_synth) {
  # Ensure same shape
  if (!all(dim(M_orig) == dim(M_synth))) {
    stop("Original and synthetic matrices must have the same dimensions")
  }
  hamming_dist <- sum(abs(M_orig - M_synth))
  similarity <- 1 - (hamming_dist / (nrow(M_orig) * ncol(M_orig)))
  return(similarity)
}

### Helper: computes dropout points ###
compute_dropout_points <- function(M) {
  apply(M, 1, function(row) {
    for (t in seq_along(row)) {
      if (all(row[t:length(row)] == 0)) {
        return(t)
      }
    }
    return(length(row) + 1)  # never dropped out (after last time)
  })
}

### Helper: calculates the dropout distribution ###
compute_dropout_distribution <- function(dropout_points, support = NULL) {
  tab <- table(factor(dropout_points, levels = support))
  dist <- as.numeric(tab) / sum(tab)
  names(dist) <- names(tab)
  return(dist)
}

### Helper: KL-divergence ###
kl_divergence <- function(P, Q, epsilon = 1e-8) {
  P <- P + epsilon
  Q <- Q + epsilon
  P <- P / sum(P)
  Q <- Q / sum(Q)
  sum(P * log(P / Q))
}

### Helper: calculates the dropout divergence ###
dropout_divergence <- function(M_orig, M_synth) {
  d_orig <- compute_dropout_points(M_orig)
  d_synth <- compute_dropout_points(M_synth)
  
  # Define common support
  all_points <- sort(unique(c(d_orig, d_synth)))
  
  P <- compute_dropout_distribution(d_orig, support = all_points)
  Q <- compute_dropout_distribution(d_synth, support = all_points)
  
  divergence <- kl_divergence(P, Q)
  return(divergence)
}

### Helper: individual trajectories ###
plot_subject_specific_baselines <- function(rank_data, n_per_group = 10, n_groups = 6, 
                                            max_time, x_interval) {
  
  # Step 1: Find baseline time per subject
  subject_baseline <- rank_data %>%
    arrange(data_type, `Subject ID`, time) %>%
    group_by(data_type, `Subject ID`) %>%
    slice(1) %>%  # Select the first observed row (baseline)
    ungroup() %>%
    select(`Subject ID`, data_type, baseline_rank = smoothed_rank)
  
  # Step 2: Sample subjects per rank group *and* data_type
  sampled_ids <- subject_baseline %>%
    group_by(data_type, baseline_rank) %>%
    sample_n(min(n(), n_per_group), replace = FALSE) %>%
    ungroup()
  
  # Step 4: Subset full data for sampled subjects
  plot_data <- rank_data %>%
    semi_join(sampled_ids, by = c("Subject ID", "data_type")) %>%
    left_join(sampled_ids, by = c("Subject ID", "data_type")) %>%
    mutate(
      baseline_rank = factor(baseline_rank, labels = paste("Group", 1:n_groups)),
      subject_id = factor(`Subject ID`),
      data_type = factor(data_type, levels = c("Original", "Synthetic"))  # Ensure consistent facet order
    )
  
  # Generate a distinct color for each subject per group
  plot_data <- plot_data %>%
    group_by(data_type, baseline_rank) %>%
    mutate(subject_plot_id = factor(as.numeric(as.factor(`Subject ID`)))) %>%
    ungroup()
  
  # Step 5: Faceted plot by rank group (columns) and data type (rows)
  
  ggplot(plot_data, aes(x = time, y = value, group = `Subject ID`, color = subject_plot_id)) +
    geom_line(size = 0.9, alpha = 0.85) +
    geom_point(size = 0.7, alpha = 0.85) +
    facet_grid(data_type ~ baseline_rank) +
    labs(
      x = "Time", y = "Value", color = "Subject"
    ) +
    theme_minimal() +
    scale_color_manual(values = hue_pal()(n_per_group)) +  # Only need 20 colors
    scale_x_continuous(breaks = seq(0, max_time, by = x_interval), limits = c(0,max_time)) +
    theme(
      strip.text = element_text(size = 12),
      legend.position = "none",
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12)
    )
  
}


### Helper: Hungarian algortihm ###
# Finds the row permutation P which minimizes the norm of X - Y[P, ]
# X = n x p matrix
# Y = n x p matrix
opt_perm <- function(X, Y){
  HungarianSolver(-1*X%*%t(Y))$pairs[, 2]
}

### Helper: easurement density relative to total measurements ###
compute_relative_density <- function(M) {
  time_counts <- colSums(M, na.rm = TRUE)
  total_counts <- sum(M, na.rm = TRUE)
  time_counts / total_counts
}

# Label definitions
eye_labels <- c(
  "no response"   = "No eye opening",
  "spontaneously" = "Spontaneous eye opening",
  "to speech"     = "Eye opening to verbal command",
  "to pain"       = "Eye opening to pain"
)

verbal_labels <- c(
  "oriented"               = "Alert and oriented",
  "confused"               = "Confused",
  "inappropriate words"    = "Inappropriate words",
  "incomprehensible sounds"= "Incomprehensible sounds",
  "no response"            = "No verbal response"
)

motor_labels <- c(
  "obeys commands"     = "Obeys commands",
  "localizes pain"     = "Localizes pain",
  "flex-withdraws"     = "Withdraws from pain",
  "abnormal flexion"   = "Abnormal flexion",
  "abnormal extension" = "Extension",
  "no response"        = "No motor response"
)

label_lookup <- list(
  "Glascow coma scale eye opening"     = eye_labels,
  "Glascow coma scale verbal response" = verbal_labels,
  "Glascow coma scale motor response"  = motor_labels
)

