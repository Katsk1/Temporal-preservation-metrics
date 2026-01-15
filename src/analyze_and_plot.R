### Summary statistics plots ###

# NOTE: categorical variables with more than two classes are assumed to be
# set as either factors or characters, binary variables can be left as numeric

analyze_and_plot_data <- function(original_data,
                                  synthetic_data,
                                  id_col = "id", # column with the subject ids
                                  time_col = "time", # column indicating the time
                                  max_time = numeric(), # maximum time after which observations are cut out, 
                                                        # if left empty then the maxium from the data is used
                                  time_unit = 1, # the time unit, i.e., precision, default 1 (hour), assumed to be same for all variables
                                  x_interval = 24, # grid space for time axis when plotting, default 24 (hours)
                                  quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9),
                                  u_grid = NULL, # lags for variogram
                                  bandwidth = 1, # bandwidth to kernel smoothing
                                  n_per_group = 20, # sample size for individual trajectories
                                  meas_iter = 10, # number of iterations to calculate measurement similarity matrices
                                  meas_n = 1000, # sample size to subset measurement similarity matrices to perform Hungarian algorithm
                                  save_path = NULL, # path to save plots, if null printed out
                                  mean_quantile = TRUE,  # form mean-quantile plots
                                  cat_prop = TRUE, # form proportions plots
                                  rank_stab = TRUE, # form rank-order plots
                                  meas_sim = TRUE) { # form measurement similarity plots
  
  
  # Get all variables except ID and time
  vars <- setdiff(names(original_data), c(id_col, time_col))
  
  # Base colors used in plotting
  base_colors <- c("#1f77b4",  # Blue
                   "#ff7f0e",  # Orange
                   "#2ca02c",  # Green
                   "#9467bd",  # Purple
                   "#17becf", # Teal
                   "#d62728",  # Red
                   "#8c564b",  # Brown
                   "#e377c2",  # Pink
                   "#7f7f7f",  # Gray
                   "#bcbd22",  # Lime
                   "#aec7e8",  # Light Blue
                   "#ffbb78",  # Light Orange
                   "#98df8a"   # Light Green
  )
  
  for (var_name in vars) {
    
    save_path_var <- paste0(save_path,"/",gsub(" ", "_", tolower(var_name)), "/")
    dir.create(save_path_var, recursive = TRUE, showWarnings = FALSE)
    
    # Check if the variable is numeric or categorical
    is_cat <- is_categorical(original_data[[var_name]])
    
    # Process the original and synthetic data
    orig_data <- process_single_dataset(original_data, var_name, "Original",
                                        id_col, time_col)
    synth_data <- process_single_dataset(synthetic_data, var_name, "Synthetic",
                                         id_col, time_col)
    
    # Combine original and synthetic datasets
    combined_data <- bind_rows(orig_data, synth_data)
    
    # Get min and max times across the combined data
    time_vals <- combined_data[[time_col]]
    time_grid <- seq(min(time_vals, na.rm = TRUE), max(time_vals, na.rm = TRUE), by = time_unit)
    
    # Set max time for plotting
    max_time <- ifelse(length(max_time) == 0, max(time_grid), max_time)
    
    ### Smoothing proportions for categorical variables ###
    if (is_cat) {
      
      if (cat_prop) {
        
        smoothed_proportions <- bind_rows(
          compute_smoothed_proportions(orig_data, time_col, bandwidth, time_grid, "Original"),
          compute_smoothed_proportions(synth_data, time_col, bandwidth, time_grid, "Synthetic")
        )
        
        # Recode class labels
        if(var_name != "Glascow coma scale total"){
          
          smoothed_proportions <- smoothed_proportions %>%
            rowwise() %>%
            mutate(
              class = recode(class, !!!label_lookup[[var_name]])
            ) %>%
            ungroup()
        }
        
        class_levels <- unique(smoothed_proportions$class)
        
        # Compute last time for each subject
        subject_end_times <- combined_data %>%
          group_by(`Subject ID`, data_type) %>%
          summarise(last_time = max(Time), .groups = "drop")
        
        # Create a subject × time grid with data_type retained
        subject_time_grid <- subject_end_times %>%
          select(`Subject ID`, data_type, last_time) %>%
          crossing(time = sort(unique(combined_data$Time))) %>%
          filter(time <= last_time)
        
        # Count subjects at each time and data_type
        subjects_remaining <- subject_time_grid %>%
          group_by(time, data_type) %>%
          summarise(subjects_in = n_distinct(`Subject ID`), .groups = "drop")
        
        # Normalize within each data_type
        followup_curve <- subjects_remaining %>%
          group_by(data_type) %>%
          mutate(prop_remaining = subjects_in / max(subjects_in)) %>%
          ungroup()
        
        # Global color mapping
        color_map <- setNames(base_colors[1:length(class_levels)], class_levels)
        
        # Proportions plot
        proportions_plot <- ggplot(smoothed_proportions, aes(x = time, y = smoothed_proportion, fill = class)) +
          geom_area(position = "stack", alpha = 0.8) +
          scale_fill_manual(values = color_map, drop = FALSE) +
          facet_wrap(~data_type) +
          labs(
            title = NULL,
            x = NULL, y = "Proportion", fill = "Class"
          ) +
          scale_x_continuous(breaks = seq(0, max_time, by = x_interval), 
                            limits = c(0, max_time), expand = c(0, 0)) +
          scale_y_continuous(expand = c(0,0)) +
          theme_minimal() +
          theme_minimal() +
          theme(
            axis.title.y = element_text(margin = margin(r = 15), size = 12),
            legend.position = "top",
            legend.box = "horizontal",
            legend.direction = "horizontal",
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 12),
            plot.margin = margin(5, 10, 5, 10),
            strip.text = element_text(size = 12),
            panel.spacing = unit(1, "lines")
          ) +
          guides(fill = guide_legend(nrow = 1))
        
        legend <- get_legend( proportions_plot)
        plot_nolegend <-  proportions_plot + theme(legend.position = "none",
                                                    plot.margin = margin(5, 10, 15, 10))
        
        x_breaks <- seq(0, max_time, by = x_interval)
        x_limits <- range(x_breaks)

        # Follow-up plot
        # Follow-up plot
        followup_plot <- ggplot(followup_curve, aes(x = time, y = prop_remaining, color = data_type)) +
          geom_line(size = 1) +
          scale_x_continuous(breaks = x_breaks, limits = x_limits, expand = c(0, 0)) +
          scale_y_continuous(
            labels = scales::percent_format(accuracy = 1),
            breaks = function(lims) {
              span <- diff(lims)
              step <- dplyr::case_when(
                span <= 0.25 ~ 0.05,  # 5% if very tight range
                span <= 0.5  ~ 0.10,  # 10%
                span <= 0.75 ~ 0.20,  # 20%
                TRUE         ~ 0.25   # 25% for wide ranges
              )
              lo <- floor(lims[1] / step) * step
              hi <- ceiling(lims[2] / step) * step
              seq(lo, hi, by = step)
            },
            expand = expansion(mult = c(0, 0.02))
          ) +
          scale_color_manual(values = c("Original" = base_colors[1], "Synthetic" = base_colors[2])) +
          labs(y = "Subjects in follow-up (%)", x = "Time (hours)") +
          facet_wrap(~ data_type, ncol = 2, scales = "free_y") +
          theme_minimal() +
          theme(
            axis.title.y = element_text(margin = margin(r = 15), size = 12),
            axis.title.x = element_text(size = 12),
            strip.text = element_blank(),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 12),
            legend.position = "none",
            plot.margin = margin(0, 9, 0, 9),
            panel.spacing = unit(1, "lines")
          )
        
        # Convert to grobs
        base_grob <- ggplotGrob(plot_nolegend)
        followup_grob <- ggplotGrob(followup_plot)
        
        # Remove right y-axis from follow-up plot (second facet)
        axis_l_indices <- which(grepl("axis-l", followup_grob$layout$name))
        right_y_axis_index <- axis_l_indices[1] # second facet's left axis
        
        followup_grob$grobs[[right_y_axis_index]] <- nullGrob()
        followup_yaxis_pos <- followup_grob$layout[followup_grob$layout$name == followup_grob$layout$name[right_y_axis_index], ]
        followup_grob$widths[followup_yaxis_pos$l] <- unit(0, "cm")
        
        # Find panel columns in base plot layout (where panels sit)
        panel_cols <- base_grob$layout[base_grob$layout$name == "panel-1-1", "l"]
        
        # All widths to left of panel are axis + labels
        yaxis_widths_base <- base_grob$widths[1:(panel_cols - 1)]
        
        # Similarly for follow-up plot, find the panel column
        panel_cols_followup <- followup_grob$layout[followup_grob$layout$name == "panel-1-1", "l"]
        
        # Replace the widths left of panel in follow-up with those from base plot
        followup_grob$widths[1:(panel_cols_followup - 1)] <- yaxis_widths_base
        
        # Align x-axis tick label widths between base and follow-up
        
        base_x_axis_indices <- which(base_grob$layout$name %in% c("axis-b-1-1", "axis-b-2-1"))
        followup_x_axis_indices <- which(followup_grob$layout$name %in% c("axis-b-1-1", "axis-b-2-1"))
        
        for (i in seq_along(base_x_axis_indices)) {
          base_axis_grob <- base_grob$grobs[[base_x_axis_indices[i]]]
          followup_axis_grob <- followup_grob$grobs[[followup_x_axis_indices[i]]]
          
          followup_axis_grob$widths <- base_axis_grob$widths
          
          followup_grob$grobs[[followup_x_axis_indices[i]]] <- followup_axis_grob
        }
        
        # Combine the plots using patchwork
        final_proportions_plot <- proportions_plot / followup_grob  +
          plot_layout(heights = c(10, 2.5)) +
          plot_annotation(theme = theme(plot.margin = margin(10, 10, 10, 10)))
        
        # Save or display
        if (!is.null(save_path)) {
          ggsave(file.path(save_path_var, paste0("proportions_", var_name, ".pdf")), 
                 final_proportions_plot, width = 16.54, height = 8.8, units = "in", dpi = 600)
        } else {
          print(proportions_plot)
        }
        
        smooth_trans_orig <- compute_smoothed_transitions(orig_data, 
                                                          id_col = id_col, 
                                                          time_col = time_col, 
                                                          value_col = "value", 
                                                          bandwidth = bandwidth, 
                                                          time_grid = time_grid)
        smooth_trans_orig$data_type <- "Original"
        smooth_trans_synth <- compute_smoothed_transitions(synth_data, 
                                                           id_col = id_col, 
                                                           time_col = time_col, 
                                                           value_col = "value", 
                                                           bandwidth = bandwidth, 
                                                           time_grid = time_grid)
        smooth_trans_synth$data_type <- "Synthetic"
        smooth_trans_all <- bind_rows(smooth_trans_orig, smooth_trans_synth)
        
        gc()
        
        # Recode 'from' and 'to'
        if(var_name != "Glascow coma scale total"){
          smooth_trans_all <- smooth_trans_all %>%
            rowwise() %>%
            mutate(
              from = recode(from, !!!label_lookup[[var_name]]),
              to   = recode(to,   !!!label_lookup[[var_name]])
            ) %>%
            ungroup()
        }
        
        # Global list of all possible `to` states
        all_tos <- sort(unique(c(smooth_trans_all$from, smooth_trans_all$to)))
        
        # Global color mapping
        color_map <- setNames(base_colors[1:length(all_tos)], all_tos)
        
        # Set correct types
        smooth_trans_all <- smooth_trans_all %>%
          mutate(
            from = factor(from, levels = all_tos),
            to = factor(to, levels = all_tos),  # <- crucial
            data_type = factor(data_type, levels = c("Original", "Synthetic"))
          )
        
        # Calculate start times and probs
        start_times <- smooth_trans_all %>%
          group_by(from, data_type) %>%
          summarise(start_time = min(time), .groups = "drop")
        
        start_probs <- smooth_trans_all %>%
          inner_join(start_times, by = c("from", "data_type")) %>%
          filter(time == start_time)
        
        # Calculate total start probs by 'from' and 'to' across data_types to get ordering
        start_ordering <- start_probs %>%
          group_by(from, to) %>%
          summarise(start_prob_sum = sum(prob), .groups = "drop") %>%
          arrange(from, desc(start_prob_sum)) %>%
          group_by(from) %>%
          summarise(to_order = list(to), .groups = "drop")  # to_order is ordered vector of to levels for each from
        
        # Create all unique time points from the plot
        all_times <- sort(unique(smooth_trans_all$time))
        
        gc()
        
        # Loop over each 'from'
        from_states <- unique(smooth_trans_all$from)
        
        walk(from_states, function(f) {
          df_f <- filter(smooth_trans_all, from == f)
          
          # Extract ordering for this 'from'
          to_levels_ordered <- start_ordering %>% 
            filter(from == f) %>% pull(to_order) %>% .[[1]]
          
          # Reorder factor levels of 'to' so highest start prob is first (bottom of stack)
          df_f <- df_f %>%
            mutate(to = factor(to, levels = rev(to_levels_ordered)))
          
          # Base probability plot
          trans_plot <- ggplot(df_f, aes(x = time, y = prob, fill = to)) +
            geom_area(position = "stack", alpha = 0.8) +
            facet_wrap(~ data_type, ncol = 2) +
            labs(y = "Probability", x = "Time (hours)", 
                 fill = paste0("From '", f, "' to")) +
            scale_fill_manual(
              values = color_map,
              drop = FALSE,
              breaks = levels(smoothed_proportions$class)) +
            scale_x_continuous(breaks = seq(0, max_time, by = x_interval), 
                               limits = c(0, max_time), expand = c(0, 0)) +
            scale_y_continuous(expand = c(0,0)) +
            theme_minimal() +
            theme(
              axis.title.y = element_text(margin = margin(r = 15), size = 12),
              axis.title.x = element_text(size = 12),
              legend.position = "top",
              legend.box = "horizontal",
              legend.direction = "horizontal",
              strip.text = element_text(size = 12),
              legend.title = element_text(size = 12),
              legend.text = element_text(size = 12),
              plot.margin = margin(5, 10, 5, 10),
              panel.spacing = unit(1, "lines")
            ) +
            guides(fill = guide_legend(nrow = 1))
          
          if (!is.null(save_path)) {
            safe_name <- gsub("[^A-Za-z0-9]", "_", f)
            ggsave(file.path(save_path_var, paste0("trans_plot__", var_name, "_from_", safe_name, ".pdf")),
                   trans_plot, width = 16.54, height = 8.8, units = "in", dpi = 600)
          } else {
            print(trans_plot)
          }
        })
        
        gc()
        
      }
      
    } else {
      
      # Compute last time for each subject
      subject_end_times <- combined_data %>%
        group_by(`Subject ID`, data_type) %>%
        summarise(last_time = max(Time), .groups = "drop")
      
      # Create a subject × time grid with data_type retained
      subject_time_grid <- subject_end_times %>%
        select(`Subject ID`, data_type, last_time) %>%
        crossing(time = sort(unique(combined_data$Time))) %>%
        filter(time <= last_time)
      
      # Count subjects at each time and data_type
      subjects_remaining <- subject_time_grid %>%
        group_by(time, data_type) %>%
        summarise(subjects_in = n_distinct(`Subject ID`), .groups = "drop")
      
      # Normalize within each data_type
      followup_curve <- subjects_remaining %>%
        group_by(data_type) %>%
        mutate(prop_remaining = subjects_in / max(subjects_in)) %>%
        ungroup()
      
      ### Metrics ###
      
      if (mean_quantile){
      # Mean, quantiles and between-individual variance
      quantile_results <- bind_rows(
        compute_smoothed_summary(orig_data, var_name, time_col, 
                                 quantiles, bandwidth, time_grid, "Original"),
        compute_smoothed_summary(synth_data, var_name, time_col, quantiles, 
                                 bandwidth, time_grid, "Synthetic")
      )
      
      ### Plots ###
      
      # Define the color map dynamically
      quantile_labels <- paste0("q", quantiles * 100)
      all_stats <- c("mean", quantile_labels)
      
      # Get min and max quantile column names dynamically
      min_quant_label <- quantile_labels[which.min(quantiles)]
      max_quant_label <- quantile_labels[which.max(quantiles)]
      
      # Extract bounds from quantile_results
      quantile_bounds <- quantile_results %>%
        select(time, data_type, !!min_quant_label, !!max_quant_label) %>%
        rename(lower = !!min_quant_label, upper = !!max_quant_label)
      
      # Join bounds to combined_data and flag outliers
      combined_with_bounds <- combined_data %>%
        rename(time = Time) %>%  # Rename for join
        left_join(quantile_bounds, by = c("time", "data_type")) %>%
        mutate(outside = value < lower | value > upper)
      
      
      # Long-format for plotting
      plot_data_long <- quantile_results %>%
        pivot_longer(cols = all_of(all_stats), names_to = "statistic", values_to = "value")
      
      plot_data_long <- plot_data_long %>%
        mutate(statistic = case_when(
          statistic == "mean" ~ "Mean",
          grepl("^q\\d{2}$", statistic) ~ toupper(statistic),
          TRUE ~ str_to_title(statistic)
        ))
      
      stat_levels <- c("Mean","Q95","Q75","Q50","Q25","Q5")  # adjust based on what's in your data
      plot_data_long <- plot_data_long %>%
        mutate(statistic = factor(statistic, levels = stat_levels))
      
      # Color assignment:
      base_colors_sub <- base_colors[1:(length(quantiles)+1)]
      names(base_colors_sub ) <- str_to_title(c("mean", "q50", setdiff(quantile_labels, "q50")))
      
      # Without outliers
      y_breaks <- scales::pretty_breaks(n = 6)(range(plot_data_long$value))
      y_limits <- range(plot_data_long$value)
      
      # Mean and Quantiles Plot
      mean_quantile_plot <- ggplot() +
        geom_line(data = plot_data_long,
                  aes(x = time, y = value, color = statistic),
                  size = 1.2) +
        facet_wrap(~data_type) +
        scale_color_manual(values = base_colors_sub[names(base_colors_sub) %in% unique(plot_data_long$statistic)]) +
        labs(
          title = NULL, x = "Time", y = var_name, color = "Statistic"
        ) +
        scale_x_continuous(breaks = seq(0, max_time, by = x_interval), 
                           limits = c(0, max_time), expand = c(0, 0)) +
        scale_y_continuous(breaks = y_breaks, limits = y_limits,
                           expand = expansion(mult = c(0, 0.07))) +
        theme_minimal()+
        theme(
          axis.title.x = element_blank(),
          axis.title.y = element_text(margin = margin(r = 15), size = 12),
          strip.text = element_text(size = 12),
          legend.position = "top",
          legend.box = "horizontal",
          legend.direction = "horizontal",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          plot.margin = margin(5, 10, 5, 10),
          panel.spacing = unit(1, "lines")
        ) +
        guides(color = guide_legend(nrow = 1))
      
      legend <- get_legend(mean_quantile_plot)
      plot_nolegend <- mean_quantile_plot + theme(legend.position = "none",
                                         plot.margin = margin(5, 10, 15, 10))
      
      x_breaks <- seq(0, max_time, by = x_interval)
      x_limits <- range(plot_data_long$time)
      
      # Follow-up plot
      followup_plot <- ggplot(followup_curve, aes(x = time, y = prop_remaining, color = data_type)) +
        geom_line(size = 1) +
        scale_x_continuous(breaks = x_breaks, limits = x_limits, expand = c(0,0)) +
        scale_y_continuous(
          labels = scales::percent_format(accuracy = 1),
          breaks = function(lims) {
            span <- diff(lims)
            step <- dplyr::case_when(
              span <= 0.25 ~ 0.05,  # 5% if very tight range
              span <= 0.5  ~ 0.10,  # 10%
              span <= 0.75 ~ 0.20,  # 20%
              TRUE         ~ 0.25   # 25% for wide ranges
            )
            lo <- floor(lims[1] / step) * step
            hi <- ceiling(lims[2] / step) * step
            seq(lo, hi, by = step)
          },
          expand = expansion(mult = c(0, 0.02))
        ) +
        scale_color_manual(values = c("Original" = base_colors[1], "Synthetic" = base_colors[2])) +
        labs(y = "Subjects in follow-up (%)", x = "Time (hours)") +
        facet_wrap(~ data_type, ncol = 2, scales = "free_y") +
        theme_minimal() +
        theme(
          axis.title.y = element_text(margin = margin(r = 15), size = 12),
          axis.title.x = element_text(size = 12),
          strip.text = element_blank(),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.position = "none",
          plot.margin = margin(0, 9, 5, 9),
          panel.spacing = unit(1, "lines")
        )
      
      # Convert to grobs
      base_grob <- ggplotGrob(plot_nolegend)
      followup_grob <- ggplotGrob(followup_plot)
      
      # Remove right y-axis from follow-up plot (second facet)
      axis_l_indices <- which(grepl("axis-l", followup_grob$layout$name))
      right_y_axis_index <- axis_l_indices[1] # second facet's left axis
      
      followup_grob$grobs[[right_y_axis_index]] <- nullGrob()
      followup_yaxis_pos <- followup_grob$layout[followup_grob$layout$name == followup_grob$layout$name[right_y_axis_index], ]
      followup_grob$widths[followup_yaxis_pos$l] <- unit(0, "cm")
      
      # Find panel columns in base plot layout (where panels sit)
      panel_cols <- base_grob$layout[base_grob$layout$name == "panel-1-1", "l"]
      
      # All widths to left of panel are axis + labels
      yaxis_widths_base <- base_grob$widths[1:(panel_cols - 1)]
      
      # Similarly for follow-up plot, find the panel column
      panel_cols_followup <- followup_grob$layout[followup_grob$layout$name == "panel-1-1", "l"]
      
      # Replace the widths left of panel in follow-up with those from base plot
      followup_grob$widths[1:(panel_cols_followup - 1)] <- yaxis_widths_base
      
      # Align x-axis tick label widths between base and follow-up
      
      base_x_axis_indices <- which(base_grob$layout$name %in% c("axis-b-1-1", "axis-b-2-1"))
      followup_x_axis_indices <- which(followup_grob$layout$name %in% c("axis-b-1-1", "axis-b-2-1"))
      
      for (i in seq_along(base_x_axis_indices)) {
        base_axis_grob <- base_grob$grobs[[base_x_axis_indices[i]]]
        followup_axis_grob <- followup_grob$grobs[[followup_x_axis_indices[i]]]
        
        followup_axis_grob$widths <- base_axis_grob$widths
        
        followup_grob$grobs[[followup_x_axis_indices[i]]] <- followup_axis_grob
      }
      
      # Combine the plots using patchwork
      final_mean_quantile_plot <- mean_quantile_plot / followup_grob  +
        plot_layout(heights = c(10, 2.5)) +
        plot_annotation(theme = theme(plot.margin = margin(10, 10, 10, 10)))
      
      # With outliers
      y_breaks <- scales::pretty_breaks(n = 6)(range(combined_with_bounds$value))
      y_limits <- range(combined_with_bounds$value)
      
      mean_quantile_plot_out <- ggplot() +
        # Gray dots for values outside quantiles
        geom_point(data = filter(combined_with_bounds, outside),
                   aes(x = time, y = value),
                   color = "gray80", alpha = 0.6, size = 0.6) +
        geom_line(data = plot_data_long,
                  aes(x = time, y = value, color = statistic),
                  size = 1.2) +
        facet_wrap(~data_type) +
        scale_color_manual(values = base_colors_sub[names(base_colors_sub) %in% unique(plot_data_long$statistic)]) +
        labs(
          title = NULL, x = "Time (hours)", y = var_name, color = "Statistic"
        ) +
        scale_x_continuous(breaks = seq(0, max_time, by = x_interval), 
                           limits = c(0, max_time), expand = c(0, 0)) +
        scale_y_continuous(breaks = y_breaks, limits = y_limits,
                           expand = expansion(mult = c(0, 0.07))) +
        theme_minimal()+
        theme(
          axis.title.x = element_blank(),
          axis.title.y = element_text(margin = margin(r = 15), size = 12),
          strip.text = element_text(size = 12),
          legend.position = "top",
          legend.box = "horizontal",
          legend.direction = "horizontal",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          plot.margin = margin(5, 10, 5, 10),
          panel.spacing = unit(1, "lines")
        ) +
        guides(color = guide_legend(nrow = 1))
      
      # Residuals
      residual_orig <- compute_residuals(orig_data, var_name, id_col, time_col, bandwidth , time_grid, "Original")
      residual_synth <- compute_residuals(synth_data, var_name, id_col, time_col, bandwidth, time_grid, "Synthetic")
      residual_combined <- bind_rows(residual_orig, residual_synth)
      gc()
      
      # Variance plot
      variance_df <- residual_combined %>%
        select(time, variable, data_type, variance) %>%
        distinct() %>%
        rename(value = variance) %>%
        mutate(statistic = "Variance")
      
      # Variogram
      
      if(is.null(u_grid)){
        u_grid <- time_grid
      }
      
      variogram_orig <- compute_variogram_fast(residual_orig, u_grid, bandwidth, "Original", var_name)
      variogram_synth <- compute_variogram_fast(residual_synth, u_grid, bandwidth, "Synthetic", var_name)
      variogram_df <- bind_rows(variogram_orig, variogram_synth)
      
      # Combine for plotting
      var_vario_df <- bind_rows(variance_df, variogram_df)
      y_breaks <- scales::pretty_breaks(n = 8)(0:ceiling(max(var_vario_df$value)))
      y_limits <- range(y_breaks)
      
      vario_plot <- ggplot(var_vario_df, aes(x = time, y = value, color = statistic)) +
        geom_line(linewidth = 1.2) +
        facet_wrap(~data_type) +
        labs(
          title = NULL, x = "Time (hours) / lag", y = "Value", color = "Statistic"
        ) +
        scale_x_continuous(breaks = seq(0, max_time, by = x_interval), limits = c(0, max_time)) +
        scale_y_continuous(breaks = y_breaks, limits = y_limits, expand = expansion(mult = c(0, 0.07))) +
        scale_color_manual(values = base_colors[1:2]) +
        theme_minimal() +
        theme(
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(margin = margin(r = 15), size = 12),
          strip.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.position = "top",
          plot.margin = margin(5, 10, 5, 10),
          panel.spacing = unit(1, "lines")
        ) +
        guides(color = guide_legend(nrow = 1))
      
        rm(residual_orig, residual_synth, residual_combined, variance_df, variogram_df, var_vario_df)

        # Save or display
        if (!is.null(save_path)) {
          ggsave(file.path(save_path_var, paste0("mean_quantiles_", tolower(var_name), ".pdf")), 
                 final_mean_quantile_plot, width = 8.27*1.5, height = 4.4*1.6, units = "in", dpi = 600)
          ggsave(file.path(save_path_var, paste0("mean_quantiles_", tolower(var_name), "_outliers.pdf")), 
                 mean_quantile_plot_out, width = 8.27*1.5, height = 4.4*1.2, units = "in", dpi = 600)
          ggsave(file.path(save_path_var, paste0("variogram_", tolower(var_name), ".pdf")), 
                 plot = vario_plot, width = 8.27*1.5, height = 4.4*1.2, units = "in", dpi = 600)
        } else {
          print(final_mean_quantile_plot)
          print(mean_quantile_plot_out)
          print(vario_plot)
        }
      }
    
      if (rank_stab){
        
        rank_orig <- compute_smoothed_ranks(orig_data, var_name, time_col, 
                                            quantiles, bandwidth, time_grid, "Original")
        
        rank_synth <- compute_smoothed_ranks(synth_data, var_name, time_col, 
                                             quantiles, bandwidth, time_grid, "Synthetic")
        
        rank_df <- bind_rows(rank_orig,rank_synth)
        
        variability_df <- rank_df %>%
          group_by(`Subject ID`, data_type) %>%
          summarise(variability = bounded_rank_variability(smoothed_rank), .groups = "drop")

        box_p <- ggplot(variability_df, aes(x = data_type, y = variability, fill = data_type)) +
          geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 2, outlier.alpha = 0.5) +
          labs(title = NULL,
               x = "Data type",
               y = "Rank variability") +
          theme_minimal() +
          scale_fill_manual(values = base_colors[1:2]) +
          scale_y_continuous(breaks = seq(-1,1,0.5), limits = c(-1,1), expand = c(0.01,0.01)) +
          theme(legend.position = "none")
        
        max_y <- variability_df %>%
          filter(!is.na(variability)) %>%
          group_by(data_type) %>%
          summarise(max_y = max(density(variability)$y)) %>%
          pull(max_y) %>%
          max()
        
        y_breaks <- pretty(c(0, ceiling(max_y)), n = 8)
        y_lims <- range(y_breaks)

        dens_p <- ggplot(variability_df, aes(x = variability, fill = data_type)) +
          geom_density(alpha = 0.7) +
          labs(title = NULL,
               y = "Density",
               x = "Rank variability") +
          theme_minimal() +
          scale_fill_manual(values = base_colors[1:2]) +
          scale_x_continuous(breaks = seq(-1,1,0.5), limits = c(-1,1), expand = c(0.01,0.01)) +
          scale_y_continuous(
          breaks = y_breaks, limits = c(0, ceiling(max_y))) +
          theme(legend.position = "none",
                axis.title.x = element_text(size = 12),
                axis.title.y = element_text(size = 12))

        # Combine side by side
        rank_plot <- box_p + dens_p
        
        indiv_plot <- plot_subject_specific_baselines(rank_df, n_per_group, 
                                                      n_groups = length(quantiles)+1,
                                                      max_time, x_interval)
        
        if (!is.null(save_path)) {
          ggsave(file.path(save_path_var, paste0("rank_order_", tolower(var_name), ".pdf")),
                 rank_plot, width = 8.27*1.5, height = 4.4*1.2, units = "in", dpi = 600)
          ggsave(file.path(save_path_var, paste0("indiv_traject_", tolower(var_name),"_n",n_per_group, ".pdf")), 
                 indiv_plot, width = 8.27*1.5, height = 4.4*1.2, units = "in", dpi = 600)
        } else {
          print(box_p)
          print(dens_p)
          print(indiv_plot)
        }
      }
    }
    
    if(meas_sim){
      
      # Get data
      orig_data <- original_data %>%
        select(all_of(id_col), all_of(time_col), all_of(var_name))
      
      synth_data <- synthetic_data %>%
        select(all_of(id_col), all_of(time_col), all_of(var_name))
      
      M_orig_full <- create_measurement_matrix_fast(orig_data, time_grid)
      M_synth_full <- create_measurement_matrix_fast(synth_data, time_grid)
      
      # Compute for original and synthetic
      density_orig <- compute_relative_density(M_orig_full)
      density_synth <- compute_relative_density(M_synth_full)
      
      # Build dataframe for plotting
      density_df <- data.frame(
        Time = seq_along(density_orig)-1,
        Original = density_orig,
        Synthetic = density_synth
      ) |>
        tidyr::pivot_longer(cols = c("Original", "Synthetic"),
                            names_to = "DataType",
                            values_to = "RelativeDensity")
      
      max_y <- max(c(density_orig, density_synth), na.rm = TRUE)
      y_breaks <- pretty(c(0, max_y), n = 6)
      y_lims <- range(pretty(c(0, max_y), n = 6))

      # Plot
      relat_meas <- ggplot(density_df, aes(x = Time, y = RelativeDensity, color = DataType)) +
        geom_line(size = 1) +
        labs(
          title = NULL,
          x = "Time (hours)",
          y = "Proportion of total measurements"
        ) +
        scale_color_manual(values = base_colors[1:2]) +
        scale_x_continuous(breaks = seq(0, max_time, by = x_interval), expand = c(0.01,0.01)) +
        scale_y_continuous(
          limits = y_lims,
          expand = c(0, 0),
          breaks = y_breaks) +
        theme_minimal() +
        theme(axis.title.x = element_text(size = 12),
              axis.title.y = element_text(margin = margin(r = 15), size = 12),
              strip.text = element_text(size = 12),
              legend.title = element_blank(),
              legend.text = element_text(size = 12),
              legend.position = "top")
      
      
      # Save or print
      if (!is.null(save_path)) {
        ggsave(file.path(save_path_var, paste0("real_meas_dens_", tolower(var_name), ".pdf")), 
               relat_meas, width = 8.27*1.5, height = 4.4*1.2, units = "in", dpi = 600)
      } else {
        print(relat_meas)
      }
      
      # Compute the statitics
      
      frobenius_norms <- numeric(meas_iter)
      frobenius_norms_ref <- numeric(meas_iter)
      
      similarity_scores <- numeric(meas_iter)
      similarity_scores_ref <- numeric(meas_iter)
      
      kl_scores <- numeric(meas_iter)
      kl_scores_ref <- numeric(meas_iter)
      
      for (i in seq_len(meas_iter)) {
        
        M_orig <- M_orig_full[sample(1:nrow(M_orig_full), meas_n),]
        M_orig_ref <- M_orig_full[sample(1:nrow(M_orig_full), meas_n),]
        M_synth <- M_synth_full[sample(1:nrow(M_synth_full), meas_n),]
        
        # Optimize and sort
        opt_order_ref <- opt_perm(M_orig, M_orig_ref)
        M_orig_ref_opt <- M_orig_ref[opt_order_ref,]
        
        opt_order_synth <- opt_perm(M_orig, M_synth)
        M_synth_opt <- M_synth[opt_order_synth,]
        
        frobenius_norms_ref[i] <- norm(M_orig - M_orig_ref_opt, type = "F")
        similarity_scores_ref[i] <- measurement_similarity(M_orig, M_orig_ref_opt)
        kl_scores_ref[i] <- dropout_divergence(M_orig, M_orig_ref_opt)
        
        frobenius_norms[i] <- norm(M_orig - M_synth_opt, type = "F")
        similarity_scores[i] <- measurement_similarity(M_orig, M_synth_opt)
        kl_scores[i] <- dropout_divergence(M_orig, M_synth_opt)
      }
      
      # Add mean row
      mean_row <- data.frame(
        Source = c("Reference", "Synthetic"),
        FrobeniusNorm = c(mean(frobenius_norms_ref, na.rm = TRUE), mean(frobenius_norms, na.rm = TRUE)),
        SimilarityScore = c(mean(similarity_scores_ref, na.rm = TRUE), mean(similarity_scores, na.rm = TRUE)),
        DropoutDivergence = c(mean(kl_scores_ref, na.rm = TRUE), mean(kl_scores, na.rm = TRUE))
      )
      
      # Write to CSV
      if(!is.null(save_path)){
        write.csv(mean_row, 
                  file = paste0(save_path_var, "meas_similarity_scores.csv"), 
                  row.names = FALSE)
      } else {
        
        print(mean_row)
      }
      
    }
  
    gc()  # Force garbage collection
    
  }
}
