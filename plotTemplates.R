# Grady Conwell 2024
# Templates
# style: theme_minimal(), viridis color palette, clean labels

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Load required libraries
if (!require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot2)
}
if (!require(dplyr)) {
  install.packages("dplyr")
  library(dplyr)
}
if (!require(viridis)) {
  install.packages("viridis")
  library(viridis)
}
if (!require(RColorBrewer)) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}

# Global plot configuration - easily change these to update all templates
plot_config <- list(
  font_family = "serif",
  color_palette = "Set3",
  # Default to some lovely pastels.
  # I did not lift these straight off of Apple's website.
  default_colors = c(
    "#5D7DAA",
    "#CF6980",
    "#817FB2",
    "#558D62",
    "#D77752",
    "#E8D25D"
  ),
  base_size = 20,
  point_size = 2,
  line_size = 1.2,
  smooth_size = 1.2,
  alpha_points = 0.3,
  alpha_bars = 0.8,
  alpha_hist = 0.7,
  alpha_box = 0.7,
  legend_position = "bottom"
)

# Template 1: Scatter Plot
create_scatter_plot <- function(data, x_col, y_col, color_col = NULL, title = "Scatter Plot", x_label = NULL, y_label = NULL, manual_colors = NULL) {
  p <- ggplot(data, aes(x = .data[[x_col]], y = .data[[y_col]], color = if(!is.null(color_col)) .data[[color_col]] else NULL)) +
    geom_point(size = plot_config$point_size, alpha = 1) +
    labs(title = title, x = if(!is.null(x_label)) x_label else x_col, y = if(!is.null(y_label)) y_label else y_col, color = NULL) +
    theme_minimal(base_size = plot_config$base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = plot_config$legend_position,
      text = element_text(family = plot_config$font_family)
    )
  
  if (!is.null(manual_colors)) {
    p <- p + scale_color_manual(values = manual_colors)
  } else {
    p <- p + scale_color_manual(values = plot_config$default_colors)
  }
  
  return(p)
}

# Example usage:
# p1 <- create_scatter_plot(data, "x", "y", "category", "My Scatter Plot")
# print(p1)

# Template 2: Bar Plot
create_bar_plot <- function(data, x_col, y_col, fill_col = NULL, title = "Bar Plot", x_label = NULL, y_label = NULL, manual_colors = NULL) {
  p <- ggplot(data, aes(x = .data[[x_col]], y = .data[[y_col]], fill = if(!is.null(fill_col)) .data[[fill_col]] else NULL)) +
    geom_bar(stat = "identity", alpha = plot_config$alpha_bars) +
    labs(title = title, x = if(!is.null(x_label)) x_label else x_col, y = if(!is.null(y_label)) y_label else y_col, fill = NULL) +
    theme_minimal(base_size = plot_config$base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = plot_config$legend_position,
      text = element_text(family = plot_config$font_family)
    )
  
  if (!is.null(manual_colors)) {
    p <- p + scale_fill_manual(values = manual_colors)
  } else {
    p <- p + scale_fill_manual(values = plot_config$default_colors)
  }
  
  return(p)
}

# Example usage:
# summary_data <- data %>% group_by(category) %>% summarise(mean_y = mean(y))
# p2 <- create_bar_plot(summary_data, "category", "mean_y", "category", "Mean Y by Category")
# print(p2)

# Template 3: Line Plot
create_line_plot <- function(data, x_col, y_col, color_col = NULL, title = "Line Plot", x_label = NULL, y_label = NULL, manual_colors = NULL) {
  p <- ggplot(data, aes(x = .data[[x_col]], y = .data[[y_col]], color = if(!is.null(color_col)) .data[[color_col]] else NULL)) +
    geom_line(size = plot_config$line_size) +
    geom_point(size = plot_config$point_size) +
    labs(title = title, x = if(!is.null(x_label)) x_label else x_col, y = if(!is.null(y_label)) y_label else y_col, color = NULL) +
    theme_minimal(base_size = plot_config$base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = plot_config$legend_position,
      text = element_text(family = plot_config$font_family)
    )
  
  if (!is.null(manual_colors)) {
    p <- p + scale_color_manual(values = manual_colors)
  } else {
    p <- p + scale_color_manual(values = plot_config$default_colors)
  }
  
  return(p)
}

# Example usage:
# time_data <- data.frame(time = 1:10, value = cumsum(rnorm(10)), group = rep(c("A", "B"), each = 5))
# p3 <- create_line_plot(time_data, "time", "value", "group", "Time Series")
# print(p3)

# Template 4: Histogram
create_histogram <- function(data, x_col, fill_col = NULL, title = "Histogram", x_label = NULL, y_label = NULL, manual_colors = NULL) {
  p <- ggplot(data, aes(x = .data[[x_col]], fill = if(!is.null(fill_col)) .data[[fill_col]] else NULL)) +
    geom_histogram(alpha = plot_config$alpha_hist, bins = 30) +
    labs(title = title, x = if(!is.null(x_label)) x_label else x_col, y = if(!is.null(y_label)) y_label else "Count", fill = NULL) +
    theme_minimal(base_size = plot_config$base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = plot_config$legend_position,
      text = element_text(family = plot_config$font_family)
    )
  
  if (!is.null(manual_colors)) {
    p <- p + scale_fill_manual(values = manual_colors)
  } else {
    p <- p + scale_fill_manual(values = plot_config$default_colors)
  }
  
  return(p)
}

# Example usage:
# p4 <- create_histogram(data, "y", "category", "Distribution of Y")
# print(p4)

# Template 5: Box Plot
create_box_plot <- function(data, x_col, y_col, fill_col = NULL, title = "Box Plot", x_label = NULL, y_label = NULL, manual_colors = NULL) {
  p <- ggplot(data, aes(x = .data[[x_col]], y = .data[[y_col]], fill = if(!is.null(fill_col)) .data[[fill_col]] else NULL)) +
    geom_boxplot(alpha = plot_config$alpha_box) +
    labs(title = title, x = if(!is.null(x_label)) x_label else x_col, y = if(!is.null(y_label)) y_label else y_col, fill = NULL) +
    theme_minimal(base_size = plot_config$base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = plot_config$legend_position,
      text = element_text(family = plot_config$font_family)
    )
  
  if (!is.null(manual_colors)) {
    p <- p + scale_fill_manual(values = manual_colors)
  } else {
    p <- p + scale_fill_manual(values = plot_config$default_colors)
  }
  
  return(p)
}

# Example usage:
# p5 <- create_box_plot(data, "category", "y", "category", "Y by Category")
# print(p5)

# Template 6: Superimposed Scatter Plot
create_superimposed_scatter <- function(data1, data2, x_col, y_col, label1 = "Dataset 1", label2 = "Dataset 2", title = "Superimposed Scatter Plot", x_label = NULL, y_label = NULL, manual_colors = NULL) {
  # Combine data with labels
  data1$dataset <- label1
  data2$dataset <- label2
  combined <- bind_rows(data1, data2)
  
  p <- ggplot(combined, aes(x = .data[[x_col]], y = .data[[y_col]], color = dataset, shape = dataset)) +
    geom_point(size = plot_config$point_size, alpha = 1) +
    geom_smooth(method = "lm", se = FALSE, size = plot_config$smooth_size) +
    scale_shape_manual(values = c(16, 17)) +  # circle and triangle
    labs(title = title, x = if(!is.null(x_label)) x_label else x_col, y = if(!is.null(y_label)) y_label else y_col, color = NULL) +
    theme_minimal(base_size = plot_config$base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = plot_config$legend_position,
      text = element_text(family = plot_config$font_family)
    )
  
  if (!is.null(manual_colors)) {
    p <- p + scale_color_manual(values = manual_colors)
  } else {
    p <- p + scale_color_manual(values = plot_config$default_colors)
  }
  
  return(p)
}

# Example usage:
# data1 <- data.frame(x = rnorm(20), y = rnorm(20))
# data2 <- data.frame(x = rnorm(20) + 1, y = rnorm(20) + 1)
# p6 <- create_superimposed_scatter(data1, data2, "x", "y", "Sensor 1", "Sensor 2", "Temperature Comparison")
# print(p6)

# Template 7: Multi-Line Plot
# data_list: list of data frames
# labels: vector of labels for each dataset
create_multi_line_plot <- function(data_list, x_col, y_col, labels, title = "Multi-Line Plot", x_label = NULL, y_label = NULL, manual_colors = NULL, show_points = TRUE, show_lines = TRUE, add_smooth = FALSE, smooth_method = "loess", smooth_se = FALSE, smooth_span = 0.75) {
  combined <- bind_rows(lapply(seq_along(data_list), function(i) {
    df <- data_list[[i]]
    df$dataset <- labels[i]
    df
  }))
  
  p <- ggplot(combined, aes(x = .data[[x_col]], y = .data[[y_col]], color = dataset))
  
  if (show_lines) {
    p <- p + geom_line(size = plot_config$line_size)
  }
  
  if (show_points) {
    p <- p + geom_point(size = plot_config$point_size, alpha = plot_config$alpha_points)
  }
  
  if (add_smooth) {
    p <- p + geom_smooth(method = smooth_method, se = smooth_se, span = smooth_span, size = plot_config$smooth_size)
  }
  
  p <- p + labs(title = title, x = if(!is.null(x_label)) x_label else x_col, y = if(!is.null(y_label)) y_label else y_col, color = NULL) +
    theme_minimal(base_size = plot_config$base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = plot_config$legend_position,
      axis.title.x = element_text(margin = margin(t = 15)),
      axis.title.y = element_text(margin = margin(r = 15)),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
      text = element_text(family = plot_config$font_family)
    )
  
  if (!is.null(manual_colors)) {
    p <- p + scale_color_manual(values = manual_colors)
  } else {
    p <- p + scale_color_manual(values = plot_config$default_colors)
  }
  
  return(p)
}

# Example usage for a pretty temp sensor plot:
#   x_label = "Time",
#   y_label = "Temperature (°C)",
#   show_points = TRUE,
#   show_lines = FALSE,
#   add_smooth = TRUE,
#   smooth_span = 0.2

# Example usage:
# sensor1 <- data.frame(time = 1:10, temp = cumsum(rnorm(10)) + 20)
# sensor2 <- data.frame(time = 1:10, temp = cumsum(rnorm(10)) + 22)
# p7 <- create_multi_line_plot(list(sensor1, sensor2), "time", "temp", c("Sensor 1", "Sensor 2"), "Temperature Over Time")
# print(p7)

# Template 8: Multi-Line Plot Example 
# # Assuming you have data frames for each sensor
# sensor1_data <- data.frame(time = your_times, temperature = sensor1_readings)
# sensor2_data <- data.frame(time = your_times, temperature = sensor2_readings)

# p <- create_multi_line_plot(
#   list(sensor1_data, sensor2_data), 
#   "time", 
#   "temperature", 
#   c("Sensor 1", "Sensor 2"), 
#   "Temperature Comparison"
# )
# print(p)

# Template 9: Faceted Line Plot
create_faceted_line_plot <- function(data_list, x_col, y_col, labels, types, title = "Faceted Plot", x_label = NULL, y_label = NULL, manual_colors = NULL, show_points = TRUE, show_lines = TRUE, add_smooth = FALSE, smooth_span = 0.75) {
  # data_list: list of data frames
  # labels: vector of labels for each dataset (for the legend/color)
  # types: vector of types for each dataset (for faceting)
  
  combined <- bind_rows(lapply(seq_along(data_list), function(i) {
    df <- data_list[[i]]
    df$dataset <- labels[i]
    df$facet_type <- types[i]
    df
  }))
  
  p <- ggplot(combined, aes(x = .data[[x_col]], y = .data[[y_col]], color = dataset))
  
  if (show_lines) {
    p <- p + geom_line(linewidth = plot_config$line_size)
  }
  
  if (show_points) {
    p <- p + geom_point(size = plot_config$point_size, alpha = plot_config$alpha_points)
  }
  
  if (add_smooth) {
    p <- p + geom_smooth(method = "loess", se = FALSE, span = smooth_span, linewidth = plot_config$smooth_size)
  }
  
  p <- p + facet_grid(rows = vars(facet_type), scales = "free_y", switch = "y") +
    labs(title = title, x = if(!is.null(x_label)) x_label else x_col, y = NULL, color = NULL) +
    theme_minimal(base_size = plot_config$base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = plot_config$legend_position,
      axis.title.x = element_text(margin = margin(t = 15)),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
      text = element_text(family = plot_config$font_family),
      strip.text = element_text(face = "bold", size = plot_config$base_size * 0.9),
      strip.placement = "outside"
    )
  
  if (!is.null(manual_colors)) {
    p <- p + scale_color_manual(values = manual_colors)
  } else {
    p <- p + scale_color_manual(values = plot_config$default_colors)
  }
  
  return(p)
}

# Example usage:
# create_faceted_line_plot(list(data1, data2, data3), "time", "value", c("Sensor 1", "Sensor 2", "Sensor 3"), c("Temp", "Fan", "Current"))