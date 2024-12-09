# =============================================================================
# Title:        Plots Script for Ribosome Footprinting data from Bohlen et al., 2020
# Author:       Lauren A. Shelby, Ziyue Cheng
# Date:         December 9, 2024
# Description:  This script has plot-specific code for ribosome footprinting 
#               data in TSV format. 
# =============================================================================

# Variables ---

sample_names = c(
  "SRR10416856",  # 80S WT
  "SRR11553468"  # 40S WT
)

elements = c("cds_end", "cds_start", "stst")

.heatmap_y_min = 15
.heatmap_y_max = 70 
.heatmap_x_min = -50
.heatmap_x_max = 50

.background = "#EBEBEB"

.distance_40to80 = -3 # convert 40S coordinate to 80S
.length_40to80 = 4    # convert 40S coordinate to 80S

sample80 <- read.csv("sample80.csv")
sample40 <- read.csv("sample40.csv")

# Functions ----

library(ggplot2)
library(dplyr)
library(DMwR2)
library(FNN)


# Find the 'box' of length and position peak ranges for heatmapping 
find_top_range = function(sample, sample_type = "CDS Start", count_fraction = 0.75) {
  # CDS Start range
  
  # get the largest point when position < 0
  top = sample %>%
    filter(type == sample_type, distance2element < 0) %>%
    arrange(desc(count_normalized)) %>%
    slice(1)
  
  # get the largest count
  max_count = top$count_normalized
  
  values = sample %>%
    filter(type == sample_type, distance2element < 0, count >= max_count * count_fraction) %>%
    arrange(desc(count_normalized))
  
  # distance range
  lower_range_distance = top$distance2element
  
  while((lower_range_distance - 1) %in% values$distance2element) {
    lower_range_distance = lower_range_distance - 1
  }
  
  higher_range_distance = top$distance2element
  
  while((higher_range_distance + 1) %in% values$distance2element) {
    higher_range_distance = higher_range_distance + 1
  }
  
  values_in_range = values %>%
    filter(distance2element %in% lower_range_distance:higher_range_distance)
  
  # length range
  lower_range_length = top$length
  
  while((lower_range_length - 1) %in% values_in_range$length) {
    lower_range_length = lower_range_length - 1
  }
  
  higher_range_length = top$length
  
  while((higher_range_length + 1) %in% values_in_range$length) {
    higher_range_length = higher_range_length + 1
  }
  
  return(c(
    lower_range_distance,
    higher_range_distance,
    lower_range_length,
    higher_range_length
  ))
}

# Defining 'box' limits
.heatmap_80S_y_higher_line = find_top_range(sample80)[4]
.heatmap_80S_y_lower_line = find_top_range(sample80)[3]
.heatmap_80S_x_left_line = find_top_range(sample80)[1]
.heatmap_80S_x_right_line = find_top_range(sample80)[2]

.heatmap_40S_y_higher_line = find_top_range(sample40)[4]
.heatmap_40S_y_lower_line = find_top_range(sample40)[3]
.heatmap_40S_x_left_line = find_top_range(sample40)[1]
.heatmap_40S_x_right_line = find_top_range(sample40)[2]

# read count table HIDDEN FUNCTION
.getCountTable = function(ribo_type, element, prime = "5'") {
  # read different data set based on ribo type
  samples = switch (
    ribo_type,
    "80S" = .sample_80S,
    "40S" = .sample_40S
  )
  
  # read all count tables and add them together, take into account if strange genes were removed
  if (prime == "5'") {
    #.heatmap_x_max = 20
    if (file.exists(paste0("transcript-data/", samples[1], "_read_to_", tolower(element), "_new.tsv"))) {
      table = read.delim(paste0("transcript-data/", samples[1], "_read_to_", tolower(element), "_new.tsv"))
      if (length(samples) != 1) {
        for (i in 2:length(samples)) {
          table = table %>%
            rbind(read.delim(paste0("transcript-data/", samples[i], "_read_to_", tolower(element), "_new.tsv")))
        }
      }
    } else {
      table = read.delim(paste0("transcript-data/", samples[1], "_read_to_", tolower(element), ".tsv"))
      if (length(samples) != 1) {
        for (i in 2:length(samples)) {
          table = table %>%
            rbind(read.delim(paste0("transcript-data/", samples[i], "_read_to_", tolower(element), ".tsv")))
        }
      } 
    }
  } else if (prime == "3'") {
    #.heatmap_x_max = 50
    if (file.exists(paste0("transcript-3prime/", samples[1], "_read_to_", tolower(element), "_new.tsv"))) {
      table = read.delim(paste0("transcript-3prime/", samples[1], "_read_to_", tolower(element), "_new.tsv"))
      if (length(samples) != 1) {
        for (i in 2:length(samples)) {
          table = table %>%
            rbind(read.delim(paste0("transcript-3prime/", samples[i], "_read_to_", tolower(element), "_new.tsv")))
        }
      }
    } else {
      table = read.delim(paste0("transcript-3prime/", samples[1], "_read_to_", tolower(element), ".tsv"))
      if (length(samples) != 1) {
        for (i in 2:length(samples)) {
          table = table %>%
            rbind(read.delim(paste0("transcript-3prime/", samples[i], "_read_to_", tolower(element), ".tsv")))
        }
      } 
    }
  }
  
  
  table = table %>%
    filter(distance2element >= .heatmap_x_min, distance2element <= .heatmap_x_max) %>%
    filter(length >= .heatmap_y_min, length <= .heatmap_y_max) %>%
    group_by(distance2element, length) %>%
    summarise_all(function(x) sum(x, na.rm = T)) %>%
    ungroup()
  
  return(table)
}

# Creating a function to add in all possible lenxpos values to the count table df for imputation
# Same as above but for smoothed heatmaps 
getCountTable2 <- function(ribo_type, element, kn=4, prime = "5'") {
  # read different data set based on ribo type
  samples = switch (
    ribo_type,
    "80S" = .sample_80S,
    "40S" = .sample_40S
  )
  
  #prime = ensym(prime)
  
  # read all count tables and add them together, take into account if strange genes were removed
  if (prime == "5'") {
    if (file.exists(paste0("transcript-data/", samples[1], "_read_to_", tolower(element), "_new.tsv"))) {
      table = read.delim(paste0("transcript-data/", samples[1], "_read_to_", tolower(element), "_new.tsv"))
      if (length(samples) != 1) {
        for (i in 2:length(samples)) {
          table = table %>%
            rbind(read.delim(paste0("transcript-data/", samples[i], "_read_to_", tolower(element), "_new.tsv")))
        }
      }
    } else {
      table = read.delim(paste0("transcript-data/", samples[1], "_read_to_", tolower(element), ".tsv"))
      if (length(samples) != 1) {
        for (i in 2:length(samples)) {
          table = table %>%
            rbind(read.delim(paste0("transcript-data/", samples[i], "_read_to_", tolower(element), ".tsv")))
        }
      } 
    }
  } else if (prime == "3'"){
    if (file.exists(paste0("transcript-3prime/", samples[1], "_read_to_", tolower(element), "_new.tsv"))) {
      table = read.delim(paste0("transcript-3prime/", samples[1], "_read_to_", tolower(element), "_new.tsv"))
      if (length(samples) != 1) {
        for (i in 2:length(samples)) {
          table = table %>%
            rbind(read.delim(paste0("transcript-3prime/", samples[i], "_read_to_", tolower(element), "_new.tsv")))
        }
      }
    } else {
      table = read.delim(paste0("transcript-3prime/", samples[1], "_read_to_", tolower(element), ".tsv"))
      if (length(samples) != 1) {
        for (i in 2:length(samples)) {
          table = table %>%
            rbind(read.delim(paste0("transcript-3prime/", samples[i], "_read_to_", tolower(element), ".tsv")))
        }
      } 
    }
  }
  
  table = table %>%
    filter(distance2element >= .heatmap_x_min, distance2element <= .heatmap_x_max) %>%
    filter(length >= .heatmap_y_min, length <= .heatmap_y_max) %>%
    group_by(distance2element, length) %>%
    summarise_all(function(x) sum(x, na.rm = T)) %>%
    ungroup()
  
  # Generate all possible combinations of length and distance2element
  all_combinations <- expand.grid(
    length = unique(table$length),
    distance2element = unique(table$distance2element))
  
  # Merge the all_combinations with the original data frame
  df_full <- merge(all_combinations, table, by = c("length", "distance2element"), all.x = TRUE)
  
  df_full <- df_full %>% 
    select(distance2element, length, count_normalized)
  
  # Impute missing count_normalized values with kNN imputation
  df_imp <- knnImputation(df_full, k=kn)
  
  # Prepare the data for k-NN
  coords <- df_imp[, c("length", "distance2element")]
  count_normalized <- df_imp$count_normalized
  
  # Apply k-NN regression to smooth the count_normalized values
  knn_result <- knn.reg(train = coords, test = coords, y = count_normalized, k = kn)
  
  # Extract the smoothed values
  df_imp$count_normalized_smoothed <- knn_result$pred
  
  return(df_imp)
}

# tweak heatmap: add lines, change axis breaks
.addHeatmapLines = function(ggplot_object, ribo_type) {
  .heatmap_y_lower_line = ifelse(ribo_type == "80S", .heatmap_80S_y_lower_line, .heatmap_40S_y_lower_line)
  .heatmap_y_higher_line  = ifelse(ribo_type == "80S", .heatmap_80S_y_higher_line, .heatmap_40S_y_higher_line)
  .heatmap_x_left_line   = ifelse(ribo_type == "80S", .heatmap_80S_x_left_line, .heatmap_40S_x_left_line)
  .heatmap_x_right_line  = ifelse(ribo_type == "80S", .heatmap_80S_x_right_line, .heatmap_40S_x_right_line)
  
  ggplot_object = ggplot_object + 
    # lower border
    geom_hline(yintercept = .heatmap_y_lower_line - 0.5, linetype="dashed", color = "black") +
    geom_text(aes(.heatmap_x_max + 2, .heatmap_y_lower_line - 2.5, label = .heatmap_y_lower_line)) + 
    # higher border
    geom_hline(yintercept = .heatmap_y_higher_line + 0.5, linetype="dashed", color = "black") +
    geom_text(aes(.heatmap_x_max + 2, .heatmap_y_higher_line + 2.5, label = .heatmap_y_higher_line)) + 
    # left border
    geom_vline(xintercept = .heatmap_x_left_line - 0.5, linetype="dashed", color = "black") +
    geom_text(aes(.heatmap_x_left_line - 2.5, .heatmap_y_max + 2, label = .heatmap_x_left_line)) + 
    # right border
    geom_vline(xintercept = .heatmap_x_right_line + 0.5, linetype="dashed", color = "black") +
    geom_text(aes(.heatmap_x_right_line + 2.5, .heatmap_y_max + 2, label = .heatmap_x_right_line))
  return(ggplot_object)
}

.changeTheme = function(ggplot_object) {
  ggplot_object = ggplot_object +
    coord_cartesian(
      ylim = c(.heatmap_y_min - 0.5, .heatmap_y_max + 3),
      xlim = c(.heatmap_x_min - 0.5, .heatmap_x_max + 3),
      expand = FALSE
    ) +
    scale_y_continuous(breaks = seq(.heatmap_y_min, .heatmap_y_max, by = 5)) +
    scale_x_continuous(breaks = seq(.heatmap_x_min, .heatmap_x_max, by = 5)) + 
    theme(
      panel.background = element_rect(fill = .background),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.x = element_line(colour = 'black', linewidth=0.5, linetype='solid'),
      axis.line.y = element_line(colour = 'black', linewidth=0.5, linetype='solid'),
      axis.text  = element_text(size=12),
      axis.title = element_text(size=14), # ,face="bold"
      plot.title = element_text(size=18) # , hjust = 0.5
    ) 
  
  return(ggplot_object)
}


# get limits of a numeric vector c(min, max)
.getHeatmapLimits = function(num_vector, log10_fill) {
  return(
    c(ifelse(log10_fill, 1, 0), ceiling(max(num_vector, na.rm = T)))
  )
}

# get 4 breaks based on that limit
.generate4Breaks = function(limit, log10_fill = F) {
  if (log10_fill) {
    return(c(
      1,
      ceiling(limit[2] / 1000),
      floor(limit[2] / 100),
      floor(limit[2] / 10),
      floor(limit[2])
    ))
  } else {
    return(c(
      limit[1],
      floor(0.75 * limit[1] + 0.25 * limit[2]),
      floor(0.5  * limit[1] + 0.5  * limit[2]),
      floor(0.25 * limit[1] + 0.75 * limit[2]),
      limit[2]
    ))
  }
}

# get either 40S or 80S heatmap to that element 
getRiboHeatmap = function(ribo_type, element, log10_fill = F, plot_limits = NULL, prime = "5'") {
  sample = .getCountTable(ribo_type, element, prime)
  
  if (is.null(plot_limits))
    plot_limits = .getHeatmapLimits(sample$count_normalized, log10_fill)
  
  plot_breaks = .generate4Breaks(plot_limits, log10_fill)
  
  if (log10_fill) { # log 10 heatmap
    this_heatmap = sample %>% ggplot + 
      geom_tile(aes(distance2element, length, fill= count_normalized)) + 
      scale_fill_gradientn(
        colours = c("#0000AA", "#0055FF", "#00AAFF", "#80FFFF", "#FFFF80", "#FFAA00", "#FF5500", "#AA0000"),
        trans = "log10", na.value = "white",
        breaks = plot_breaks, limits = plot_limits
      )
  } else {
    this_heatmap = sample %>% ggplot + 
      geom_tile(aes(distance2element, length, fill= count_normalized)) + 
      scale_fill_gradient(
        low = "white", high = "red4", na.value = .background,
        breaks = plot_breaks, limits = plot_limits
      )
  }
  
  # add labels
  this_heatmap = this_heatmap +
    labs(
      title = paste0(ribo_type, " around ", gsub("new", "", gsub("_", " ", toupper(element)))),
      x = paste0("Footprint distance to ", gsub("new", "", gsub("_", " ", toupper(element)))),
      y = "Length (nt)",
      fill = "Normalized\ncount"
    )
  # add lines
  this_heatmap %>% .addHeatmapLines(ribo_type) %>% .changeTheme %>% return
}

# Get either 40S or 80S heatmap to that element smoothed on a specified kNN
getRiboHeatmapSmooth = function(ribo_type, element, log10_fill = F, kn = 4, plot_limits = NULL, prime = "5'") {
  sample = getCountTable2(ribo_type, element, kn, prime)
  
  if (is.null(plot_limits))
    plot_limits = .getHeatmapLimits(sample$count_normalized_smoothed, log10_fill)
  
  plot_breaks = .generate4Breaks(plot_limits, log10_fill)
  
  if (log10_fill) { # log 10 heatmap
    this_heatmap = sample %>% ggplot + 
      geom_tile(aes(distance2element, length, fill= count_normalized_smoothed)) + 
      scale_fill_gradientn(
        colours = c("#0000AA", "#0055FF", "#00AAFF", "#80FFFF", "#FFFF80", "#FFAA00", "#FF5500", "#AA0000"),
        trans = "log10", na.value = "white",
        breaks = plot_breaks, limits = plot_limits
      )
  } else {
    this_heatmap = sample %>% ggplot + 
      geom_tile(aes(distance2element, length, fill= count_normalized_smoothed)) + 
      scale_fill_gradient(
        low = "white", high = "red4", na.value = .background,
        breaks = plot_breaks, limits = plot_limits
      )
  }
  
  # add labels
  this_heatmap = this_heatmap +
    labs(
      title = paste0(ribo_type, " around ", gsub("new", "", gsub("_", " ", toupper(element))), " (with smoothing)"),
      x = paste0("Footprint distance to ", gsub("new", "", gsub("_", " ",toupper(element)))),
      y = "Length (nt)",
      fill = "Normalized\ncount"
    )
  # add lines
  this_heatmap %>% .addHeatmapLines(ribo_type) %>% .changeTheme %>% return
}

# Creating a plot that displays the 5' and 3' position densities for a given element on the same density plot 
getElementDensityOverlayPos <- function(ribo_type, sample, colors, x_limits = c(-100, 100), y_limits = c(0, 0.1)) {
  ribotype <- ensym(ribo_type)
  
  # Filter and prepare data for 5' plot
  sample_5 <- sample %>%
    filter(ribo == ribotype, POS5 < .heatmap_x_max & POS5 > .heatmap_x_min) %>%
    group_by(element)
  
  # Filter and prepare data for 3' plot
  sample_3 <- sample %>%
    filter(ribo == ribotype, POS3 < .heatmap_x_max & POS3 > .heatmap_x_min) %>%
    group_by(element)
  
  # Combine 5' and 3' data for overlay plotting
  plot <- ggplot() +
    # Add 5' density layer (below)
    geom_density(data = sample_5, aes(x = POS5, y = after_stat(density), color = element), size = 1, linetype = "solid") +
    # Add 3' density layer (on top)
    geom_density(data = sample_3, aes(x = POS3, y = after_stat(density), color = element), size = 1, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
    ggtitle("5' and 3' end position density distribution for all elements") +
    ylab("FP end frequency") +
    xlab("Distance from start codon (nt)") +
    scale_color_manual(values = colors) +
    #scale_x_continuous(limits = x_limits) +
    #scale_y_continuous(limits = y_limits, breaks = seq(0, y_limits[2], by = 0.025)) +
    theme_minimal() +
    theme(
      legend.position = "right",
      axis.text.y = element_text(size = 8)
    )
  
  return(plot)
}

#Same as above but length density
getElementDensityOverlayLen <- function(ribo_type, sample, colors, x_limits = c(-100, 100), y_limits = c(0, 0.1)) {
  ribotype <- ensym(ribo_type)
  
  # Filter and prepare data for 5' plot
  sample_5 <- sample %>%
    filter(ribo == ribotype, LEN < .heatmap_y_max & LEN > .heatmap_y_min) %>%
    group_by(element)
  
  # Combine 5' and 3' data for overlay plotting
  plot <- ggplot() +
    # Add 5' density layer (below)
    geom_density(data = sample_5, aes(x = LEN, y = after_stat(density), color = element), size = 1, linetype = "solid") +
    #geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
    ggtitle("Footprint length density distribution for all elements") +
    ylab("Frequency") +
    xlab("Footprint Length") +
    scale_color_manual(values = colors) +
    theme_minimal() +
    theme(
      legend.position = "right",
      axis.text.y = element_text(size = 8)
    )
  
  return(plot)
}

