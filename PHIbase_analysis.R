# Script to look at the PHIbase minimap data 

library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(dplyr)
library(lubridate)
library(scales)
library(patchwork)

setwd("~/OneDrive - Norwich BioScience Institutes/Air_Samples/Tiptree")


# Script to create heatmaps from Tiptree data ---

# Reading in the data --------
raw_f_a <- read.csv("Tiptree_taxa_counts_initial_analysis/Q20_Tiptree_combined_data.csv", header=TRUE)
filter_f_a <- raw_f_a %>% 
  select(taxaID, read_count, species, barcode, total_reads, Sample, Date.collected, repeat., location)

# Dec data
raw_dec_data <- read.csv("Tiptree_taxa_counts_initial_analysis/Tiptree_Dec_unfiltered.csv", header = TRUE)

# Process the data ---
#Want to join the two datasets
all_heatmap_data <- rbind(filter_f_a , raw_dec_data)

# Process the date column
all_heatmap_data$Date.collected = as.Date(all_heatmap_data$Date.collected, format="%d/%m/%Y")

#Adding a col for month
all_heatmap_data$month <- 
  format(as.Date(all_heatmap_data$Date.collected, format="%d/%m/%Y"),"%b")

#reorder month cols
all_heatmap_data$month <- 
  factor(all_heatmap_data$month,
         levels=c("Dec", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug"))

# Change labels and order
all_heatmap_data <- all_heatmap_data %>%
  mutate(location = factor(location,
                           levels = c("F", "RG", "LG"), 
                           labels = c("Field", "Greenhouse 1", "Greenhouse 2")))

# Need to normalise the reads for each location, e.g. hits per 10,000 reads
# Therefore we need the total number of reads for each sample, not just the classified reads
all_heatmap_filtered <- 
  all_heatmap_data %>% 
  filter(read_count >=10 ) %>%                      #filter out any species with less than 10 reads
  mutate(
    hits_per_100.000 = (read_count * 100000)/(total_reads)) %>% 
  mutate(
    log_hits = log10(hits_per_100.000)) %>%     #log scale-Mutate the data to add a new row
  mutate(Sample = as.character(Sample))



# Plotting ---
# Heatmap function:
plot_species_heatmap <- function(data,
                                 x = "Sample",
                                 y = "species",
                                 xlab = x,
                                 ylab = y,
                                 fill = "log_hits", # Can change this if don't want to plot logs
                                 flip_axes = FALSE,
                                 title = "Species Abundance Heatmap",
                                 save_path = NULL,
                                 width = 12,
                                 height = 18) {
  
  # Compute abundance
  species_abundance <- data %>%
    group_by(species) %>%
    summarise(total_abundance = sum(log_hits, na.rm = TRUE)) %>%
    arrange(desc(total_abundance))
  
  # Apply the order to the main data (most abundant species at the bottom)
  data <- data %>%
    mutate(species = factor(species, levels = species_abundance$species))
  
  # Flip x and y if needed
  if (flip_axes) {
    temp <- x
    x <- y
    y <- temp
  }
  
  # Dynamically map aesthetics
  aes_mapping <- aes_string(x = x, y = y, fill = fill)
  
  # Create plot
  p <- ggplot(data, mapping = aes_mapping) +
    geom_tile() +
    
    # scale_fill_distiller(
    #   name = "Hits per 100,000",
    #   palette = "YlOrRd",
    #   labels = function(x) scales::comma(10^x)
    #   ) +
    
    scale_fill_viridis_c(
      name = "Hits per 100,000",
      labels = function(x) round(10^x),
      option = "B",  # Options: "D", "C", "B", etc.
      direction = 1
    ) +

    labs(title = title, x = xlab, y = ylab) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = if (flip_axes) 70 else 0, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      legend.key.size = unit(1.5, "cm"),
      legend.key.width = unit(0.4, "cm"),
      axis.line = element_line(color = "black", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank()
    )
  
  # Save plot if path provided
  if (!is.null(save_path)) {
    ggsave(filename = save_path, plot = p, width = width, height = height)
  }
  
  return(p)
}

#Preparing filtered data ----

#filtering out taxa which are only present in <10 samples - less stringent when grouped
all_heatmap_filtered_10taxa <-
  all_heatmap_filtered %>% 
  group_by(species) %>% 
  filter(n()>=10)

# filtering out taxa which are only present in <15 samples - for all the samples to be plotted seperate
all_heatmap_filtered_15taxa <- all_heatmap_filtered %>% 
  group_by(species) %>% 
  filter(n()>=15)


# Call the function ---
# All samples-
plot_species_heatmap(
  data = all_heatmap_filtered_15taxa,
  x = "Sample",
  y = "species",
  ylab="Species",
  title = "",
  save_path = "Graphs/Heatmaps/heatmap_by_sample.pdf",
  width = 18
)

# Group by month -
plot_species_heatmap(
  data = all_heatmap_filtered_10taxa,
  x = "month",
  y = "species",
  xlab="Month",
  ylab="Species",
  title = "",
  save_path = "Graphs/Heatmaps/heatmap_by_month.pdf"
)

#Month rotated
plot_species_heatmap(
  data = all_heatmap_filtered_10taxa,
  x = "month",
  y = "species",
  ylab="Month",
  xlab="Species",
  title = "",
  flip_axes = TRUE,
  save_path = "Graphs/Heatmaps/heatmap_by_month_rot.pdf",
  width = 18,
  height = 12
)

# Plot by location -
plot_species_heatmap(
  data = all_heatmap_filtered_10taxa,
  x = "location",
  y = "species",
  xlab = "Location",
  ylab="Species",
  title = "",
  save_path = "Graphs/Heatmaps/heatmap_by_location.pdf"
)

# Rotated plot
plot_species_heatmap(
  data = all_heatmap_filtered_10taxa,
  x = "location",
  y = "species",
  ylab = "Location",
  xlab="Species",
  title = "",
  flip_axes = TRUE,
  save_path = "Graphs/Heatmaps/heatmap_by_location_rot.pdf",
  width = 18,
  height = 12
)


# Top 20 Species ----------------------------------------------------------



# Bar plot of the most abundant species --------

abundance_summary <- all_heatmap_filtered %>%
  group_by(species) %>%
  summarise(total_hits = sum(hits_per_100.000, na.rm = TRUE)) %>%
  arrange(desc(total_hits)) %>%
  slice_max(total_hits, n = 20)  # top 20 species

# Reorder data
abundance_summary <- abundance_summary %>%
  mutate(species = factor(species, levels = rev(species)))  # bar chart top-down

most_abundant <- ggplot(abundance_summary, aes(x = species, y = total_hits)) +
  scale_y_continuous(
    labels = scales::label_comma(accuracy = 1)
  ) +
  geom_col(fill = "#4FBDC5") +
  coord_flip() +
  labs(
    x = "Species",
    y = "Total Hits per 100,000"
  ) + theme_minimal(base_size = 14) +
  theme(
    plot.margin = margin(t = 10, r = 20, b = 10, l = 30),
    axis.line = element_line(color = "black", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

most_abundant

ggsave(filename = "Graphs/PHIbase_Top20_Bar.pdf", plot = most_abundant, width = 10, height = 6)

# Same top 20 but this time don't summarise ----
top_20_species_names <- abundance_summary$species
# Dataset with just the top 20 
top_20_sp <- all_heatmap_filtered %>%
  filter(species %in% top_20_species_names)

# Get Mean & SD
top_20_summary <- top_20_sp %>%
  group_by(species) %>%
  summarise(
    mean_hits = mean(hits_per_100.000, na.rm = TRUE),
    sd_hits = sd(hits_per_100.000, na.rm = TRUE),
    n = n(),
    se_hits = sd_hits / sqrt(n)
  ) %>%
  arrange(desc(mean_hits))

#Nice order
top_20_summary$species <- factor(top_20_summary$species, levels = rev(top_20_summary$species))


# Plot with error bars
most_abundant_with_error <- ggplot(top_20_summary, aes(x = species, y = mean_hits)) +
  geom_col(fill = "#4FBDC5") +
  geom_errorbar(aes(ymin = mean_hits - se_hits, ymax = mean_hits + se_hits), width = 0.4) +
  coord_flip() +
  scale_y_continuous(labels = scales::label_comma(accuracy = 1)) +
  labs(
    x = "Species",
    y = "Mean Hits per 100,000"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.margin = margin(t = 10, r = 20, b = 10, l = 30),
    axis.line = element_line(color = "black", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

most_abundant_with_error

ggsave(filename = "Graphs/PHIbase_Top20_Bar_ErrorBars.pdf", plot = most_abundant_with_error, width = 10, height = 6)



# Top 10 species stacked bar  ----
top_10_species <- all_heatmap_filtered %>%
  group_by(species) %>%
  summarise(total_hits = sum(hits_per_100.000, na.rm = TRUE)) %>%
  arrange(desc(total_hits)) %>%
  slice_max(total_hits, n = 10) %>% 
  pull(species)

top_species_data <- all_heatmap_filtered %>%
  filter(species %in% top_10_species)

top_species_data$species <- 
  factor(top_species_data$species, levels = top_10_species)


# Specific colours
species_colors <- setNames(
  RColorBrewer::brewer.pal(10, "Set3")[1:length(top_10_species)],
  top_10_species
)

# Function to plot stacked top 10 species, with a grouping variable 
plot_top_species_stacked <- function(data, 
                                     group_var, 
                                     x_label = "Group",
                                     y_label = "Total Hits per 100,000",
                                     save_path = NULL,
                                     width = 10,
                                     height = 8) {
  group_var <- rlang::sym(group_var)
  
  # Summarise hits by group (month or location)
  summary_data <- data %>%
    group_by(!!group_var, species) %>%
    summarise(total_hits = sum(hits_per_100.000, na.rm = TRUE), .groups = "drop")
  
  # Order species by total abundance for stacked bar order
  species_order <- summary_data %>%
    group_by(species) %>%
    summarise(total = sum(total_hits)) %>%
    arrange(desc(total)) %>%
    pull(species)
  
  summary_data$species <- factor(summary_data$species, levels = species_order)
  
  # Plot
 p <-  ggplot(summary_data, aes(x = !!group_var, y = total_hits, fill = species)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = species_colors, name = "Species") +
    scale_y_continuous(
     labels = scales::label_comma(accuracy = 1)  # ensures no decimals
    ) +
    labs(
      x = x_label,
      y = y_label
    ) +
    theme_minimal(base_size = 14) +
    theme(
     plot.title = element_text(hjust = 0.5),
     axis.line = element_line(color = "black", linewidth = 0.3),
     panel.grid.minor = element_blank(),
     panel.grid.major = element_blank()
   )
 
 # Save plot if path provided
 if (!is.null(save_path)) {
   ggsave(filename = save_path, plot = p, width = width, height = height)
 }
 
 return(p)
 
}

#By Month 
#reorder month cols
top_species_data$month <- 
  factor(top_species_data$month,
         levels=c("Dec", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug"))

plot_top_species_stacked(
  data = top_species_data,
  group_var = "month",
  x_label = "Month",
  save_path = "Graphs/stacked_bar_by_month.pdf",
)

#By location
plot_top_species_stacked(
  data = top_species_data,
  group_var = "location",
  x_label = "Location",
  save_path = "Graphs/stacked_bar_by_location.pdf",
)


