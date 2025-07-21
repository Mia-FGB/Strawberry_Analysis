# Script to plot heatmaps of the PHIbase alignment data from Tiptree

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
all_heatmap_data <- rbind(filter_f_a , raw_dec_data) %>% 
  rename(Species = species) 

# Process the date column
all_heatmap_data$Date.collected = as.Date(all_heatmap_data$Date.collected, format="%d/%m/%Y")

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

#Preparing filtered data ----

#Adding a col for month
all_heatmap_filtered$month <- 
  format(as.Date(all_heatmap_filtered$Date.collected, format="%d/%m/%Y"),"%b")

#filtering out taxa which are only present in <10 samples - less stringent when grouped
all_heatmap_filtered_10taxa <-
  all_heatmap_filtered %>% 
  group_by(Species) %>% 
  filter(n()>=10)

# filtering out taxa which are only present in <15 samples - for all the samples to be plotted seperate
all_heatmap_filtered_15taxa <- all_heatmap_filtered %>% 
  group_by(Species) %>% 
  filter(n()>=15)

#reorder month cols
all_heatmap_filtered_10taxa$month <- 
  factor(all_heatmap_filtered_10taxa$month,
         # levels=c("Dec", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug"))
         levels=c("Aug", "Jul", "Jun", "May", "Apr", "Mar", "Feb", "Dec"))

# Change labels and order
all_heatmap_filtered_10taxa <- all_heatmap_filtered_10taxa %>%
  mutate(location = factor(location,
                           levels = c("F", "RG", "LG"), 
                           labels = c("Field", "Greenhouse 1", "Greenhouse 2")))



# Plotting ---
# Heatmap function:
plot_species_heatmap <- function(data,
                                 x = "Sample",
                                 y = "Species",
                                 xlab = x,
                                 fill = "log_hits",
                                 title = "Species Abundance Heatmap",
                                 save_path = NULL,
                                 width = 12,
                                 height = 18,
                                 theme_extra = NULL,
                                 flip = FALSE) {
  
  # Compute abundance
  species_abundance <- data %>%
    group_by(Species) %>%
    summarise(total_abundance = sum(log_hits, na.rm = TRUE)) %>%
    arrange(desc(total_abundance))
  
  # Apply the order to the main data (most abundant species at the bottom)
  data <- data %>%
    mutate(Species = factor(Species, levels = species_abundance$Species))
  
  # Swap x/y if flip = TRUE
  if (flip) {
    aes_mapping <- aes_string(x = y, y = x, fill = fill)
  } else {
    aes_mapping <- aes_string(x = x, y = y, fill = fill)
  }
  
  # Create plot
  p <- ggplot(data, mapping = aes_mapping) +
    geom_tile() +
    
    scale_fill_viridis_c(
      name = "Hits per 100,000",
      labels = function(x) round(10^x),
      option = "C",  # Options: "D", "C", "B", etc.
      direction = 1
    ) +
  
    theme_bw() +
    labs(
      title = title,
      x = ifelse(flip, y, xlab),
      y = ifelse(flip, xlab, y)
    ) +
    theme(
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 16),
      plot.title = element_text(size = 18, hjust = 0.5),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.key.size = unit(1.5, "cm"),
      legend.key.width = unit(0.4, "cm")
    ) +
    theme_extra
  
  # Save plot if path provided
  if (!is.null(save_path)) {
    ggsave(filename = save_path, plot = p, width = width, height = height)
  }
  
  return(p)
}


# Call the function ---
# All samples-
plot_species_heatmap(
  data = all_heatmap_filtered_15taxa,
  x = "Sample",
  y = "Species",
  title = "",
  save_path = "Graphs/Heatmaps/heatmap_by_sample.pdf",
  width = 18
)


# Group by month - Vertical plot
plot_species_heatmap(
  data = all_heatmap_filtered_10taxa,
  x = "month",
  y = "Species",
  xlab="Month",
  title = "",
  save_path = "Graphs/Heatmaps/heatmap_by_month.pdf",
  height = 18,
  width = 8
)

# Group by month - Horizonta plot
plot_species_heatmap(
  data = all_heatmap_filtered_10taxa,
  x = "month",
  y = "Species",
  xlab="Month",
  title = "",
  save_path = "Graphs/Heatmaps/heatmap_by_month_rot.pdf",
  flip = TRUE,
  height = 6,
  width = 18,
  theme_extra = theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
)

# Plot by location - vertical 
plot_species_heatmap(
  data = all_heatmap_filtered_10taxa,
  x = "location",
  y = "Species",
  xlab = "Location",
  title = "",
  save_path = "Graphs/Heatmaps/heatmap_by_location.pdf",
  height = 18,
  width = 12
)

# Plot by location - horizontal  
plot_species_heatmap(
  data = all_heatmap_filtered_10taxa,
  x = "location",
  y = "Species",
  xlab = "Location",
  title = "",
  save_path = "Graphs/Heatmaps/heatmap_by_location_rot.pdf",
  flip = TRUE,
  height = 5,
  width = 18,
  theme_extra = theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
)

 