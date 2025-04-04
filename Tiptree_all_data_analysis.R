#This script brings together my Dec sequencing and other sequencing data 
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(dplyr)
library(lubridate)
library(scales)
library(patchwork)

setwd("~/OneDrive - Norwich BioScience Institutes/Air_Samples/Tiptree")

#Prepping the data---------------------

#Summarised Feb - August data from PHIbase mapping
feb_aug_data <- read.csv("Tiptree_taxa_counts_initial_analysis/Tiptree_Feb_Aug_3_genus.csv", header=TRUE)
feb_aug_data$Date.collected = as.Date(feb_aug_data$Date.collected, format="%d/%m/%Y")

#Dec data, needs to be processed
raw_dec_data <- read.csv("Tiptree_taxa_counts_initial_analysis/Tiptree_Dec_unfiltered.csv", header = TRUE)

#Process the December data
dec_data <- raw_dec_data %>% 
  filter(grepl('Podosphaera|Botrytis|Phytophthora', species)) %>% 
  select(Sample, Date.collected, location, total_reads, read_count, species) %>% 
  #splitting group and species
  separate(species, c("Genus", "species"), sep = " ") %>% 
  #want to group by genus, location & date 
  group_by(Date.collected, location, Genus, Sample) %>% 
  #sum up the count data
  mutate(genus_read = sum(read_count)) %>% 
  #only need one row per genus for each sample
  distinct(Sample, .keep_all=TRUE) %>% 
  #no longer want the species col 
  subset(select = -(species)) %>% 
  #adding some normalised values
  mutate(hits_per_1000 = (genus_read * 1000)/(total_reads)) %>% 
  #log scale-Mutate the data to add a new row
  mutate(log_hits = log10(hits_per_1000)) %>% 
  #adding a row for percentage 
  mutate(percent = (genus_read/total_reads)*100) %>% 
  #no longer want the read_count as this was specific to the species
  subset(select = -(read_count)) 

# Change date formats
dec_data$Date.collected <- dmy(dec_data$Date.collected)

# Combine the datasets
all_data <- rbind(feb_aug_data, dec_data)

# Change the locations
standardise_location <- function(df) {
  # Check which column exists: "location" or "Location"
  loc_col <- if ("location" %in% names(df)) {
    "location"
  } else if ("Location" %in% names(df)) {
    "Location"
  } else {
    stop("No 'location' or 'Location' column found.")
  }
  
  # Rename to a working temp name, mutate, then rename back
  df <- df %>%
    rename(temp_loc = !!loc_col) %>%
    mutate(temp_loc = case_when(
      temp_loc == "F"  ~ "Field",
      temp_loc == "LG" ~ "Greenhouse 2",
      temp_loc == "RG" ~ "Greenhouse 1",
      TRUE ~ temp_loc
    )) %>%
    rename(!!loc_col := temp_loc)
  
  return(df)
}


all_data <- standardise_location(all_data)

# Setting up constants for all graphs------

# To keep consistent limits -
date_range <- range(all_data$Date.collected, na.rm = TRUE)
# Can cal this in all graphs
month_scale <- scale_x_date(
  date_breaks = "1 month",
  date_labels = "%b",
  limits = date_range
)

custom_theme <- theme_minimal(base_size = 14) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.3),
    panel.grid.major = element_blank()
  )

# Considering combining the plots with patchwork -----
genus_plots <- list()
disease_plots <- list()
temp_plots <- list()
humidity_plots <- list()
combined_plots <- list()


# Sequence data plots  ---------------------------------------------

#  - Script updated April 2025

#mean, min & max data for each location & genus for "errorbar"
data_min_max <- all_data %>%
  group_by(Genus, location, Date.collected) %>%
  summarise(
    mean = mean(log_hits, na.rm = TRUE),
    min = min(log_hits, na.rm = TRUE),
    max = max(log_hits, na.rm = TRUE),
    .groups = "drop"
  )

# Function to create plots for each genus - facetted by location 
create_genus_plot <- function(genus_name) {
  genus_data <- filter(data_min_max, Genus == genus_name)
  
  ggplot(genus_data, aes(x = Date.collected, y = mean)) +
    geom_point(size = 3) +
    geom_line(size = 0.8) +
    geom_errorbar(aes(ymin = min, ymax = max), width = 7, color = "grey50") +
    facet_wrap(~location, ncol = 3) +
    month_scale +
    labs(
      x = " ", # Want it empty so i can stack
      y = "Hits per 1000 (log)"
    ) +
    custom_theme +
    scale_y_continuous(expand = c(0, 0))
}

# To create seperate plots for each Genus
genus_list <- unique(data_min_max$Genus)
for (genus in genus_list) {
  plot <- create_genus_plot(genus)
  print(plot)
  filename <- paste0("Graphs/Sequence_Data/Line_0425_", genus, ".pdf")
  ggsave(filename, plot = plot, width = 24, height = 5, units = "in")
}

# Disease score data -------------
disease_score <- read.csv("Metadata/tiptree_disease_data.csv", header=TRUE)

#Change location 
disease_score <- standardise_location(disease_score)

# filter out N/A
filter_disease_score <- disease_score %>% 
  filter(score != 'N/A') %>% 
  subset(select = -c(total_reads, read_count))


# make the correct format
filter_disease_score$date = as.Date(filter_disease_score$date, format="%d/%m/%Y")
filter_disease_score$score = as.integer(filter_disease_score$score)


# Function to create and save plot for a specific species
disease_score_plot <- function(species_name, save_path = "Graphs/Disease_Score/", width = 24, height = 5) {
  # Filter data
  species_data <- filter(filter_disease_score, species == species_name)
  
  # Skip if no data
  if (nrow(species_data) == 0) return(NULL)
  
  # Create plot
  p <- ggplot(species_data, aes(x = date, y = score, group = location)) +
    geom_point(size = 3) +
    facet_grid(cols = vars(location)) +
    month_scale +
    labs(
      x = " ",
      y = "Disease Score",
    ) +
    custom_theme +
    scale_y_continuous(limits = c(-0.5,5))
  
  # Save plot
  filename <- paste0(save_path, "disease_score_0425_", species_name, ".pdf")
  ggsave(filename, plot = p, width = width, height = height, units = "in")
  
  return(p)
}

# Run function for the 3 species - to create seperate plots
for (genus in genus_list) {
  disease_score_plot(genus)
}


# Environmental Data --------------


Temp_Hum <- read.csv("Metadata/Tiptree_Temp_Hum.csv", header=TRUE, na.strings=c(""," "))

# Change locations 
Temp_Hum <- Temp_Hum %>%
  mutate(Location = case_when(
    Location == "NGS1" ~ "Greenhouse 1",
    Location == "NGS2" ~ "Greenhouse 2",
    Location == "Field" ~ "Field"
  )) %>%
  filter(!is.na(Location))  # Remove all other values (turned into NA)

# Use a class to work with the datetime easily - not working yet with the diff date formats
Temp_Hum$X.Date <- as.POSIXct(Temp_Hum$X.Date,format = "%d-%m-%Y %H:%M:%OS")
Temp_Hum$Humidity <- as.numeric(Temp_Hum$Humidity)
Temp_Hum$Temp <- as.numeric(Temp_Hum$Temp)

# Group & calculate daily averages 
Temp_Hum_Avg <- Temp_Hum %>%
  mutate(date_only = as.Date(X.Date)) %>%
  group_by(date_only, Location) %>%
  summarise(
    mean_temp = mean(Temp, na.rm = TRUE),
    sd_temp   = sd(Temp, na.rm = TRUE),
    mean_hum  = mean(Humidity, na.rm = TRUE),
    sd_hum    = sd(Humidity, na.rm = TRUE),
    .groups = "drop"
  )

# Temp Graph --
Temp_all_graph <-  ggplot(Temp_Hum_Avg, aes(x = date_only, y = mean_temp)) +
  geom_ribbon(aes(ymin = mean_temp - sd_temp, ymax = mean_temp + sd_temp),
              fill = "grey80") +
  geom_line(colour = "black")+
  facet_grid(~factor(Location, levels=c('Field', 'Greenhouse 1', 'Greenhouse 2'))) +
  xlab(label = "Month") + 
  ylab(label = "Temperature (°C)") + 
  month_scale +
  custom_theme

Temp_all_graph

ggsave(filename= "Graphs/Weather_Data/all_temp.pdf", plot = Temp_all_graph, width=24, height=5)

plot_env_variable <- function(data, x_var, y_mean, y_sd, y_label, output_file = NULL) {
  p <- ggplot(data, aes(x = {{ x_var }}, y = {{ y_mean }})) +
    geom_ribbon(aes(ymin = {{ y_mean }} - {{ y_sd }}, ymax = {{ y_mean }} + {{ y_sd }}),
                fill = "grey80") +
    geom_line(color = "black") +
    facet_grid(~factor(Location, levels = c("Field", "Greenhouse 1", "Greenhouse 2"))) +
    xlab("Month") +
    ylab(y_label) +
    month_scale +
    custom_theme
  
  if (!is.null(output_file)) {
    ggsave(filename = output_file, plot = p, width = 24, height = 5, units = "in")
  }
  
  return(p)
}

# Temperature plot
plot_env_variable(
  data = Temp_Hum_Avg,
  x_var = date_only,
  y_mean = mean_temp,
  y_sd = sd_temp,
  y_label = "Temperature (°C)",
  output_file = "Graphs/Weather_Data/all_temp.pdf"
)

# Humidity plot
plot_env_variable(
  data = Temp_Hum_Avg,
  x_var = date_only,
  y_mean = mean_hum,
  y_sd = sd_hum,
  y_label = "Humidity (%)",
  output_file = "Graphs/Weather_Data/all_humidity.pdf"
)

# Attempting to combine the plots into one rather than individually ---

for (genus in genus_list) {
  # Sequence plot (keep Y label, remove X)
  genus_plot <- create_genus_plot(genus) +
    xlab(NULL)
  
  # Disease score plot (remove title and X label)
  disease_plot <- disease_score_plot(genus, save_path = NULL) +
    labs(title = NULL) +
    xlab(NULL) +
    theme(strip.text = element_blank())
  
  # Temp plot (remove X label only)
  temp_plot <- plot_env_variable(
    Temp_Hum_Avg, date_only, mean_temp, sd_temp,
    "Temperature (°C)", output_file = NULL
  ) +
    xlab(NULL) +
    theme(strip.text = element_blank())
  
  # Humidity plot (keep X label as final one)
  humidity_plot <- plot_env_variable(
    Temp_Hum_Avg, date_only, mean_hum, sd_hum,
    "Humidity (%)", output_file = NULL
  ) +
    theme(strip.text = element_blank())
  
  # Combine all
  combined_plot <- genus_plot /
    disease_plot /
    temp_plot /
    humidity_plot +
    plot_annotation(title = paste("Combined Data for Genus:", genus))
  
  # Save it
  ggsave(
    filename = paste0("Graphs/Combined_Genus/combined_", genus, ".pdf"),
    plot = combined_plot,
    width = 14,
    height = 16,
    units = "in"
  )
}


# Fungicide Spraying  ------------------

spraying <- read.csv("Metadata/Tiptree_spraying.csv", header=TRUE)

# Need the date to be ordered
spraying$Date = as.Date(spraying$Date, format="%d/%m/%Y")

spraying <- standardise_location(spraying)

#Remove NA
na_remove_spraying <- spraying %>% 
  filter(protection != 'N/A')

# Plot to show the dates of the fungicide applications in different locations 
spraying_dates <-  
  ggplot(na_remove_spraying, aes(x = Date, y = protection, colour = protection)) +
  geom_point(shape = 4, size = 6) +
  facet_grid(cols= vars(Location)) +
  xlab(label = "Month") + 
  ylab(label = "Fungicide") + 
  month_scale +
  ggtitle("Spraying dates") +
  custom_theme

spraying_dates

ggsave(filename= "Graphs/Fungicide/application_dates.pdf", plot = spraying_dates, width=24, height=5)

