#This script brings together my Dec sequencing and other sequencing data 
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(dplyr)
library(lubridate)
library(scales)
library(patchwork)
library(stringr)

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
seq_data <- rbind(feb_aug_data, dec_data)

# Adding a hits per 100k column
seq_data <- seq_data %>%
  mutate(hits_per_100000 = hits_per_1000 * 100) %>%
  mutate(log_hits_100k = log10(hits_per_100000))

# Change Location col
seq_data <- seq_data %>%
  rename(Location = location) %>%  # Rename the column
  mutate(Location = case_when(
    Location == "RG" ~ "Greenhouse 1",
    Location == "LG" ~ "Greenhouse 2",
    Location == "F"  ~ "Field",
  ))

# Setting up constants for all graphs------

# To keep consistent limits -
date_range <- range(seq_data$Date.collected, na.rm = TRUE)
date_range <- c(date_range[1] - days(7), date_range[2] + days(7)) # To keep eror bar caps

#Consistent breaks 
month_scale <- scale_x_date(
  date_breaks = "1 month",
  date_labels = "%b",
  limits = date_range
)

# Same theme 
custom_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )


# Sequence data plots  ---------------------------------------------

#mean, min & max data for each location & genus for "errorbar"
data_min_max <- seq_data %>%
  group_by(Genus, Location, Date.collected) %>%
  summarise(
    mean = mean(log_hits_100k, na.rm = TRUE),
    min = min(log_hits_100k, na.rm = TRUE),
    max = max(log_hits_100k, na.rm = TRUE),
    .groups = "drop"
  )

# Function to create plots for each genus - faceted by location 
create_genus_plot <- function(genus_name) {
  genus_data <- filter(data_min_max, Genus == genus_name)
  
  ggplot(genus_data, aes(x = Date.collected, y = mean)) +
    # geom_point(size = 3) +
    geom_line(size = 0.8) +
    geom_errorbar(aes(ymin = min, ymax = max), linewidth = 0.6, width = 5, color = "grey50") +
    scale_y_continuous(
      labels = function(x) scales::comma(10^x), # So they are actual values but plotted on log scale
      expand = c(0, 0)
    ) +
    facet_wrap(~Location, ncol = 3) +
    labs(
      x = " ", # Want it empty so i can stack
      y = "Hits per 100,000" # As the axis labels are converted
    ) +
    month_scale +
    custom_theme
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
disease_score <- read.csv("Metadata/disease_score.csv", header=TRUE)

# Convert na values to 0
disease_score$Score[is.na(disease_score$Score)] <- 0

# make the correct format
disease_score$Date = as.Date(disease_score$Date, format="%d/%m/%Y")

# Function to create and save plot for a specific species
disease_score_plot <- function(species_name, save_path = "Graphs/Disease_Score/", width = 24, height = 5) {
  # Filter data
  species_data <- filter(disease_score, Disease == species_name)
  
  # Skip if no data
  if (nrow(species_data) == 0) return(NULL)
  
  # Create plot
  p <- ggplot(species_data, aes(x = Date, y = Score, group = Location)) +
    geom_point(size = 3) +
    facet_grid(cols = vars(Location)) +
    month_scale +
    labs(
      x = " ",
      y = "Disease Score",
    ) +
    custom_theme +
    scale_y_continuous(limits = c(-0.5,5))
  
  print(p)
  
  # Save plot
  filename <- paste0(save_path, "disease_score_0425_", species_name, ".pdf")
  ggsave(filename, plot = p, width = width, height = height, units = "in")
  
  return(p)
}

# Run function for the 3 species - to create separate plots
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

#Function for environmental data plotting 

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


# Fungicide Spraying  ------------------

# Data Prep 
spraying <-  read.csv("Metadata/spraying.csv") %>% 
  mutate(Fungicide = str_remove(Sprayed, "\\s*\\(.*\\)")) %>% 
  select(where(~ !all(is.na(.))))     # drop NA-only columns

#Fungicide info
frac <- read.csv("Metadata/fung_frac_risk.csv") #Sheet of Fungicide info

#Merge on frac to get the protection info from fungicides
spraying <- spraying %>%
  left_join(frac, by = "Fungicide") 

# Remove NA and Insecticide
spraying <- spraying %>%
  filter(!is.na(Target) & !(Target %in% c("n/a", "Insecticide"))) 

spray_deduplicated <- spraying %>%
  select(Fungicide, Location, Date, Target) %>% 
  distinct()  # Remove duplicated rows

spray_deduplicated$Date <- as.Date(spray_deduplicated$Date, format="%d/%m/%Y")

# Plot Fungicide just to see -------
# spraying_dates <-  
#   ggplot(spray_deduplicated, aes(x = Date, y = Target, colour = Target)) +
#   geom_point(shape = 4, size = 6) +
#   facet_grid(cols= vars(Location)) +
#   xlab(label = "Month") + 
#   ylab(label = "Fungicide") + 
#   ggtitle("Spraying dates") +
#   month_scale +
#   custom_theme
# 
# spraying_dates
# 
# ggsave(filename= "Graphs/Fungicide/application_dates.pdf", plot = spraying_dates, width=24, height=5)

# Setting up the lines to use in combined plots ----
spray_deduplicated$Target <- factor(spray_deduplicated$Target,
                               levels = c("Botrytis & Podosphaera", "Botrytis", "Podosphaera"))

spray_colours <- c(
  "Botrytis & Podosphaera" = "#E32E60",
  "Botrytis" = "#305182",
  "Podosphaera" = "#4FBDC5"
)

# Spray lines with color mapping - change the aesthetics here
spray_lines <- list(
  geom_vline(
    data = spray_deduplicated,
    aes(xintercept = Date, color = Target),
   # linetype = "dashed",
     linewidth = 0.6,
     alpha = 0.6
  ),
  scale_color_manual(values = spray_colours, name = "Fungicide Target") # to add a legend 
)

# Combined Plots -------------


genus_list <- unique(data_min_max$Genus)

for (genus in genus_list) {
  # Sequence plot (keep Y label, remove X)
  genus_plot <- create_genus_plot(genus) +
    xlab(NULL) +
    spray_lines +
    theme(legend.key.height = unit(2, "lines")) +
    guides(color = guide_legend(override.aes = list(linewidth = 2)))
  
  # Disease score plot (remove title and X label)
  disease_plot <- disease_score_plot(genus, save_path = NULL) +
    labs(title = NULL) +
    xlab(NULL) +
    spray_lines +
    theme(legend.position = "none") + # So the legend isn't repeated
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
  
  # Adding "No data" to the Field Humidity plot
  humidity_plot <- humidity_plot + 
    geom_text(
      data = Temp_Hum_Avg %>% filter(Location == "Field") %>% slice(1),  # one row per facet
      aes(x = mean(date_range), y = 65),  # adjust position
      label = "No data", color = "grey60", size = 8
    )
  
  # Combine all plots
  combined_plot <- (genus_plot /
    disease_plot /
    temp_plot /
    humidity_plot )

  filename <- paste0("Graphs/Combined_Genus/combined_no_title_", genus, ".pdf")
  
  ggsave(
    filename = filename,
    plot = combined_plot,
    width = 14,
    height = 16,
    units = "in"
  )
}





