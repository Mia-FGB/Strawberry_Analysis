# Looking at the fungicide data
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(lubridate)

setwd("~/OneDrive - Norwich BioScience Institutes/Air_Samples/Tiptree")


# Data 
data <-  read.csv("Metadata/Tiptree_spraying.csv")
frac <- read.csv("Metadata/fung_frac_risk.csv") #Sheet of Fungicide info

# Change Location column 
data <- data %>%
  mutate(Location = case_when(
    Location == "F" ~ "Field",
    Location == "LG" ~ "Greenhouse 2",
    Location == "RG" ~ "Greenhouse 1",
    TRUE ~ Location  # keep original if it doesn't match any of the above
  ))

# Looking at Protection type ------

# Prepare the data
spray_summary <- data %>%
  mutate(
    Sprayed = na_if(Sprayed, "N/A"),   # Convert "N/A" string to actual NA
    Date = dmy(Date)                   # Convert Date string to Date object
  ) %>%
  filter(!is.na(Sprayed)) %>%          # Filter out the NA Sprayed rows
  mutate(
    Month = floor_date(Date, "month")  # Group by month (1st day of each month)
  ) %>%
  group_by(Location, Month, protection) %>%
  summarise(Spray_Count = n(), .groups = "drop")

# Plotting --

plot <- ggplot(spray_summary, aes(x = Month, y = Spray_Count, fill = protection)) +
  geom_bar(stat = "identity",
           position = "stack",
          # position = "dodge". # To have bars side by side
           ) +
  facet_wrap(~Location, ncol = 3) +  # Each location gets its own panel
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  scale_fill_manual(
    values = c(
      "Botrytis" = "#E32E60",
      "Broad" = "#305182",
      "Podosphaera" = "#4FBDC5"
    )) +
  labs(
    # title = "Fungicide Sprays per Month by Location",
    x = "Month",
    y = "Number of Applications",
    fill = "Protection Type"
  ) + 
  theme_minimal(base_size=14) +
  theme( axis.line = element_line(color = "black", linewidth = 0.3),
         panel.grid.major = element_blank())  +
  scale_y_continuous(limits = c(0, 8), expand = c(0, 0))

plot

ggsave("Graphs/Fungicide/fungicide_sprays_protection.pdf", plot = plot, width = 15, height = 5)

# Looking at FRAC risk -----
# Join datasets (assuming `spray_data` and `frac_data`)
spray_risk_data <- data %>%
  mutate(
    Fungicide = Sprayed,
    Sprayed = na_if(Sprayed, "N/A"),
    Date = dmy(Date)
  ) %>%
  filter(!is.na(Sprayed)) %>%
  left_join(frac, by = "Fungicide") # # this expands sprays by chemical as some fungicides have multiple chemicals

# Define order
risk_levels <- c("N/A", "Low - Medium", "Medium", "Medium - High", "High")

spray_risk_data <- spray_risk_data %>%
  mutate(Resistance_Risk = factor(Resistance_Risk, levels = risk_levels))

# define colours - own palette
# risk_colors <- c(
#  "N/A" = "grey80",
#  "Low - Medium" = "#58A65C",         # green
#  "Medium" = "#EFC94C",      # yellow-ish
#  "Medium - High" = "#F29E4C", # orange
#  "High" = "#D64550"         # red
# )

# Get 4 Brewer colors 
risk_colors <- brewer.pal(4, "YlOrRd")
# Add grey for N/A
risk_colors <- c("N/A" = "grey70", 
                  setNames(risk_colors, risk_levels[-1]))  # assign brewer colors to remaining risk levels

# Plot 1 - Total sprays per FRAC risk
spray_frac <- spray_risk_data %>%
  count(Location, Resistance_Risk) %>%
  ggplot(aes(x = Resistance_Risk, y = n, fill = Resistance_Risk)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Location) +
  scale_fill_manual(values = risk_colors, na.translate = FALSE) +
  labs(
    x = "FRAC Resistance Risk",
    y = "Number of Applications",
    fill = "Resistance Risk",
   # title = "Fungicide Use by FRAC Risk and Location"
  ) +
  theme_minimal(base_size=14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color = "black", linewidth = 0.3),
        panel.grid.major = element_blank())  +
  scale_y_continuous(expand = c(0, 0))

spray_frac

ggsave("Graphs/Fungicide/fungicide_sprays_frac.pdf", plot = spray_frac, width = 14, height = 6)

# Plot 2 - Stacked (think I prefer this)

# Reverse order
spray_risk_data <- spray_risk_data %>%
  mutate(Resistance_Risk = factor(Resistance_Risk, levels = rev(risk_levels)))

stacked_plot <- spray_risk_data %>%
  count(Location, Resistance_Risk) %>%
  ggplot(aes(x = Location, y = n, fill = Resistance_Risk)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = risk_colors, na.translate = FALSE) +
  labs(
    x = "Location",
    y = "Number of Applications",
    fill = "FRAC Risk",
    # title = "Fungicide Sprays by Location and FRAC Risk"
  )  +
  theme_minimal(base_size=14) +
  theme(axis.line = element_line(color = "black", linewidth = 0.3),
        panel.grid.major = element_blank())  +
  scale_y_continuous(expand = c(0, 0))

ggsave("Graphs/Fungicide/fungicide_sprays_frac_stacked.pdf", plot = stacked_plot, width = 8, height = 5)
