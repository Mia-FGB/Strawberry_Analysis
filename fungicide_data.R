# Looking at the fungicide data
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(lubridate)
library(stringr)

setwd("~/OneDrive - Norwich BioScience Institutes/Air_Samples/Tiptree")


# Data 
data <-  read.csv("Scripts/Supplementary_Data/spraying.csv") %>% 
  mutate(Fungicide = str_remove(Sprayed, "\\s*\\(.*\\)")) %>% 
  select(where(~ !all(is.na(.))))     # drop NA-only columns

#Fungicide info
frac <- read.csv("Scripts/Supplementary_Data/fung_frac_risk.csv") #Sheet of Fungicide info

#Merge on frac to get the protection info from fungicides
merged <- data %>%
  left_join(frac, by = "Fungicide") 

# Remove NA and Insecticide
merged <- merged %>%
  filter(!is.na(Target) & !(Target %in% c("n/a", "Insecticide"))) 


# Looking at Protection type ------

#The merge led to extra rows per fungicide as some have multiple active chemicals
spray_deduplicated <- merged %>%
  select(Fungicide, Location, Date, Target) %>% 
  distinct()  # Remove duplicated rows

# Prepare the data
spray_summary <- spray_deduplicated %>%
  mutate(Date = dmy(Date)) %>%              # Convert Date string to Date object 
  mutate(
    Month = floor_date(Date, "month")        # Group by month (1st day of each month)
  ) %>%
  group_by(Location, Month, Target) %>%
  summarise(Spray_Count = n(), .groups = "drop")

# Plotting --

spray_summary$Target <- factor(spray_summary$Target,
                               levels = c("Botrytis & Podosphaera", "Botrytis", "Podosphaera"))


plot <- ggplot(spray_summary, aes(x = Month, y = Spray_Count, fill = Target)) +
  geom_bar(stat = "identity",
           position = "stack",
          # position = "dodge". # To have bars side by side
           ) +
  facet_wrap(~Location, ncol = 3) +  # Each location gets its own panel
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  scale_fill_manual(
    values = c(
      "Botrytis & Podosphaera" = "#E32E60",
      "Botrytis" = "#305182",
      "Podosphaera" = "#4FBDC5"
    )) +
  labs(
    # title = "Fungicide Sprays per Month by Location",
    x = "Month",
    y = "Number of Applications",
    fill = "Disease Target"
  ) + 
  theme_minimal(base_size=14) +
  theme( axis.line = element_line(color = "black", linewidth = 0.3),
         panel.grid.major = element_blank())  +
  scale_y_continuous(limits = c(0, 8), expand = c(0, 0))

plot

ggsave("Graphs/Fungicide/fungicide_sprays_protection.pdf", plot = plot, width = 15, height = 5)

# Looking at FRAC risk -----------

merged <- merged %>%
  mutate(Date = dmy(Date))

# Define order
risk_levels <- c("N/A", "Low - Medium", "Medium", "Medium - High", "High")

merged <- merged %>%
  mutate(resistance_risk = factor(resistance_risk, levels = risk_levels))

# colours --- 4 Brewer colors 
risk_colors <- brewer.pal(4, "YlOrRd")
# Add grey for N/A
risk_colors <- c("N/A" = "grey70", 
                  setNames(risk_colors, risk_levels[-1]))  # assign brewer colors to remaining risk levels

# Total sprays per FRAC risk plot --------
spray_frac <- merged %>%
  count(Location, resistance_risk) %>%
  ggplot(aes(x = resistance_risk, y = n, fill = resistance_risk)) +
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
  scale_y_continuous(expand = c(0, 0), breaks = scales::pretty_breaks(n = 5))

spray_frac

ggsave("Graphs/Fungicide/fungicide_sprays_frac.pdf", plot = spray_frac, width = 14, height = 6)

# Stacked FRAC resistance plot ----------
# (think I prefer this)

# Reverse order
merged <- merged %>%
  mutate(resistance_risk = factor(resistance_risk, levels = rev(risk_levels)))

stacked_plot <- merged %>%
  count(Location, resistance_risk) %>%
  ggplot(aes(x = Location, y = n, fill = resistance_risk)) +
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

stacked_plot

ggsave("Graphs/Fungicide/fungicide_sprays_frac_stacked.pdf", plot = stacked_plot, width = 8, height = 5)
