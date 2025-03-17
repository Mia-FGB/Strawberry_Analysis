#This script brings together my Dec sequencing and other sequencing data 
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(dplyr)
library(lubridate)
library(scales)
library(patchwork)

#Prepping the data---------------------
#Summarised Feb - August data from PHIbase mapping
feb_aug_data <- read.csv("Tiptree_Feb_Aug_3_genus.csv", header=TRUE)

#Dec data, needs to be processed
raw_dec_data <- read.csv("Tiptree_Dec_unfiltered.csv", header = TRUE)

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

#Want to join the two datasets
all_data <- rbind(feb_aug_data, dec_data)

#need the date to be ordered
all_data$Date.collected = as.Date(all_data$Date.collected, format="%d/%m/%Y")
feb_aug_data$Date.collected = as.Date(feb_aug_data$Date.collected, format="%d/%m/%Y")

#Saved as all_sequence_data.csv

# Read number table ---------------------------
loc_read_num <- all_data %>% 
  group_by(location, Genus) %>% 
  summarise(mean_read_num = mean(genus_read), mean_percentage_reads = mean(percent))

mon_read_num <- all_data %>% 
  group_by(Date.collected, Genus) %>% 
  summarise(mean_read_num = mean(genus_read), mean_percentage_reads = mean(percent))

#Graphs ---------------------------------------------
#mean, min & max data for each location & genus for errorbar
data_min_max <- feb_aug_data %>% #Change depending in data being used
  group_by(Genus, location, Date.collected) %>% 
  summarise(mean = mean(log_hits), min = min(log_hits), max = max(log_hits))

#Specific for just Feb - Aug data
break.vec <- c(as.Date("2022-02-01"), seq(from=as.Date("2022-02-01"), to=as.Date("2022-08-31"),by="2 months")) 

#Sequencing data----------
# Creating the plot with consistent x-axes and labels on all facets
# genus_normalised_hits_error <- ggplot(data_min_max, aes(x = Date.collected, y = mean, group = location)) +
#   geom_line() + 
#   geom_point() +
#   geom_errorbar(aes(ymin = min, ymax = max), width = 7, colour = "grey39") +
#   facet_grid(facets = vars(Genus), cols = vars(location), scales ="free_x") + 
#   # Removing gridlines
#   theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
#                      strip.background = element_blank(),
#                      strip.placement = "outside",
#                      panel.spacing = unit(1, "lines")) +
#   xlab(label = "Date collected (Month)") + ylab(label = "Hits per 1000 (log)") +
#   # Using abbreviated month names
#   scale_x_date(breaks = break.vec, date_labels = "%b", limits = range(break.vec)) +
#   ggtitle("Hits per 1000 from sequencing")
# 
# # Print the plot
# print(genus_normalised_hits_error)

#ggsave(filename= "Images_Tiptree/genus_normalised_hits_error.pdf", plot = genus_normalised_hits_error, width=7, height=4)

#Making sep genus plots instead---

# Function to create plot for a specific genus
create_genus_plot <- function(genus_name) {
  p <- ggplot(data_min_max[data_min_max$Genus == genus_name, ], aes(x = Date.collected, y = mean, group = location)) +
      facet_grid(cols = vars(location), scales ="free_x") + 
      geom_line(size = 1.2) + 
      geom_point(size = 4) +
      geom_errorbar(aes(ymin = min, ymax = max), width = 10, colour = "grey39") +
      theme_bw() + 
      theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
      xlab("Date collected (Month)") + 
      ylab("Hits per 1000 (log)") +
      scale_x_date(breaks = break.vec, date_labels = "%b", limits = range(break.vec)) +
      ggtitle(paste("Hits per 1000 from sequencing - Genus:", genus_name))
  return(p)
}

# Creating plots for each genus
genus_list <- unique(data_min_max$Genus)

for (genus in genus_list) {
  plot <- create_genus_plot(genus)
  filename <- paste0("Images_Tiptree/Not_Dec_Data/hits_plot_0524_", genus, ".pdf")
  ggsave(filename, plot = plot, width = 24, height = 5, units = "in")
}

#Disease score data --------
disease_score <- read.csv("tiptree_disease_data.csv", header=TRUE)

#need the date to be ordered
disease_score$date = as.Date(disease_score$date, format="%d/%m/%Y")
#filter out N/A
filter_disease_score <- disease_score %>% 
  filter(score != 'N/A') %>% 
  subset(select = -c(total_reads, read_count))

# genus_disease_score <- ggplot(filter_disease_score, aes(x = date, y = score, group = location)) +
#   geom_point()+
#   facet_grid(facets =  vars(species), cols= vars(location)) +
#   #Removing gridlines
#   theme_bw() + theme( panel.grid.minor = element_blank()) +
#   xlab(label = "Month") + ylab(label = "Disease score")  +
#   scale_x_date(breaks = break.vec, date_labels = "%m", limits=range(break.vec)) +
#   ggtitle("Disease score")
# 
# ggsave(filename= "Images_Tiptree/disease_score_tiptree.pdf", plot = genus_disease_score, width=7, height=4)

#Individual Genus plots ----

# Function to create plot for a specific species
create_diseasesc_plot <- function(species_name) {
  # Filter data for the given species
  species_data <- filter(filter_disease_score, species == species_name)
  
  p <- ggplot(species_data, aes(x = date, y = score, group = location)) +
    facet_grid(cols = vars(location)) + 
    geom_line(size = 1.2) + 
    geom_point(size = 4) +
    theme_bw() + 
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    xlab("Date collected (Month)") + 
    ylab("Disease Score") +
    scale_x_date(breaks = break.vec, date_labels = "%b", limits = range(break.vec)) +
    ggtitle(paste("Disease Score - Species:", species_name))
  return(p)
}


# Creating and saving plots for each species
species_list <- unique(filter_disease_score$species)

for (species in species_list) {
  plot <- create_diseasesc_plot(species)
  filename <- paste0("Images_Tiptree/Not_Dec_Data/disease_score_0524", species, ".pdf")
  ggsave(filename, plot = plot, width = 24, height = 5, units = "in")
}

#Spraying ------------------
spraying <- read.csv("Tiptree_spraying.csv", header=TRUE)

#need the date to be ordered
spraying$Date = as.Date(spraying$Date, format="%d/%m/%Y")
na_remove_spraying <- spraying %>% 
  filter(protection != 'N/A')

#3_genus_scaled_spraying
spraying_dates <-  ggplot(na_remove_spraying, aes(x = Date, y = protection, colour = protection)) +
  geom_point(shape = 4, size = 6)+
  facet_grid(cols= vars(Location)) + theme_bw() +
  xlab(label = "Month") + ylab(label = "Fungicide") + 
  #to make it a square and match the others
  scale_x_date(breaks = break.vec, date_labels = "%b", limits=range(break.vec)) +
  ggtitle("Spraying dates")

ggsave(filename= "Images_Tiptree/Not_Dec_Data/spraying_dates.pdf", plot = spraying_dates, width=24, height=5)

#Temperature and hum--------------
Temp_Hum <- read.csv("Tiptree_Temp_Hum.csv", header=TRUE, na.strings=c(""," "))

#Use a class to work with the datetime easily - not working yet with the diff date formats
Temp_Hum$X.Date <- as.POSIXct(Temp_Hum$X.Date,format = "%d-%m-%Y %H:%M:%OS")
Temp_Hum$Humidity <- as.numeric(Temp_Hum$Humidity)
Temp_Hum$Temp <- as.numeric(Temp_Hum$Temp)
#removing some weird rows
Temp_Hum <- Temp_Hum %>% filter(!row_number() %in% c(60465, 26327, 13096, 47234))

Temp_Hum_Avg <- Temp_Hum %>% 
  mutate (date_only = date(X.Date)) %>% 
  group_by(date_only, Location) %>% 
  summarise(mean_temp = mean(Temp,na.rm=TRUE), sd_temp = sd(Temp, na.rm = TRUE),
            mean_hum = mean(Humidity, na.rm=TRUE), sd_hum = sd(Humidity, na.rm = TRUE))

Temp_all_graph <-  ggplot(Temp_Hum_Avg, aes(x = date_only, y = mean_temp)) +
  geom_ribbon(aes(ymin = mean_temp - sd_temp, ymax = mean_temp + sd_temp), fill = "grey70") +
  geom_line(colour = "black")+
  #facet_grid(cols= vars(Location)) +
  facet_grid(~factor(Location, levels=c('Field', 'NGS2', 'NGS1'))) +
  xlab(label = "Month") + ylab(label = "Temperature (°C)") + 
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  #to make it a square and match the others
  scale_x_date(breaks = break.vec, date_labels = "%b", limits=range(break.vec)) 

#ggsave(filename= "Images_Tiptree/all_temp.pdf", plot = Temp_all_graph, width=7, height=2)
ggsave(filename= "../Presntation\ graphics/PPT_Not_Dec_Data/all_temp.pdf", plot = Temp_all_graph, width=24, height=5)

hum_all_graph <-  ggplot(Temp_Hum_Avg, aes(x = date_only, y = mean_hum)) +
  geom_line(colour = "black")+
  geom_errorbar(aes(ymin = mean_hum - sd_hum, ymax = mean_hum + sd_hum), width = 10, colour = "grey39") +
  #facet_grid(cols= vars(Location)) +
  facet_grid(~factor(Location, levels=c('Field', 'NGS2', 'NGS1'))) +
  xlab(label = "Month") + ylab(label = "Humidity (%)") + 
  theme_bw() + theme( panel.grid.minor = element_blank()) +
  #to make it a square and match the others
  scale_x_date(breaks = break.vec, date_labels = "%m", limits=range(break.vec)) 

ggsave(filename= "Images_Tiptree/all_hum.pdf", plot = hum_all_graph, width=7, height=2)

#Heatmaps------------------------------------------------------------
#Graph Splitting all the samples to own column
#reading in the data
raw_f_a <- read.csv("Q20_Tiptree_combined_data.csv", header=TRUE)
filter_f_a <- raw_f_a %>% 
  select(taxaID, read_count, species, barcode, total_reads, Sample, Date.collected, repeat., location)

#Want to join the two datasets
all_heatmap_data <- rbind(filter_f_a , raw_dec_data)

#need the date to be ordered
all_heatmap_data$Date.collected = as.Date(all_heatmap_data$Date.collected, format="%d/%m/%Y")

#Need to normalise the reads for each location, e.g. hits per 10,000 reads
#Therefore we need the total number of reads for each sample, not just the classified reads
all_heatmap_filtered <- 
  all_heatmap_data %>% 
  #filter out any species with less than 10 reads
  filter(read_count >=10 ) %>% 
  mutate(hits_per_100.000 = (read_count * 100000)/(total_reads)) %>% 
  #log scale-Mutate the data to add a new row
  mutate(log_hits = log10(hits_per_100.000)) %>% 
  #don't want sample to be number but chr 
  mutate(Sample = as.character(Sample))

#filtering out taxa which are only present in <15 samples
all_heatmap_filtered_15taxa <- all_heatmap_filtered %>% 
  group_by(species) %>% 
  filter(n()>=15)

#Graph
#filtered out any species with less than 10 reads per sample
#filtered out species which were only present in 15 samples or less
all_heatmap <- ggplot(data = all_heatmap_filtered_15taxa, mapping = aes(x = Sample,
                                                       y = species,
                                                       fill = log_hits)) +
  geom_tile() +
  xlab(label = "Sample") + ylab(label = "Species") +
  scale_fill_distiller(name = "Hits per 100,000", palette = "RdBu") +
  theme_bw() +
  ggtitle(label = "Tiptree PHIbase Species Abundance (log hits per 100,000 > 10)") +
  #Editing axis title size
  theme(axis.text.y = element_text(vjust = 0.5, hjust=1, size = 22,)) +
  theme(axis.title.y = element_text(size = 22, vjust = 0.5)) +
  theme(axis.title.x = element_text(size = 22, hjust = 0.5)) +
  theme(axis.text.x = element_text(vjust=0.5, size = 15,)) +
  # Change legend key size and key width
  theme(legend.key.size = unit(6, "cm"),legend.key.width = unit(1.5,"cm"))+
  #Change legend font size
  theme(legend.title = element_text(size = 15)) +
  theme(legend.text = element_text(size =15)) +
  #Changing the title size
  theme(plot.title = element_text(size = 25, hjust = 0.5))

all_heatmap
#Saving the plot to right dimensions, much taller y axis
ggsave(filename= "Images_Tiptree/Tiptree_heatmap_sample_all_data_2504.png", plot =all_heatmap, width=30, height=40)

#Date heatmap-------------
#lots of empty spaces so going to group into months
#Adding a col for month
all_heatmap_filtered$month <- format(as.Date(all_heatmap_filtered$Date.collected, format="%d/%m/%Y"),"%b")

#filtering out taxa which are only present in <10 samples
all_heatmap_filtered_10taxa <- all_heatmap_filtered %>% 
  group_by(species) %>% 
  filter(n()>=10)

#reorder cols
all_heatmap_filtered_10taxa$month <- factor(all_heatmap_filtered_10taxa$month, levels=c("Aug", "Jul", "Jun", "May", "Apr", "Mar", "Feb", "Dec"))

all_date_heatmap <- ggplot(data = all_heatmap_filtered_10taxa, mapping = aes(x = species,
                                                            y = month,
                                                            fill = log_hits)) +
  geom_tile() +
  xlab(label = "Species") + ylab(label = "Month") +
  scale_fill_distiller(name = "Hits per 100,000", palette = "RdBu") +
  theme_bw() +
  ggtitle(label = "Tiptree PHIbase Species Abundance grouped by collection month (log hits per 100,000 > 10)") +
  #Editing axis title size
  theme(axis.text.y = element_text(vjust = 0.5, hjust=1, size = 20)) +
  theme(axis.title.y = element_text(size = 20, vjust = 0.5)) +
  theme(axis.title.x = element_text(size = 20, hjust = 0.5)) +
  theme(axis.text.x = element_text(vjust=0.5, hjust = 1, size = 20, angle = 90)) +
  # Change legend key size and key width
  theme(legend.key.size = unit(6, "cm"),legend.key.width = unit(1.5,"cm"))+
  #Change legend font size
  theme(legend.title = element_text(size = 15)) +
  theme(legend.text = element_text(size =15)) +
  #Changing the title size
  theme(plot.title = element_text(size = 25, hjust = 0.5))

all_date_heatmap
#Saving the plot to right dimensions, much taller y axis
ggsave(filename= "Images_Tiptree/Tiptree_heatmap_month_all_data_2504.png", plot = all_date_heatmap, width=45, height=20)

#Old----------
A_NGS1$hum <- as.numeric(A_NGS1$hum)
#Setting lims for the graph
lims <- as.POSIXct(strptime(c("2022-12-01 00:00:00","2022-09-01 00:00:00"), format = "%Y-%m-%d %H:%M"))

#Daily Average
A_NGS1_avg <- A_NGS1 %>% 
  #using lubridate package to extract just the date (ignoring the time)
  mutate (date_only = date(Date)) %>% 
  group_by(date_only) %>% 
  summarise(mean_temp = mean(Temp,na.rm=TRUE), sd_temp = sd(Temp, na.rm = TRUE),
            mean_hum = mean(hum, na.rm=TRUE), sd_hum = sd(hum, na.rm = TRUE))

#Avg Temp NGS1
A_NGS1_Avg_Temp <- ggplot(A_NGS1_avg, aes(x = date_only, y = mean_temp)) + geom_line(stat = "identity") +
  geom_errorbar(aes(ymin = mean_temp - sd_temp, ymax = mean_temp + sd_temp), color = "grey 30", alpha = .5) +
  xlab(label = "Month") + ylab(label = "Temperature") + 
  #to make it a square and match the others
  scale_x_date(breaks = break.vec, date_labels = "%m", limits=range(break.vec)) +
  theme_bw() +  theme( panel.grid.minor = element_blank()) +
  ylim(-10,40) +
  ggtitle("NGS1")

#Avg hum NGS1
A_NGS1_Avg_Hum <- ggplot(A_NGS1_avg, aes(x = date_only, y = mean_hum)) + geom_line(stat = "identity") +
  geom_errorbar(aes(ymin = mean_hum - sd_hum, ymax = mean_hum + sd_hum), color = "grey 30", alpha = .5) +
  xlab(label = "Month") + ylab(label = "hum (%)") + 
  #to make it a square and match the others
  scale_x_date(breaks = break.vec, date_labels = "%m", limits=range(break.vec)) +
  theme_bw() + theme( panel.grid.minor = element_blank())  +
  ylim(30,100) 




