library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(dplyr)
library(lubridate)
library(scales)
library(patchwork)

setwd("~/OneDrive - Norwich BioScience Institutes/Air_Samples/Tiptree")


#Heatmaps------------------------------------------------------------


#Graph Splitting all the samples to own column
#reading in the data
raw_f_a <- read.csv("Tiptree_taxa_counts_initial_analysis/Q20_Tiptree_combined_data.csv", header=TRUE)
filter_f_a <- raw_f_a %>% 
  select(taxaID, read_count, species, barcode, total_reads, Sample, Date.collected, repeat., location)

#Dec data
raw_dec_data <- read.csv("Tiptree_taxa_counts_initial_analysis/Tiptree_Dec_unfiltered.csv", header = TRUE)

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

# Same as above but ordering species by abundance and rotating 02.04.25 -------------
# Calculate total abundance per species across all months
species_abundance <- all_heatmap_filtered_10taxa %>%
  group_by(species) %>%
  summarise(total_abundance = sum(log_hits))

# Sort species by total abundance in descending order
species_abundance <- species_abundance %>%
  arrange((total_abundance))

# Reorder species in the heatmap plot
all_heatmap_filtered_10taxa$species <- factor(all_heatmap_filtered_10taxa$species, 
                                              levels = species_abundance$species)

#reorder month cols
all_heatmap_filtered_10taxa$month <- factor(all_heatmap_filtered_10taxa$month, levels=c("Dec", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug"))

# Now, plot the heatmap with species ordered by abundance and on the y
all_date_heatmap <- ggplot(data = all_heatmap_filtered_10taxa, mapping = aes(x = month, 
                                                                             y = species, 
                                                                             fill = log_hits)) +
  geom_tile() +
  xlab(label = "Month") + 
  ylab(label = "Species") + 
  scale_fill_distiller(name = "Hits per 100,000", palette = "RdBu") +
  theme_bw() + 
  #ggtitle(label = "Tiptree PHIbase Species Abundance grouped by collection month (log hits per 100,000 > 10)") +
  theme(
    axis.title.y = element_text(size = 20, vjust = 0.5),
    axis.text.y  = element_text(size = 12, vjust = 0.5, hjust = 1),
    axis.title.x = element_text(size = 20, hjust = 0.5),
    axis.text.x  = element_text(size = 15, vjust = 0.5, hjust = 1),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.key.size = unit(5, "cm"), legend.key.width = unit(0.5, "cm"))

all_date_heatmap

# Saving the plot to the right dimensions
ggsave(filename= "../Q20_phibase_mapping/Images_Tiptree/Tiptree_heatmap_month_all_data_reordered.pdf", 
       plot = all_date_heatmap, width=12, height=18)
