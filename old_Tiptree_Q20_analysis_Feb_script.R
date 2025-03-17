#Tiptree Q20 - reads - heatmap - read_count - disease score - temperature - humidity 
#Setting the working directory & adding packages
setwd("~/OneDrive - Norwich BioScience Institutes/Air_Samples/Tiptree/Q20_phibase_mapping/")
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(stringi)
library(tidyr)
library(dplyr)
library(lubridate)
library(scales)
library(patchwork)
 

#Prepping the data---------------------
#reading in the data
phi <- read.csv("Q20_Tiptree_combined_data.csv", header=TRUE)

#Need to normalise the reads for each location, e.g. hits per 10,000 reads
#Therefore we need the total number of reads for each sample, not just the classified reads
phi_filter <- 
  phi %>% 
  #filter out any species with less than 10 reads
  #lots of blank space so raising to 20
  filter(read_count >=10 ) %>% 
  mutate(hits_per_10.000 = (read_count * 10000)/(total_reads)) %>% 
  #log scale-Mutate the data to add a new row
  mutate(log_hits = log10(hits_per_10.000)) %>% 
  #adding a row for percentage 
  mutate(percent = (read_count/total_reads)*100) %>% 
  #don't want sample to be number but chr 
  mutate(Sample = as.character(Sample))

#Adding a col for month
phi_filter$month <- format(as.Date(phi_filter$Date.collected, format="%d/%m/%Y"),"%m")


#--------------------------------------------------------------------------
#Tracking species of interest  --------------------------------------
#--------------------------------------------------------------------------

#Want all the graphs to have the same x axis, setting it here
break.vec <- c(as.Date("2022-02-03"),
               seq(from=as.Date("2022-02-03"), to=as.Date("2022-08-20"),by="month")) 


#Sample name, Date, Location, Pass reads,
#Reads matching Botrytis, Phytophthora, Podosphera

phi_summarised <- phi %>% 
  filter(grepl('Podosphaera|Botrytis|Phytophthora', species)) %>% 
  select(Sample, Date.collected, location, total_reads, read_count, species) %>% 
  #Only keeping reads from samples with enough DNA as the others are outliers
  #will remove 5 samples
  filter(total_reads >= 3200) %>% 
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
#Writing out to a csv so I can load it to my newer script
write.csv(phi_summarised, "Tiptree_Feb_Aug_3_genus.csv", row.names = FALSE)

#need the date to be ordered
phi_summarised$Date.collected = as.Date(phi_summarised$Date.collected, format="%d/%m/%Y")

#3_genus_scaled_sequence 
genus_scaled_sequence <- ggplot(phi_summarised, aes(x = Date.collected, y = log_hits, group = location, colour = location)) +
  geom_point() + 
  facet_grid(facets =  vars(Genus), cols = vars(location)) + theme_bw() +
  xlab(label = "Date collected") + ylab(label = "Hits per 1000 (log)") +
  scale_x_date(breaks = break.vec, date_labels = "%m", limits=range(break.vec)) +
  ggtitle("Hits per 1000 from sequencing")

#Trying to get mean, min & max data for each location & genus
#Can use this to plot errorbars
phi_min_max <- phi %>% 
  filter(grepl('Podosphaera|Botrytis|Phytophthora', species)) %>% 
  select(Sample, Date.collected, location, total_reads, read_count, species) %>%
  #Only keeping reads from samples with enough DNA as the others are outliers, will remove 5 samples
  filter(total_reads >= 3200) %>% 
  #splitting group and species
  separate(species, c("Genus", "species"), sep = " ") %>% 
  #no longer want the species col 
  subset(select = -(species)) %>% 
  group_by(Genus,Sample, location, Date.collected, total_reads) %>% 
  #Grouping all the species into genus counts
  summarise(genus_count = (sum(read_count))) %>% 
  #adding some normalised values
  mutate(log_hits = (log10((genus_count * 1000)/(total_reads)))) %>% 
  #need to change grouping
  ungroup() %>% 
  group_by(Genus, location, Date.collected) %>% 
  summarise(mean = mean(log_hits), min = min(log_hits), max = max(log_hits))

phi_min_max$Date.collected = as.Date(phi_min_max$Date.collected, format="%d/%m/%Y")

genus_scaled_sequence_error <- ggplot(phi_min_max, aes(x = Date.collected, y = mean, group = location)) +
  geom_line() + 
  geom_point()+
  geom_errorbar(aes(ymin = min, ymax = max),  width=7, colour = "grey 39") +
  facet_grid(facets =  vars(Genus), cols = vars(location)) + theme_bw() +
  xlab(label = "Date collected") + ylab(label = "Hits per 1000 (log)") +
  scale_x_date(breaks = break.vec, date_labels = "%m", limits=range(break.vec)) +
  ggtitle("Hits per 1000 from sequencing")

genus_scaled_sequence_error

#--------------------------------------------------------------------------
# Disease scoring  --------------------------------------
#--------------------------------------------------------------------------

disease_score <- read.csv("tiptree_disease_data.csv", header=TRUE)

#need the date to be ordered
disease_score$date = as.Date(disease_score$date, format="%d/%m/%Y")

#filter out N/A
filter_disease_score <- disease_score %>% 
  filter(score != 'N/A') %>% 
  filter(date != "2022-12-26")

#3_genus_scaled_score
genus_scaled_score <- ggplot(filter_disease_score, aes(x = date, y = score, group = location)) +
  geom_point()+
  facet_grid(facets =  vars(species), cols= vars(location)) + theme_bw() +
  xlab(label = "Month") + ylab(label = "Disease score")  +
  scale_x_date(breaks = break.vec, date_labels = "%m", limits=range(break.vec)) +
  ggtitle("Disease score")

#--------------------------------------------------------------------------
#Spraying data ------------------------------
#--------------------------------------------------------------------------

spraying <- read.csv("Tiptree_spraying.csv", header=TRUE)

#need the date to be ordered
spraying$Date = as.Date(spraying$Date, format="%d/%m/%Y")
na_remove_spraying <- spraying %>% 
  filter(protection != 'N/A')

#3_genus_scaled_spraying
genus_scaled_spraying <-  ggplot(na_remove_spraying, aes(x = Date, y = protection)) +
  geom_point()+
  facet_grid(cols= vars(Location)) + theme_bw() +
  xlab(label = "Month") + ylab(label = "Fungicide") + 
  #to make it a square and match the others
  scale_x_date(breaks = break.vec, date_labels = "%m", limits=range(break.vec)) +
  ggtitle("Spraying dates")


#Trying to get all my plots the same size so going to try and save them with set dimensions 
#Then rearrange with ppt/illustrator
ggsave(filename= "genus_scaled_sequence_error_set.png", plot = genus_scaled_sequence_error, width=10, height=6)
ggsave(filename= "genus_scaled_score_set.png", plot = genus_scaled_score, width=10, height=6)
ggsave(filename= "genus_scaled_spraying_set.png", plot = genus_scaled_spraying, width=9, height=2)

#Older non error bar graph
ggsave(filename= "genus_scaled_sequence_set.png", plot = genus_scaled_sequence, width=10, height=6)

#------------------------------------------------------------
#Temp and humidity ------------------------------------------
#------------------------------------------------------------
#Need all graphs to have same limit on the y axis
#e.g. -40 to 40

#RG = right greenhouse, empty, roundA/NGS1-------------------
#Turning any blanks/spaces into NA's to not affect mean calc
NGS1 <- read.csv("NGS1_Temp_Hum_all.csv", header=TRUE, na.strings=c(""," "))

#Use a class to work with the datetime easily
NGS1$Date <- as.POSIXct(NGS1$Date,format="%d-%m-%Y %H:%M:%OS")
NGS1$Humidity <- as.numeric(NGS1$Humidity)
#Setting lims for the graph
lims <- as.POSIXct(strptime(c("2022-02-03 00:00:00","2022-08-20 00:00:00"), format = "%Y-%m-%d %H:%M"))

#Temp--
#plotting everything -> quite difficult to see changes 
ggplot(NGS1, aes(x = Date, y = Temp)) + geom_line() +
  xlab(label = "Month") + ylab(label = "Temperature") + 
  #to make it a square and match the others
  scale_x_datetime(labels=date_format("%m"),
                   breaks = date_breaks("month"),
                   limits = lims) +
  theme_bw() +
  ggtitle("Temperature NGS1")

#Getting average temp & humidity for each day 
A_NGS1_avg <- NGS1 %>% 
  #using lubridate package to extract just the date (ignoring the time)
  mutate (date_only = date(Date)) %>% 
  group_by(date_only) %>% 
  summarise(mean_temp = mean(Temp,na.rm=TRUE), sd_temp = sd(Temp, na.rm = TRUE),
            mean_hum = mean(Humidity, na.rm=TRUE), sd_hum = sd(Humidity, na.rm = TRUE))

#Avg Temp NGS1
A_NGS1_Avg_Temp <- ggplot(A_NGS1_avg, aes(x = date_only, y = mean_temp)) + geom_line(stat = "identity") +
  geom_errorbar(aes(ymin = mean_temp - sd_temp, ymax = mean_temp + sd_temp), color = "grey 30", alpha = .5) +
  xlab(label = "Month") + ylab(label = "Temperature") + 
  #to make it a square and match the others
  scale_x_date(breaks = break.vec, date_labels = "%m", limits=range(break.vec)) +
  theme_bw() +  theme( panel.grid.minor = element_blank()) +
  ylim(-10,40) +
  ggtitle("NGS1")

#Avg humidity NGS1
A_NGS1_Avg_Hum <- ggplot(A_NGS1_avg, aes(x = date_only, y = mean_hum)) + geom_line(stat = "identity") +
  geom_errorbar(aes(ymin = mean_hum - sd_hum, ymax = mean_hum + sd_hum), color = "grey 30", alpha = .5) +
  xlab(label = "Month") + ylab(label = "Humidity (%)") + 
  #to make it a square and match the others
  scale_x_date(breaks = break.vec, date_labels = "%m", limits=range(break.vec)) +
  theme_bw() + theme( panel.grid.minor = element_blank())  +
  ylim(30,100)

#LG = left greenhouse, full, round B/NGS2---------------
NGS2 <- read.csv("NGS2_Temp_Hum_all.csv", header=TRUE, na.strings=c(""," "))

#Use a class to work with the datetime easily
NGS2$Date <- as.POSIXct(NGS2$Date,format="%d-%m-%Y %H:%M:%OS")
NGS2$Humidity <- as.numeric(NGS2$Humidity)
#Setting lims for the graph
lims <- as.POSIXct(strptime(c("2022-02-03 00:00:00","2022-08-20 00:00:00"), format = "%Y-%m-%d %H:%M"))

#Calculating averages
NGS2_avg <- NGS2 %>% 
  #using lubridate package to extract just the date (ignoring the time)
  mutate (date_only = date(Date)) %>% 
  group_by(date_only) %>% 
  summarise(mean_temp = mean(Temp,na.rm=TRUE), sd_temp = sd(Temp, na.rm = TRUE),
            mean_hum = mean(Humidity, na.rm=TRUE), sd_hum = sd(Humidity, na.rm = TRUE))

#Plotting the average Temp NGS2
NGS2_Avg_Temp <- ggplot(NGS2_avg, aes(x = date_only, y = mean_temp)) + geom_line(stat = "identity") +
  geom_errorbar(aes(ymin = mean_temp - sd_temp, ymax = mean_temp + sd_temp), color = "grey 30", alpha = .5) +
  xlab(label = "Month") + ylab(label = "Temperature") + 
  #to make it a square and match the others
  scale_x_date(breaks = break.vec, date_labels = "%m", limits=range(break.vec)) +
  ylim(-10,40) +
  theme_bw() + theme( panel.grid.minor = element_blank()) +
  ggtitle("NGS2")

#Plotting average humidity NGS2
NGS2_Avg_Hum <- ggplot(NGS2_avg, aes(x = date_only, y = mean_hum)) + geom_line(stat = "identity") +
  geom_errorbar(aes(ymin = mean_hum - sd_hum, ymax = mean_hum + sd_hum), color = "grey 30", alpha = .5) +
  xlab(label = "Month") + ylab(label = "Humidity (%)") + 
  #to make it a square and match the others
  scale_x_date(breaks = break.vec, date_labels = "%m", limits=range(break.vec)) +
  theme_bw() + theme( panel.grid.minor = element_blank()) +
  ylim(30,100) 

#Field Temp---------------------------
Field <- read.csv("Field_Temp_01-01-22_24-08-22.xls.csv", header=TRUE, na.strings=c(""," "))

#Use a class to work with the datetime easily
Field$Date <- as.POSIXct(Field$Date,format="%d/%m/%Y %H:%M")

#Getting average temp & humidity for each day 
Field_avg <- Field %>% 
  #using lubridate package to extract just the date (ignoring the time)
  mutate (date_only = date(Date)) %>% 
  group_by(date_only) %>% 
  summarise(mean_temp = mean(Temp,na.rm=TRUE), sd_temp = sd(Temp, na.rm = TRUE))

#Plotting the average Temp Field
#Plotting the average Temp Field
Field_Avg_Temp <- ggplot(Field_avg, aes(x = date_only, y = mean_temp)) + geom_line(stat = "identity") +
  geom_errorbar(aes(ymin = mean_temp - sd_temp, ymax = mean_temp + sd_temp), color = "grey 30", alpha = .5) +
  xlab(label = "Month") + ylab(label = "Temperature") + 
  #to make it a square and match the others
  scale_x_date(breaks = break.vec, date_labels = "%m", limits=range(break.vec)) +
  theme_bw() + theme( panel.grid.minor = element_blank()) +
  ylim(-10,40) +
  ggtitle("Field")

#Putting these graphs into a grid 
#Uses patchwork package
Temp_Hum_Graphs <- (Field_Avg_Temp + NGS2_Avg_Temp + A_NGS1_Avg_Temp) / 
  (plot_spacer() + NGS2_Avg_Hum + A_NGS1_Avg_Hum) +
  plot_annotation(title = 'Average Daily Temperature and Humidity')

ggsave(filename= "Temp_Hum_Graphs.png", plot = Temp_Hum_Graphs, width=10, height=5)

#Old stuff I didn't use
#--------------------------------------------------------------------------

#Heatmaps
#------------------------------------------------------------
#Graph Splitting all the samples to own column---------------
#------------------------------------------------------------

phi.heatmap <- ggplot(data = phi_filter, mapping = aes(x = Sample,
                                                       y = species,
                                                       fill = log_hits)) +
  geom_tile() +
  xlab(label = "Sample") + ylab(label = "Species") +
  scale_fill_distiller(name = "Hits per 10,000", palette = "RdBu") +
  theme_bw() +
  ggtitle(label = "Tiptree Air-Seq Species Abundance (log hits per 10,000 > 10)") +
  #Editing axis title size
  theme(axis.text.y = element_text(vjust = 0.5, hjust=1, size = 20,)) +
  theme(axis.title.y = element_text(size = 20, vjust = 0.5)) +
  theme(axis.title.x = element_text(size = 20, hjust = 0.5)) +
  theme(axis.text.x = element_text(vjust=0.5, size = 15,)) +
  # Change legend key size and key width
  theme(legend.key.size = unit(6, "cm"),legend.key.width = unit(1.5,"cm"))+
  #Change legend font size
  theme(legend.title = element_text(size = 15)) +
  theme(legend.text = element_text(size =15)) +
  #Changing the title size
  theme(plot.title = element_text(size = 20, hjust = 0.5))

phi.heatmap
#Saving the plot to right dimensions, much taller y axis
ggsave(filename= "Tiptree_heatmap_by_sample.png", plot = phi.heatmap, width=30, height=40)



#------------------------------------------------------------
#Graph Splitting all the Locations to own column---------------
#------------------------------------------------------------

phi.loc.heatmap <- ggplot(data = phi_filter, mapping = aes(x = location,
                                                           y = species,
                                                           fill = log_hits)) +
  geom_tile() +
  xlab(label = "Location") + ylab(label = "Species") +
  scale_fill_distiller(name = "Hits per 10,000", palette = "RdBu") +
  theme_bw() +
  ggtitle(label = "Tiptree Air-Seq Species Abundance by location (log hits per 10,000 > 10)") +
  #Editing axis title size
  theme(axis.text.y = element_text(vjust = 0.5, hjust=1, size = 20,)) +
  theme(axis.title.y = element_text(size = 20, vjust = 0.5)) +
  theme(axis.title.x = element_text(size = 20, hjust = 0.5)) +
  theme(axis.text.x = element_text(vjust=0.5, size = 15,)) +
  # Change legend key size and key width
  theme(legend.key.size = unit(6, "cm"),legend.key.width = unit(1.5,"cm"))+
  #Change legend font size
  theme(legend.title = element_text(size = 15)) +
  theme(legend.text = element_text(size =15)) +
  #Changing the title size
  theme(plot.title = element_text(size = 20, hjust = 0.5))

phi.loc.heatmap
#Saving the plot to right dimensions, much taller y axis
ggsave(filename= "Tiptree_heatmap_by_location.png", plot = phi.loc.heatmap, width=30, height=40)


#------------------------------------------------------------
#Graph Splitting all the Dates to own column---------------
#------------------------------------------------------------



phi.date.heatmap <- ggplot(data = phi_filter, mapping = aes(x = month,
                                                            y = species,
                                                            fill = log_hits)) +
  geom_tile() +
  xlab(label = "Month") + ylab(label = "Species") +
  scale_fill_distiller(name = "Hits per 10,000", palette = "RdBu") +
  theme_bw() +
  ggtitle(label = "Tiptree Air-Seq Species Abundance by Month of collection (log hits per 10,000 > 10)") +
  #Editing axis title size
  theme(axis.text.y = element_text(vjust = 0.5, hjust=1, size = 20,)) +
  theme(axis.title.y = element_text(size = 20, vjust = 0.5)) +
  theme(axis.title.x = element_text(size = 20, hjust = 0.5)) +
  theme(axis.text.x = element_text(vjust=0.5, size = 15,)) +
  # Change legend key size and key width
  theme(legend.key.size = unit(6, "cm"),legend.key.width = unit(1.5,"cm"))+
  #Change legend font size
  theme(legend.title = element_text(size = 15)) +
  theme(legend.text = element_text(size =15)) +
  #Changing the title size
  theme(plot.title = element_text(size = 20, hjust = 0.5))

phi.date.heatmap
#Saving the plot to right dimensions, much taller y axis
ggsave(filename= "Tiptree_heatmap_by_date.png", plot = phi.date.heatmap, width=30, height=40)
#Think I need a new dataframe where I input 0's where there is no data 

#Graph of podosphera through months --------------------------------------

phi_paphanis <- phi_filter %>%  
  #only include the p. aphanis reads
  filter (species== 'Podosphaera aphanis') %>% 
  group_by(month, location) %>% 
  summarise(avg_norm_read = mean(hits_per_10.000))
#There are some missing rows as P aphanis was not found in all places every month

#phi_paphanis <- phi_filter %>%  
#only include the p. aphanis reads
#  filter (species== 'Podosphaera aphanis') %>% 
#  group_by(Date.collected, location) %>% 
#  summarise(avg_norm_read = mean(hits_per_10.000))


ggplot(phi_paphanis, aes(x=month, y = avg_norm_read, colour = location, group = location)) +
  geom_point() +geom_line() + theme_bw() +
  ggtitle(label = "Average normalised P.apahnis reads over time")

#Graph of botrytis through months --------------------------------------


phi_bcinerea <- phi_filter %>%  
  #only include the p. aphanis reads
  filter (species== 'Botrytis cinerea') %>% 
  group_by(month, location) %>% 
  summarise(avg_norm_read = mean(hits_per_10.000))

#phi_bcinerea <- phi_filter %>%  
#only include the p. aphanis reads
#  filter (species== 'Botrytis cinerea') %>% 
#  group_by(Date.collected, location) %>% 
#  summarise(avg_norm_read = mean(hits_per_10.000))

ggplot(phi_bcinerea, aes(x=month, y = avg_norm_read, colour = location, group = location )) +
  geom_point(size = 0.5) + geom_line() + theme_bw() +
  ggtitle(label = "Average normalised B.cinerea reads over time")

