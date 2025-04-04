#Comparing total read count and DNA yield of sample 

#Setting the working directory & adding packages
setwd("~/OneDrive - Norwich BioScience Institutes/Air_Samples/Tiptree")
library(dplyr)
library(ggplot2)

#Reading in the data file
data <- read.csv("DNA_yield_read.csv", header=TRUE)

#log the data so it's better spread out
#adding a column to the dataset for logged reads
data <- 
  data %>% 
  #log scale-Mutate the data to add a new row
  mutate(log_no_reads = log10(no_reads)) %>% 
  mutate(log_extraction_yield = log10(extraction_yield)) %>% 
  mutate(log_T7_yield = log10(T7_yield))

#Graphs--------------
#want to plot a line graph with DNA yield on x and read count on y
#graph for extraction & T7 yields seperate (otherwise would have to re-do my table)

#Extraction-----
ggplot(data, aes(x= extraction_yield, y = no_reads)) +
  geom_point(col = "grey", size = 1, alpha = 0.8) +
  #adding in a linear model regression line 
  geom_smooth(method = "lm", se = FALSE, col = "firebrick") +
  xlab(label = "DNA Extraction yield (ng/ul) ") + ylab(label = "Number of basecalled reads") +
  theme_bw()

#Logged extraction--------
ggplot(data, aes(x= log_extraction_yield, y = log_no_reads)) +
  geom_point(col = "grey50", size = 1) +
  #adding in a linear model regression line 
  geom_smooth(method = "lm", se = FALSE, col = "chartreuse3") +
  xlab(label = "DNA Extraction yield (ng/ul) (Log10) ") + ylab(label = "Number of basecalled reads (Log10)") +
  ggtitle("Extraction DNA yield and number of basecalled reads after seqeuncing (Log 10)") +
  theme_bw() +
  #changing the title text appearance 
  theme(plot.title = element_text(face = "bold", size = 10))


#T7 yield------
#making T7 a integer column 
data$T7_yield <- as.numeric(data$T7_yield)
#graph
ggplot(data, aes(x= T7_yield, y = no_reads)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)

#T7 logged 
ggplot(data, aes(x= log_T7_yield, y = log_no_reads)) +
  geom_point(col = "grey50", size = 1) +
  #adding in a linear model regression line 
  geom_smooth(method = "lm", se = FALSE, col = "chartreuse3") +
  xlab(label = "T7 DNA yield (ng/ul) (Log10) ") + ylab(label = "Number of basecalled reads (Log10)") +
  ggtitle("DNA yield after T7 step and number of basecalled reads after seqeuncing (Log 10)") +
  theme_bw() +
  #changing the title text appearance 
  theme(plot.title = element_text(face = "bold", size = 10))

