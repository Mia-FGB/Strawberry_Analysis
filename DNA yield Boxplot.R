#Tiptree data 
#Boxplot of DNA yield by location 

#Setting the working directory 
setwd("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/Tiptree/Q20_phibase_mapping")

#packages
library(dplyr)
library(ggplot2)
library(lubridate)
library(forcats)

#Importing the data
yield_data <- read.csv('Tiptree_DNA_yield.csv.csv', header = TRUE)

#cleaning the data
#filtering out the too low variable
clean_yield_data <- filter(yield_data, first_extraction !='<0.05', location != "")
#want to read 1st extraction as an integer 
yield_data$first_extraction <- as.numeric(yield_data$first_extraction)
yield_data$After_WGA <- as.numeric(yield_data$After_WGA)

#need the date to be ordered
#clean_yield_data$Date = as.Date(clean_yield_data$Date, format="%d/%m/%Y")

#manually ordering
yield_data <- yield_data %>%
  mutate(Date = fct_relevel(Date, 
                            "10/02/2022", "09/03/2022", "05/04/2022", "11/05/2022", 
                            "07/06/2022", "13/07/2022", "11/08/2022"))

#first extraction
ggplot(yield_data, aes(x= Date, y = first_extraction, group = Date))+
  geom_boxplot()+
  ylab(label = "DNA First Extraction yield (ng/ul) ") +
  xlab(label = "Date of collection") +
  theme_bw()

#After WGA
ggplot(yield_data, aes(x= Date, y = After_WGA, group = Date))+
  geom_boxplot()+
  ylab(label = "After WGA DNA yield (ng/ul) ") +
  xlab(label = "Date of collection") +
  theme_bw()

# plot both boxplots in a single graph
ggplot() +
  geom_boxplot(data = yield_data, aes(x = Date, y = After_WGA),
               colour = "#ff7f00") +
  geom_boxplot(data = yield_data, aes(x = Date, y = first_extraction),
               colour = "#377eb8") +
  ylab(label = "DNA Yield (ng/ul)") +
  xlab(label = "Date of collection") +
  theme_bw()

