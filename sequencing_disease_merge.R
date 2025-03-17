#--------------------------
#Statistical tests 26.02.24 --------
#--------------------------
#Script to try and understand the relationship between the disease score & sequencing data and apply statistical tests
setwd("/Users/berelsom/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/Tiptree/Q20_phibase_mapping")
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(dplyr)
library(broom)
library(zoo)
library(lubridate)
library(astsa)

#Prepping the data---------------------
setwd("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/Tiptree/Q20_phibase_mapping")
#Sequencing data
seq_data = read.csv("all_sequence_data.csv")
seq_data$Date.collected <- as.Date(seq_data$Date.collected, format="%d/%m/%Y")

#Disease score data 
disease_score <- read.csv("tiptree_disease_data.csv", header=TRUE)
#need the date to be ordered
disease_score$date = as.Date(disease_score$date, format="%d/%m/%Y")
#filter out N/A
disease_score <- disease_score %>% 
  filter(score != 'N/A') %>% 
  subset(select = -c(total_reads, read_count))
disease_score$score <- as.numeric(disease_score$score)

# Merge --------
#Merging on location, date & species (same as genus in these datasets)
merged_data <- merge(
  seq_data, disease_score,
  by.x = c("location", "Date.collected", "Genus"),
  by.y = c("location", "date", "species"),
  all.x = TRUE,
  all.y = TRUE
)
#Does leave me with lots of half blank rows where the dates aren't the same

# Round dates to the beginning of each month
merged_data$Month <- floor_date(merged_data$Date.collected, unit = "month")

# Group data by month and calculate mean for both normalised_hits and score
#Ignoring the different locations and just using them as repeats 
monthly_data <- merged_data %>%
  group_by(Month, Genus) %>%
  summarise(mean_norm_hits = mean(hits_per_1000, na.rm = TRUE),
            sd_norm_hits = sd(hits_per_1000, na.rm = TRUE),
            mean_log_hits = mean(log_hits, na.rm = TRUE),
            sd_log_hits = sd(log_hits, na.rm = TRUE),
            mean_score = mean(score, na.rm = TRUE),
            sd_score = sd(score, na.rm = TRUE))

#Exploratory data analysis------
#Summary statistics ------------


n_distinct(merged_data$location) #Expect and got 3 unique locations
n_distinct(merged_data$Sample)   #52 unique samples
n_distinct(merged_data$score)    #Expected and got 6 unique scores (NA 0 1 2 3 4)

# Compute summary statistics 
summary(merged_data$total_reads)
summary(merged_data$genus_read)

#Number of samples in each location 
location_count <- table(merged_data$location) #- nearly even but not quite the same in each, due to disease score data?
# Create a subset of the data with non-missing values in the total_reads column
data_with_reads <- merged_data[!is.na(merged_data$total_reads), ]
# Count the number of samples with data in the total_reads column per location
samples_with_reads_per_location <- table(data_with_reads$location) #Still uneven number of replicates in each location

#Range of dates collected over
summarise(merged_data, FirstDate=first(Date.collected),LastDate=last(Date.collected))

#Histogram plots--------
# Plot a histogram of the frequency distribution of samples collected on different dates
hist(merged_data$Date.collected, 
     main = "Frequency Distribution of Samples Collected on Different Dates",
     xlab = "Date Collected",
     ylab = "Frequency",
     breaks = 20,
     xlim = as.Date(c("2021-12-01", "2022-12-31"))
     )

# Plot a histogram of the frequency distribution of scores - high skew towards the lower scores
hist(merged_data$score, 
     main = "Frequency Distribution of scores",
     xlab = "Disease score",
     ylab = "Frequency",
     breaks = 5,
)

#Relationship between read count (normalised) & disease scores over time------
#Scatter plot time on x , read count on left y, disease score on right y
# Plot the relationship between read counts and disease scores over time

break.vec <- c(as.Date("2021-12-01"), seq(from=as.Date("2021-12-01"), to=as.Date("2022-09-01"),by="month")) 

#Not split by Genus - no clear relationship between the two
ggplot(merged_data, aes(x = Date.collected)) +
  geom_point(aes(y = log_hits), shape = 17, size = 2, colour = "red") +
  geom_smooth(aes(y = log_hits), method = "lm", se = FALSE, colour = "red") +  # Add line of best fit for log_hits
  labs(
    x = "Date Collected",
    y = "Normalised Read Count (log)",
    title = "Relationship between Read Counts and Disease Scores over Time") +
  theme_minimal() +
  theme(legend.position = "right") +
  # Add a second y-axis for disease scores
  geom_point(aes(y = score), shape = 15, size = 2, colour = "blue") + 
  geom_smooth(aes(y = score), method = "lm", se = FALSE, colour = "blue") +  # Add line of best fit for disease score
  # Set y-axis label and title for the right y-axis
  scale_y_continuous(sec.axis = sec_axis(~., name = "Disease Scores")) +
  scale_x_date(date_labels = "%Y-%m-%d")

#Doesn't work all on one axis - potentially could have 6 colours but would be v confusing
ggplot(merged_data, aes(x = Date.collected)) +
  geom_point(aes(y = log_hits, color = Genus), shape = 17, size = 2) +
  geom_smooth(aes(y = log_hits, color = Genus), method = "lm", se = FALSE) +  # Add line of best fit for log_hits
  labs(
    x = "Date Collected",
    y = "Normalised Read Count (log)",
    title = "Relationship between Read Counts and Disease Scores over Time For Each Genera") +
  theme_minimal() +
  theme(legend.position = "right") +
  # Add a second y-axis for disease scores
  geom_point(aes(y = score, color = Genus), shape = 15, size = 2) + 
  geom_smooth(aes(y = score, color = Genus), method = "lm", se = FALSE) +  # Add line of best fit for disease score
  # Set y-axis label and title for the right y-axis
  scale_y_continuous(sec.axis = sec_axis(~., name = "Disease Scores")) +
  scale_x_date(date_labels = "%Y-%m-%d")


# Instead facet so each Genus is on it's own graph
genus_rc_ds <- ggplot(merged_data, aes(x = Date.collected, y = log_hits)) +
  geom_point(colour = "#9DD0BD", shape = 17) + #triangle
  geom_smooth(aes(y = log_hits), method = "lm", se = FALSE, colour = "#009E73") + #line of best fit
  labs(
    x = "Date Collected (Month)",
    y = "Normalised Read Count (log)",
    title = "Relationship between Read Counts and Disease Scores over Time"
    ) +
  facet_wrap(~ Genus, scales = "free_y") +
  theme_minimal() +
  # Add a second y-axis for disease scores
  geom_point(aes(y = score), color = "#E7B293", shape = 15) + #square
  geom_smooth(aes(y = score), method = "lm", se = FALSE, colour = "#D15D0D") +
  # Set y-axis label and title for the right y-axis
  scale_y_continuous(sec.axis = sec_axis(~., name = "Disease Scores")) +
  scale_x_date(breaks = break.vec, date_labels = "%m", limits=range(break.vec)) 

ggsave(
  filename= "Images_Tiptree/genus_read_count_disease_score.pdf",
  plot = genus_rc_ds,
  width=9,
  height=3
  )

#Plotting the mean of each group per date -----

#mean, min & max data for each location & genus for errorbar 
#can't integrate disease score into these graphs 
data_min_max <- merged_data %>% 
  drop_na(hits_per_1000) %>%  # Drop rows with missing values in log_hits
  group_by(Genus, Date.collected) %>% 
  summarise(mean = mean(hits_per_1000, na.rm = TRUE), 
            min = min(hits_per_1000, na.rm = TRUE), 
            max = max(hits_per_1000, na.rm = TRUE))

#these aren't errorbars!! just min and max values as there are only 2 replicates from each
genus_normalised_hits_error <- 
  ggplot(data_min_max,
         aes(x = Date.collected,
             y = mean
             )) +
  geom_line() + 
  geom_point()+
  geom_errorbar(aes(ymin = min, ymax = max),  width=7, colour = "grey 39") +
  facet_wrap(~ Genus, scales = "free_y") +
  #Removing gridlines
  theme_bw() + theme( panel.grid.minor = element_blank()) +
  xlab(label = "Date collected") + ylab(label = "Hits per 1000 (log)") +
  scale_x_date(breaks = break.vec, date_labels = "%m", limits=range(break.vec)) +
  ggtitle("Hits per 1000 from sequencing")

#Try location instead of Genus - no visible relationship as expected ----
location_rc_ds <- ggplot(merged_data, aes(x = Date.collected, y = log_hits)) +
  geom_point(colour = "#9DD0BD", shape = 17) + #triangle
  geom_smooth(aes(y = log_hits), method = "lm", se = FALSE, colour = "#009E73") +
  labs(
    x = "Date Collected (Month)",
    y = "Normalised Read Count (log)",
    title = "Relationship between Read Counts and Disease Scores over Time"
  ) +
  facet_wrap(~ location, scales = "free_y") +
  theme_minimal() +
  # Add a second y-axis for disease scores
  geom_point(aes(y = score), color = "#E7B293", shape = 15) + #square
  geom_smooth(aes(y = score), method = "lm", se = FALSE, colour = "#D15D0D") +
  # Set y-axis label and title for the right y-axis
  scale_y_continuous(sec.axis = sec_axis(~., name = "Disease Scores")) +
  scale_x_date(breaks = break.vec, date_labels = "%m", limits=range(break.vec)) 

#Plotting monthly------

# Plot the aggregated data
month_rc_ds <- ggplot(monthly_data, aes(x = Month)) +
  
  #Unlogged data - difficult to interpret the relationship
  #geom_line(aes(y = mean_norm_hits), colour = "#009E73") +
  #geom_errorbar(aes(
    #ymin=mean_norm_hits-sd_norm_hits,
    #ymax=mean_norm_hits+sd_norm_hits),
    #width=.2,
    #colour = "#9DD0BD"
    #) +
  
  #Logged data
  geom_line(aes(y = mean_log_hits), colour = "#009E73") +
  geom_errorbar(aes(
   ymin=mean_log_hits-sd_log_hits,
   ymax=mean_log_hits+sd_log_hits),
   width=.2,
   colour = "#9DD0BD"
   ) +
  geom_line(aes(y = mean_score), colour = "#D15D0D") +
  geom_errorbar(aes(
    ymin=mean_score-sd_score,
    ymax=mean_score+sd_score),
    width=.2,
    colour = "#E7B293"
  ) +
  facet_wrap(~ Genus, scales = "free_y") +
  labs(x = "Month", y = "Normalised Read Count (log)") +
  theme_minimal() +
  ggtitle("Aggregated Data by Month") +
  scale_x_date(breaks = break.vec, date_labels = "%m", limits=range(break.vec)) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "Disease Score")) 

ggsave(
  filename= "Images_Tiptree/monthly_aggregate_read_count_disease_score.pdf",
  plot = month_rc_ds,
  width=9,
  height=3
)


# -------------------------------------------------------------------------
#Pearsons correlation --------
# Group the data by Genus
grouped_gen_data <- merged_data %>%
  group_by(Genus)
  

group_loc_data <- merged_data %>% #no correlation with location
  group_by(location)

# Calculate Pearson correlation coefficient within each group
correlation_grouped <- summarise(
  group_loc_data,
  correlation = cor(hits_per_1000, score, use = "pairwise.complete.obs"))

# Calculate Pearson correlation coefficient and its p-value
# Perform correlation test within each group
correlation_test <- group_loc_data %>%
  do(tidy(cor.test(.$hits_per_1000, .$score, method = "pearson")))


#Cross correlation analysis --------
#I don't think I have enough data for this to work
# Remove rows with missing values - both need to be the same length
whole_data <- na.omit(merged_data) #only 36 values so low power

# Create time series objects for disease score and read count
score_ts <- ts(whole_data$score, start = min(whole_data$Date.collected), frequency = 1)
read_count_ts <- ts(whole_data$hits_per_1000, start = min(whole_data$Date.collected), frequency = 1)

# Perform cross-correlation analysis
cross_correlation <- ccf(score_ts, read_count_ts, lag.max = 20)  # Adjust lag.max as needed

# Plot cross-correlation
plot(cross_correlation, main = "Cross-Correlation Plot: Disease Score vs. Read Count")

#For each Genus --
# 1. Group the data by a categorical variable (e.g., 'Genus')
grouped_data <- split(whole_data, whole_data$Genus)

# 2. Compute cross-correlation for each group
cross_correlation_results <- lapply(grouped_data, function(group) {
  score_ts <- ts(group$score, start = min(group$Date.collected), frequency = 1)
  read_count_ts <- ts(group$total_reads, start = min(group$Date.collected), frequency = 1)
  ccf_result <- ccf(score_ts, read_count_ts, lag.max = 10)
  return(ccf_result)
})

# 3. Visualize the cross-correlation results for each group
par(mfrow=c(2, 2))  # Set up a multi-panel plot
for (i in 1:length(cross_correlation_results)) {
  plot(cross_correlation_results[[i]], main = paste("Cross-Correlation Plot for", names(cross_correlation_results)[i]))
}

#Trying a different package --
#Not realy showing anything
astsa::lag2.plot(read_count_ts, score_ts, 6)


# -------------------------------------------------------------------------


#Comparing disease score & read count from the month prior -- 28.02.24 -----
# Use table() to count the number of entries for each date
table(seq_data$Date.collected)
# those with 0 aren't in the analysis - ignore for now and see how it works
#ultimately will need to add them back in :(


#Need to lag first and then merge on disease score data
#Probably want to get average of the repeats and use this in the lag analysis
group_seq_data <- seq_data %>%
  group_by(Date.collected, location, Genus) %>% 
  filter(total_reads>=1000) %>%  #Remove rows where the read count is low, think these will are outliers due to technical issue
  summarise(mean_norm_hits = mean(hits_per_1000, na.rm = TRUE), #Calculate mean of the groups
            sd_norm_hits   = sd(hits_per_1000,   na.rm = TRUE),
            ) %>% 
  mutate(mean_log_hits = log(mean_norm_hits)) %>% #better to log the mean values
  arrange(Date.collected, location, Genus) # Sort the data by location, genus, and Date.collected

# Group by location and genus, then lag the hits_per_1000 column
lag_data <- group_seq_data  %>%
  group_by(location, Genus) %>%
  mutate(prior_norm_hits    = lag(mean_norm_hits, 2),
         prior_sd_norm_hits = lag(sd_norm_hits, 2),
         prior_log_hits    = lag(mean_log_hits, 2),
         )

#Merging on location, date & species (same as genus in these datasets)
merg_lag_data <- merge(
  lag_data, disease_score,
  by.x = c("location", "Date.collected", "Genus"),
  by.y = c("location", "date", "species"),
  all.x = TRUE,
  all.y = TRUE
)
#Make score numeric
merg_lag_data$score <- as.numeric(merg_lag_data$score)

#Plot the lag data
# Instead facet so each Genus is on it's own graph
prior_plot <- ggplot(merg_lag_data, aes(x = Date.collected, y = mean_log_hits)) +
 
   #Normal data
  geom_point(colour = "#90F5D2", shape = 17, alpha = .5) + #triangle
  geom_smooth(aes(y = mean_log_hits), method = "lm", se = FALSE, colour = "#90F5D2") + 
  
  #Month Prior data
  geom_point(colour = "#3DBCF5", shape = 19, alpha = .5) + #circle
  geom_smooth(aes(y = prior_log_hits), method = "lm", se = FALSE, colour = "#3DBCF5") + 
  labs(
    x = "Date Collected (Month)",
    y = "Normalised Read Count (Log)",
    title = "Relationship between Disease Scores and Read counts from Two Time Points"
  ) +
  facet_wrap(~ Genus, scales = "free_y") +
  theme_minimal() +
  # Add a second y-axis for disease scores
  geom_point(aes(y = score), color = "#F4925B", shape = 15, alpha = .5) + #square
  geom_smooth(aes(y = score), method = "lm", se = FALSE, colour = "#F4925B") +
  # Set y-axis label and title for the right y-axis
  scale_y_continuous(sec.axis = sec_axis(~., name = "Disease Scores")) +
  scale_x_date(breaks = break.vec, date_labels = "%m", limits=range(break.vec)) 

ggsave(
  filename= "Images_Tiptree/prior_read_count_disease_score.pdf",
  plot = prior_plot,
  width=9,
  height=3
)

#Pearson's correlation on prior data
# Group the data by Genus
grouped_lag_data <- merg_lag_data %>%
  group_by(Genus)

# Calculate Pearson correlation coefficient within each group
pr_corr_grouped <- summarise(
  grouped_lag_data,
  correlation = cor(prior_norm_hits, score, use = "pairwise.complete.obs"))

# Calculate Pearson correlation coefficient within each group
corr_grouped <- summarise(
  grouped_lag_data,
  correlation = cor(mean_norm_hits, score, use = "pairwise.complete.obs"))

# Calculate Pearson correlation coefficient and its p-value
# Perform correlation test within each group
pr_corr_test <- grouped_lag_data %>%
  do(tidy(cor.test(.$prior_norm_hits, .$score, method = "pearson")))

corr_test <- grouped_lag_data %>%
  do(tidy(cor.test(.$mean_norm_hits, .$score, method = "pearson")))
