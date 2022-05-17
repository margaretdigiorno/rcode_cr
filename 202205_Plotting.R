library(dplyr)
library(lubridate)
library(ggplot2)
library(gridExtra)
library(grid)

#### Load Data ####
slopefiles = list.files(path = "Data/WQ_Slopes", pattern = "V2")
slopefiles = paste0("Data/WQ_Slopes/", slopefiles)
slopes = lapply(slopefiles, read.csv)

apr = rbind(slopes[[1]], slopes[[2]], slopes[[3]])
feb = rbind(slopes[[4]], slopes[[5]], slopes[[6]], slopes[[7]])
jul = rbind(slopes[[8]], slopes[[9]], slopes[[10]], slopes[[11]])

# Store the slope/difference thresholds
tempslope = 0.5
tempdiff = 0.1
ecslope = 3
ecdiff = 0.6

# Make the dfs into a list to loop through
slopedfs = list(feb, jul, apr)

#### ID Anomalies ####

# Use mutate to create anomaly columns (0 : background, 1: anomaly, NA : no slope
for (i in 1:length(slopedfs)){
    slopedfs[[i]]$TIMESTAMP <- as_datetime(slopedfs[[i]]$TIMESTAMP, tz = "America/Los_Angeles")
    slopedfs[[i]] <- slopedfs[[i]] %>% mutate(anom_dtemp = as.factor(case_when(is.na(temp_deepslope) ~ "NA",
                                                          abs(temp_deepslope) > tempslope & abs(temp_deepdiff) > tempdiff ~ "1",
                                                          TRUE ~ "0")))
    
    slopedfs[[i]] <- slopedfs[[i]] %>% mutate(anom_stemp = as.factor(case_when(is.na(temp_surfslope) ~ "NA",
                                                          abs(temp_surfslope) > tempslope & abs(temp_surfdiff) > tempdiff ~ "1",
                                                          TRUE ~ "0")))
    
    slopedfs[[i]] <- slopedfs[[i]] %>% mutate(anom_dec = as.factor(case_when(is.na(ec_deepslope) ~ "NA",
                                                          abs(ec_deepslope) > ecslope & abs(ec_deepdiff) > ecdiff ~ "1",
                                                          TRUE ~ "0")))
    
    slopedfs[[i]] <- slopedfs[[i]] %>% mutate(anom_sec = as.factor(case_when(is.na(ec_surfslope) ~ "NA",
                                                          abs(ec_surfslope) > ecslope & abs(ec_surfdiff) > ecdiff ~ "1",
                                                          TRUE ~ "0")))
}

# Put data back into individual dataframes
feb = slopedfs[[1]]
jul = slopedfs[[2]]
apr = slopedfs[[3]]

#### Plot anomalies to check the thresholds ####
anoms_timeseries <- function(df, month){
  for (i in unique(day(df$TIMESTAMP))){
    temp = df[day(df$TIMESTAMP) == i,]
    plotted <- grid.arrange(
      (ggplot(temp)+geom_point(aes(x = TIMESTAMP, y = surfcond, color = anom_sec))),
      (ggplot(temp)+geom_point(aes(x = TIMESTAMP, y = deepcond, color = anom_dec))),
      (ggplot(temp)+geom_point(aes(x = TIMESTAMP, y = surftemp, color = anom_stemp))),
      (ggplot(temp)+geom_point(aes(x = TIMESTAMP, y = deeptemp, color = anom_dtemp))),
      top=textGrob(paste0(month, i)),
      ncol = 1
    )
    ggsave(paste0(month, '_', i,'.png'), plot = plotted, width = 10.5, height = 8, units = 'in')
  }
}

anoms_timeseries(feb, 'February')
anoms_timeseries(jul, 'July')
anoms_timeseries(apr, 'April')

#### Average Uncorrected Values for Neil ####
# Calculate average uncorrected background values for FloaTEM processing
mean(jul$surf_Cuncorr_Avg[jul$anom_sec == 0])
mean(jul$deep_Cuncorr_Avg[jul$anom_dec == 0])
mean(apr$surf_Cuncorr_Avg[apr$anom_sec == 0])
mean(apr$deep_Cuncorr_Avg[apr$anom_dec == 0])

#### Left/Right Profiles ####

