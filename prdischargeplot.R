library(dplyr)
library(ggplot2)

#### Load the discharge data ####
setwd("~/columbiariver")
data <- read.csv("public-data/02-raw-data/pr_discharge_samplingwindow.txt", comment.char = "#", sep = "\t")
data$datetime <- as.POSIXct(data$datetime, tz = "America/Los_Angeles")

ggplot(data) + geom_rect(aes(xmin = as.POSIXct("2021-02-16"), xmax = as.POSIXct("2021-02-20"), ymin = -Inf, ymax = +Inf), fill = 'red', alpha = 0.5)+
  geom_rect(aes(xmin = as.POSIXct("2022-04-18"), xmax = as.POSIXct("2022-04-23"), ymin = -Inf, ymax = +Inf), fill = 'red', alpha = 0.5)+
  geom_rect(aes(xmin = as.POSIXct("2021-07-19"), xmax = as.POSIXct("2021-07-24"), ymin = -Inf, ymax = +Inf), fill = 'red', alpha = 0.5)+
  geom_line(aes(x = datetime, y = X151851_00060)) + labs(x ="", y = "Discharge at Priest Rapids (cfs)")
  