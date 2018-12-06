library(tidyverse)
library(lubridate)
setwd("/Users/walkerharrison/563_Project")
logs <- read_csv('logs.csv')
names(logs) <- c("pid", "date", "mem", "cpu")

logs.mem <- logs %>% select(-cpu)

logs.mem.time <- logs.mem %>%
  group_by(pid, time = paste(hour(date), minute(date))) %>%
  summarize(mem = mean(mem, na.rm = T)) %>%
  spread(pid, mem) %>%
  select(-time)
  
X <- unname(as.matrix(logs.mem.time[, !is.na(colSums(logs.mem.time))]))


