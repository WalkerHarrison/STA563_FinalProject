#### mider

setwd("/Users/walkerharrison/563_Project")
file <- 'b1.csv'
source('estimateH2.R') 
source('estimateH3.R')
source('estimateH4.R')
source('mider.R')
source('plot_mider.R')

X <- read.csv(file, header = FALSE)

dat <- mider(X, 10, 3)

plot_mider(dat$dist, dat$con_array, dat$T, 1:10)



logs <- read_csv('logs_big.csv')

logs.cpu.time <- logs %>%
  group_by(pid, time = paste(hour(time), minute(time))) %>%
  summarize(cpu = mean(cpu, na.rm = T)) %>%
  spread(pid, cpu) %>%
  select(-time)

X <- unname(as.matrix(logs.cpu.time[, !is.na(colSums(logs.cpu.time))]))

dat <- mider(X, 10, 3)
attach(dat)

pids <- names(logs.cpu.time)[!is.na(colSums(logs.cpu.time))]
procname = c("WindowServer",
             "mds_stores",
             "Finder",
             "Spotlight",
             "DropBox",
             "CloudApp",
             "reversetemplated",
             "Google Chrome",
             "Google Chrome Helper",
             "RStudio",
             "WebKit",
             "WebKit",
             "Messages",
             "WebKit")

plot_mider(dist, con_array, T, 1:16)

Xt <- as.tibble(X); names(Xt) <- c('Citr-M',	'AMP-M',	'Pi',	'F26BP',	'F16BP',	'DHAP',	'F6P',	'G6P',	'Citr-I',	'AMP-I')
X
