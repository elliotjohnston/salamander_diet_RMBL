
######### Historical FS Density #########

rm(list = ls())

library(lubridate)
library(dplyr)
library(ggplot2)
library(scales)
library(visreg)
library(emmeans)
library(car)

##import and format data
fs_density_historical <- read.csv(file = "Data/Historical_FS_Data.csv", stringsAsFactors = FALSE)

fs_density_historical$X <- NULL
fs_density_historical$Date = mdy(fs_density_historical$Date)
fs_density_historical$Pond = as.factor(fs_density_historical$Pond)
fs_density_historical$Survey_Type = as.factor(fs_density_historical$Survey_Type)
fs_density_historical$Pond_Section = as.factor(fs_density_historical$Pond_Section)

fs_density_historical_subset = subset(fs_density_historical, Survey_Type == "Fairy_Shrimp")
fs_density_historical_subset$Notes <- NULL

#calculate density (FS/m^3)
fs_density_historical_subset$FS_Density_Transect = fs_density_historical_subset$Number_FS / fs_density_historical_subset$Sample_Volume_m3

#average transect densities for a given pond & date
fs_density_pond = fs_density_historical_subset %>% group_by(Date, Pond) %>% dplyr::summarize(FS_Density_Pond = mean(FS_Density_Transect))
fs_density_historical_subset = merge(x = fs_density_historical_subset, y = fs_density_pond, by = c("Date", "Pond"))
rm(fs_density_pond)

#take log of densities 
fs_density_historical_subset$log_FS_Density_Transect = log(fs_density_historical_subset$FS_Density_Transect + 1)

fs_density_historical_subset$log_FS_Density_Pond = log(fs_density_historical_subset$FS_Density_Pond + 1)


#standard error approach throughout this code: going to calculate uncertainty based on variability at
#the level of pond transects
standard_errors_density = fs_density_historical_subset %>% group_by(Date, Pond) %>% dplyr::summarize(SE_Density = sd(log_FS_Density_Transect)/sqrt(n()))


#plot change in density in ponds over time
fs_density_pond = fs_density_historical_subset[c("Date", "Pond", "FS_Density_Pond", "log_FS_Density_Pond")]
fs_density_pond = fs_density_pond[!duplicated(fs_density_pond), ]

fs_density_pond = merge(fs_density_pond, standard_errors_density, by = c("Date", "Pond"))
fs_density_pond$lower_density = fs_density_pond$log_FS_Density_Pond - fs_density_pond$SE_Density
fs_density_pond$upper_density = fs_density_pond$log_FS_Density_Pond + fs_density_pond$SE_Density

#subset to just ponds and dates of interest - excluding dates later in the season because FS
#decline dramatically by then
fs_density_pond_subset = subset(fs_density_pond, (Pond == 6 | Pond == 8 | Pond == 10 | Pond == 11 | 
                                Pond == 13 | Pond == 15 | Pond == 51 | Pond == 52)
                                & (Date != "2004-07-20" & Date != "2016-07-25"))

fs_plot_density_historical =
  ggplot(fs_density_pond_subset, aes(x=Date, y=log_FS_Density_Pond, fill = Pond)) +
  geom_errorbar(aes(ymin=lower_density, ymax=upper_density), color = "black", width = 0, lwd = 1) +
  geom_line(aes(group = Pond, color = Pond), lwd = 2) +
  geom_point(aes(color = Pond), size = 10, alpha = 0.8) +
  #scale_x_date(limits = as.Date(c('2021-06-15','2021-08-07'))) +
  scale_y_continuous(limits = c(0, 9.21034), breaks=c(log(1),log(10),log(100),log(1000),log(10000)), labels=c("0", "10", "100", "1,000", "10,000")) +
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x)) +
  ylab(expression(paste("Mean Fairy Shrimp Density Pond", " ", "(count/", m^{3}, ") ± SE"))) +
  theme_bw(33) +
  theme(axis.title.x = element_blank(),
        panel.grid.minor = element_blank())

fs_plot_density_historical 

### add in current data to historical data and plot

fs_density_current <- read.csv(file = "Data/fs_pond.csv", stringsAsFactors = FALSE)
fs_density_current$X <- NULL

fs_density_current = subset(fs_density_current, (Pond == 6 | Pond == 8 | Pond == 10 | Pond == 11 | 
                                                Pond == 13 | Pond == 15 | Pond == 51 | Pond == 52) 
                                                & (Date == "2021-06-30" | Date == "2021-06-29" |
                                                   Date == "2022-06-27"))

fs_density_current = fs_density_current[,c("Date", "Pond", "FS_Density_Pond", "log_FS_Density_Pond", "SE_Density", "lower_density", "upper_density")]

fs_density_pond_subset = rbind(fs_density_pond_subset, fs_density_current)


xlimits <- as.POSIXct(strptime(c("2002-01-01", "2023-01-01"), 
                               format = "%Y-%m-%d"))

xbreaks <- as.POSIXct(strptime(c("2002-01-01", "2006-01-01", "2010-01-01",
                                 "2014-01-01", "2018-01-01", "2022-01-01"), 
                               format = "%Y-%m-%d"))

fs_density_pond_subset$Date = as.POSIXct(fs_density_pond_subset$Date, format="%Y-%m-%d")

fs_plot_density_historical_v2 =
  ggplot(fs_density_pond_subset, aes(x=Date, y=log_FS_Density_Pond, fill = Pond)) +
  geom_errorbar(aes(ymin=lower_density, ymax=upper_density), color = "black", width = 0, lwd = 1) +
  geom_line(aes(group = Pond, color = Pond), lwd = 2) +
  geom_point(aes(color = Pond), size = 10, alpha = 0.8) +
  scale_x_datetime(breaks = xbreaks,
                   labels = c("2002", "2006", "2010", "2014", "2018", "2022"),
                   limits = xlimits) +
  scale_y_continuous(limits = c(0, 9.21034), breaks=c(log(1),log(10),log(100),log(1000),log(10000)), labels=c("0", "10", "100", "1,000", "10,000")) +
  ylab(expression(paste("Mean Fairy Shrimp Density Pond", " ", "(count/", m^{3}, ") ± SE"))) +
  theme_bw(33) +
  theme(axis.title.x = element_blank(),
        panel.grid.minor = element_blank())

fs_plot_density_historical_v2

#ggsave(fs_plot_density_historical_v2, filename = "Outputs/fs_density_plot_historical.png",  width = 18, height = 12, dpi = 100)


## Low FS Ponds
low_fs = subset(fs_density_pond_subset, Pond == "6" | Pond == "10" | Pond == "11")

hex_codes1 <- hue_pal()(8) #get hex codes used in above plot                           
hex_codes1
show_col(hex_codes1)

fs_plot_density_historical_low =
  ggplot(low_fs, aes(x=Date, y=log_FS_Density_Pond, fill = Pond)) +
  scale_fill_manual(values=c("#F8766D", "#7CAE00", "#00BE67"), aesthetics = c("colour", "fill")) +
  geom_errorbar(aes(ymin=lower_density, ymax=upper_density), color = "black", width = 0, lwd = 1) +
  geom_line(aes(group = Pond, color = Pond), lwd = 2) +
  geom_point(aes(color = Pond), size = 10, alpha = 0.8) +
  scale_x_datetime(breaks = xbreaks,
                   labels = c("2002", "2006", "2010", "2014", "2018", "2022"),
                   limits = xlimits) +
  scale_y_continuous(limits = c(0, 9.21034), breaks=c(log(1),log(10),log(100),log(1000),log(10000)), labels=c("0", "10", "100", "1,000", "10,000")) +
  ylab(expression(paste("Mean Fairy Shrimp Density Pond", " ", "(count/", m^{3}, ") ± SE"))) +
  theme_bw(33) +
  theme(axis.title.x = element_blank(),
        panel.grid.minor = element_blank())

fs_plot_density_historical_low

#ggsave(fs_plot_density_historical_low, filename = "Outputs/fs_plot_density_historical_low.png",  width = 18, height = 12, dpi = 100)


## High FS Ponds
high_fs = subset(fs_density_pond_subset, Pond == "13" | Pond == "15" | Pond == "52")

hex_codes1 <- hue_pal()(8)                       
hex_codes1
show_col(hex_codes1)

fs_plot_density_historical_high =
  ggplot(high_fs, aes(x=Date, y=log_FS_Density_Pond, fill = Pond)) +
  scale_fill_manual(values=c("#00BFC4", "#00A9FF", "#FF61CC"), aesthetics = c("colour", "fill")) +
  geom_errorbar(aes(ymin=lower_density, ymax=upper_density), color = "black", width = 0, lwd = 1) +
  geom_line(aes(group = Pond, color = Pond), lwd = 2) +
  geom_point(aes(color = Pond), size = 10, alpha = 0.8) +
  scale_x_datetime(breaks = xbreaks,
                   labels = c("2002", "2006", "2010", "2014", "2018", "2022"),
                   limits = xlimits) +
  scale_y_continuous(limits = c(0, 9.21034), breaks=c(log(1),log(10),log(100),log(1000),log(10000)), labels=c("0", "10", "100", "1,000", "10,000")) +
  ylab(expression(paste("Mean Fairy Shrimp Density Pond", " ", "(count/", m^{3}, ") ± SE"))) +
  theme_bw(33) +
  theme(axis.title.x = element_blank(),
        panel.grid.minor = element_blank())

fs_plot_density_historical_high

#ggsave(fs_plot_density_historical_high, filename = "Outputs/fs_plot_density_historical_high.png",  width = 18, height = 12, dpi = 100)


###### Correlation Across Yrs ######
fs_density_pond_subset$Year = year(fs_density_pond_subset$Date)
fs_density_pond_subset$Year = as.factor(fs_density_pond_subset$Year)
fs_density_pond_subset = fs_density_pond_subset %>% group_by(Year) %>% mutate(FS_Rank = order(order(FS_Density_Pond, decreasing=TRUE)))

fs_yrs_density = glm(FS_Rank ~ Pond, family = "poisson", data = fs_density_pond_subset)
plot(fs_yrs_density)
car::Anova(fs_yrs_density)


###### Mean Densities ######

fs_pond_rank = fs_density_pond_subset %>% group_by(Pond) %>% summarise(Mean_FS_Rank = mean(FS_Rank))

#not using below code anymore#
fs_density_pond_subset = fs_density_pond_subset %>% group_by(Pond) %>% mutate(log_FS_Density_Pond_Mean = mean(log_FS_Density_Pond))
fs_density_pond_subset = fs_density_pond_subset %>% group_by(Pond) %>% mutate(log_FS_Density_Pond_SE = sd(log_FS_Density_Pond)/sqrt(n()))
fs_density_pond_subset$log_FS_Density_Pond_Mean_lwr = fs_density_pond_subset$log_FS_Density_Pond_Mean - fs_density_pond_subset$log_FS_Density_Pond_SE*1.96
fs_density_pond_subset$log_FS_Density_Pond_Mean_upr = fs_density_pond_subset$log_FS_Density_Pond_Mean + fs_density_pond_subset$log_FS_Density_Pond_SE*1.96
fs_density_pond_subset$log_FS_Density_Pond_Mean_lwr[fs_density_pond_subset$log_FS_Density_Pond_Mean_lwr < 0] <- 0

fs_historical_mean = unique(fs_density_pond_subset[,c("Pond", "log_FS_Density_Pond_Mean", "log_FS_Density_Pond_SE", "log_FS_Density_Pond_Mean_lwr", "log_FS_Density_Pond_Mean_upr")])

sd(fs_density_pond_subset$FS_Density_Pond)
