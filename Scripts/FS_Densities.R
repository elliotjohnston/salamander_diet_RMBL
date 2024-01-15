
######### FS Density #########

rm(list = ls())
library(lubridate)
library(glmmTMB)
library(dplyr)
library(ggplot2)
library(scales)
library(visreg)
library(emmeans)
library(car)
library(lme4)
library(DHARMa)

##import and format data
fs_density_transect <- read.csv(file = "Data/FS_Count.csv", stringsAsFactors = FALSE)
fs_biomass <- read.csv(file = "Data/FS_Biomass.csv", stringsAsFactors = FALSE)

fs_density_transect$X <- NULL
fs_density_transect$Pond = as.factor(fs_density_transect$Pond)
fs_density_transect$Hydroperiod = as.factor(fs_density_transect$Hydroperiod)
fs_biomass$Pond = as.factor(fs_biomass$Pond)
fs_density_transect$Pond_Section = as.factor(fs_density_transect$Pond_Section)
fs_biomass$Pond_Section = as.factor(fs_biomass$Pond_Section)

#convert FS length from cm to mm
fs_biomass$Length_cm = fs_biomass$Length_cm*10
names(fs_biomass)[names(fs_biomass) == 'Length_cm'] <- 'Length_mm'

##add in Biomass values to fs_biomass df based on LW Regression formula
#Formula: log(y) = 2.22706*log(x) - 3.69561
#y = Biomass, x = Length
fs_biomass$Biomass_mg = 2.22706*(log(fs_biomass$Length_mm)) - 3.69561
fs_biomass$Biomass_mg = exp(fs_biomass$Biomass_mg)

summary(fs_biomass$Biomass_mg) #values look sensible

#any Notes that require attention before working with data? Notes the same between two dataframes
fs_density_transect$Notes = as.factor(fs_density_transect$Notes)
levels(fs_density_transect$Notes) #no concerns

#drop 6/17/21 surveys. Day of piloting methods. Also didn't take photos of FS.
fs_density_transect$Date = mdy(fs_density_transect$Date)
fs_biomass$Date = mdy(fs_biomass$Date)

fs_density_transect = subset(fs_density_transect, Date != "2021-06-17")
fs_biomass = subset(fs_biomass, Date != "2021-06-17")

#drop ponds that only sampled once when figuring out which ponds to survey. Drop ponds 2,3,4,7,9,44,56.
fs_density_transect = subset(fs_density_transect, Pond == 6 | Pond == 8 | Pond == 10 | Pond == 11 | Pond == 13 | Pond == 15 | Pond == 51 | Pond == 52 | Pond == 55)
fs_biomass = subset(fs_biomass, Pond == 6 | Pond == 8 | Pond == 10 | Pond == 11 | Pond == 13 | Pond == 15 | Pond == 51 | Pond == 52 | Pond == 55)

#need to remove one transect from dataframes that was not sampled. Pond 8 South, 6/18/2021. Seems OK to calculate pond 
#average from 7 rather than 8 transects for this one date.
fs_density_transect = subset(fs_density_transect, Number_FS != "Not sampled")
fs_biomass = subset(fs_biomass, Number_FS != "Not sampled")

#merge count data from biomass df into count df
fs_count = fs_biomass[,c("Date", "Pond", "Pond_Section", "FS_ID")]
fs_count$FS_ID[is.na(fs_count$FS_ID)] <- 0 #set NAs to 0
fs_count = fs_count %>% group_by(Date, Pond, Pond_Section) %>% summarise(Number_FS = max(FS_ID)) #extract max ID number for each group, which equals count

fs_density_transect = merge(x = fs_density_transect, y = fs_count, by = c("Date", "Pond", "Pond_Section"))
fs_density_transect$Number_FS.x <- NULL
names(fs_density_transect)[names(fs_density_transect) == "Number_FS.y"] <- "Number_FS"
fs_density_transect = fs_density_transect %>% relocate(Number_FS, .before = Picture)

#set aside opportunistic Max transects
fs_max = subset(fs_density_transect, Pond_Section == "East(Max)")
fs_density_transect = subset(fs_density_transect, Pond_Section != "East(Max)")

#calculate density (FS/m^3)
fs_density_transect$Sample_Volume = as.numeric(fs_density_transect$Sample_Volume)
fs_density_transect$FS_Density_Tran = fs_density_transect$Number_FS / fs_density_transect$Sample_Volume
fs_density_transect = fs_density_transect %>% relocate(FS_Density_Tran, .before = Picture)

#average transect densities for a given pond & date
fs_density_pond = fs_density_transect %>% group_by(Date, Pond) %>% dplyr::summarize(FS_Density_Pond = mean(FS_Density_Tran))
fs_density_transect = merge(x = fs_density_transect, y = fs_density_pond, by = c("Date", "Pond"))
fs_density_transect = fs_density_transect %>% relocate(FS_Density_Pond, .after = FS_Density_Tran)
rm(fs_density_pond)

#take log of densities 
fs_density_transect$FS_Density_Tran = fs_density_transect$FS_Density_Tran + 1 #add 1 for log transformation
fs_density_transect$log_FS_Density_Tran = log(fs_density_transect$FS_Density_Tran)
fs_density_transect$FS_Density_Tran = fs_density_transect$FS_Density_Tran - 1 #undo 
fs_density_transect = fs_density_transect %>% relocate(log_FS_Density_Tran, .after = FS_Density_Tran)

fs_density_transect$FS_Density_Pond = fs_density_transect$FS_Density_Pond + 1
fs_density_transect$log_FS_Density_Pond = log(fs_density_transect$FS_Density_Pond)
fs_density_transect$FS_Density_Pond = fs_density_transect$FS_Density_Pond - 1
fs_density_transect = fs_density_transect %>% relocate(log_FS_Density_Pond, .after = FS_Density_Pond)

#standard error approach throughout this code: going to calculate uncertainty based on variability at
#the level of pond transects. Finer-scale measurements (e.g., biomass of each FS at a given transect)
#will be averaged to compute an transect-level estimate.
standard_errors_density = fs_density_transect %>% group_by(Date, Pond) %>% dplyr::summarize(SE_Density = sd(log_FS_Density_Tran)/sqrt(n()))


#plot change in density in ponds over time
fs_density_pond = fs_density_transect[c("Date", "Pond", "Hydroperiod", "FS_Density_Pond", "log_FS_Density_Pond")]
fs_density_pond = fs_density_pond[!duplicated(fs_density_pond), ]

fs_density_pond = merge(fs_density_pond, standard_errors_density, by = c("Date", "Pond"))
fs_density_pond$lower_density = fs_density_pond$log_FS_Density_Pond - fs_density_pond$SE_Density*1.96
fs_density_pond$upper_density = fs_density_pond$log_FS_Density_Pond + fs_density_pond$SE_Density*1.96

fs_density_pond_temporary = subset(fs_density_pond, Hydroperiod == "Temporary")

fs_density_pond$Week = week(fs_density_pond$Date)
fs_density_pond = fs_density_pond %>% relocate(Week, .after = Date)
fs_density_pond$Week = as.factor(fs_density_pond$Week)

fs_density_pond = fs_density_pond %>% mutate(Ordinal_Day = lubridate::yday(Date))
fs_density_pond = fs_density_pond %>% relocate(Ordinal_Day, .after = Date)

fs_density_pond$Year = year(fs_density_pond$Date)
fs_density_pond = fs_density_pond %>% relocate(Year, .before = Date)
fs_density_pond$Year = as.factor(fs_density_pond$Year)


fs_plot_density =
  ggplot(fs_density_pond, aes(x=Ordinal_Day, y=log_FS_Density_Pond)) +
  facet_wrap(~Pond)+
  geom_errorbar(aes(ymin=lower_density, ymax=upper_density), color = "black", width = 0, lwd = 2) +
  geom_line(aes(group = Year, color = Year), lwd = 3) +
  geom_point(aes(color = Year), size = 10, alpha = 0.8) +
  scale_y_continuous(breaks=c(log(1),log(10),log(100),log(1000),log(10000), log(100000)), labels=c("0", "10", "100", "1,000", "10,000", "100,000")) +
  scale_x_continuous(limits = c(165,217), breaks=c(166, 182, 196, 213), labels=c("Jun 15", "Jul 1", "Jul 15", "Aug 1")) +
  ylab(expression(paste("Mean Fairy Shrimp Density Pond", " ", "(count/", m^{3}, ") ± 95% CI"))) +
  scale_colour_manual(values = c("orange", "steelblue"))  + 
  theme_bw(30) +
  theme(axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(2, "lines"))

fs_plot_density

##stats
fs_density_transect$Year = year(fs_density_transect$Date)
fs_density_transect = fs_density_transect %>% relocate(Year, .before = Date)
fs_density_transect$Year = as.factor(fs_density_transect$Year)

fs_density_transect$Week = week(fs_density_transect$Date)
fs_density_transect = fs_density_transect %>% relocate(Week, .after = Date)
fs_density_transect$Week = as.numeric(as.character(fs_density_transect$Week))

fs_density_transect$Transect_ID = interaction(fs_density_transect$Pond, fs_density_transect$Pond_Section, sep="_")
fs_density_transect = fs_density_transect %>% relocate(Transect_ID, .after = Pond_Section)

fs_density_transect$FS_Density_Tran_PresAbs <- ifelse(fs_density_transect[,c("FS_Density_Tran")] == 0, 0, 1)

m1_density = glmmTMB(FS_Density_Tran_PresAbs ~ Year + Pond + I(Week^2) + (1|Transect_ID), family = binomial, data = fs_density_transect)

m1_density_simres <- simulateResiduals(fittedModel = m1_density)
plot(m1_density_simres) #looks good
hist(m1_density_simres) #looks good

summary(m1_density)

car::Anova(m1_density)
emmeans(m1_density, c("Pond"), type = "response")
emmeans(m1_density, c("Year"), type = "response")

#conditional model
fs_density_transect_no_zeros = subset(fs_density_transect, log_FS_Density_Tran > 0)

m2_density = lmer(log(FS_Density_Tran + 1) ~ Pond + Week + Year + I(Week^2) + Pond:Week + Pond:I(Week^2) 
                  + (1|Transect_ID), data = fs_density_transect_no_zeros)


plot(m2_density)

summary(m2_density)

car::Anova(m2_density)

m2_emmeans1 = emmeans(m2_density, pairwise ~ Pond | Week, type = "response")
contrast(m2_emmeans1) #same as below

emmeans(m2_density, c("Year"), type = "response")
emmeans(m2_density, c("Pond"), type = "response")
emmeans(m2_density, c("Pond", "Week"), type = "response")

#hard to say if FS abundances are actually higher in ponds with the highest densities (e.g., 15). Could
#be a function of ponds with really small volumes. BUT even if this skews estimates a bit, density still
#gives us a good idea of salamander encounter rate 


######### FS Biomass #########

#need to correct biomass estimates for sampling intensity -> biomass/m^3
fs_biomass_transect = fs_biomass %>% group_by(Date, Pond, Pond_Section) %>% dplyr::summarize(FS_Biomass_Tran = sum(Biomass_mg))

fs_biomass_transect = merge(fs_biomass_transect, fs_density_transect[,c("Date", "Pond", "Pond_Section", "Sample_Volume")], by = c("Date", "Pond", "Pond_Section"))

#calculate density (Biomass/m^3)
fs_biomass_transect$FS_Biomass_Tran_m3 = fs_biomass_transect$FS_Biomass_Tran / fs_biomass_transect$Sample_Volume

#average transect densities for a given pond & date
fs_biomass_pond = fs_biomass_transect %>% group_by(Date, Pond) %>% dplyr::summarize(FS_Biomass_Pond_m3 = mean(FS_Biomass_Tran_m3))
fs_biomass_transect = merge(x = fs_biomass_transect, y = fs_biomass_pond, by = c("Date", "Pond"))

rm(fs_biomass_pond)

#take log of biomass densities 
fs_biomass_transect$FS_Biomass_Tran_m3 = fs_biomass_transect$FS_Biomass_Tran_m3 + 1
fs_biomass_transect$log_FS_Biomass_Tran_m3 = log(fs_biomass_transect$FS_Biomass_Tran_m3)
fs_biomass_transect$FS_Biomass_Tran_m3 = fs_biomass_transect$FS_Biomass_Tran_m3 - 1
fs_biomass_transect = fs_biomass_transect %>% relocate(log_FS_Biomass_Tran_m3, .after = FS_Biomass_Tran_m3)

fs_biomass_transect$FS_Biomass_Pond_m3 = fs_biomass_transect$FS_Biomass_Pond_m3 + 1
fs_biomass_transect$log_FS_Biomass_Pond_m3 = log(fs_biomass_transect$FS_Biomass_Pond_m3)
fs_biomass_transect$FS_Biomass_Pond_m3 = fs_biomass_transect$FS_Biomass_Pond_m3 - 1
fs_biomass_transect = fs_biomass_transect %>% relocate(log_FS_Biomass_Pond_m3, .after = FS_Biomass_Pond_m3)

standard_errors_biomass = fs_biomass_transect %>% group_by(Date, Pond) %>% dplyr::summarize(SE_Biomass = sd(log_FS_Biomass_Tran_m3)/sqrt(n()))


#plot change in density in ponds over time
fs_biomass_pond = fs_biomass_transect[c("Date", "Pond", "FS_Biomass_Pond_m3", "log_FS_Biomass_Pond_m3")]
fs_biomass_pond = fs_biomass_pond[!duplicated(fs_biomass_pond), ]

fs_biomass_pond = merge(fs_biomass_pond, standard_errors_biomass, by = c("Date", "Pond"))
fs_biomass_pond$lower_biomass = fs_biomass_pond$log_FS_Biomass_Pond_m3 - fs_biomass_pond$SE_Biomass*1.96
fs_biomass_pond$upper_biomass = fs_biomass_pond$log_FS_Biomass_Pond_m3 + fs_biomass_pond$SE_Biomass*1.96

fs_biomass_pond$Year = year(fs_biomass_pond$Date)
fs_biomass_pond = fs_biomass_pond %>% relocate(Year, .before = Date)
fs_biomass_pond$Year = as.factor(fs_biomass_pond$Year)

fs_biomass_pond = fs_biomass_pond %>% mutate(Ordinal_Day = lubridate::yday(Date))
fs_biomass_pond = fs_biomass_pond %>% relocate(Ordinal_Day, .after = Date)


fs_plot_biomass =
  ggplot(fs_biomass_pond, aes(x=Ordinal_Day, y=log_FS_Biomass_Pond_m3)) +
  facet_wrap(~Pond)+
  geom_errorbar(aes(ymin=lower_biomass, ymax=upper_biomass), color = "black", width = 0, lwd = 1) +
  geom_line(aes(group = Year, color = Year), lwd = 1.5) +
  geom_point(aes(color = Year), size = 6, alpha = 0.8) +
  scale_y_continuous(breaks=c(log(1),log(10),log(100),log(1000),log(10000), log(100000)), labels=c("0", "10", "100", "1,000", "10,000", "100,000")) +
  scale_x_continuous(limits = c(165,217), breaks=c(166, 182, 196, 213), labels=c("Jun 15", "Jul 1", "Jul 15", "Aug 1")) +
  ylab(expression(paste("Mean Fairy Shrimp Biomass Pond", " ", "(mg/", m^{3}, ") ± 95% CI"))) +
  scale_colour_manual(values = c("blue", "red"))  + 
  theme_bw(30) +
  theme(axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(2, "lines"))

fs_plot_biomass

ggsave(fs_plot_biomass, filename = "Outputs/fs_plot_biomass.png",  width = 18, height = 12, dpi = 100)


##stats
fs_biomass_transect$Year = year(fs_biomass_transect$Date)
fs_biomass_transect = fs_biomass_transect %>% relocate(Year, .before = Date)
fs_biomass_transect$Year = as.factor(fs_biomass_transect$Year)

fs_biomass_transect$Week = week(fs_biomass_transect$Date)
fs_biomass_transect = fs_biomass_transect %>% relocate(Week, .after = Date)
fs_biomass_transect$Week = as.numeric(as.character(fs_biomass_transect$Week))

fs_biomass_transect$Transect_ID = interaction(fs_biomass_transect$Pond, fs_biomass_transect$Pond_Section, sep="_")
fs_biomass_transect = fs_biomass_transect %>% relocate(Transect_ID, .after = Pond_Section)

#I have semi-continuous zero-inflated data. 
#could eventually go down the hurdle model route (https://drizopoulos.github.io/GLMMadaptive/articles/ZeroInflated_and_TwoPart_Models.html)
#but for now going to analyze data separately as binary (all) and conditional (positive values only)
fs_biomass_transect$FS_Biomass_Tran_PresAbs <- ifelse(fs_biomass_transect[,c("FS_Biomass_Tran_m3")] == 0, 0, 1)

m1_biomass = glmmTMB(FS_Biomass_Tran_PresAbs ~ Year + Pond + I(Week^2) + (1|Transect_ID), family = binomial, data = fs_biomass_transect)
#model won't converge with both linear and quadratic Week terms, so just including quadratic

m1_biomass_simres <- simulateResiduals(fittedModel = m1_biomass)
plot(m1_biomass_simres) #looks good
hist(m1_biomass_simres) #looks good

summary(m1_biomass)

car::Anova(m1_biomass)
emmeans(m1_biomass, c("Pond"), type = "response")
emmeans(m1_biomass, c("Year"), type = "response")

#conditional model
fs_biomass_transect_no_zeros = subset(fs_biomass_transect, log_FS_Biomass_Tran_m3 > 0)

m2_biomass = lmer(log(FS_Biomass_Tran_m3 + 1) ~ Pond + Week + Year + I(Week^2) + Pond:Week + Pond:I(Week^2) 
                  + (1|Transect_ID), data = fs_biomass_transect_no_zeros)


plot(m2_biomass)

summary(m2_biomass)

car::Anova(m2_biomass)

m2_emmeans1 = emmeans(m2_biomass, pairwise ~ Pond | Week, type = "response")
contrast(m2_emmeans1) #same as below

emmeans(m2_biomass, c("Year"), type = "response")
emmeans(m2_biomass, c("Pond"), type = "response")
emmeans(m2_biomass, c("Pond", "Week"), type = "response")

#good resource: https://stats.oarc.ucla.edu/r/seminars/interactions-r

#write.csv(fs_biomass_transect, "Data/fs_biomass_transect.csv")


#visreg plots contains a prediction line, confidence band, and partial residuals (?)
visreg(m1_biomass, xvar = "Pond", whitespace = 0.4,
       points.par = list(cex = 1.1, col = "red"))

#emmeans(m1_biomass, "Pond", data = fs_biomass_pond)



## **code for analyses below has not been checked for compatability with 2022 data** ##


######### FS Biomass Per Individual #########

#calculate biomass/individual at the transect level, and then calculate SEs

fs_biomass_transect = merge(fs_biomass_transect, fs_count, by = c("Date", "Pond", "Pond_Section"))
names(fs_biomass_transect)[names(fs_biomass_transect) == "Number_FS"] <- "Number_FS_Tran"
fs_biomass_transect$Biomass_Per_FS_Tran = fs_biomass_transect$FS_Biomass_Tran / fs_biomass_transect$Number_FS_Tran

fs_biomass_transect$Biomass_Per_FS_Tran[is.nan(fs_biomass_transect$Biomass_Per_FS_Tran)] <- 0 #set NaN to 0

standard_errors_indiv = fs_biomass_transect %>% group_by(Date, Pond) %>% dplyr::summarize(SE_Indiv = sd(Biomass_Per_FS_Tran)/sqrt(n()))

#average transects for a given pond & date
fs_biomass_scratch = fs_biomass_transect %>% group_by(Date, Pond) %>% dplyr::summarize(Biomass_Per_FS_Pond = mean(Biomass_Per_FS_Tran))
fs_biomass_transect = merge(x = fs_biomass_transect, y = fs_biomass_scratch, by = c("Date", "Pond"))

rm(fs_biomass_scratch)

#plot change in density in ponds over time
fs_biomass_per_indiv_pond = fs_biomass_transect[c("Date", "Pond", "Biomass_Per_FS_Pond")]
fs_biomass_per_indiv_pond = fs_biomass_per_indiv_pond[!duplicated(fs_biomass_per_indiv_pond), ]


fs_biomass_per_indiv_pond = merge(fs_biomass_per_indiv_pond, standard_errors_indiv, by = c("Date", "Pond"))
fs_biomass_per_indiv_pond$lower_indiv = fs_biomass_per_indiv_pond$Biomass_Per_FS_Pond - fs_biomass_per_indiv_pond$SE_Indiv
fs_biomass_per_indiv_pond$upper_indiv = fs_biomass_per_indiv_pond$Biomass_Per_FS_Pond + fs_biomass_per_indiv_pond$SE_Indiv

fs_biomass_per_FS =
  ggplot(fs_biomass_per_indiv_pond, aes(x=Date, y=Biomass_Per_FS_Pond, fill = Pond)) +
  geom_errorbar(aes(ymin=lower_indiv, ymax=upper_indiv), color = "black", width = 0, lwd = 1) +
  geom_line(aes(group = Pond, color = Pond), lwd = 2) +
  geom_point(aes(color = Pond), size = 10, alpha = 0.8) +
  scale_x_date(limits = as.Date(c('2021-06-15','2021-08-07'))) +
  ylim(0,6) +
  ylab("Mean Fairy Shrimp Biomass/Individual in Pond (mg) ± SE") +
  theme_bw(30) +
  theme(axis.title.x = element_blank())

fs_biomass_per_FS 

#looks like increasing variability over time in a few ponds. Interesting. Could either be earliest hatching
#FS getting big and new, small FS hatching later throughout summer. Or, FS could all hatch around the same
#time and there are winners and losers in terms of growth.

ggsave(fs_biomass_per_FS, filename = "Outputs/fs_biomass_per_FS.png",  width = 18, height = 12, dpi = 100)


#stats
m1_biomass_per_FS = lm(Biomass_Per_FS_Pond ~ Date*Pond, data = fs_biomass_per_indiv_pond)
#plot(m1_biomass) #looks good

car::Anova(m1_biomass_per_FS, type = 3) #Pond and Date:Pond significant
#Date might be significant if took out zero points?

#visreg plots contains a prediction line, confidence band, and partial residuals (?)
visreg(m1_biomass_per_FS, xvar = "Date", by = "Pond", whitespace = 0.4,
       points.par = list(cex = 1.1, col = "red"))

visreg(m1_biomass_per_FS, xvar = "Pond", whitespace = 0.4,
       points.par = list(cex = 1.1, col = "red"))

emmeans(m1_biomass_per_FS, "Pond", data = fs_density_pond)


######### FS BPI ~ FS Density #########

#make fs pond and fs transect dataframes
fs_pond = merge(fs_biomass_pond, fs_biomass_per_indiv_pond, by = c("Date", "Pond"))
fs_pond = merge(fs_pond, fs_density_pond, by = c("Date", "Pond"))
fs_pond = fs_pond %>% relocate(c("FS_Density_Pond", "log_FS_Density_Pond", "SE_Density", "lower_density", "upper_density"), .after = Pond)

fs_transect = merge(x = fs_biomass_transect[,c("Date", "Pond", "Pond_Section", "FS_Biomass_Tran_m3", "log_FS_Biomass_Tran_m3", "Number_FS_Tran", "Biomass_Per_FS_Tran")],
                    y = fs_density_transect[,c("Date", "Pond", "Pond_Section", "FS_Density_Tran", "log_FS_Density_Tran")],
                    by = c("Date", "Pond", "Pond_Section"))

fs_transect = fs_transect[,c("Date", "Pond", "Pond_Section", "Number_FS_Tran", "FS_Density_Tran", "log_FS_Density_Tran", "FS_Biomass_Tran_m3", "log_FS_Biomass_Tran_m3", "Biomass_Per_FS_Tran")]


fs_meansize_v_density =
  ggplot(fs_pond, aes(x=log_FS_Density_Pond, y=Biomass_Per_FS_Pond, fill = Pond)) +
  geom_errorbar(aes(ymin=lower_indiv, ymax=upper_indiv), color = "black", width = 0, lwd = 1) +
  geom_errorbarh(aes(xmin=lower_density, xmax=upper_density), color = "black", height = 0, size = 1) +
  geom_point(aes(color = Pond), size=10, alpha = 0.8) +
  scale_x_continuous(breaks=c(log(1),log(10),log(100),log(1000),log(10000)), labels=c("0", "10", "100", "1,000", "10,000")) +
  xlab(expression(paste("Mean Fairy Shrimp Density Pond", " ", "(count/", m^{3}, ") ± SE"))) +
  ylab("Mean Fairy Shrimp Biomass/Individual in Pond (mg) ± SE") +
  ylim(0,6) +
  theme_bw(42) +
  theme(axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        panel.grid.minor = element_blank())

print(fs_meansize_v_density)

ggsave(fs_meansize_v_density, filename = "Outputs/fs_meansize_v_density.png",  width = 18, height = 12, dpi = 100)

#mean FS size is more variable at lower densities! Biological interpretation? How much of this is driven
#by low samples sizes (e.g, pond 55)? Put sample sizes next to pond in legend when present this.

#stats
#what to run? Quantile regression?

######### FS Max Size ~ FS Density #########

fs_biomass_max_pond = fs_biomass %>% group_by(Date, Pond) %>% summarise(FS_Biomass_Max_Pond = max(Biomass_mg))
fs_biomass_max_tran = fs_biomass %>% group_by(Date, Pond, Pond_Section) %>% summarise(FS_Biomass_Max_Tran = max(Biomass_mg))

fs_pond = merge(fs_pond, fs_biomass_max_pond, by = c("Date", "Pond"))
fs_transect = merge(fs_transect, fs_biomass_max_tran, by = c("Date", "Pond", "Pond_Section"))

fs_maxsize_v_density =
  ggplot(fs_transect, aes(x=log_FS_Density_Tran, y=FS_Biomass_Max_Tran, fill = Pond)) +
  #geom_errorbarh(aes(xmin=lower_density, xmax=upper_density), color = "black", height = 0, size = 1) +
  geom_point(aes(color = Pond), size=10, alpha = 0.8) +
  scale_x_continuous(breaks=c(log(1),log(10),log(100),log(1000),log(10000)), labels=c("0", "10", "100", "1,000", "10,000")) +
  xlab(expression(paste("Fairy Shrimp Density Transect", " ", "(count/", m^{3}, ")"))) +
  ylab("Max Fairy Shrimp Biomass Transect (mg)") +
  theme_bw(42) +
  theme(axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        panel.grid.minor = element_blank())

print(fs_maxsize_v_density)

ggsave(fs_maxsize_v_density, filename = "Outputs/fs_maxsize_v_density.png",  width = 18, height = 12, dpi = 100)

#stats
#what to run? Quantile regression?

######### FS Mode Size ~ FS Density #########

# given that size is a quantitative trait and there probably
#aren't that many individuals per size, I am going to round to the nearest 0.1 decimal place

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

fs_biomass_mode_pond = fs_biomass %>% group_by(Date, Pond) %>% summarise(FS_Biomass_Mode_Pond = getmode(round(Biomass_mg, digits = 1)))
fs_biomass_mode_tran = fs_biomass %>% group_by(Date, Pond, Pond_Section) %>% summarise(FS_Biomass_Mode_Tran = getmode(round(Biomass_mg, digits = 1)))

fs_pond = merge(fs_pond, fs_biomass_mode_pond, by = c("Date", "Pond"))
fs_transect = merge(fs_transect, fs_biomass_mode_tran, by = c("Date", "Pond", "Pond_Section"))

fs_modesize_v_density =
  ggplot(fs_transect, aes(x=log_FS_Density_Tran, y=FS_Biomass_Mode_Tran, fill = Pond)) +
  #geom_errorbarh(aes(xmin=lower_density, xmax=upper_density), color = "black", height = 0, size = 1) +
  geom_point(aes(color = Pond), size=10, alpha = 0.8) +
  scale_x_continuous(breaks=c(log(1),log(10),log(100),log(1000),log(10000)), labels=c("0", "10", "100", "1,000", "10,000")) +
  xlab(expression(paste("Fairy Shrimp Density Transect", " ", "(count/", m^{3}, ")"))) +
  ylab("Mode Fairy Shrimp Biomass Transect (mg)") +
  theme_bw(42) +
  theme(axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22))

print(fs_modesize_v_density)

ggsave(fs_modesize_v_density, filename = "Outputs/fs_modesize_v_density.png",  width = 18, height = 12, dpi = 100)

#not much of a different story than mean and max plots

#write file with FS info for use in subsequent R scripts
#write.csv(fs_pond, "Data/fs_pond.csv")
#write.csv(fs_transect, "Data/fs_transect.csv")


