
##### Salamander Movement ######

rm(list = ls())

library(lubridate)
library(dplyr)
library(ggplot2)
library(ajfhelpR)
library(plotly)
library(tidyr)
library(fishmethods)
library(FSA)
library(nlstools)
library(lme4)
library(emmeans)
library(paletteer)
library(MASS)
library(psych)
library(Hmisc)
detach(package:Hmisc)
library(ggeffects)
library(sjPlot)

salamander <- read.csv(file = "Data/Salamander_Master.csv", stringsAsFactors = FALSE)

#assign column classes
salamander$Field1 <- NULL
salamander$ID <- as.factor(salamander$ID)
salamander$INDIV_ID <- as.factor(salamander$INDIV_ID)
salamander$DATEOFYEAR = dmy(salamander$DATEOFYEAR)
salamander$RECAPTURE <- as.factor(salamander$RECAPTURE)
salamander$MORPH <- as.factor(salamander$MORPH)
salamander$SEX <- as.factor(salamander$SEX)
salamander$CLIP <- as.factor(salamander$CLIP)
salamander$CONDITION <- as.factor(salamander$CONDITION)
salamander$CLOACA <- as.factor(salamander$CLOACA)
salamander$GILL <- as.factor(salamander$GILL)

#align ponds names with with other data sets. Just use pond number rather than L and U for lower and upper
salamander <- salamander %>% 
  mutate(POND = case_when(
    POND == 'U06' ~ '56', POND == 'U05' ~ '55', POND == 'U02' ~ '52', POND == 'U2' ~ '52', POND == 'U01' ~ '51', 
    POND == 'U1' ~ '51', POND == 'L52' ~ '52', POND == 'L41' ~ '41', POND == 'L39' ~ '39', POND == 'L38' ~ '38', 
    POND == 'L37' ~ '37', POND == 'L36' ~ '36', POND == 'L35' ~ '35', POND == 'L34' ~ '34', POND == 'L27' ~ '27', 
    POND == 'L18B' ~ '18B', POND == 'L18' ~ '18', POND == 'L17' ~ '17', POND == 'L16' ~ '16', POND == 'L15' ~ '15', 
    POND == 'L14' ~ '14', POND == 'L13' ~ '13', POND == 'L12' ~ '12', POND == 'L012' ~ '12', POND == 'L11' ~ '11', 
    POND == 'L10' ~ '10', POND == 'L09' ~ '9', POND == 'L8' ~ '8', POND == 'L08' ~ '8', POND == 'L07' ~ '7', 
    POND == 'L06' ~ '6', POND == 'L05A' ~ '5A', POND == 'L05' ~ '5', POND == 'L04' ~ '4', POND == 'L03' ~ '3', 
    POND == 'L02' ~ '2', POND == 'L01A' ~ '1A', POND == 'L01' ~ '1', POND == '1L01' ~ '1', POND == '' ~ 'DELETE', 
    POND == 'L02/L03' ~ 'DELETE', POND == 'Terrestrial' ~ 'DELETE', POND == 'Unknown' ~ 'DELETE',
    TRUE ~ POND))

salamander$POND = as.factor(salamander$POND)

salamander = subset(salamander, INDIV_ID != "" & !is.na(POND) & !(MORPH == "Larvae" | MORPH == "Hatchling" | MORPH == "") & 
                    !(SEX == "Immature" | SEX == "2014" | SEX == "") & CONDITION != "dead" & POND != "DELETE")

salamander$SEX = car::recode(salamander$SEX,"'female'='Female'")
salamander$SEX = car::recode(salamander$SEX,"'Female?'='Female'")
salamander$SEX = car::recode(salamander$SEX,"'female?'='Female'")
salamander$SEX = car::recode(salamander$SEX,"'male'='Male'")
salamander$SEX = car::recode(salamander$SEX,"'male?'='Male'")
salamander$SEX = car::recode(salamander$SEX,"'Male?'='Male'")

salamander$MORPH = car::recode(salamander$MORPH,"'meta'='Meta'")
salamander$MORPH = car::recode(salamander$MORPH,"'paedo'='Paedo'")
salamander$MORPH = car::recode(salamander$MORPH,"'Paedo?'='Paedo'")

#need to assess how many salamanders IDs have multiple recorded sexes
sex_check =  salamander %>% 
  group_by(INDIV_ID) %>%
  filter(n_distinct(SEX)>1)

sex_check$SEX = as.character(sex_check$SEX)

mode <- function(x) {
  ux <- unique(na.omit(x))
  tx <- tabulate(match(x, ux))
  if(length(ux) != 1 & sum(max(tx) == tx) > 1) {
    if (is.character(ux)) return(NA_character_) else return(NA_real_)
  }
  max_tx <- tx == max(tx)
  return(ux[max_tx])
}

#for manuscript, need to make sure that there aren't cases when a salamander was 'Female?' 
#for a bunch of initial captures and then confirmed 'Male' only for the most recent capture or two. With
#the SEX level reassignment that I did above, the mode function would select 'Female' for this 
#hypothetical individual

sex_check = sex_check %>% group_by(INDIV_ID) %>% summarise(SEX = mode(SEX))
sex_check$SEX = as.factor(sex_check$SEX)
INDIV_sex_change = subset(sex_check, !is.na(SEX))
INDIV_sex_drop = subset(sex_check, is.na(SEX))
INDIV_sex_drop$SEX = car::recode(INDIV_sex_drop$SEX,"NA='DROP'")

salamander = merge(salamander, INDIV_sex_change, by = "INDIV_ID", all.x = T)
salamander$SEX.y[ is.na(salamander$SEX.y) ] <- salamander$SEX.x[ is.na(salamander$SEX.y) ]
salamander$SEX.x <- NULL
names(salamander)[names(salamander) == 'SEX.y'] <- 'SEX'

salamander = merge(salamander, INDIV_sex_drop, by = "INDIV_ID", all.x = T)
salamander$SEX.y <- factor(salamander$SEX.y, levels = c("DROP", "KEEP"))
salamander$SEX.y[is.na(salamander$SEX.y)] <- "KEEP"
salamander = subset(salamander, SEX.y == "KEEP")
salamander$SEX.y <- NULL
names(salamander)[names(salamander) == 'SEX.x'] <- 'SEX'
salamander = salamander %>% relocate(SEX, .after = COHORT)
#end of sex check#

#extreme ends of WEIGHT category data entry errors? Seem like it. Check in with Howard.
salamander = subset(salamander, SVL < 120 & SVL > 59)
salamander = subset(salamander, WEIGHT < 101 & WEIGHT > 6)

salamander = salamander[ , which(names(salamander) %in% c("INDIV_ID", "DATEOFYEAR", "POND", "MORPH",
                                                           "COHORT", "SEX", "SVL", "WEIGHT"))]
salamander = droplevels(salamander)
summary(salamander)

##subset by morph 
paedomorphs = subset(salamander, MORPH == "Paedo")
salamander = subset(salamander, MORPH == "Meta")

salamander_np = subset(salamander, POND == 6 | POND == 8 | POND == 10 | POND == 11 | POND == 13 | POND == 15 | POND == 51 | POND == 52)
salamander_np = droplevels(salamander_np)


### xMeta & FS '21####

#how many salamander were captured at each pond on each sampling occasion? 
#do salamander leave ponds when FS densities drop?

salamander_np_2021 = subset(salamander_np, DATEOFYEAR > "2021-01-01" & DATEOFYEAR < "2022-01-01")

salamander_np_2021_counts = salamander_np_2021 %>% group_by(DATEOFYEAR, POND) %>% summarise(Meta_Count = n())


meta_pond_counts =
  ggplot(salamander_np_2021_counts, aes(x=DATEOFYEAR, y=Meta_Count, fill = POND)) +
  geom_line(aes(group = POND, color = POND), lwd = 2) +
  geom_point(aes(color = POND), size = 10, alpha = 0.8) +
  scale_x_date(limits = as.Date(c('2021-06-14','2021-07-27'))) +
  ylim(0,30) +
  theme_bw(33) +
  theme(axis.title.x = element_blank(),
        panel.grid.minor = element_blank())

meta_pond_counts 

#for clarity, make individuals plots for each pond with two axes (one for FS densities 
#and one for Meta counts).

fs_pond = read.csv(file = "Data/fs_pond.csv", stringsAsFactors = FALSE)
fs_pond$X <- NULL
fs_pond$Date = ymd(fs_pond$Date)
fs_pond$Pond = as.factor(fs_pond$Pond)

#pond 8

meta_pond8 = subset(salamander_np_2021_counts, POND == "8")
colnames(meta_pond8) = c("Date", "Pond", "Value")
meta_pond8$Dataset = "Metamorph Abundance"
fs_pond8 = subset(fs_pond, Pond == "8" & Date < "2022-01-01")
fs_pond8 = fs_pond8[,c("Date", "Pond", "FS_Density_Pond")]
colnames(fs_pond8) = c("Date", "Pond", "Value")
fs_pond8$Dataset = "Fairy Shrimp Density"
fs_pond8$Date = ymd(fs_pond8$Date)
fs_pond8$Pond = as.factor(fs_pond8$Pond)
fs_pond8$Value[fs_pond8$Value == 1] <- 0 
plotdata8 = rbind(meta_pond8, fs_pond8)

mycolors <- c("Metamorph Abundance"="blue", "Fairy Shrimp Density"="red")

pond8_meta_fs = 
  ggplot(plotdata8, aes(x=Date, y=Value, group=Dataset, color=Dataset)) +
  geom_path(lwd = 5) +
  geom_point(size = 20) +
  scale_y_continuous(limits = c(0,40)) +
  scale_color_manual(name="Dataset", values = mycolors) +
  scale_x_date(limits = as.Date(c('2021-06-14','2021-08-07'))) +
  theme(legend.position = "none",
    axis.title = element_blank(),
    axis.text.y = element_text(margin=margin(r=10)),
    axis.text = element_text(size = 50),
    axis.ticks = element_line(size = 2),
    axis.ticks.length = unit(.25, "cm"),
    axis.text.x=element_text(margin=margin(t=10)))

pond8_meta_fs

#ggsave(pond8_meta_fs, filename = "Outputs/pond8_meta_fs.png",  width = 18, height = 12, dpi = 100)


### different approach ###

##pond 8

meta_pond8 = subset(salamander_np_2021_counts, POND == "8")
colnames(meta_pond8) = c("Date", "Pond", "Count")

fs_pond8 = subset(fs_pond, Pond == "8" & Date < "2022-01-01")
fs_pond8 = fs_pond8[,c("Date", "Pond", "log_FS_Biomass_Pond_m3", "SE_Biomass")]
fs_pond8$Date = ymd(fs_pond8$Date)
fs_pond8$Pond = as.factor(fs_pond8$Pond)

#fs
pond8_fs = 
  ggplot(data = fs_pond8, aes(x=Date, y=log_FS_Biomass_Pond_m3)) +
  geom_point(data = fs_pond8, aes(x=Date, y=log_FS_Biomass_Pond_m3), size = 15) +
  geom_errorbar(data = fs_pond8, aes(ymin=log_FS_Biomass_Pond_m3-SE_Biomass*1.96, ymax=log_FS_Biomass_Pond_m3+SE_Biomass*1.96), color = "black", width = 0, lwd = 2) +
  geom_path(data = fs_pond8, aes(x=Date, y=log_FS_Biomass_Pond_m3), lwd = 3) +
  scale_y_continuous(limits = c(log(1), log(100)), breaks=c(log(1),log(10),log(100)), labels=c("0", "10", "100")) +
  scale_x_date(limits = as.Date(c('2021-06-14','2021-08-07'))) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_text(margin=margin(r=10)),
        axis.text = element_text(size = 50),
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x=element_text(margin=margin(t=10)))

pond8_fs

#meta
pond8_meta = 
  ggplot(data = meta_pond8, aes(x=Date, y=Count)) +
  geom_point(data = meta_pond8, aes(x=Date, y=Count), size = 15) +
  geom_path(data = meta_pond8, aes(x=Date, y=Count), lwd = 3) +
  scale_y_continuous(limits = c(0, 30)) +
  scale_x_date(limits = as.Date(c('2021-06-14','2021-08-07'))) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_text(margin=margin(r=10)),
        axis.text = element_text(size = 50),
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x=element_text(margin=margin(t=10)))

pond8_meta

#body condition
#choose dates with enough metas captured to equal a meaningful average
salamander_np_2021_pond8 = subset(salamander_np_2021, POND == "8" &
                                    (DATEOFYEAR == "2021-06-18" | DATEOFYEAR == "2021-06-24" | 
                                       DATEOFYEAR == "2021-06-30" | DATEOFYEAR == "2021-07-05" |
                                       DATEOFYEAR == "2021-07-12" | DATEOFYEAR == "2021-07-27") &
                                    !is.na(SVL))

salamander_np_2021_pond8$BODY_CONDITION = salamander_np_2021_pond8$WEIGHT / salamander_np_2021_pond8$SVL

salamander_np_2021_pond8 = salamander_np_2021_pond8 %>% group_by(DATEOFYEAR) %>% mutate(Mean_Condition = mean(BODY_CONDITION))
salamander_np_2021_pond8 = salamander_np_2021_pond8 %>% group_by(DATEOFYEAR) %>% mutate(SE_Condition = sd(BODY_CONDITION)/sqrt(n()))
salamander_condition_pond8 = unique(salamander_np_2021_pond8[,c("DATEOFYEAR", "POND", "Mean_Condition", "SE_Condition")])

pond8_meta_condition = 
  ggplot(salamander_condition_pond8, aes(x=DATEOFYEAR, y=Mean_Condition)) +
  geom_errorbar(aes(ymin=Mean_Condition-SE_Condition*1.96, ymax=Mean_Condition+SE_Condition*1.96), color = "black", width = 0, lwd = 2) +
  geom_path(lwd = 2) +
  geom_point(size = 10) +
  scale_y_continuous(limits = c(0.1,0.4)) +
  scale_x_date(limits = as.Date(c('2021-06-14','2021-08-07'))) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_text(margin=margin(r=10)),
        axis.text = element_text(size = 50),
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x=element_text(margin=margin(t=10)))

pond8_meta_condition

### end different approach ###

#pond 10

meta_pond10 = subset(salamander_np_2021_counts, POND == "10")
colnames(meta_pond10) = c("Date", "Pond", "Value")
meta_pond10$Dataset = "Metamorph Abundance"
fs_pond10 = subset(fs_pond, Pond == "10" & Date < "2022-01-01")
fs_pond10 = fs_pond10[,c("Date", "Pond", "FS_Density_Pond")]
colnames(fs_pond10) = c("Date", "Pond", "Value")
fs_pond10$Dataset = "Fairy Shrimp Density"
fs_pond10$Date = ymd(fs_pond10$Date)
fs_pond10$Pond = as.factor(fs_pond10$Pond)
fs_pond10$Value[fs_pond10$Value == 1] <- 0 
plotdata10 = rbind(meta_pond10, fs_pond10)

pond10_meta_fs = 
  ggplot(plotdata10, aes(x=Date, y=Value, group=Dataset, color=Dataset)) +
  geom_path(lwd = 5) +
  geom_point(size = 20) +
  scale_y_continuous(limits = c(0,25)) +
  scale_color_manual(name="Dataset", values = mycolors) +
  scale_x_date(limits = as.Date(c('2021-06-14','2021-08-07'))) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_text(margin=margin(r=10)),
        axis.text = element_text(size = 50),
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x=element_text(margin=margin(t=10)))

pond10_meta_fs

#ggsave(pond10_meta_fs, filename = "Outputs/pond10_meta_fs.png",  width = 18, height = 12, dpi = 100)


#pond 11

meta_pond11 = subset(salamander_np_2021_counts, POND == "11")
colnames(meta_pond11) = c("Date", "Pond", "Value")
meta_pond11$Dataset = "Metamorph Abundance"
fs_pond11 = subset(fs_pond, Pond == "11" & Date < "2022-01-01")
fs_pond11 = fs_pond11[,c("Date", "Pond", "FS_Density_Pond")]
colnames(fs_pond11) = c("Date", "Pond", "Value")
fs_pond11$Dataset = "Fairy Shrimp Density"
fs_pond11$Date = ymd(fs_pond11$Date)
fs_pond11$Pond = as.factor(fs_pond11$Pond)
fs_pond11$Value[fs_pond11$Value == 1] <- 0 
fs_pond11_sec_axis <- within(fs_pond11, { Value = Value/10 })
plotdata11 = rbind(meta_pond11, fs_pond11_sec_axis)

pond11_meta_fs = 
  ggplot(plotdata11, aes(x=Date, y=Value, group=Dataset, color=Dataset)) +
  geom_path(lwd = 5) +
  geom_point(size = 20) +
  scale_y_continuous(limits = c(0, 15), name="Metamorph Abundance", sec.axis = sec_axis(~ 10*., name="Fairy Shrimp Density")) +
  scale_color_manual(name="Dataset", values = mycolors) +
  scale_x_date(limits = as.Date(c('2021-06-14','2021-08-07'))) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_text(color = mycolors["Metamorph Abundance"], margin=margin(r=10)),
        axis.text.y.right = element_text(color = mycolors["Fairy Shrimp Density"], margin=margin(l=10)),
        axis.text = element_text(size = 50),
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x=element_text(margin=margin(t=10)))

pond11_meta_fs

#ggsave(pond11_meta_fs, filename = "Outputs/pond11_meta_fs.png",  width = 18, height = 12, dpi = 100)


#pond 6

meta_pond6 = subset(salamander_np_2021_counts, POND == "6")
colnames(meta_pond6) = c("Date", "Pond", "Value")
meta_pond6$Dataset = "Metamorph Abundance"
fs_pond6 = subset(fs_pond, Pond == "6" & Date < "2022-01-01")
fs_pond6 = fs_pond6[,c("Date", "Pond", "FS_Density_Pond")]
colnames(fs_pond6) = c("Date", "Pond", "Value")
fs_pond6$Dataset = "Fairy Shrimp Density"
fs_pond6$Date = ymd(fs_pond6$Date)
fs_pond6$Pond = as.factor(fs_pond6$Pond)
fs_pond6$Value[fs_pond6$Value == 1] <- 0 
plotdata6 = rbind(meta_pond6, fs_pond6)

pond6_meta_fs = 
  ggplot(plotdata6, aes(x=Date, y=Value, group=Dataset, color=Dataset)) +
  geom_path(lwd = 5) +
  geom_point(size = 20) +
  scale_y_continuous(limits = c(0,15)) +
  scale_color_manual(name="Dataset", values = mycolors) +
  scale_x_date(limits = as.Date(c('2021-06-14','2021-08-07'))) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_text(margin=margin(r=10)),
        axis.text = element_text(size = 50),
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x=element_text(margin=margin(t=10)))

pond6_meta_fs

#ggsave(pond6_meta_fs, filename = "Outputs/pond6_meta_fs.png",  width = 18, height = 12, dpi = 100)


### xMeta & FS '22####

salamander_np_2022 = subset(salamander_np, DATEOFYEAR > "2022-01-01")

salamander_np_2022_counts = salamander_np_2022 %>% group_by(DATEOFYEAR, POND) %>% summarise(Meta_Count = n())

#pond 8

meta_pond8 = subset(salamander_np_2022_counts, POND == "8")
colnames(meta_pond8) = c("Date", "Pond", "Value")
meta_pond8$Dataset = "Metamorph Abundance"
fs_pond8 = subset(fs_pond, Pond == "8" & Date > "2022-01-01")
fs_pond8 = fs_pond8[,c("Date", "Pond", "FS_Density_Pond")]
colnames(fs_pond8) = c("Date", "Pond", "Value")
fs_pond8$Dataset = "Fairy Shrimp Density"
fs_pond8$Date = ymd(fs_pond8$Date)
fs_pond8$Pond = as.factor(fs_pond8$Pond)
fs_pond8$Value[fs_pond8$Value == 1] <- 0 
plotdata8 = rbind(meta_pond8, fs_pond8)

mycolors <- c("Metamorph Abundance"="blue", "Fairy Shrimp Density"="red")

pond8_meta_fs_2022 = 
  ggplot(plotdata8, aes(x=Date, y=Value, group=Dataset, color=Dataset)) +
  geom_path(lwd = 5) +
  geom_point(size = 20) +
  scale_y_continuous(limits = c(0,40)) +
  scale_color_manual(name="Dataset", values = mycolors) +
  scale_x_date(limits = as.Date(c('2022-06-14','2022-08-07'))) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_text(margin=margin(r=10)),
        axis.text = element_text(size = 50),
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x=element_text(margin=margin(t=10)))

pond8_meta_fs_2022

#ggsave(pond8_meta_fs_2022, filename = "Outputs/pond8_meta_fs_2022.png",  width = 18, height = 12, dpi = 100)


#pond 6

meta_pond6 = subset(salamander_np_2022_counts, POND == "6")
colnames(meta_pond6) = c("Date", "Pond", "Value")
meta_pond6$Dataset = "Metamorph Abundance"
fs_pond6 = subset(fs_pond, Pond == "6" & Date > "2022-01-01")
fs_pond6 = fs_pond6[,c("Date", "Pond", "FS_Density_Pond")]
colnames(fs_pond6) = c("Date", "Pond", "Value")
fs_pond6$Dataset = "Fairy Shrimp Density"
fs_pond6$Date = ymd(fs_pond6$Date)
fs_pond6$Pond = as.factor(fs_pond6$Pond)
fs_pond6$Value[fs_pond6$Value == 1] <- 0 
plotdata6 = rbind(meta_pond6, fs_pond6)

mycolors <- c("Metamorph Abundance"="blue", "Fairy Shrimp Density"="red")

pond6_meta_fs_2022 = 
  ggplot(plotdata6, aes(x=Date, y=Value, group=Dataset, color=Dataset)) +
  geom_path(lwd = 5) +
  geom_point(size = 20) +
  #scale_y_continuous(limits = c(0,40)) +
  scale_color_manual(name="Dataset", values = mycolors) +
  scale_x_date(limits = as.Date(c('2022-06-14','2022-08-07'))) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_text(margin=margin(r=10)),
        axis.text = element_text(size = 50),
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x=element_text(margin=margin(t=10)))

pond6_meta_fs_2022

#ggsave(pond6_meta_fs_2022, filename = "Outputs/pond6_meta_fs_2022.png",  width = 18, height = 12, dpi = 100)


#pond 10

meta_pond10 = subset(salamander_np_2022_counts, POND == "10")
colnames(meta_pond10) = c("Date", "Pond", "Value")
meta_pond10$Dataset = "Metamorph Abundance"
fs_pond10 = subset(fs_pond, Pond == "10" & Date > "2022-01-01")
fs_pond10 = fs_pond10[,c("Date", "Pond", "FS_Density_Pond")]
colnames(fs_pond10) = c("Date", "Pond", "Value")
fs_pond10$Dataset = "Fairy Shrimp Density"
fs_pond10$Date = ymd(fs_pond10$Date)
fs_pond10$Pond = as.factor(fs_pond10$Pond)
fs_pond10$Value[fs_pond10$Value == 1] <- 0 
plotdata10 = rbind(meta_pond10, fs_pond10)

mycolors <- c("Metamorph Abundance"="blue", "Fairy Shrimp Density"="red")

pond10_meta_fs_2022 = 
  ggplot(plotdata10, aes(x=Date, y=Value, group=Dataset, color=Dataset)) +
  geom_path(lwd = 5) +
  geom_point(size = 20) +
  #scale_y_continuous(limits = c(0,40)) +
  scale_color_manual(name="Dataset", values = mycolors) +
  scale_x_date(limits = as.Date(c('2022-06-14','2022-08-07'))) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_text(margin=margin(r=10)),
        axis.text = element_text(size = 50),
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x=element_text(margin=margin(t=10)))

pond10_meta_fs_2022

#ggsave(pond10_meta_fs_2022, filename = "Outputs/pond10_meta_fs_2022.png",  width = 18, height = 12, dpi = 100)


#pond 11

meta_pond11 = subset(salamander_np_2022_counts, POND == "11")
colnames(meta_pond11) = c("Date", "Pond", "Value")
meta_pond11$Dataset = "Metamorph Abundance"
fs_pond11 = subset(fs_pond, Pond == "11" & Date > "2022-01-01")
fs_pond11 = fs_pond11[,c("Date", "Pond", "FS_Density_Pond")]
colnames(fs_pond11) = c("Date", "Pond", "Value")
fs_pond11$Dataset = "Fairy Shrimp Density"
fs_pond11$Date = ymd(fs_pond11$Date)
fs_pond11$Pond = as.factor(fs_pond11$Pond)
fs_pond11$Value[fs_pond11$Value == 1] <- 0 
plotdata11 = rbind(meta_pond11, fs_pond11)

mycolors <- c("Metamorph Abundance"="blue", "Fairy Shrimp Density"="red")

pond11_meta_fs_2022 = 
  ggplot(plotdata11, aes(x=Date, y=Value, group=Dataset, color=Dataset)) +
  geom_path(lwd = 5) +
  geom_point(size = 20) +
  #scale_y_continuous(limits = c(0,40)) +
  scale_color_manual(name="Dataset", values = mycolors) +
  scale_x_date(limits = as.Date(c('2022-06-14','2022-08-07'))) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_text(margin=margin(r=10)),
        axis.text = element_text(size = 50),
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x=element_text(margin=margin(t=10)))

pond11_meta_fs_2022

#ggsave(pond11_meta_fs_2022, filename = "Outputs/pond11_meta_fs_2022.png",  width = 18, height = 12, dpi = 100)


#pond 52
meta_pond52 = subset(salamander_np_2022_counts, POND == "52")
colnames(meta_pond52) = c("Date", "Pond", "Value")
meta_pond52$Dataset = "Metamorph Abundance"
fs_pond52 = subset(fs_pond, Pond == "52" & Date > "2022-01-01")
fs_pond52 = fs_pond52[,c("Date", "Pond", "FS_Density_Pond")]
colnames(fs_pond52) = c("Date", "Pond", "Value")
fs_pond52$Dataset = "Fairy Shrimp Density"
fs_pond52$Date = ymd(fs_pond52$Date)
fs_pond52$Pond = as.factor(fs_pond52$Pond)
fs_pond52$Value[fs_pond52$Value == 1] <- 0 
plotdata52 = rbind(meta_pond52, fs_pond52)

mycolors <- c("Metamorph Abundance"="blue", "Fairy Shrimp Density"="red")

pond52_meta_fs_2022 = 
  ggplot(plotdata52, aes(x=Date, y=Value, group=Dataset, color=Dataset)) +
  geom_path(lwd = 5) +
  geom_point(size = 20) +
  #scale_y_continuous(limits = c(0,40)) +
  scale_color_manual(name="Dataset", values = mycolors) +
  scale_x_date(limits = as.Date(c('2022-06-14','2022-08-07'))) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_text(margin=margin(r=10)),
        axis.text = element_text(size = 50),
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x=element_text(margin=margin(t=10)))

pond52_meta_fs_2022

#ggsave(pond52_meta_fs_2022, filename = "Outputs/pond52_meta_fs_2022.png",  width = 18, height = 12, dpi = 100)


#pond 15
meta_pond15 = subset(salamander_np_2022_counts, POND == "15")
colnames(meta_pond15) = c("Date", "Pond", "Value")
meta_pond15$Dataset = "Metamorph Abundance"
fs_pond15 = subset(fs_pond, Pond == "15" & Date > "2022-01-01")
fs_pond15 = fs_pond15[,c("Date", "Pond", "FS_Density_Pond")]
colnames(fs_pond15) = c("Date", "Pond", "Value")
fs_pond15$Dataset = "Fairy Shrimp Density"
fs_pond15$Date = ymd(fs_pond15$Date)
fs_pond15$Pond = as.factor(fs_pond15$Pond)
fs_pond15$Value[fs_pond15$Value == 1] <- 0 
plotdata15 = rbind(meta_pond15, fs_pond15)

mycolors <- c("Metamorph Abundance"="blue", "Fairy Shrimp Density"="red")

pond15_meta_fs_2022 = 
  ggplot(plotdata15, aes(x=Date, y=Value, group=Dataset, color=Dataset)) +
  geom_path(lwd = 5) +
  geom_point(size = 20) +
  #scale_y_continuous(limits = c(0,40)) +
  scale_color_manual(name="Dataset", values = mycolors) +
  scale_x_date(limits = as.Date(c('2022-06-14','2022-08-07'))) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_text(margin=margin(r=10)),
        axis.text = element_text(size = 50),
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x=element_text(margin=margin(t=10)))

pond15_meta_fs_2022

#ggsave(pond15_meta_fs_2022, filename = "Outputs/pond15_meta_fs_2022.png",  width = 18, height = 12, dpi = 100)



#### x21/22 Pond Mvmt #####
#look at meta movement between non-permanent ponds. Do they always move to a pond with a higher fs density?

#ideas:
#(1) could do map with arrows that has green arrows for moving to pond with higher FS density and red
#arrows for moving to pond with lower FS density, and then scale arrow size by how many salamanders made 
#that category of movement 

#(2) Logistic regression. Seen in any pond again (Y/N)? ~ % of peak FS density when meta left pond


#(1)
#Figure out how many qualifying meta movements are in the dataset. Use 2021 and 2022 data 
#combined or analyze each year separately?

#2021
salamander_2021 = subset(salamander, DATEOFYEAR > "2021-01-01" & DATEOFYEAR < "2022-01-01" & !is.na(POND))

#select IDs that have at least two unique values for POND
salamander_2021_mvmt = 
  salamander_2021 %>% 
  group_by(INDIV_ID) %>%
  filter(n_distinct(POND)>1)


#add in pond hydroperiod column
salamander_2021_mvmt$Hydroperiod <- ifelse(salamander_2021_mvmt$POND == "1A" |
                                           salamander_2021_mvmt$POND == "10" | 
                                           salamander_2021_mvmt$POND == "11" |
                                           salamander_2021_mvmt$POND == "13" |
                                           salamander_2021_mvmt$POND == "15" |
                                           salamander_2021_mvmt$POND == "38" |
                                           salamander_2021_mvmt$POND == "51" |
                                           salamander_2021_mvmt$POND == "52" |
                                           salamander_2021_mvmt$POND == "6" |
                                           salamander_2021_mvmt$POND == "8", "Non-permanent", "Permanent")

salamander_2021_mvmt = salamander_2021_mvmt %>% relocate(Hydroperiod, .after = POND)

#there are 30 unique metas with qualifying movement (between 14 ponds)
levels(droplevels(salamander_2021_mvmt$INDIV_ID))
levels(droplevels(salamander_2021_mvmt$POND))

#make schematic for ponds 1, 12, 5, 9 (permanent) and 1A, 10, 11, 13, 15, 51, 52, 6, 8 (non-permanent)
#one meta: 12 -> 11 -> 10 (all others [29/30] only moved between two ponds)


##2022
salamander_2022 = subset(salamander, DATEOFYEAR > "2022-01-01" & !is.na(POND))

#select IDs that have at least two unique values for POND
salamander_2022_mvmt = 
  salamander_2022 %>% 
  group_by(INDIV_ID) %>%
  filter(n_distinct(POND)>1)

#add in pond hydroperiod column
salamander_2022_mvmt$Hydroperiod <- ifelse(salamander_2022_mvmt$POND == "10" | 
                                             salamander_2022_mvmt$POND == "11" |
                                             salamander_2022_mvmt$POND == "1A" |
                                             salamander_2022_mvmt$POND == "15" |
                                             salamander_2022_mvmt$POND == "6" |
                                             salamander_2022_mvmt$POND == "8", "Non-permanent", "Permanent")

salamander_2022_mvmt = salamander_2022_mvmt %>% relocate(Hydroperiod, .after = POND)

#there are 12 unique metas with qualifying movement (between 9 ponds)
levels(droplevels(salamander_2022_mvmt$INDIV_ID))
levels(droplevels(salamander_2022_mvmt$POND))

#one meta: 12 -> 10 -> 11 (all others [11/12] only moved between two ponds)

##in both 2021 and 2022, there were 6 movements between non-permanent ponds (with majority btwn. 10 & 11)


#2021 - merge in closest fs data for each date-pond so can see if metas are moving from low to high fs density
salamander_2021_mvmt_np = subset(salamander_2021_mvmt, Hydroperiod == "Non-permanent" & POND != "1A" & POND != "38")
salamander_2021_mvmt_np = as.data.frame(salamander_2021_mvmt_np) #ungroups df
salamander_2021_mvmt_np[,c("POND")] <- factor(salamander_2021_mvmt_np[,c("POND")], levels=c("6", "8", "10", "11", "13", "15", "51", "52", "55"))

temp <- list()

for (i in 1:nrow(salamander_2021_mvmt_np)) {       
  fs_pond_subset <- subset(fs_pond, Pond==(salamander_2021_mvmt_np[i, 6]))
  temp[[i]] <- ajfhelpR::date_near(fs_pond_subset$Date, salamander_2021_mvmt_np[i, 3], sidepref = 'l')
}

temp

dates = temp %>% purrr::reduce(c) #make list into vector
rm(temp)

salamander_2021_mvmt_np$Date_Sample_FS <- dates
salamander_2021_mvmt_np$Date_Sample_FS = ymd(salamander_2021_mvmt_np$Date_Sample_FS)
abs(salamander_2021_mvmt_np$Date_Sample_FS - salamander_2021_mvmt_np$DATEOFYEAR) #time difference btwn GF and FS sampling

#merge in fs biomass data
names(fs_pond)[names(fs_pond) == 'Date'] <- 'Date_Sample_FS'
names(salamander_2021_mvmt_np)[names(salamander_2021_mvmt_np) == 'POND'] <- 'Pond'

salamander_2021_mvmt_np = merge(salamander_2021_mvmt_np, fs_pond[,c("Date_Sample_FS", "Pond", "log_FS_Biomass_Pond_m3",
                                                              "SE_Biomass", "lower_biomass", "upper_biomass")], 
                                                               by = c("Date_Sample_FS", "Pond"), all.x = T)

salamander_2021_mvmt_np = salamander_2021_mvmt_np %>% relocate(Date_Sample_FS, .after = WEIGHT)
salamander_2021_mvmt_np = salamander_2021_mvmt_np %>% relocate(Pond, .after = MORPH)
salamander_2021_mvmt_np = salamander_2021_mvmt_np[,c("INDIV_ID", "DATEOFYEAR", "Pond", "COHORT", "SEX", "SVL", "WEIGHT", "log_FS_Biomass_Pond_m3", "SE_Biomass")]

salamander_2021_mvmt_np = 
  salamander_2021_mvmt_np %>% 
  group_by(INDIV_ID) %>%
  filter(n_distinct(Pond)>1)

#do salamanders typically move to ponds with higher FS densities? Movement btwn 10 and 11 may not
#follow pattern bc they're so close?


#Q: how many metas are not showing up in the inter-pond movement data? I.e., how many were captured in just
#one permanent/non-permanent pond during the season???

salamander_2021_no_mvmt = 
  salamander_2021 %>% 
  group_by(INDIV_ID) %>%
  filter(n_distinct(POND)<2)

levels(droplevels(salamander_2021$INDIV_ID)) #167 metas captured in 2021
levels(droplevels(salamander_2021_no_mvmt$INDIV_ID)) #137 metas (82%) made no inter-pond movements
levels(droplevels(salamander_2021_mvmt$INDIV_ID)) #30 metas (18%) moved btween ponds


salamander_2021_no_mvmt$Hydroperiod <- ifelse(salamander_2021_no_mvmt$POND == "1" | 
                                                salamander_2021_no_mvmt$POND == "12" |
                                                salamander_2021_no_mvmt$POND == "5" |
                                                salamander_2021_no_mvmt$POND == "9" |
                                                salamander_2021_no_mvmt$POND == "18", "Permanent", "Non-permanent")

salamander_2021_no_mvmt = salamander_2021_no_mvmt %>% relocate(Hydroperiod, .after = POND)

count = salamander_2021_no_mvmt[,c("INDIV_ID", "POND", "Hydroperiod")]
count = unique(count)
count = droplevels(count)
summary(as.factor(count$Hydroperiod)) 
summary(count$POND) 

#76/137 metas (55%) that made no inter-pond movements were captured in non-permanent ponds
#61/137 metas (45%) that made no inter-pond movements were captured in permanent ponds


#2022 - merge in closest fs data for each date-pond so can see if metas are moving from low to high fs density
salamander_2022_mvmt_np = subset(salamander_2022_mvmt, Hydroperiod == "Non-permanent" & POND != "1A")
salamander_2022_mvmt_np = as.data.frame(salamander_2022_mvmt_np) #ungroups df
salamander_2022_mvmt_np[,c("POND")] <- factor(salamander_2022_mvmt_np[,c("POND")], levels=c("6", "8", "10", "11", "13", "15", "51", "52", "55"))

temp <- list()

for (i in 1:nrow(salamander_2022_mvmt_np)) {       
  fs_pond_subset <- subset(fs_pond, Pond==(salamander_2022_mvmt_np[i, 6]))
  temp[[i]] <- ajfhelpR::date_near(fs_pond_subset$Date, salamander_2022_mvmt_np[i, 3], sidepref = 'l')
}

temp

dates = temp %>% purrr::reduce(c) #make list into vector
rm(temp)

salamander_2022_mvmt_np$Date_Sample_FS <- dates
salamander_2022_mvmt_np$Date_Sample_FS = ymd(salamander_2022_mvmt_np$Date_Sample_FS)
abs(salamander_2022_mvmt_np$Date_Sample_FS - salamander_2022_mvmt_np$DATEOFYEAR) #time difference btwn GF and FS sampling

#merge in fs biomass data
names(fs_pond)[names(fs_pond) == 'Date'] <- 'Date_Sample_FS'
names(salamander_2022_mvmt_np)[names(salamander_2022_mvmt_np) == 'POND'] <- 'Pond'

salamander_2022_mvmt_np = merge(salamander_2022_mvmt_np, fs_pond[,c("Date_Sample_FS", "Pond", "log_FS_Biomass_Pond_m3",
                                                                    "SE_Biomass", "lower_biomass", "upper_biomass")], 
                                by = c("Date_Sample_FS", "Pond"), all.x = T)

salamander_2022_mvmt_np = salamander_2022_mvmt_np %>% relocate(Date_Sample_FS, .after = WEIGHT)
salamander_2022_mvmt_np = salamander_2022_mvmt_np %>% relocate(Pond, .after = MORPH)
salamander_2022_mvmt_np = salamander_2022_mvmt_np[,c("INDIV_ID", "DATEOFYEAR", "Pond", "COHORT", "SEX", "SVL", "WEIGHT", "log_FS_Biomass_Pond_m3", "SE_Biomass")]

salamander_2022_mvmt_np = 
  salamander_2022_mvmt_np %>% 
  group_by(INDIV_ID) %>%
  filter(n_distinct(Pond)>1)


#Q: how many metas are not showing up in the inter-pond movement data? I.e., how many were captured in just
#one permanent/non-permanent pond during the season???

salamander_2022_no_mvmt = 
  salamander_2022 %>% 
  group_by(INDIV_ID) %>%
  filter(n_distinct(POND)<2)

levels(droplevels(salamander_2022$INDIV_ID)) #147 metas captured in 2022
levels(droplevels(salamander_2022_no_mvmt$INDIV_ID)) #135 metas (92%) did not move between ponds
levels(droplevels(salamander_2022_mvmt$INDIV_ID)) #12 metas (8%) moved between ponds

salamander_2022_no_mvmt$Hydroperiod <- ifelse(salamander_2022_no_mvmt$POND == "1" | 
                                                salamander_2022_no_mvmt$POND == "12" |
                                                salamander_2022_no_mvmt$POND == "5" |
                                                salamander_2022_no_mvmt$POND == "9" |
                                                salamander_2022_no_mvmt$POND == "14" |
                                                salamander_2022_no_mvmt$POND == "3" |
                                                salamander_2022_no_mvmt$POND == "53" |
                                                salamander_2022_no_mvmt$POND == "54" |
                                                salamander_2022_no_mvmt$POND == "18", "Permanent", "Non-permanent")


salamander_2022_no_mvmt = salamander_2022_no_mvmt %>% relocate(Hydroperiod, .after = POND)

count = salamander_2022_no_mvmt[,c("INDIV_ID", "POND", "Hydroperiod")]
count = unique(count)
summary(as.factor(count$Hydroperiod)) 
#96/135 metas (71%) that made no inter-pond movements were captured in non-permanent ponds
#39/135 metas (29%) that made no inter-pond movements were captured in permanent ponds


#(2) Logistic regression. Seen in any pond again (Y/N)? ~ % of peak FS density when meta left pond
#not enough data to answer this question? Only have 12 data points (metas leaving
#a non-permanent pond)


##### Decision Tree #####
#Mvmt Patterns Long-Term
#can divide annual meta movements into four broad categories? Mine long-term dataset

#(1) don't return to pond system ---> avg, median, max years btwn. captures?
salamander$Year = year(salamander$DATEOFYEAR)
salamander = salamander %>% relocate(Year, .before = DATEOFYEAR)
salamander$INDIV_ID_Year = interaction(salamander$INDIV_ID, salamander$Year, sep="_")
salamander = salamander %>% relocate(INDIV_ID_Year, .after = Year)
salamander$INDIV_ID_Year = as.factor(salamander$INDIV_ID_Year)
salamander = droplevels(salamander)

nlevels(salamander$INDIV_ID) #1827 unique metas
nrow(salamander) #6915 total captures
range(salamander$Year) #1990-2022
nlevels(salamander$INDIV_ID_Year) #4145 meta-years

salamander_recap = 
  salamander %>% 
  group_by(INDIV_ID) %>%
  filter(n_distinct(Year)>1)

salamander_recap_yrs = salamander_recap[,c("INDIV_ID", "Year", "POND")]

salamander_recap_yrs_return = 
  salamander_recap_yrs %>% 
  group_by(INDIV_ID) %>% 
  mutate(Years_To_Return = Year - lead(Year))

salamander_recap_yrs_return = subset(salamander_recap_yrs_return, Years_To_Return > 0)

mean(salamander_recap_yrs_return$Years_To_Return, na.rm = TRUE) #2.1 years
median(salamander_recap_yrs_return$Years_To_Return, na.rm = TRUE) #1 year
max(salamander_recap_yrs_return$Years_To_Return, na.rm = TRUE) #18 years
hist(salamander_recap_yrs_return$Years_To_Return,
     main = "",
     xlab = "Years Between Metamorph Captures",
     ylab = "Frequency",
     xlim = range(0, 18),
     ylim = c(0, 1500))

salamander_recap_yrs_return = droplevels(salamander_recap_yrs_return)
nlevels(salamander_recap_yrs_return$INDIV_ID) #729 metas captured in at least two years
nrow(salamander_recap_yrs_return) #2200 capture intervals


#(2) return to one pond in a given year
salamander_mvmt = salamander[,c("INDIV_ID", "DATEOFYEAR", "Year", "INDIV_ID_Year", "POND")]
salamander_mvmt = droplevels(salamander_mvmt)
levels(salamander_mvmt$POND)

salamander_mvmt$Hydroperiod <- ifelse(salamander_mvmt$POND == "1" | 
                                        salamander_mvmt$POND == "12" |
                                        salamander_mvmt$POND == "5" |
                                        salamander_mvmt$POND == "9" |
                                        salamander_mvmt$POND == "14" |
                                        salamander_mvmt$POND == "3" |
                                        salamander_mvmt$POND == "18", "Permanent", "Non-permanent")

salamander_mvmt$Hydroperiod = as.factor(salamander_mvmt$Hydroperiod)
salamander_mvmt$Year = as.factor(salamander_mvmt$Year)

salamander_mvmt_one_pond = 
  salamander_mvmt %>% 
  group_by(INDIV_ID, Year) %>%
  filter(n_distinct(POND)==1)

summary(salamander_mvmt_one_pond) 
salamander_mvmt_one_pond = droplevels(salamander_mvmt_one_pond)
nlevels(salamander_mvmt_one_pond$INDIV_ID_Year) #3747 meta-years

temp = salamander_mvmt_one_pond[,c("INDIV_ID_Year", "Hydroperiod")]
temp = unique(temp)
summary(temp$Hydroperiod)
#2215/3773 (59%) of one pond meta-years were in non-permanent ponds
#1532/3773 (41%) of one pond meta-years were in permanent ponds

##what is the average percentage of an individuals movements over its lifetime to one pond (vs. visiting multiple)
#analyze only for metas that have been captured in at least two years

num_ponds_visited = 
  salamander_recap_yrs %>% 
  dplyr::group_by(INDIV_ID, Year) %>%
  summarize(Ponds_Visited = n_distinct(POND))

num_ponds_visited$Ponds_Visited = as.factor(num_ponds_visited$Ponds_Visited)

num_ponds_visited = 
  num_ponds_visited %>% 
  group_by(INDIV_ID) %>%
  count(INDIV_ID, Ponds_Visited)

num_ponds_visited = 
  num_ponds_visited %>% 
  group_by(INDIV_ID) %>%
  mutate(Percent = n/sum(n))

sum(num_ponds_visited$n) #3104 annual movements

num_ponds_visited = droplevels(num_ponds_visited)
nlevels(num_ponds_visited$INDIV_ID) #786 unique metas

num_ponds_visited = subset(num_ponds_visited, Ponds_Visited == "1")
sum(num_ponds_visited$n)

mean(num_ponds_visited$Percent) 


#(3) return to multiple ponds in a given year

salamander_mvmt_mult_pond = 
  salamander_mvmt %>% 
  group_by(INDIV_ID, Year) %>%
  filter(n_distinct(POND)>1)

salamander_mvmt_mult_pond = unique(salamander_mvmt_mult_pond)
salamander_mvmt_mult_pond = droplevels(salamander_mvmt_mult_pond)
nlevels(salamander_mvmt_mult_pond$INDIV_ID_Year) #398 meta-years

mult_2021 = subset(salamander_mvmt_mult_pond, Year == "2021")
mult_2021 = mult_2021[,c("INDIV_ID", "Year", "POND", "Hydroperiod")]
mult_2021 = unique(mult_2021)
mult_2022 = subset(salamander_mvmt_mult_pond, Year == "2022")
mult_2022 = mult_2022[,c("INDIV_ID", "Year", "POND", "Hydroperiod")]
mult_2022 = unique(mult_2022)

##how many of these are movements utilize permanent and non-permanent ponds?
salamander_mvmt_mult_pond_hydro = 
  salamander_mvmt_mult_pond %>% 
  group_by(INDIV_ID, Year) %>%
  filter(n_distinct(Hydroperiod)>1)

num_mvmts_hydro = unique(salamander_mvmt_mult_pond_hydro[,c("INDIV_ID", "Year", "INDIV_ID_Year")]) 
num_mvmts_hydro = droplevels(num_mvmts_hydro)
nlevels(num_mvmts_hydro$INDIV_ID_Year) #283 meta-years


#how do FS densities influence any of these movement types, if at all?
#do any of these movement patterns vary consistently by age or sex??

  
#Q: do those metas that had a lot of FS in their guts and only visited np ponds in 2021 show up again
#in 2022, or is that their 'quota' for a few years?


##### Meta Abund. Long-Term #####
meta_abund_year = salamander[,c("INDIV_ID", "Year")] 
meta_abund_year = unique(meta_abund_year)

meta_abund_year = meta_abund_year %>% group_by(Year) %>% count()

#plot total metas captured among all ponds during a field season 1990-2022
meta_abund_LT = 
  ggplot(meta_abund_year, aes(x=Year, y=n)) +
  geom_path(lwd = 5) +
  geom_point(size = 20) +
  #scale_y_continuous(limits = c(0,40)) +
  #scale_color_manual(name="Dataset", values = mycolors) +
  #scale_x_date(limits = as.Date(c('2022-06-14','2022-08-07'))) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_text(margin=margin(r=10)),
        axis.text = element_text(size = 50),
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x=element_text(margin=margin(t=10)))

meta_abund_LT

#ggsave(meta_abund_LT, filename = "Outputs/meta_abund_LT.png",  width = 18, height = 12, dpi = 100)


##plot number of metas in a given pond during a given field season
salamander_pond_counts = salamander %>% group_by(Year, POND) %>% count()

meta_abund_2010_pond8 = subset(salamander, Year == 2010 & POND == "8")

meta_abund_2010_pond8_counts = meta_abund_2010_pond8 %>% group_by(DATEOFYEAR, POND) %>% summarise(Meta_Count = n())


meta_2010_pond8 =
  ggplot(meta_abund_2010_pond8_counts, aes(x=DATEOFYEAR, y=Meta_Count)) +
  geom_line(aes(group = POND), lwd = 5) +
  geom_point(size = 20, alpha = 0.8) +
  #scale_x_date(limits = as.Date(c('2021-06-14','2021-07-27'))) +
  #ylim(0,30) +
  labs(y = "") +
  theme_bw(50) +
  theme(axis.title.x = element_blank(),
        panel.grid.minor = element_blank())

meta_2010_pond8
  

#ggsave(meta_2010_pond8, filename = "Outputs/metapond8_2005.png",  width = 18, height = 12, dpi = 100)

#### . ####

#### Num. Metas Into Ponds ####
#make heat map
#fill in missing cohort values for individuals that have a known cohort recorded at some point
salamander = salamander %>% group_by(INDIV_ID) %>% fill(COHORT)

#subset to include only metas with age and SVL info. Even though hurts sample size a lot (all due to lack
#of cohort data), need age and svl data to relate individual info to population growth curve
salamander = subset(salamander, !is.na(COHORT) & !is.na(SVL))
salamander = droplevels(salamander)

#need to subset movement dataframe to eliminate consecutive captures of an individual in the same
#pond in the same year
salamander = merge(salamander, unique(salamander_mvmt[,c("POND", "Hydroperiod")]), all.x = T)
salamander = salamander %>% relocate(Hydroperiod, .after = POND)
salamander = unique(salamander)

salamander_heatmap = subset(salamander, POND == "1" | POND == "10" | POND == "11" | POND == "15" |
                           POND == "5" | POND == "51" | POND == "52" | POND == "6" |
                           POND == "8" | POND == "9" | POND == "12" | POND == "13")

salamander_heatmap = salamander_heatmap %>% distinct(INDIV_ID_Year, POND, .keep_all = T)

#summarize
salamander_heatmap_sum = salamander_heatmap %>% group_by(POND, Year) %>% mutate(Unique_Meta_Visits = length(POND))
salamander_heatmap_sum = salamander_heatmap_sum %>% group_by(POND, Year) %>% mutate(Total_Surveys = n_distinct(DATEOFYEAR))
salamander_heatmap_sum = unique(salamander_heatmap_sum[,c("POND", "Year", "Unique_Meta_Visits", "Total_Surveys")])
salamander_heatmap_sum$Unique_Metas_Per_Survey = salamander_heatmap_sum$Unique_Meta_Visits / salamander_heatmap_sum$Total_Surveys

salamander_heatmap_sum = droplevels(salamander_heatmap_sum)
salamander_heatmap_sum$POND <- factor(salamander_heatmap_sum$POND,levels = c("15", "13",
                                                                        "52", "51", "11", "10",
                                                                       "8", "6", "12",
                                                                       "9", "5", "1"))
heatmap_plot = 
  ggplot(salamander_heatmap_sum, aes(x = Year, y = POND, fill = Unique_Metas_Per_Survey)) +
  geom_raster() +
  scale_x_continuous(breaks = seq(1990, 2020, by = 5), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradientn(colours = paletteer_dynamic("cartography::wine.pal", 20), na.value = "grey50") +
  theme_bw(35) + 
  labs(y = "Pond", fill = "Mean Unique\n Metamorphs Per\n Survey Occassion") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0)),
        legend.key.size = unit(2, 'cm'))

heatmap_plot

#ggsave(heatmap_plot, filename = "Outputs/heatmap_plot.png",  width = 20, height = 12, dpi = 100)
  

#stats. include an offset, which allows you to use rate data with a poisson distribution, which requires
#integer data. In other words, will model meta counts in a pond-year corrected for the total number of
#surveys in a pond-year
m1_movement = glmer(Unique_Meta_Visits ~ POND + (1|Year), offset = log(Total_Surveys), family = poisson, data = salamander_heatmap_sum)

plot(m1_movement)
car::Anova(m1_movement)
summary(m1_movement)

pond_visitation_means = emmeans(m1_movement, c("POND"), type = "response")
posthoc = as.data.frame(contrast(pond_visitation_means, interaction = "pairwise", adjust = "Bonferroni"))
pond_visitation_means_df = as.data.frame.list(pond_visitation_means)

#the Poisson distribution has only one parameter which is the rate parameter, λ (lambda).
#this is the mean number of events.

#### Landscape Factors ####

POND = c(1,5,9,12,6,8,10,11,13,15,51,52)
Hydroperiod = c(rep("Permanent", 4), rep("Semipermanent", 4), rep("Temporary", 2), rep("Semipermanent", 2))
Size = c(2237,462,254,715,91,414,139,96,22,26,124,346)
Dist_Avg = c(73.6, 127.1, 94.5, 80.8, 107.5, 79.2, 101.1, 107.8, 98.1, 141.6, 205.4, 207.5)
Dist_Avg_Perms = c(54, 116.3, 71.7, 54, 70.8, 40.3, 60.3, 67.8, 62.3, 105.8, 198, 207.3)
Elevation_Avg = c(15.5, 18.8, 15.7, 14.8, 17.5, 15.3, 16.0, 16.3, 16.3, 22, 62.7, 69.8)
Elevation_Avg_Perms = c(3.9, 5.6, 4, 3, 3.5, 2.2, 3.8, 3.7, 2.5, 8.5, 66.3, 75.3)

landscape = as.data.frame(cbind(POND, Hydroperiod, Size, Dist_Avg, Dist_Avg_Perms, Elevation_Avg, Elevation_Avg_Perms))
landscape[,1:2] <- lapply(landscape[,1:2], as.factor)
landscape[,3:7] <- lapply(landscape[,3:7], as.numeric)

#### Pond Quality Factors ####
#paedo densities
paedomorphs$Year = year(paedomorphs$DATEOFYEAR)
paedomorphs = paedomorphs %>% relocate(Year, .before = DATEOFYEAR)
paedomorphs$INDIV_ID_Year = interaction(paedomorphs$INDIV_ID, paedomorphs$Year, sep="_")
paedomorphs = paedomorphs %>% relocate(INDIV_ID_Year, .after = Year)
paedomorphs$INDIV_ID_Year = as.factor(paedomorphs$INDIV_ID_Year)

paedomorphs = subset(paedomorphs, POND == "1" | POND == "5" | POND == "12") #not including Pond 9 bc it has only one obs
paedomorphs = droplevels(paedomorphs)
hist(paedomorphs$SVL)

#keep only one obs for a unique paedo in a pond for a given year -- trying to calculate the max paedo individuals
#that were seen in a pond in a given year
paedomorphs = paedomorphs %>% distinct(INDIV_ID_Year, POND, .keep_all = T)

paedo_count_total = paedomorphs %>% group_by(Year, POND) %>% summarise(Total_Unique_Paedos = sum(n())) #number of paedos
paedo_svl_total = paedomorphs %>% group_by(Year, POND) %>% summarise(Total_Unique_Paedos_SVL = sum(SVL)) #total paedo SVL
paedo_count_svl = merge(paedo_count_total, paedo_svl_total, by = c("Year", "POND"))

plot(Total_Unique_Paedos ~ Total_Unique_Paedos_SVL, data = paedo_count_svl) #two variables are highly correlated

paedo_count_svl = merge(paedo_count_svl, landscape[,c("POND", "Size")], by = "POND")

#density of paedo SVL
paedo_count_svl$SVL_Density = paedo_count_svl$Total_Unique_Paedos_SVL / paedo_count_svl$Size


paedo_count_svl = 
  paedo_count_svl %>%
  tidyr::complete(POND, Year, fill = list(Total_Unique_Paedos = NA, Total_Unique_Paedos_SVL = NA, 
                                          Size = NA, SVL_Density = NA))


##fairy shrimp densities
POND_FS = c(6,8,10,11,13,15,51,52)
Mean_Historical_Density = c(7,4,6,4.2,2,2.2,4.7,3)

FS_Historical = as.data.frame(cbind(POND_FS, Mean_Historical_Density))
colnames(FS_Historical) <- c("POND", "Mean_Historical_Density")

FS_Historical$Mean_Historical_Density = ordered(FS_Historical$Mean_Historical_Density, levels = c(0,2,2.2,3,4,4.2,4.7,6,7))

#merge landscape and pond quality metrics into main dataframe
meta_mvmt_factors = merge(salamander_heatmap_sum, FS_Historical, all.x = TRUE)
meta_mvmt_factors$Mean_Historical_Density[is.na(meta_mvmt_factors$Mean_Historical_Density)] <- 0

meta_mvmt_factors =  merge(meta_mvmt_factors, paedo_count_svl, all = TRUE)
meta_mvmt_factors = meta_mvmt_factors[!is.na(meta_mvmt_factors$Unique_Meta_Visits),]

meta_mvmt_factors = meta_mvmt_factors[ , -which(names(meta_mvmt_factors) %in% c("Total_Unique_Paedos", "Total_Unique_Paedos_SVL", "Size"))]
names(meta_mvmt_factors)[names(meta_mvmt_factors) == 'SVL_Density'] <- 'Paedo_SVL_Density'
names(meta_mvmt_factors)[names(meta_mvmt_factors) == 'Mean_Historical_Density'] <- 'Historical_FS_Density'

meta_mvmt_factors$Paedo_SVL_Density[is.na(meta_mvmt_factors$Paedo_SVL_Density)] <- 0

meta_mvmt_factors$Paedo_SVL_Density[meta_mvmt_factors$POND == "9"] <- NA 
meta_mvmt_factors$Paedo_SVL_Density[meta_mvmt_factors$POND == "12" & meta_mvmt_factors$Year == "1996"] <- NA 
meta_mvmt_factors$Paedo_SVL_Density[meta_mvmt_factors$POND == "12" & meta_mvmt_factors$Year == "1999"] <- NA 
meta_mvmt_factors$Paedo_SVL_Density[meta_mvmt_factors$POND == "12" & meta_mvmt_factors$Year == "2000"] <- NA 

meta_mvmt_factors = merge(meta_mvmt_factors, landscape, all.x = TRUE)

#meta_mvmt_factors$Historical_FS_Density = as.factor(meta_mvmt_factors$Historical_FS_Density)
#meta_mvmt_factors$Historical_FS_Density = as.numeric(meta_mvmt_factors$Historical_FS_Density)


#divide analyses for this part into permanent and nonpermanent ponds
meta_mvmt_factors_perms = subset(meta_mvmt_factors, Hydroperiod == "Permanent")
meta_mvmt_factors_perms = na.omit(meta_mvmt_factors_perms) #to be able to compare landscape and biotic factors
#with AIC need to have same sized dataframe. Remove NAs due to no paedo surveys in pond 9 (and a few years in pond 12)

meta_mvmt_factors_nonperms = subset(meta_mvmt_factors, Hydroperiod != "Permanent")
#cor.test(meta_mvmt_factors_nonperms$Historical_FS_Density, 
#         meta_mvmt_factors_nonperms$Size, method = "spearman")

##model perms
#make dataframe with correlations
corr_var_perms = meta_mvmt_factors_perms[,c("Size", "Dist_Avg", "Elevation_Avg", "Paedo_SVL_Density")]

corr.test(corr_var_perms)
library(Hmisc)

#make function
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

perms_corr<-rcorr(as.matrix(corr_var_perms))
perms_corr = flattenCorrMatrix(perms_corr$r, perms_corr$P)
perms_corr$cor = abs(perms_corr$cor)
detach(package:Hmisc)

#model
m1_perms_movement = glmer(Unique_Meta_Visits ~ scale(Size) + scale(Elevation_Avg) + (1|Year), 
                        offset = log(Total_Surveys), family = poisson, data = meta_mvmt_factors_perms)

plot(m1_perms_movement)
car::Anova(m1_perms_movement)
summary(m1_perms_movement)
plot(Unique_Metas_Per_Survey ~ scale(Elevation_Avg), data = meta_mvmt_factors_perms)
boxplot(Unique_Metas_Per_Survey ~ scale(Elevation_Avg), data = meta_mvmt_factors_perms)


perm_predict1 <- ggpredict(m1_perms_movement, "Size [all]")
#plot(perm_predict1, add.data = TRUE)
plot(perm_predict1, residuals = TRUE)

perm_predict2 <- ggpredict(m1_perms_movement, "Elevation_Avg [all]")
#plot(perm_predict2, add.data = TRUE)
plot(perm_predict2, residuals = TRUE)

m2_perms_movement = glmer(Unique_Meta_Visits ~ scale(Paedo_SVL_Density) + (1|Year), 
                          offset = log(Total_Surveys), family = poisson, data = meta_mvmt_factors_perms)

AIC(m1_perms_movement, m2_perms_movement)


###model nonperms
#make dataframe with correlations
corr_var_nonperms = meta_mvmt_factors_nonperms[,c("Size", "Dist_Avg", "Elevation_Avg")]

corr.test(corr_var_nonperms)

library(Hmisc)
#make function
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

nonperms_corr<-rcorr(as.matrix(corr_var_nonperms))
nonperms_corr = flattenCorrMatrix(nonperms_corr$r, nonperms_corr$P)
nonperms_corr$cor = abs(nonperms_corr$cor)
detach(package:Hmisc)

#model
m1_nonperms_fs = glmer(Unique_Meta_Visits ~ Historical_FS_Density + (1|Year), 
           offset = log(Total_Surveys), family = poisson, data = meta_mvmt_factors_nonperms)

m2_nonperms_landscape = glmer(Unique_Meta_Visits ~ scale(Size) + scale(Dist_Avg) + (1|Year), 
           offset = log(Total_Surveys), family = poisson, data = meta_mvmt_factors_nonperms)

m3_nonperms_fs_landscape = glmer(Unique_Meta_Visits ~ scale(Size) + Historical_FS_Density + (1|Year), 
                              offset = log(Total_Surveys), family = poisson, data = meta_mvmt_factors_nonperms)

car::Anova(m3_nonperms_fs_landscape)
summary(m3_nonperms_fs_landscape)
tab_model(m3_nonperms_fs_landscape)
m3_nonperms_fs_landscape_predict <- ggpredict(m3_nonperms_fs_landscape, terms = "Historical_FS_Density")
plot(m3_nonperms_fs_landscape_predict, residuals = T)

AIC(m1_nonperms_fs, m2_nonperms_landscape, m3_nonperms_fs_landscape)

#landscape stats and plot
car::Anova(m2_nonperms_landscape)
m2_nonperms_landscape_predict <- ggpredict(m2_nonperms_landscape, terms = "Size [all]")
plot(m2_nonperms_landscape_predict, residuals = T)
nonperms_pondsize_plot = plot(m2_nonperms_landscape_predict, residuals = F)

nonperms_pondsize_plot2 = nonperms_pondsize_plot + 
  geom_line(size = 3) +
  ylab("Metamorph ID Visits") +
  xlab(expression(paste("Pond Size", " ", "(", m^{2}, ")"))) + 
  coord_cartesian(xlim = c(0,400), ylim = c(0, 4)) +
  scale_y_continuous(breaks = c(0,1,2,3,4)) +
  scale_x_continuous(breaks = c(0,100,200,300,400)) +
  theme_bw(50) +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 40),
        axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)))

nonperms_pondsize_plot2

#ggsave(nonperms_pondsize_plot2, filename = "Outputs/meta_movement_landscape_pond_size.png",  width = 13, height = 12, dpi = 100)

#fs stats and plot
plot(m1_nonperms_fs)
car::Anova(m1_nonperms_fs)
summary(m1_nonperms_fs)
pond_visitation_fs_means = emmeans(m1_nonperms_fs, c("Historical_FS_Density"), type = "response")
posthoc2 = as.data.frame(contrast(pond_visitation_fs_means, interaction = "pairwise", adjust = "Bonferroni"))


m1_nonperms_fs_predict <- ggpredict(m1_nonperms_fs, terms = "Historical_FS_Density")
plot(m1_nonperms_fs_predict, add.data = TRUE)
plot(m1_nonperms_fs_predict)
nonperms_fs_plot = plot(m1_nonperms_fs_predict, residuals = F, dot.size = 10, line.size = 4)

nonperms_fs_plot2 = nonperms_fs_plot + 
  geom_point(size = 12) +
  ylab("Metamorph ID Visits") +
  xlab("Mean Historical Fairy Shrimp Density Rank") + 
  coord_cartesian(ylim = c(0, 4)) +
  scale_y_continuous(breaks = c(0,1,2,3,4)) +
  scale_x_reverse(breaks = c(1,2,3,4,5,6,7,8)) +
  theme_bw(50) +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 40),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)))

nonperms_fs_plot2

#ggsave(nonperms_fs_plot2, filename = "Outputs/meta_movement_FS_plot.png",  width = 13, height = 12, dpi = 100)



#### Which Metas Into Ponds ####

#make growth curve for all metas
salamander_growth = salamander
  
salamander_growth$AGE = salamander_growth$Year - salamander_growth$COHORT
salamander_growth = subset(salamander_growth, AGE > 0)
salamander_growth_male = subset(salamander_growth, SEX == "Male")
salamander_growth_female = subset(salamander_growth, SEX == "Female")

plot(SVL ~ AGE, data = salamander_growth)
plot(SVL ~ AGE, data = salamander_growth_male)
plot(SVL ~ AGE, data = salamander_growth_female)

fishmethods::growth(size=salamander_growth$SVL,age=salamander_growth$AGE,
       error=1, Sinf=101,K=0.1,t0=-12) #good starting values

#resource for below code: https://danstich.github.io/we-r-nycafs/fishStats.html

#define von Bertalanffy growth function
vbmod <- SVL ~ Linf * (1 - exp(-K * (AGE - t0)))


#fit the von Bertalanffy growth function using nonlinear least squares (nls) optimization
growth_mod <- nls(vbmod, data = salamander_growth, start = list(Linf = 101, K = 0.1, t0 = -12))
summary(growth_mod)

#predict
growth_predict = predict(growth_mod)

##bootstrap 95% CI
#get the desired growth function from a list of those that are available in FSA
vbO <- vbFuns("typical")

#fit the model to the data using nls, like we did before
vb_fit <- nls(SVL~vbO(AGE,Linf,K, t0), data = salamander_growth, start = list(Linf = 101, K = 0.1, t0 = -12))

#now, bootstrap the model fitting process
boot_fit <- nlsBoot(vb_fit, niter = 9999)

#predict length at age from the model (t is age). Here, we tell R to predict length at 
#each unique age in our original data, and calculate some bootstrapped confidence
#intervals
boot_preds <- data.frame(predict(boot_fit, vbO, t = sort(unique(salamander_growth$AGE))))
names(boot_preds) <- c("AGE", "fit", "lwr", "upr")

growth_preds <- merge(salamander_growth, boot_preds, by = "AGE")

#plot
growth_curve_plot = 
  ggplot(growth_preds, aes(x = AGE, y = SVL)) +
  geom_jitter(width = 0.1, alpha = 0.15, size = 5) +
  geom_line(aes(y = fit), linewidth = 2) +
  geom_ribbon(
    aes(x = AGE, ymin = lwr, ymax = upr, color = NULL), alpha = 0.3) +
  scale_y_continuous(limits = c(59,120)) +
  xlab("Age (years)") +
  ylab("Snout-Vent Length (mm)") +
  theme_bw(40)

growth_curve_plot

#ggsave(growth_curve_plot, filename = "Outputs/growth_curve_plot.png",  width = 14, height = 12, dpi = 100)


#calculate residuals 
growth_preds = growth_preds[,c("DATEOFYEAR", "Year", "INDIV_ID", "INDIV_ID_Year",
                               "MORPH", "POND", "COHORT", "SEX", "WEIGHT", "SVL",
                               "fit", "lwr", "upr")]

names(growth_preds)[names(growth_preds) == 'fit'] <- 'SVL_fit'
names(growth_preds)[names(growth_preds) == 'lwr'] <- 'SVL_fit_lwr'
names(growth_preds)[names(growth_preds) == 'upr'] <- 'SVL_fit_upr'

growth_preds$SVL_residual = growth_preds$SVL - growth_preds$SVL_fit 
hist(growth_preds$SVL_residual)

growth_preds = droplevels(growth_preds)
nlevels(growth_preds$INDIV_ID)

num_captures = growth_preds %>% group_by(INDIV_ID) %>% summarise(N = n())
range(num_captures$N)
mean(num_captures$N)
sd(num_captures$N)

#create df for SVL residual data for metamorphs that were captured at least twice. In other words, eliminate
#metas that were captured in one pond in one year only. These individuals don't have the opportunity to learn
#about pond resources, and they're more likely to end up in ponds based on landscape factors (e.g., biggest pond)
growth_preds_multiple_obs = growth_preds %>% 
  add_count(INDIV_ID, name = "Num_Captures") %>% 
  filter(Num_Captures > 1)

#add in mean SVL residual data for each meta over their lifetime
growth_preds_avg = growth_preds_multiple_obs %>% group_by(INDIV_ID) %>% summarise(SVL_residual_mean = mean(SVL_residual))
hist(growth_preds_avg$SVL_residual_mean)
growth_preds = inner_join(growth_preds, growth_preds_avg, by = "INDIV_ID")  
growth_preds = droplevels(growth_preds)

#subset to just ponds we care about
growth_preds_subset = subset(growth_preds, POND == "1" | POND == "10" | POND == "11" | POND == "15" |
                              POND == "5" | POND == "51" | POND == "52" | POND == "6" |
                              POND == "8" | POND == "9" | POND == "12" | POND == "13")

#model
growth_preds_subset_reorder <- growth_preds_subset %>% 
  mutate(POND = case_when(
    POND == '1' ~ '1', POND == '5' ~ '2', POND == '9' ~ '3', POND == '12' ~ '4', POND == '6' ~ '5', 
    POND == '8' ~ '6', POND == '10' ~ '7', POND == '11' ~ '8', POND == '13' ~ '11', POND == '15' ~ '12', 
    POND == '51' ~ '9', POND == '52' ~ '10'))

growth_preds_subset_reorder$POND = as.factor(growth_preds_subset_reorder$POND)

m1_svl_movement = lmer(SVL_residual ~ POND + (1|INDIV_ID) + (1|Year), data = growth_preds_subset_reorder)

plot(m1_svl_movement)
car::Anova(m1_svl_movement)
summary(m1_svl_movement)
pond_which_meta_means = emmeans(m1_svl_movement, c("POND"), type = "response")
posthoc3 = as.data.frame(contrast(pond_which_meta_means, interaction = "pairwise", adjust = "Bonferroni"))


growth_residuals_predict <- ggpredict(m1_svl_movement, "POND")
which_metas_plot_base = plot(growth_residuals_predict, line.size = 3.5)

which_metas_plot_base2 =
  which_metas_plot_base + 
  geom_hline(yintercept=0, linetype="dashed", color = "dodgerblue", size=3) +
  geom_point(size = 10) +
  ylab("Metamorph SVL Residual") +
  xlab("Pond Number") + 
  theme_bw(50) +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 40),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)))

#ggsave(which_metas_plot_base2, filename = "Outputs/which_metas_plot_base.png",  width = 14, height = 12, dpi = 100)


which_metas_plot = plot(growth_residuals_predict, add.data = TRUE, dot.size = 8, line.size = 3.5)

which_metas_plot2 = which_metas_plot + 
  geom_hline(yintercept=0, linetype="dashed", color = "dodgerblue", size=3) +
  geom_point(size = 10) +
  ylab("Metamorph SVL Residual") +
  xlab("Pond Number") + 
  #scale_y_continuous(limits = c(-15,15), breaks = c(-15,-10,-5,0,5,10,15)) +
  theme_bw(50) +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 40),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)))

which_metas_plot2

#ggsave(which_metas_plot2, filename = "Outputs/which_metas_plot_base.png",  width = 14, height = 12, dpi = 100)



#there is some shrinkage going on. More extreme INDIV_ID clusters brought back toward the all cluster (population)
#average. Those clusters with more observations will generally exhibit less shrinkage, and those with 
#fewer observations the opposite. Many of the metas have a handful of captures only. 
#Some good resources:
#https://m-clark.github.io/posts/2019-05-14-shrinkage-in-mixed-models/
#https://www.tjmahr.com/plotting-partial-pooling-in-mixed-effects-models/

#interpretation of emmeans from m1_svl_movement model: no ponds regularly draw metas that are bigger for
#their size, but some seem to really draw metas that are small for their size.

#merge in landscape and biotic factors
names(meta_mvmt_factors)[names(meta_mvmt_factors) == 'Pond'] <- 'POND'
svl_movement = merge(growth_preds_subset, meta_mvmt_factors[,c("POND", "Year", "Historical_FS_Density", "Paedo_SVL_Density",
                            "Hydroperiod", "Size", "Dist_Avg", "Dist_Avg_Perms", "Elevation_Avg", "Elevation_Avg_Perms")],
                            by = c("POND", "Year"), all.x = T)

svl_movement_male_perm = subset(svl_movement, SEX == "Male" & Hydroperiod == "Permanent")
svl_movement_female_perm = subset(svl_movement, SEX == "Female" & Hydroperiod == "Permanent")
svl_movement_male_nonperm = subset(svl_movement, SEX == "Male" & Hydroperiod != "Permanent")
svl_movement_female_nonperm = subset(svl_movement, SEX == "Female" & Hydroperiod != "Permanent")
svl_movement_nonperm = subset(svl_movement, Hydroperiod != "Permanent")

svl_movement_male_perm = na.omit(svl_movement_male_perm)
svl_movement_female_perm = na.omit(svl_movement_female_perm)


#male perms
m1_svl_movement_male_perm = lmer(SVL_residual ~ POND + (1|INDIV_ID) + (1|Year), data = svl_movement_male_perm)
plot(m1_svl_movement_male_perm)
car::Anova(m1_svl_movement_male_perm)
summary(m1_svl_movement_male_perm)
residuals_svl_movement_male_perm = emmeans(m1_svl_movement_male_perm, c("POND"), type = "response")
residuals_svl_movement_male_perm = as.data.frame.list(residuals_svl_movement_male_perm)

m2_svl_movement_male_perm = lmer(SVL_residual ~ scale(Size) + scale(Elevation_Avg) + scale(Paedo_SVL_Density) + 
                                (1|INDIV_ID) + (1|Year), data = svl_movement_male_perm)
plot(m2_svl_movement_male_perm)
car::Anova(m2_svl_movement_male_perm)
summary(m2_svl_movement_male_perm)

svl_movement_male_perm_predict1 <- ggpredict(m2_svl_movement_male_perm, "Size [all]")
plot(svl_movement_male_perm_predict1, add.data = TRUE)
plot(svl_movement_male_perm_predict1, residuals = TRUE)

#female perms
m1_svl_movement_female_perm = lmer(SVL_residual ~ POND + (1|INDIV_ID) + (1|Year), data = svl_movement_female_perm)
plot(m1_svl_movement_female_perm)
car::Anova(m1_svl_movement_female_perm)

m2_svl_movement_female_perm = lmer(SVL_residual ~ scale(Size) + scale(Elevation_Avg) + scale(Paedo_SVL_Density) + 
                                   (1|INDIV_ID) + (1|Year), data = svl_movement_female_perm)
plot(m2_svl_movement_female_perm)
car::Anova(m2_svl_movement_female_perm)

#male nonperms
m1_svl_movement_male_nonperm = lmer(SVL_residual ~ POND + (1|INDIV_ID) + (1|Year), data = svl_movement_male_nonperm)
plot(m1_svl_movement_male_nonperm)
car::Anova(m1_svl_movement_male_nonperm)
summary(m1_svl_movement_male_nonperm)

svl_movement_male_nonperm_predict1 <- ggpredict(m1_svl_movement_male_nonperm, "POND")
plot(svl_movement_male_nonperm_predict1, add.data = TRUE)
plot(svl_movement_male_nonperm_predict1, residuals = TRUE)

m2_svl_movement_male_nonperm = lmer(SVL_residual ~ Historical_FS_Density + 
                                     (1|INDIV_ID) + (1|Year), data = svl_movement_male_nonperm)
plot(m2_svl_movement_male_nonperm)
car::Anova(m2_svl_movement_male_nonperm)

m2_svl_movement_male_nonperm_predict1 <- ggpredict(m2_svl_movement_male_nonperm, "Historical_FS_Density")
plot(m2_svl_movement_male_nonperm_predict1, add.data = TRUE)
plot(m2_svl_movement_male_nonperm_predict1, residuals = TRUE)


#female nonperms
m1_svl_movement_female_nonperm = lmer(SVL_residual ~ POND + (1|INDIV_ID) + (1|Year), data = svl_movement_female_nonperm)
plot(m1_svl_movement_female_nonperm)
car::Anova(m1_svl_movement_female_nonperm)
summary(m1_svl_movement_female_nonperm)

svl_movement_female_nonperm_predict1 <- ggpredict(m1_svl_movement_female_nonperm, "POND")
plot(svl_movement_female_nonperm_predict1, add.data = TRUE)
plot(svl_movement_female_nonperm_predict1, residuals = TRUE)

m2_svl_movement_female_nonperm = lmer(SVL_residual ~ Historical_FS_Density + 
                                      (1|INDIV_ID) + (1|Year), data = svl_movement_female_nonperm)
plot(m2_svl_movement_female_nonperm)
car::Anova(m2_svl_movement_female_nonperm)

m2_svl_movement_female_nonperm_predict1 <- ggpredict(m2_svl_movement_female_nonperm, "Historical_FS_Density")
plot(m2_svl_movement_female_nonperm_predict1, add.data = TRUE)
plot(m2_svl_movement_female_nonperm_predict1, residuals = TRUE)


#both sexes, nonperms
m2_svl_movement_nonperms = lmer(SVL_residual ~ Historical_FS_Density + 
                         (1|INDIV_ID) + (1|Year), data = svl_movement_nonperm)
plot(m2_svl_movement_nonperms)
car::Anova(m2_svl_movement_nonperms)
pond_which_meta_fs_means = emmeans(m2_svl_movement_nonperms, c("Historical_FS_Density"), type = "response")
posthoc4 = as.data.frame(contrast(pond_which_meta_fs_means, interaction = "pairwise", adjust = "Bonferroni"))


m3_svl_movement_nonperms = lmer(SVL_residual ~ scale(Size) + scale(Elevation_Avg) + (1|INDIV_ID) + (1|Year), data = svl_movement_nonperm)
plot(m3_svl_movement_nonperms)
car::Anova(m3_svl_movement_nonperms)
summary(m3_svl_movement_nonperms)

m3_svl_movement_nonperms_predict1 <- ggpredict(m3_svl_movement_nonperms, "Size")
plot(m3_svl_movement_nonperms_predict1, add.data = T)

AIC(m2_svl_movement_nonperms, m3_svl_movement_nonperms)

#write.csv(growth_preds, "Data/meta_growth_curve_residuals.csv")

#plot fs model
m2_svl_movement_nonperms_predict1 <- ggpredict(m2_svl_movement_nonperms, "Historical_FS_Density")
plot(m2_svl_movement_nonperms_predict1)
nonperms_svl_fs_plot_base = plot(m2_svl_movement_nonperms_predict1, line.size = 3.5)

nonperms_svl_fs_plot_base2 =
  nonperms_svl_fs_plot_base + 
  geom_hline(yintercept=0, linetype="dashed", color = "dodgerblue", size=3) +
  geom_point(size = 10) +
  scale_y_continuous(limits = c(-6.5, 3), breaks = c(-6,-3,0,3)) +
  scale_x_reverse(breaks = c(1,2,3,4,5,6,7,8)) +
  ylab("Metamorph SVL Residual") +
  xlab("Rank of Mean Historical Fairy Shrimp Densities") + 
  theme_bw(50) +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 40),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)))

nonperms_svl_fs_plot_base2

#ggsave(nonperms_svl_fs_plot_base2, filename = "Outputs/which_metas_fs_plot_base.png",  width = 14, height = 12, dpi = 100)


nonperms_svl_fs_plot = plot(m2_svl_movement_nonperms_predict1, add.data = TRUE, dot.size = 8, line.size = 3.5)

nonperms_svl_fs_plot2 = nonperms_svl_fs_plot + 
  geom_hline(yintercept=0, linetype="dashed", color = "dodgerblue", size=3) +
  geom_point(size = 10) +
  ylab("Metamorph SVL Residual") +
  xlab("Rank of Mean Historical Fairy Shrimp Densities") + 
  scale_x_reverse(breaks = c(1,2,3,4,5,6,7,8)) +
  theme_bw(50) +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 40),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)))

nonperms_svl_fs_plot2

#ggsave(nonperms_svl_fs_plot2, filename = "Outputs/which_metas_fs_plot.png",  width = 14, height = 12, dpi = 100)



#### SVL Resid ~ Mvmt Pat. ####
#work fromthis df: growth_preds_subset
#look at mean svl_residuals for individuals and relate it to:
#(1) how many over total capture history
#(2) how many ponds found in per year on average
#(3) number of hydroperiod switches
#(4) how long between recaptures

growth_svl_means = unique(growth_preds_subset[,c("INDIV_ID", "MORPH", "SEX", "SVL_residual_mean")])
growth_svl_means = droplevels(growth_svl_means)

#make df for metas for which we have svl residual data
capture_history_svl_means = inner_join(salamander, growth_svl_means[,c("INDIV_ID", "SVL_residual_mean")], by = "INDIV_ID")

#(1) how many over total capture history
metas_distinct_ponds_lifetime =  capture_history_svl_means %>% 
  group_by(INDIV_ID) %>%
  summarize(Total_Ponds = n_distinct(POND))

metas_distinct_years_seen =  capture_history_svl_means %>% 
  group_by(INDIV_ID) %>%
  summarize(Years_Seen = n_distinct(Year))

metas_distinct_ponds_lifetime = merge(metas_distinct_ponds_lifetime, metas_distinct_years_seen, by = "INDIV_ID")

metas_distinct_ponds_lifetime = merge(metas_distinct_ponds_lifetime, growth_svl_means, by = "INDIV_ID")
metas_distinct_ponds_lifetime$Total_Ponds_Per_Year_Seen = metas_distinct_ponds_lifetime$Total_Ponds / metas_distinct_ponds_lifetime$Years_Seen
#metas_distinct_ponds_lifetime = subset(metas_distinct_ponds_lifetime, SEX == "Female")

m1_movement_pattern = lm(SVL_residual_mean ~ Total_Ponds_Per_Year_Seen, data = metas_distinct_ponds_lifetime)

#plot(m1_movement_pattern)
car::Anova(m1_movement_pattern)
summary(m1_movement_pattern)

movement_pattern_predict1 <- ggpredict(m1_movement_pattern, "Total_Ponds_Per_Year_Seen [all]")
total_ponds_plot = plot(movement_pattern_predict1, add.data = TRUE, dot.size = 5)

total_ponds_plot2 = 
  total_ponds_plot +
  geom_line(size = 3) +
  ylab("Lifetime Mean SVL Residual") +
  xlab("Lifetime Pond Visitation") + 
  theme_bw(50) +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 40),
        axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)))

total_ponds_plot2

#Pond Exploration Index = Distinct Ponds Over Lifetime / Number of Years Seen
#so in other words, how many total ponds did you visit, corrected for the total number of years you were seen

#ggsave(total_ponds_plot2, filename = "Outputs/lifetime_exploration_plot.png",  width = 14, height = 12, dpi = 100)


#(2) how many ponds found in per year on average -- gets at switching tendency
metas_distinct_ponds_year_avg =  capture_history_svl_means %>% 
  group_by(INDIV_ID, Year) %>%
  summarise(Total_Ponds_Year = n_distinct(POND)) %>%
  summarise(Total_Ponds_Year_Avg = mean(Total_Ponds_Year)) 

metas_distinct_ponds_year_avg = merge(metas_distinct_ponds_year_avg, growth_svl_means, by = "INDIV_ID")

m2_movement_pattern = lm(SVL_residual_mean ~ Total_Ponds_Year_Avg, data = metas_distinct_ponds_year_avg)

#plot(m2_movement_pattern)
car::Anova(m2_movement_pattern)
summary(m2_movement_pattern)

movement_pattern_predict2 <- ggpredict(m2_movement_pattern, "Total_Ponds_Year_Avg")
plot(movement_pattern_predict2, add.data = TRUE)
plot(SVL_residual_mean ~ Total_Ponds_Year_Avg, data = metas_distinct_ponds_year_avg)

annual_exploration_plot = 
  ggplot(data = metas_distinct_ponds_year_avg, aes(x = Total_Ponds_Year_Avg, y = SVL_residual_mean)) +
  geom_point(size = 5, color = "grey40", alpha = 0.7)+
  geom_line(data = movement_pattern_predict2, aes(x= x, y = predicted), size=3) +
  geom_line(data = movement_pattern_predict2, aes(x = x , y = conf.low), linetype = "dashed", size = 2) +
  geom_line(data = movement_pattern_predict2, aes(x = x , y = conf.high), linetype = "dashed", size = 2) +
  ylab("Lifetime Mean SVL Residual") +
  xlab("Annual Pond Visitation") + 
  theme_bw(50) +
  theme(plot.title = element_blank(),
          axis.title = element_text(size = 40),
          axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)))

annual_exploration_plot

#ggsave(annual_exploration_plot, filename = "Outputs/annual_exploration_plot.png",  width = 14, height = 12, dpi = 100)
  

#(3) how many hydroperiod switches per year on average
metas_distinct_hydroperiod_year_avg =  capture_history_svl_means %>% 
  group_by(INDIV_ID, Year) %>%
  summarise(Total_Hydroperiods_Year = n_distinct(Hydroperiod)) %>%
  summarise(Total_Hydroperiods_Year_Avg = mean(Total_Hydroperiods_Year)) 

#(4) how long between recaptures
recap_interval = growth_preds_subset[,c("Year", "INDIV_ID", "SVL_residual")]

recap_interval_year_diff = unique(recap_interval[,c("INDIV_ID", "Year")])

recap_interval_year_diff = 
  recap_interval_year_diff %>%
  group_by(INDIV_ID) %>%
  arrange(INDIV_ID, desc(Year)) %>%
  mutate(Year_Interval = lag(Year, default = first(Year)) - Year)

recap_interval_svl_resid_diff = unique(recap_interval[,c("Year", "INDIV_ID", "SVL_residual")])

recap_interval_svl_resid_diff = recap_interval_svl_resid_diff %>% group_by(INDIV_ID, Year) %>% summarise(SVL_residual = mean(SVL_residual))

recap_interval_svl_resid_diff = 
  recap_interval_svl_resid_diff %>%
  group_by(INDIV_ID) %>%
  arrange(INDIV_ID, desc(Year)) %>%
  mutate(SVL_Resid_Change = lag(SVL_residual, default = first(SVL_residual)) - SVL_residual)

svl_delta_v_recap_interval = merge(recap_interval_year_diff, recap_interval_svl_resid_diff, 
                                   by = c("INDIV_ID", "Year"))

svl_delta_v_recap_interval = subset(svl_delta_v_recap_interval, Year_Interval > 0)
svl_delta_v_recap_interval$Year <- NULL
svl_delta_v_recap_interval$SVL_residual <- NULL
plot(SVL_Resid_Change ~ Year_Interval, data = svl_delta_v_recap_interval)

svl_delta_v_recap_interval2 = subset(svl_delta_v_recap_interval, Year_Interval < 9)

m1_svl_recap = lmer(SVL_Resid_Change ~ Year_Interval + (1|INDIV_ID), data = svl_delta_v_recap_interval)

plot(m1_svl_recap)
car::Anova(m1_svl_recap)
summary(m1_svl_recap)

m1_svl_recap_predict1 <- ggpredict(m1_svl_recap, "Year_Interval")
plot(m1_svl_recap_predict1, add.data = TRUE)

plot(SVL_Resid_Change ~ Year_Interval, data = svl_delta_v_recap_interval)

#need to work on this one. Is the data wonky and needs to be transformed further before plotting or
#is it just a function of having less data as you progress along the x-axis and the smaller chunks
#of data each are center around zero, creating this cone shape


#Delta SVL ~ FS ####
delta_svl = growth_preds_subset %>% group_by(INDIV_ID) %>% filter(length(INDIV_ID) >= 2)
quantile(growth_preds_subset$SVL_residual, probs = 0.5)
delta_svl_smallest = subset(growth_preds_subset, SVL_residual <= 0.4107126)
delta_svl_smallest = droplevels(delta_svl_smallest)

delta_svl = delta_svl %>% group_by(INDIV_ID, POND) %>% summarise(Count = n())

delta_svl = delta_svl %>% group_by(INDIV_ID) %>% mutate(Percent = Count/sum(Count))

delta_svl = subset(delta_svl, Percent > 0.65)
names(delta_svl)[names(delta_svl) == 'POND'] <- 'POND_Favored'

delta_svl = inner_join(salamander, delta_svl[,c("INDIV_ID", "POND_Favored")], by = "INDIV_ID")  

delta_svl = delta_svl %>%
  group_by(INDIV_ID) %>%
  arrange(DATEOFYEAR, .by_group = TRUE) %>%
  mutate(Year_Diff = max(Year) - min(Year)) %>%
  mutate(SVL_Diff = max(SVL) - min(SVL))

delta_svl$SVL_Growth_Rate = delta_svl$SVL_Diff / delta_svl$Year_Diff

delta_svl = subset(delta_svl, !is.infinite(SVL_Growth_Rate) & !is.nan(SVL_Growth_Rate)) #first and last capture date in same year
delta_svl = unique(delta_svl[,c("INDIV_ID", "POND_Favored", "SVL_Growth_Rate")])
delta_svl = droplevels(delta_svl)

delta_svl$POND_Favored <- factor(delta_svl$POND_Favored,levels = c("1", "5", "9", "12", "6", "8", "10", "11",
                                                                   "13", "15", "51", "52"))

delta_svl = subset(delta_svl, SVL_Growth_Rate < 10) #one outlier value that is probs a data entry error

salamander_svl_rates_plot =
  ggplot(delta_svl, aes(x=POND_Favored, y=SVL_Growth_Rate)) + 
  geom_boxplot(color="black", fill="white", alpha=08, lwd = 2, outlier.shape = NA) + 
  geom_jitter(color="black", size=8, alpha=0.7, width = 0.1, height = 0) +
  #scale_y_continuous(breaks=c(-5,0,5,10), limits=c(-5,10)) +
  labs(x = "Favored Pond", y = "SVL Growth Rate (mm/yr)")+
  theme_bw(50)

salamander_svl_rates_plot

#ggsave(salamander_svl_rates_plot, filename = "Outputs/salamander_svl_rates_plot.png",  width = 18, height = 12, dpi = 100)


m1_delta_svl = lm(SVL_Growth_Rate ~ POND_Favored, data = delta_svl)
plot(m1_delta_svl)
car::Anova(m1_delta_svl)


##support for this graph

#body condition in ponds
salamander_heatmap$Condition = salamander_heatmap$WEIGHT / salamander_heatmap$SVL

salamander_heatmap$POND <- factor(salamander_heatmap$POND,levels = c("1", "5", "9", "12", "6", "8", "10", "11",
                                                                    "51", "52", "13", "15"))

salamander_condition_pond =
  ggplot(salamander_heatmap, aes(x=POND, y=Condition)) + 
  geom_boxplot(color="black", fill="white", alpha=08, lwd = 2, outlier.shape = NA) + 
  geom_jitter(color="black", size=8, alpha=0.7, width = 0.1, height = 0) +
  #scale_y_continuous(breaks=c(-5,0,5,10), limits=c(-5,10)) +
  labs(x = "Pond", y = "Body Condition")+
  theme_bw(50)

salamander_condition_pond

#ggsave(salamander_condition_pond, filename = "Outputs/salamander_condition_pond.png",  width = 18, height = 12, dpi = 100)


#CHANGE in body condition in ponds
delta_condition = salamander %>% group_by(INDIV_ID, Year, POND) %>% filter(n() > 1)

delta_condition$Condition = delta_condition$WEIGHT / delta_condition$SVL

delta_condition = delta_condition %>%
  group_by(INDIV_ID, Year, POND) %>%
  arrange(DATEOFYEAR, .by_group = TRUE) %>%
  mutate(Condition_Change = Condition - lag(Condition, default = first(Condition))) %>%
  mutate(Date_Diff = DATEOFYEAR - lag(DATEOFYEAR, default = first(DATEOFYEAR)))

##
delta_condition = semi_join(delta_condition, delta_svl_smallest, by = "INDIV_ID")  
##

delta_condition = delta_condition[,c("POND", "INDIV_ID", "Condition_Change", "Date_Diff")]

delta_condition = subset(delta_condition, Condition_Change != 0 & POND != "1A" & !is.na(POND))
delta_condition = droplevels(delta_condition)
delta_condition = subset(delta_condition, Condition_Change < 0.6 & Condition_Change > -0.5) #outliers?

delta_condition$Date_Diff = as.numeric(delta_condition$Date_Diff)
delta_condition$Condition_Change_Rate = delta_condition$Condition_Change / delta_condition$Date_Diff
delta_condition = subset(delta_condition, !is.infinite(Condition_Change_Rate))

delta_condition$POND <- factor(delta_condition$POND,levels = c("1", "5", "9", "12", "6", "8", "10", "11",
                                                                     "13", "15", "51", "52"))

salamander_delta_condition_pond =
  ggplot(delta_condition, aes(x=POND, y=Condition_Change_Rate)) + 
  #geom_jitter(color="black", size=8, alpha=0.2, width = 0.1, height = 0) +
  geom_boxplot(color="black", fill="white", alpha=0.8, lwd = 2, outlier.shape = NA) + 
  #scale_y_continuous(limits=c(-0.15,0.15)) +
  labs(x = "Pond", y = "Change in Body Condition")+
  theme_bw(50)

salamander_delta_condition_pond

delta_condition_reorder <- delta_condition %>% 
  mutate(POND = case_when(
    POND == '1' ~ '1', POND == '5' ~ '2', POND == '9' ~ '3', POND == '12' ~ '4', POND == '6' ~ '5', 
    POND == '8' ~ '6', POND == '10' ~ '7', POND == '11' ~ '8', POND == '13' ~ '9', POND == '15' ~ '10', 
    POND == '52' ~ '11'))

delta_condition_reorder$POND = as.factor(delta_condition_reorder$POND)

m1_delta_condition = lmer(Condition_Change_Rate ~ POND + (1|INDIV_ID), data = delta_condition_reorder)
plot(m1_delta_condition)
car::Anova(m1_delta_condition)
delta_condition_means = emmeans(m1_delta_condition, c("POND"), type = "response")
posthoc5 = as.data.frame(contrast(delta_condition_means, interaction = "pairwise", adjust = "Bonferroni"))


m1_delta_condition_predict <- ggpredict(m1_delta_condition, "POND")
plot(m1_delta_condition_predict)
plot(m1_delta_condition_predict, add.data = TRUE)
delta_condition_plot = plot(m1_delta_condition_predict, residuals = TRUE, dot.size = 4, line.size = 3,
                            ci = TRUE)

delta_condition_plot2 = delta_condition_plot + 
  geom_hline(yintercept=0, linetype="dashed", color = "dodgerblue", size=2) +
  geom_point(size = 8, color = "red") +
  ylab("Change in Metamorph Body Condition") +
  xlab("Pond Number") + 
  scale_y_continuous(limits = c(-0.3, 0.45), breaks = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4)) +
  theme_bw(50) +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 40),
        axis.text.x = element_blank(),
        axis.title.x = element_text(margin = margin(t = 60, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)))

delta_condition_plot2

#ggsave(delta_condition_plot2, filename = "Outputs/delta_condition_plot.png",  width = 13, height = 12, dpi = 100)

#do metas know where they are going? have they been there before?
salamander_pond_knowledge = unique(salamander[,c("POND", "INDIV_ID", "Year")])
  
salamander_pond_knowledge = salamander_pond_knowledge %>% group_by(INDIV_ID) %>% count(POND, sort = TRUE)

salamander_pond_knowledge = subset(salamander_pond_knowledge, POND == 1 | POND == 5 | POND == 9 | POND == 12 | POND == 6 | 
                                     POND == 8 | POND == 10 | POND == 11 | POND == 13 | POND == 15 | POND == 51 | POND == 52)

salamander_pond_knowledge = droplevels(salamander_pond_knowledge)
salamander_pond_knowledge$POND <- factor(salamander_pond_knowledge$POND,levels = c("1", "5", "9", "12", "6", "8", "10", "11",
                                                               "13", "15", "51", "52"))

salamander_pond_knowledge <- salamander_pond_knowledge %>% 
  mutate(POND = case_when(
    POND == '1' ~ '1', POND == '5' ~ '2', POND == '9' ~ '3', POND == '12' ~ '4', POND == '6' ~ '5', 
    POND == '8' ~ '6', POND == '10' ~ '7', POND == '11' ~ '8', POND == '13' ~ '9', POND == '15' ~ '10', 
    POND == '51' ~ '11', POND == '52' ~ '12'))

salamander_pond_knowledge$POND = as.factor(salamander_pond_knowledge$POND)

m1_pond_knowledge = glm.nb(n ~ POND, data = salamander_pond_knowledge, link = "log")
#plot(m1_pond_knowledge)
car::Anova(m1_pond_knowledge)

m1_pond_knowledge_predict <- ggpredict(m1_pond_knowledge, "POND")
plot(m1_pond_knowledge_predict)
plot(m1_pond_knowledge_predict, add.data = TRUE)

pond15 = subset(salamander, POND == "15")
summary(pond15$SEX)
#Female   Male 
# 22       6   =  79% females

pond52 = subset(salamander, POND == "52")
summary(pond52$SEX)
#Female   Male 
# 41       16  =  72% females


#Pond Sex Ratios####
salamander_2021_ratio = salamander_2021 %>% group_by(DATEOFYEAR, POND) %>% count(SEX)
salamander_2021_ratio = salamander_2021_ratio %>% group_by(DATEOFYEAR, POND) %>% complete(SEX, fill = list(n = 0)) 
salamander_2021_ratio = salamander_2021_ratio %>% group_by(DATEOFYEAR, POND) %>% mutate(Proportion = n / sum(n))
salamander_2021_ratio = subset(salamander_2021_ratio, SEX == "Female")
salamander_2021_ratio = subset(salamander_2021_ratio, POND == "10" | POND == "11" | POND == "15" |
         POND == "51" | POND == "52" | POND == "6" | POND == "8" | POND == "13")
salamander_2021_ratio = salamander_2021_ratio %>% arrange(POND, DATEOFYEAR)
salamander_2021_ratio = salamander_2021_ratio %>% group_by(POND) %>% mutate(Sampling_Occasion = row_number())

salamander_2021_ratio$POND <- factor(salamander_2021_ratio$POND,levels = c("52", "51", "15", "13", "11", "10",
                                                                            "8", "6"))
sexratio_heatmap_plot_2021 = 
  ggplot(salamander_2021_ratio, aes(x = Sampling_Occasion, y = POND, fill = Proportion)) +
  geom_raster() +
  scale_x_continuous(limits = c(0,12), breaks = seq(1,11, by = 2), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradientn(colours = paletteer_dynamic("cartography::wine.pal", 20), na.value = "grey50") +
  theme_bw(45) + 
  labs(y = "Pond", x = "Sampling Occasion", fill = "Percent Females") +
  theme(axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        legend.key.size = unit(2, 'cm'))

sexratio_heatmap_plot_2021

ggsave(sexratio_heatmap_plot_2021, filename = "Outputs/sexratio_heatmap_plot_2021.png",  width = 20, height = 12, dpi = 100)


salamander_2022_ratio = salamander_2022 %>% group_by(DATEOFYEAR, POND) %>% count(SEX)
salamander_2022_ratio = salamander_2022_ratio %>% group_by(DATEOFYEAR, POND) %>% complete(SEX, fill = list(n = 0)) 
salamander_2022_ratio = salamander_2022_ratio %>% group_by(DATEOFYEAR, POND) %>% mutate(Proportion = n / sum(n))
salamander_2022_ratio = subset(salamander_2022_ratio, SEX == "Female")
salamander_2022_ratio = subset(salamander_2022_ratio, POND == "10" | POND == "11" | POND == "15" |
                                 POND == "51" | POND == "52" | POND == "6" | POND == "8" | POND == "13")
salamander_2022_ratio = salamander_2022_ratio %>% arrange(POND, DATEOFYEAR)
salamander_2022_ratio = salamander_2022_ratio %>% group_by(POND) %>% mutate(Sampling_Occasion = row_number())


sexratio_heatmap_plot_2022 = 
  ggplot(salamander_2022_ratio, aes(x = Sampling_Occasion, y = POND, fill = Proportion)) +
  geom_raster() +
  #scale_x_continuous(breaks = seq(1990, 2020, by = 5), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradientn(colours = paletteer_dynamic("cartography::wine.pal", 20), na.value = "grey50") +
  theme_bw(35) + 
  labs(y = "Pond", fill = "Percent Females") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0)),
        legend.key.size = unit(2, 'cm'))

sexratio_heatmap_plot_2022


