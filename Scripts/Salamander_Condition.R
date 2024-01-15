
##### Salamander Condition #####

rm(list = ls())

#library(remotes)
#remotes::install_github("Ajfrick/ajfhelpR", force = T)
library(ajfhelpR)
library(lubridate)
library(dplyr)
library(MASS)
library(lme4)
library(lattice)
library(lmerTest)
library(ggplot2)
library(visreg)
library(viridis)
library(sjPlot)


salamander <- read.csv(file = "Data/Salamander_Master.csv", stringsAsFactors = FALSE)
meta_percent_FS_diet <- read.csv(file = "Data/Meta_Percent_FS_Diet.csv", stringsAsFactors = FALSE)
  
salamander$Field1 <- NULL
salamander$WEIGHT <- salamander$WEIGHT * 1000 #convert from g to mg to make on same scale as FS

meta_percent_FS_diet$X <- NULL


#subset data to 2021 and 2022 metas
salamander$DATEOFYEAR = dmy(salamander$DATEOFYEAR)

#align ponds names with with other data sets. Just use pond number rather than L and U for lower and upper
salamander <- salamander %>% 
  mutate(POND = case_when(
    POND == 'U05' ~ '55', POND == 'U02' ~ '52', POND == 'U01' ~ '51', POND == 'L52' ~ '52', POND == 'L38' ~ '38',
    POND == 'L18' ~ '18', POND == 'L16' ~ '16', POND == 'L15' ~ '15', POND == 'L14' ~ '14', POND == 'L13' ~ '13',
    POND == 'L12' ~ '12', POND == 'L11' ~ '11', POND == 'L10' ~ '10', POND == 'L09' ~ '9', POND == 'L08' ~ '8',
    POND == 'L06' ~ '6', POND == 'L05' ~ '5', POND == 'L01A' ~ '1A', POND == 'L01' ~ '1',
    TRUE ~ as.character(POND) #keeps pond numbers that aren't called by case_when
  )) #pond 52 listed as both U02 and L52 in database 

salamander$POND = as.factor(salamander$POND)

#subset dates and ponds
salamander_all = salamander
salamander = subset(salamander, DATEOFYEAR > "2021-01-01" & MORPH == "Meta" & CONDITION != "dead" &
                      !is.na(SVL) & !is.na(WEIGHT))
salamander = subset(salamander, POND == 6 | POND == 8 | POND == 10 | POND == 11 | POND == 13 | POND == 15 | POND == 51 | POND == 52 | POND == 55)

#write.csv(salamander, "Data/Salamander_Master_cleaned.csv")



###look at relationship between salamander body condition and FS density 

##find nearest dates of FS pond sampling to salamander capture date, for each pond

#format
fs_pond <- read.csv(file = "Data/fs_pond.csv", stringsAsFactors = FALSE)
fs_pond$X <- NULL
fs_pond$Date = ymd(fs_pond$Date)
fs_pond$Pond = as.factor(fs_pond$Pond)
names(fs_pond)[names(fs_pond) == 'Date'] <- 'Date_Sample_FS'
names(salamander)[names(salamander) == 'POND'] <- 'Pond'
salamander[,c("Pond")] <- factor(salamander[,c("Pond")], levels=c("6", "8", "10", "11", "13", "15", "51", "52", "55"))

#for loop
temp <- list()

for (i in 1:nrow(salamander)) {       
  fs_pond_subset <- subset(fs_pond, Pond==(salamander[i, 6]))
  temp[[i]] <- ajfhelpR::date_near(fs_pond_subset$Date_Sample_FS, salamander[i, 3], sidepref = 'l')
}

temp


dates = temp %>% purrr::reduce(c) #make list into vector
rm(temp)

salamander$Date_Sample_FS <- dates
salamander$Date_Sample_FS = ymd(salamander$Date_Sample_FS)

abs(salamander$Date_Sample_FS - salamander$DATEOFYEAR) #time difference
range(abs(salamander$Date_Sample_FS - salamander$DATEOFYEAR))
mean(abs(salamander$Date_Sample_FS - salamander$DATEOFYEAR))

salamander = merge(salamander, fs_pond, by = c("Date_Sample_FS", "Pond"), all.x = T)
salamander = salamander %>% relocate(c(Date_Sample_FS), .after = COMMENT)
salamander = salamander %>% relocate(Pond, .after = RECAPTURE)

##model

#4/4/22 meeting with Zach about calculating condition metric. He advocated for not attempting to calculate a body condition 
#value for each salamander capture and then using those values in further regressions. Instead, he suggested that you use one
#regression to understand effect of FS densities on salamander capture weight. In other words, when controlling for SVL and
#variation due to individual ID & Pond ID, what is the effect of FS densities on weight?

#if someone insists on having a condition value for each salamander capture, you could take the residuals from the model 
#lmer(log(WEIGHT) ~ log(SVL) + (1|POND). But you can't have the INDIV_ID random effect because then the residuals quantify
#the variation in individuals' weight (grouped by individual ID), not variation around a larger group average. Not sure how
#you would control for non-independence of repeated individuals.

nrow(salamander) / length(unique(salamander$INDIV_ID)) #an average of 2.83 obs per salamander individual

meta_condition = lmer(log(WEIGHT) ~ log_FS_Biomass_Pond_m3 + log(SVL) + (1|Pond) + (1|INDIV_ID), data = salamander)

plot(meta_condition, type=c("p","smooth"), col.line=1) #residuals look fairly homoscedastic. Apparent funnel shape may be due
#to few data points on left side, rather than true funneling of residuals 
qqmath(meta_condition) #residuals look normal

summary(meta_condition)

#use Kenward-Roger approximations for degrees of freedom in order to calculate P values (Type III SS), 
#using the lmerTest R package. Need to load lmerTest before running model.

anova(meta_condition, type = 2, ddf = "Kenward-Roger")

#great resource on visualizing multiple regressions:
#https://www.zoology.ubc.ca/~schluter/R/Model.html
#visreg FAQs: https://rdrr.io/cran/visreg/man/visreg-faq.html


visreg(meta_condition, xvar = "log_FS_Biomass_Pond_m3", type = "conditional",
       whitespace = 0.4, points.par = list(cex = 1.1, col = "red"),
       ylab="log(Salamander Weight)", xlab = "log(FS Biomass Pond)")

visreg(meta_condition, xvar = "SVL", type = "conditional", xtrans = log,
       whitespace = 0.4, points.par = list(cex = 1.1, col = "red"),
       ylab="log(Salamander mass)", xlab = "Salamander SVL")


##subset salamander dataset to just high FS ponds
salamander_highFSponds = subset(salamander, Pond == 13 | Pond == 15 | Pond == 51 | Pond == 52)

m_salamander_highFSponds = lmer(log(WEIGHT) ~ log_FS_Biomass_Pond_m3 + log(SVL) + (1|Pond) + (1|INDIV_ID), data = salamander_highFSponds)

plot(m_salamander_highFSponds, type=c("p","smooth"), col.line=1)
qqmath(m_salamander_highFSponds)

summary(m_salamander_highFSponds)

anova(m_salamander_highFSponds, type = 2, ddf = "Kenward-Roger")

visreg(m_salamander_highFSponds, xvar = "log_FS_Biomass_Pond_m3", type = "conditional",
       whitespace = 0.4, points.par = list(cex = 1.1, col = "red"),
       ylab="log(Salamander Weight)", xlab = "log(FS Biomass Pond)")

visreg(m_salamander_highFSponds, xvar = "SVL", type = "conditional",
       whitespace = 0.4, points.par = list(cex = 1.1, col = "red"),
       ylab="log(Salamander mass)", xlab = "Salamander SVL")


### look at meta weight ~ FS biomass in gut ###
fs_meta_gut <- read.csv(file = "Data/fs_meta_gut.csv", stringsAsFactors = FALSE)
fs_meta_gut$X <- NULL
names(fs_meta_gut)[names(fs_meta_gut) == 'Date_Sample_GF'] <- 'DATEOFYEAR'
names(fs_meta_gut)[names(fs_meta_gut) == 'Salamander_ID'] <- 'INDIV_ID'

fs_meta_gut$INDIV_ID = as.factor(fs_meta_gut$INDIV_ID)
fs_meta_gut$Sample_ID = as.factor(fs_meta_gut$Sample_ID)
fs_meta_gut$DATEOFYEAR = ymd(fs_meta_gut$DATEOFYEAR)

salamander_FSingut = merge(x = fs_meta_gut, y = salamander[,c("INDIV_ID", "DATEOFYEAR", "Pond", "COHORT", "SEX", "SVL", "WEIGHT", "log_FS_Biomass_Pond_m3")],
                   by = c("INDIV_ID", "DATEOFYEAR", "Pond"), all.x = T)


#model
names(salamander_FSingut)[names(salamander_FSingut) == 'log_Boat_Net'] <- 'log_FS_Biomass_Gut'

salamander_FSingut$Condition = (salamander_FSingut$WEIGHT / salamander_FSingut$SVL) / 1000
salamander_FSingut$SEX = as.factor(salamander_FSingut$SEX)
salamander_FSingut$SEX = car::recode(salamander_FSingut$SEX,"'Female?'='Female'")
salamander_FSingut = subset(salamander_FSingut, !is.na(SEX))

salamander_FSingut_positive = subset(salamander_FSingut, log_FS_Biomass_Gut > 0)


condition_sex_plot = 
  ggplot(salamander_FSingut, aes(x = log_FS_Biomass_Gut, y = Condition)) +
  geom_point(aes(colour = SEX), size = 6, alpha = 0.8) +
  geom_smooth(data = salamander_FSingut_positive, method="lm", formula= (y ~ exp(x)), se=TRUE, linetype = 1, color = "black", lwd = 2) +
  scale_colour_manual(values = c("orange", "steelblue")) + 
  scale_x_continuous(limits = c(log(1), log(1000)), breaks=c(log(1),log(10),log(100),log(1000)), labels=c("0", "10", "100", "1,000")) +
  labs(colour = "Sex", x = "Fairy Shrimp Biomass in Stomach (mg)", y = "Body Condition") +
  theme_bw(35)

condition_sex_plot

#ggsave(condition_sex_plot, filename = "Outputs/condition_sex_plot.png",  width = 16, height = 12, dpi = 100)

salamander_FSingut = merge(salamander_FSingut, meta_percent_FS_diet[,c("Sample_ID_Year", "Percent_Diet_FS")], by = "Sample_ID_Year", all.x = T)
meta_high_condition = subset(salamander_FSingut, Condition > 0.399)
mean(meta_high_condition$Percent_Diet_FS)
sd(meta_high_condition$Percent_Diet_FS)

#stats

m1_condition = lmer(Condition ~ exp(log_FS_Biomass_Gut)*SEX + (1|Pond), data = salamander_FSingut_positive)
car::Anova(m1_condition)
summary(m1_condition)
tab_model(m1_condition) #to get R-squared for a LMM


m_salamander_FSingut = lmer(log(WEIGHT) ~ log_FS_Biomass_Gut + log(SVL) + (1|Pond), data = salamander_FSingut)

plot(m_salamander_FSingut, type=c("p","smooth"), col.line=1)
qqmath(m_salamander_FSingut)

summary(m_salamander_FSingut)

anova(m_salamander_FSingut, type = 2, ddf = "Kenward-Roger")

visreg(m_salamander_FSingut, xvar = "log_FS_Biomass_Gut", type = "conditional",
       whitespace = 0.4, points.par = list(cex = 1.1, col = "red"),
       ylab="log(Salamander Weight)", xlab = "log(Fairy Shrimp Biomass in Gut)")

visreg(m_salamander_FSingut, xvar = "SVL", type = "conditional", xtrans = log,
       whitespace = 0.4, points.par = list(cex = 1.1, col = "red"),
       ylab="log(Salamander mass)", xlab = "Salamander SVL")


#### SVL residuals vs. condition ####
meta_svl_residuals <- read.csv(file = "Data/meta_growth_curve_residuals.csv", stringsAsFactors = FALSE)

svl_v_condition = merge(salamander_FSingut, unique(meta_svl_residuals[,c("Year", "INDIV_ID", "SVL_residual_mean")]), by = c("Year", "INDIV_ID"), all.x = T)
plot(SVL_residual_mean ~ Condition, data = svl_v_condition)


#### Paedo Morph Test  ##### 

#test whether there is evidence for best of a bad lot and padeomorph advantage morphs
#if so, should see bimodal distribution 

#paedos
salamander_full_dataset = paedos

paedos = subset(paedos, MORPH == "Paedo" & INDIV_ID != "" & CONDITION != "dead" & !is.na(SVL) & SVL < 200)

paedos_SVL_max = paedos %>% group_by(INDIV_ID) %>% slice(which.max(SVL))

paedos_SVL_max = subset(paedos_SVL_max, SVL > 59 & SVL < 130) #remove extreme values

hist(paedos_SVL_max$SVL) #no evidence of bimodal distribution
mean(paedos_SVL_max$SVL)

#metas
metas_SVL_max = salamander_all %>% group_by(INDIV_ID) %>% slice(which.max(SVL))

hist(metas_SVL_max$SVL)
mean(metas_SVL_max$SVL)


#compare
salamander_full_dataset = subset(salamander_full_dataset, (MORPH == "Meta" | MORPH == "Paedo") &  
                                   (SEX == "Male" | SEX == "Female") &
                                   INDIV_ID != "" & CONDITION != "dead" & !is.na(SVL) & SVL < 200)

all_SVL_max = salamander_full_dataset %>% group_by(INDIV_ID) %>% slice(which.max(SVL))


morph_max_size = 
  all_SVL_max %>%
  mutate(MORPH = factor(MORPH, levels=c("Paedo", "Meta"))) %>%
  ggplot(aes(fill = MORPH, y=SVL, x=MORPH)) + 
  geom_jitter(color="black", size=0.5, alpha=0.7) +
  geom_violin(alpha=0.5, lwd = 1) +
  geom_boxplot(alpha=0.2, width=0.2, color="black", lwd = 1, outlier.shape = NA) +
  scale_fill_viridis(discrete=T, name="") +
  xlab("") +
  ylab("") +
  ylim(0, 200) +
  theme(axis.text = element_text(size = 40),
        #legend.key.size = unit(5, 'cm'),
        #legend.position = "none",
        axis.ticks = element_line(size = 4.0),
        axis.ticks.length=unit(.300, "cm"),
        axis.text.x = element_blank())

morph_max_size

#what happens if we remove two small outlier paedo points
all_SVL_max_v2 = subset(all_SVL_max, SVL > 50)

morph_max_size_v2 = 
  all_SVL_max_v2 %>%
  mutate(MORPH = factor(MORPH, levels=c("Paedo", "Meta"))) %>%
  ggplot(aes(fill = SEX, y=SVL, x=MORPH)) + 
  geom_jitter(color="black", size=0.5, alpha=0.7) +
  geom_violin(alpha=0.5, lwd = 1) +
  geom_boxplot(alpha=0.2, width=0.2, color="black", lwd = 1, outlier.shape = NA) +
  scale_fill_viridis(discrete=T, name="") +
  xlab("") +
  ylab("") +
  theme(axis.text = element_text(size = 40),
        #legend.key.size = unit(5, 'cm'),
        #legend.position = "none",
        axis.ticks = element_line(size = 4.0),
        axis.ticks.length=unit(.300, "cm"))

morph_max_size_v2


biggest = subset(all_SVL_max_v2, SVL > 110)

biggest_size = 
  biggest %>%
  mutate(MORPH = factor(MORPH, levels=c("Paedo", "Meta"))) %>%
  ggplot(aes(fill = MORPH, y=SVL, x=MORPH)) + 
  geom_jitter(color="black", size=0.5, alpha=0.7) +
  geom_violin(alpha=0.5, lwd = 1) +
  geom_boxplot(alpha=0.2, width=0.2, color="black", lwd = 1, outlier.shape = NA) +
  scale_fill_viridis(discrete=T, name="") +
  xlab("") +
  ylab("") +
  theme(axis.text = element_text(size = 40),
        #legend.key.size = unit(5, 'cm'),
        #legend.position = "none",
        axis.ticks = element_line(size = 4.0),
        axis.ticks.length=unit(.300, "cm"))

biggest_size



### test variation in SVL among cohort classes ###
salamander_FSingut$COHORT[is.na(salamander_FSingut$COHORT)] <- "Unknown"
salamander_FSingut$COHORT = as.factor(salamander_FSingut$COHORT)

salamander_FSingut$Cohort_Class = with(salamander_FSingut,
                                      ifelse(COHORT %in% c("1996", "1998", "1999", "2000", "2001"), "Oldest", 
                                      ifelse(COHORT %in% c("2007", "2008", "2009", "2010", "2011"), "Middle", 
                                      ifelse(COHORT %in% c("2012", "2013", "2014", "2015", "2020"), "Youngest",
                                      ifelse(COHORT %in% c("Unknown"), "Unknown",
                                      "ERROR")))))

salamander_FSingut$Cohort_Class = as.factor(salamander_FSingut$Cohort_Class)
salamander_FSingut[,c("Cohort_Class")] <- factor(salamander_FSingut[,c("Cohort_Class")], levels=c("Youngest", "Middle", "Oldest", "Unknown"))
salamander_FSingut = salamander_FSingut %>% relocate(Cohort_Class, .after = COHORT)

salamander_FSingut %>% group_by(Cohort_Class) %>% summarise(Mean_SVL = mean(SVL))
salamander_FSingut %>% group_by(Cohort_Class) %>% summarise(SD_SVL = sd(SVL))



