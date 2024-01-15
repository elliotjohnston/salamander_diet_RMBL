
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
paedomorphs = droplevels(paedomorphs)
metamorphs = subset(salamander, MORPH == "Meta")
metamorphs = droplevels(metamorphs)


##### Decision Tree #####

#(1) metamorphs vs. paedomorphs
nlevels(paedomorphs$INDIV_ID) #5928 paedos (76%)
nlevels(metamorphs$INDIV_ID) #1827 metas (24%)


#(2) don't return to pond system
metamorphs$Year = year(metamorphs$DATEOFYEAR)
metamorphs = metamorphs %>% relocate(Year, .before = DATEOFYEAR)
metamorphs$INDIV_ID_Year = interaction(metamorphs$INDIV_ID, metamorphs$Year, sep="_")
metamorphs = metamorphs %>% relocate(INDIV_ID_Year, .after = Year)
metamorphs$INDIV_ID_Year = as.factor(metamorphs$INDIV_ID_Year)
metamorphs = droplevels(metamorphs)

nlevels(metamorphs$INDIV_ID) #1827 unique metas
nrow(metamorphs) #6915 total captures
range(metamorphs$Year) #1990-2022
nlevels(metamorphs$INDIV_ID_Year) #4145 meta-years

metamorphs_recap = 
  metamorphs %>% 
  group_by(INDIV_ID) %>%
  filter(n_distinct(Year)>1)

metamorphs_recap_yrs = metamorphs_recap[,c("INDIV_ID", "Year", "POND")]

metamorphs_recap_yrs_return = 
  metamorphs_recap_yrs %>% 
  group_by(INDIV_ID) %>% 
  mutate(Years_To_Return = Year - lead(Year))

metamorphs_recap_yrs_return = subset(metamorphs_recap_yrs_return, Years_To_Return > 0)

mean(metamorphs_recap_yrs_return$Years_To_Return, na.rm = TRUE) #3.2 years
median(metamorphs_recap_yrs_return$Years_To_Return, na.rm = TRUE) #2 year
max(metamorphs_recap_yrs_return$Years_To_Return, na.rm = TRUE) #19 years

hist(metamorphs_recap_yrs_return$Years_To_Return,
     main = "",
     xlab = "Years Between Metamorph Captures",
     ylab = "Frequency",
     xlim = range(0, 18),
     ylim = c(0, 1500))

metamorphs_recap_yrs_return = droplevels(metamorphs_recap_yrs_return)
nlevels(metamorphs_recap_yrs_return$INDIV_ID) #729 metas captured in at least two years
nrow(metamorphs_recap_yrs_return) #2200 capture intervals


#(3) return to one pond in a given year
metamorphs_mvmt = metamorphs[,c("INDIV_ID", "DATEOFYEAR", "Year", "INDIV_ID_Year", "POND")]
metamorphs_mvmt = droplevels(metamorphs_mvmt)
levels(metamorphs_mvmt$POND)

metamorphs_mvmt$Hydroperiod <- ifelse(metamorphs_mvmt$POND == "1" | 
                                        metamorphs_mvmt$POND == "12" |
                                        metamorphs_mvmt$POND == "5" |
                                        metamorphs_mvmt$POND == "9" |
                                        metamorphs_mvmt$POND == "14" |
                                        metamorphs_mvmt$POND == "3" |
                                        metamorphs_mvmt$POND == "18", "Permanent", "Non-permanent")

metamorphs_mvmt$Hydroperiod = as.factor(metamorphs_mvmt$Hydroperiod)
metamorphs_mvmt$Year = as.factor(metamorphs_mvmt$Year)

metamorphs_mvmt_one_pond = 
  metamorphs_mvmt %>% 
  group_by(INDIV_ID, Year) %>%
  filter(n_distinct(POND)==1)

metamorphs_mvmt_one_pond = droplevels(metamorphs_mvmt_one_pond)
nlevels(metamorphs_mvmt_one_pond$INDIV_ID_Year) #3747 meta-years

temp = metamorphs_mvmt_one_pond[,c("INDIV_ID_Year", "Hydroperiod")]
temp = unique(temp)
summary(temp$Hydroperiod)
#2215/3773 (59%) of one pond meta-years were in non-permanent ponds
#1532/3773 (41%) of one pond meta-years were in permanent ponds

##what is the average percentage of an individuals movements over its lifetime to one pond (vs. visiting multiple)
#analyze only for metas that have been captured in at least two years

num_ponds_visited = 
  metamorphs_recap_yrs %>% 
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

#(4) return to one nonpermanent pond
metamorphs_mvmt_one_pond_np = subset(metamorphs_mvmt_one_pond, Hydroperiod == "Non-permanent")
metamorphs_mvmt_one_pond_np$Hydroperiod_NP <- ifelse(metamorphs_mvmt_one_pond_np$POND == "13" | 
                                                       metamorphs_mvmt_one_pond_np$POND == "15" |
                                                       metamorphs_mvmt_one_pond_np$POND == "16" |
                                                       metamorphs_mvmt_one_pond_np$POND == "17" |
                                                       metamorphs_mvmt_one_pond_np$POND == "27" |
                                                       metamorphs_mvmt_one_pond_np$POND == "37" |
                                                       metamorphs_mvmt_one_pond_np$POND == "38" |
                                                       metamorphs_mvmt_one_pond_np$POND == "41" | 
                                                       metamorphs_mvmt_one_pond_np$POND == "7", "Temporary", "Semi-permanent")

metamorphs_mvmt_one_pond_np$Hydroperiod_NP = as.factor(metamorphs_mvmt_one_pond_np$Hydroperiod_NP)
temp_np = metamorphs_mvmt_one_pond_np[,c("INDIV_ID_Year", "Hydroperiod_NP")]
temp_np = unique(temp_np)
summary(temp_np$Hydroperiod_NP)
#2042/2215 (92%) of one pond non-permanent meta-years were in semi-permanent ponds
#173/2215 (8%) of one pond non-permanent meta-years were in temporary ponds



#(5) return to multiple ponds in a given year

metamorphs_mvmt_mult_pond = 
  metamorphs_mvmt %>% 
  group_by(INDIV_ID, Year) %>%
  filter(n_distinct(POND)>1)

metamorphs_mvmt_mult_pond = unique(metamorphs_mvmt_mult_pond)
metamorphs_mvmt_mult_pond = droplevels(metamorphs_mvmt_mult_pond)
nlevels(metamorphs_mvmt_mult_pond$INDIV_ID_Year) #398 meta-years

mult_2021 = subset(metamorphs_mvmt_mult_pond, Year == "2021")
mult_2021 = mult_2021[,c("INDIV_ID", "Year", "POND", "Hydroperiod")]
mult_2021 = unique(mult_2021)
mult_2022 = subset(metamorphs_mvmt_mult_pond, Year == "2022")
mult_2022 = mult_2022[,c("INDIV_ID", "Year", "POND", "Hydroperiod")]
mult_2022 = unique(mult_2022)

##how many of these are movements utilize permanent and non-permanent ponds?
metamorphs_mvmt_mult_pond_hydro = 
  metamorphs_mvmt_mult_pond %>% 
  group_by(INDIV_ID, Year) %>%
  filter(n_distinct(Hydroperiod)>1)

num_mvmts_hydro = unique(metamorphs_mvmt_mult_pond_hydro[,c("INDIV_ID", "Year", "INDIV_ID_Year")]) 
num_mvmts_hydro = droplevels(num_mvmts_hydro)
nlevels(num_mvmts_hydro$INDIV_ID_Year) #282 meta-years


#### Pond Immigration Rates ####

### landscape factors
POND = c(6,8,10,11,13,15,51,52)
Hydroperiod = c(rep("Semipermanent", 4), rep("Temporary", 2), rep("Semipermanent", 2))
Size = c(91,414,139,96,22,26,124,346)
Dist_Avg = c(107.5, 79.2, 101.1, 107.8, 98.1, 141.6, 205.4, 207.5)
Elevation_Avg = c(17.5, 15.3, 16.0, 16.3, 16.3, 22, 62.7, 69.8)

landscape = as.data.frame(cbind(POND, Hydroperiod, Size, Dist_Avg, Elevation_Avg))
landscape[,1:2] <- lapply(landscape[,1:2], as.factor)
landscape[,3:5] <- lapply(landscape[,3:5], as.numeric)

### pond quality factors

#fairy shrimp densities
POND_FS = c(6,8,10,11,13,15,51,52)
Mean_Historical_FS_Rank = c(7,4,6,4.2,2,2.2,4.7,3)

FS_Historical = as.data.frame(cbind(POND_FS, Mean_Historical_FS_Rank))
colnames(FS_Historical) <- c("POND", "Mean_Historical_FS_Rank")

FS_Historical$Mean_Historical_FS_Rank = ordered(FS_Historical$Mean_Historical_FS_Rank, levels = c(2,2.2,3,4,4.2,4.7,6,7))

#merge landscape and pond quality factors
meta_mvmt_factors = merge(landscape, FS_Historical, all.x = TRUE)

#meta_mvmt_factors$Mean_Historical_FS_Rank = as.factor(meta_mvmt_factors$Mean_Historical_FS_Rank)
#meta_mvmt_factors$Mean_Historical_FS_Rank = as.numeric(meta_mvmt_factors$Mean_Historical_FS_Rank)


#check correlations
corr_var_nonperms = meta_mvmt_factors[,c("Size", "Dist_Avg", "Elevation_Avg")]

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


##calculate pond visitation
#fill in missing cohort values for individuals that have a known cohort recorded at some point
metamorphs = metamorphs %>% group_by(INDIV_ID) %>% fill(COHORT)

#subset to include only metas with age and SVL info. Even though hurts sample size a lot (all due to lack
#of cohort data), need age and svl data to relate individual info to population growth curve
metamorphs = subset(metamorphs, !is.na(COHORT) & !is.na(SVL))
metamorphs = droplevels(metamorphs)

#need to subset movement dataframe to eliminate consecutive captures of an individual in the same
#pond in the same year
metamorphs = merge(metamorphs, unique(metamorphs_mvmt[,c("POND", "Hydroperiod")]), all.x = T)
metamorphs = metamorphs %>% relocate(Hydroperiod, .after = POND)
metamorphs = unique(metamorphs)

metamorphs_immigration = subset(metamorphs, POND == "10" | POND == "11" | POND == "15" |
                              POND == "51" | POND == "52" | POND == "6" | POND == "8" |  
                              POND == "13")

metamorphs_immigration = metamorphs_immigration %>% distinct(INDIV_ID_Year, POND, .keep_all = T)

#summarize
metamorphs_immigration_sum = metamorphs_immigration %>% group_by(POND, Year) %>% mutate(Unique_Meta_Visits = length(POND))
metamorphs_immigration_sum = metamorphs_immigration_sum %>% group_by(POND, Year) %>% mutate(Total_Surveys = n_distinct(DATEOFYEAR))
metamorphs_immigration_sum = unique(metamorphs_immigration_sum[,c("POND", "Year", "Unique_Meta_Visits", "Total_Surveys")])
metamorphs_immigration_sum$Unique_Metas_Per_Survey = metamorphs_immigration_sum$Unique_Meta_Visits / metamorphs_immigration_sum$Total_Surveys

metamorphs_immigration_sum = droplevels(metamorphs_immigration_sum)
metamorphs_immigration_sum$POND <- factor(metamorphs_immigration_sum$POND,levels = c("15", "13",
                                                                             "52", "51", "11", "10",
                                                                             "8", "6"))

metamorphs_immigration_factors = merge(metamorphs_immigration_sum, meta_mvmt_factors, by = "POND", all.x = T)



#model
m1_fs_rank = glmer(Unique_Meta_Visits ~ Mean_Historical_FS_Rank + (1|Year), 
                  offset = log(Total_Surveys), family = poisson, data = metamorphs_immigration_factors)

m2_landscape = glmer(Unique_Meta_Visits ~ scale(Size) + scale(Dist_Avg) + (1|Year), 
                      offset = log(Total_Surveys), family = poisson, data = metamorphs_immigration_factors)

m3_hydroperiod = glmer(Unique_Meta_Visits ~ Hydroperiod + (1|Year), 
                     offset = log(Total_Surveys), family = poisson, data = metamorphs_immigration_factors)


AIC(m1_fs_rank, m2_landscape, m3_hydroperiod)
car::Anova(m1_fs_rank)
car::Anova(m2_landscape)
car::Anova(m3_hydroperiod)

m3_fs_landscape = glmer(Unique_Meta_Visits ~ scale(Size) + Mean_Historical_FS_Rank + (1|Year), 
                     offset = log(Total_Surveys), family = poisson, data = metamorphs_immigration_factors)
#rank deficient 

m4_landscape_size_only = glmer(Unique_Meta_Visits ~ scale(Size) + (1|Year), 
                        offset = log(Total_Surveys), family = poisson, data = metamorphs_immigration_factors)

AIC(m1_fs_rank, m4_landscape_size_only, m3_hydroperiod)
car::Anova(m4_landscape_size_only)



m1_fs_rank_predict <- ggpredict(m1_fs_rank, terms = "Mean_Historical_FS_Rank")
plot(m1_fs_rank_predict)
plot(m1_fs_rank_predict, rawdata = T)

m2_landscape_predict <- ggpredict(m2_landscape, terms = "Size")
plot(m2_landscape_predict)
plot(m2_landscape_predict, rawdata = T)

m3_hydroperiod_predict <- ggpredict(m3_hydroperiod, terms = "Hydroperiod")
plot(m3_hydroperiod_predict)
plot(m3_hydroperiod_predict, rawdata = T)


m_v_f_svl = lm(SVL ~ SEX, data = salamander)
summary(m_v_f_svl)

m_v_f_svl_plot = ggplot(data = salamander, aes(x = SEX, y = SVL)) +
  geom_boxplot(lwd = 2, outlier.size = 3) +
  theme_bw(45)

#ggsave(m_v_f_svl_plot, filename = "Outputs/m_v_f_svl_plot.png",  width = 12, height = 12, dpi = 100)


m_v_f_weight = lm(WEIGHT ~ SEX, data = salamander)
summary(m_v_f_weight)

m_v_f_weight_plot = ggplot(data = salamander, aes(x = SEX, y = WEIGHT)) +
  geom_boxplot(lwd = 2, outlier.size = 3) +
  theme_bw(45)

#ggsave(m_v_f_weight_plot, filename = "Outputs/m_v_f_weight_plot.png",  width = 12, height = 12, dpi = 100)


test = salamander
test$CONDITION = test$WEIGHT / test$SVL

m_v_f_condition = lm(CONDITION ~ SEX, data = test)
summary(m_v_f_condition)

m_v_f_condition_plot = ggplot(data = test, aes(x = SEX, y = CONDITION)) +
  geom_boxplot(lwd = 2, outlier.size = 3) +
  theme_bw(45)

#ggsave(m_v_f_condition_plot, filename = "Outputs/m_v_f_condition_plot.png",  width = 12, height = 12, dpi = 100)


#### SVL Resid. + Immigration ####

metamorphs_growth = metamorphs_immigration

metamorphs_growth$AGE = metamorphs_growth$Year - metamorphs_growth$COHORT
metamorphs_growth = subset(metamorphs_growth, AGE > 0)
metamorphs_growth_male = subset(metamorphs_growth, SEX == "Male")
metamorphs_growth_female = subset(metamorphs_growth, SEX == "Female")

#### (i) all metas #### 
plot(SVL ~ AGE, data = metamorphs_growth)

fishmethods::growth(size=metamorphs_growth$SVL,age=metamorphs_growth$AGE,
                    error=1, Sinf=101,K=0.1,t0=-12) #good starting values

#resource for below code: https://danstich.github.io/we-r-nycafs/fishStats.html

#define von Bertalanffy growth function
vbmod <- SVL ~ Linf * (1 - exp(-K * (AGE - t0)))


#fit the von Bertalanffy growth function using nonlinear least squares (nls) optimization
growth_mod <- nls(vbmod, data = metamorphs_growth, start = list(Linf = 101, K = 0.1, t0 = -12))
summary(growth_mod)

#predict
growth_predict = predict(growth_mod)

##bootstrap 95% CI
#get the desired growth function from a list of those that are available in FSA
vbO <- vbFuns("typical")

#fit the model to the data using nls, like we did before
vb_fit <- nls(SVL~vbO(AGE,Linf,K, t0), data = metamorphs_growth, start = list(Linf = 101, K = 0.1, t0 = -12))

#now, bootstrap the model fitting process
boot_fit <- nlsBoot(vb_fit, niter = 9999)

#predict length at age from the model (t is age). Here, we tell R to predict length at 
#each unique age in our original data, and calculate some bootstrapped confidence
#intervals
boot_preds <- data.frame(predict(boot_fit, vbO, t = sort(unique(metamorphs_growth$AGE))))
names(boot_preds) <- c("AGE", "fit", "lwr", "upr")

growth_preds <- merge(metamorphs_growth, boot_preds, by = "AGE")

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


metamorphs_SVL_immigration_factors = merge(growth_preds, meta_mvmt_factors, by = "POND", all.x = T)

#model
m1_SVL_fs_rank = lmer(SVL_residual ~ Mean_Historical_FS_Rank + (1|Year) + (1|INDIV_ID), 
                  data = metamorphs_SVL_immigration_factors)

m2_SVL_landscape = lmer(SVL_residual ~ scale(Size) + scale(Dist_Avg) + (1|Year) + (1|INDIV_ID), 
                    data = metamorphs_SVL_immigration_factors)

m2_SVL_landscape_size_only = lmer(SVL_residual ~ scale(Size) + (1|Year) + (1|INDIV_ID), 
                        data = metamorphs_SVL_immigration_factors)

m3_SVL_hydroperiod = lmer(SVL_residual ~ Hydroperiod + (1|Year) + (1|INDIV_ID), 
                      data = metamorphs_SVL_immigration_factors)


AIC(m1_SVL_fs_rank, m3_SVL_hydroperiod, m2_SVL_landscape, m2_SVL_landscape_size_only)
car::Anova(m1_SVL_fs_rank)
car::Anova(m2_SVL_landscape)
car::Anova(m3_SVL_hydroperiod)


m1_SVL_fs_rank_predict <- ggpredict(m1_SVL_fs_rank, terms = "Mean_Historical_FS_Rank")
plot(m1_SVL_fs_rank_predict)
plot(m1_SVL_fs_rank_predict, rawdata = T)

m2_SVL_landscape_predict <- ggpredict(m2_SVL_landscape_size_only, terms = "Size")
plot(m2_SVL_landscape_predict)
plot(m2_SVL_landscape_predict, rawdata = T)

m3_SVL_hydroperiod_predict <- ggpredict(m3_SVL_hydroperiod, terms = "Hydroperiod")
plot(m3_SVL_hydroperiod_predict)
plot(m3_SVL_hydroperiod_predict, rawdata = T)


#### (ii) females #### 
plot(SVL ~ AGE, data = metamorphs_growth_female)

fishmethods::growth(size=metamorphs_growth_female$SVL,age=metamorphs_growth_female$AGE,
                    error=1, Sinf=101,K=0.1,t0=-12) #good starting values

#fit the von Bertalanffy growth function using nonlinear least squares (nls) optimization
growth_mod_female <- nls(vbmod, data = metamorphs_growth_female, start = list(Linf = 101, K = 0.1, t0 = -12))
summary(growth_mod_female)

#predict
growth_predict_female = predict(growth_mod_female)

##bootstrap 95% CI

#fit the model to the data using nls, like we did before
vb_fit_female <- nls(SVL~vbO(AGE,Linf,K, t0), data = metamorphs_growth_female, start = list(Linf = 101, K = 0.1, t0 = -12))

#now, bootstrap the model fitting process
boot_fit_female <- nlsBoot(vb_fit_female, niter = 9999)

#predict length at age from the model (t is age). Here, we tell R to predict length at 
#each unique age in our original data, and calculate some bootstrapped confidence
#intervals
boot_preds_female <- data.frame(predict(boot_fit_female, vbO, t = sort(unique(metamorphs_growth_female$AGE))))
names(boot_preds_female) <- c("AGE", "fit", "lwr", "upr")

growth_preds_female <- merge(metamorphs_growth_female, boot_preds_female, by = "AGE")

#plot
growth_curve_plot_female = 
  ggplot(growth_preds_female, aes(x = AGE, y = SVL)) +
  geom_jitter(width = 0.1, alpha = 0.15, size = 5) +
  geom_line(aes(y = fit), linewidth = 2) +
  geom_ribbon(
    aes(x = AGE, ymin = lwr, ymax = upr, color = NULL), alpha = 0.3) +
  scale_y_continuous(limits = c(59,120)) +
  xlab("Age (years)") +
  ylab("Snout-Vent Length (mm)") +
  theme_bw(40)

growth_curve_plot_female

#ggsave(growth_curve_plot_female, filename = "Outputs/growth_curve_plot_female.png",  width = 14, height = 12, dpi = 100)


#calculate residuals 
growth_preds_female = growth_preds_female[,c("DATEOFYEAR", "Year", "INDIV_ID", "INDIV_ID_Year",
                               "MORPH", "POND", "COHORT", "SEX", "WEIGHT", "SVL",
                               "fit", "lwr", "upr")]

names(growth_preds_female)[names(growth_preds_female) == 'fit'] <- 'SVL_fit'
names(growth_preds_female)[names(growth_preds_female) == 'lwr'] <- 'SVL_fit_lwr'
names(growth_preds_female)[names(growth_preds_female) == 'upr'] <- 'SVL_fit_upr'

growth_preds_female$SVL_residual = growth_preds_female$SVL - growth_preds_female$SVL_fit 
hist(growth_preds_female$SVL_residual)

growth_preds_female = droplevels(growth_preds_female)
nlevels(growth_preds_female$INDIV_ID)

num_captures = growth_preds_female %>% group_by(INDIV_ID) %>% summarise(N = n())
range(num_captures$N)
mean(num_captures$N)
sd(num_captures$N)

#create df for SVL residual data for metamorphs that were captured at least twice. In other words, eliminate
#metas that were captured in one pond in one year only. These individuals don't have the opportunity to learn
#about pond resources, and they're more likely to end up in ponds based on landscape factors (e.g., biggest pond)
growth_preds_female_multiple_obs = growth_preds_female %>% 
  add_count(INDIV_ID, name = "Num_Captures") %>% 
  filter(Num_Captures > 1)

#add in mean SVL residual data for each meta over their lifetime
growth_preds_female_avg = growth_preds_female_multiple_obs %>% group_by(INDIV_ID) %>% summarise(SVL_residual_mean = mean(SVL_residual))
hist(growth_preds_female_avg$SVL_residual_mean)
growth_preds_female = inner_join(growth_preds_female, growth_preds_female_avg, by = "INDIV_ID")  
growth_preds_female = droplevels(growth_preds_female)

metamorphs_SVL_immigration_factors_female = merge(growth_preds_female, meta_mvmt_factors,
                                                  by = "POND", all.x = T)

#model
m1_SVL_fs_rank_F = lmer(SVL_residual ~ Mean_Historical_FS_Rank + (1|Year) + (1|INDIV_ID), 
                      data = metamorphs_SVL_immigration_factors_female)

m2_SVL_landscape_F = lmer(SVL_residual ~ scale(Size) + scale(Dist_Avg) + (1|Year) + (1|INDIV_ID), 
                        data = metamorphs_SVL_immigration_factors_female)

m3_SVL_hydroperiod_F = lmer(SVL_residual ~ Hydroperiod + (1|Year) + (1|INDIV_ID), 
                          data = metamorphs_SVL_immigration_factors_female)


AIC(m1_SVL_fs_rank_F, m3_SVL_hydroperiod_F, m2_SVL_landscape_F)
car::Anova(m1_SVL_fs_rank_F)
car::Anova(m2_SVL_landscape_F)
car::Anova(m3_SVL_hydroperiod_F)


m1_SVL_fs_rank_predict_F <- ggpredict(m1_SVL_fs_rank_F, terms = "Mean_Historical_FS_Rank")
plot(m1_SVL_fs_rank_predict_F)
plot(m1_SVL_fs_rank_predict_F, rawdata = T)

m2_SVL_landscape_predict_F <- ggpredict(m2_SVL_landscape_F, terms = "Size")
plot(m2_SVL_landscape_predict_F)
plot(m2_SVL_landscape_predict_F, rawdata = T)

m3_SVL_hydroperiod_predict_F <- ggpredict(m3_SVL_hydroperiod_F, terms = "Hydroperiod")
plot(m3_SVL_hydroperiod_predict_F)
plot(m3_SVL_hydroperiod_predict_F, rawdata = T)


#### (iii) males #### 
plot(SVL ~ AGE, data = metamorphs_growth_male)

fishmethods::growth(size=metamorphs_growth_male$SVL,age=metamorphs_growth_male$AGE,
                    error=1, Sinf=101,K=0.1,t0=-12) #good starting values

#fit the von Bertalanffy growth function using nonlinear least squares (nls) optimization
growth_mod_male <- nls(vbmod, data = metamorphs_growth_male, start = list(Linf = 101, K = 0.1, t0 = -12))
summary(growth_mod_male)

#predict
growth_predict_male = predict(growth_mod_male)

##bootstrap 95% CI

#fit the model to the data using nls, like we did before
vb_fit_male <- nls(SVL~vbO(AGE,Linf,K, t0), data = metamorphs_growth_male, start = list(Linf = 101, K = 0.1, t0 = -12))

#now, bootstrap the model fitting process
boot_fit_male <- nlsBoot(vb_fit_male, niter = 9999)

#predict length at age from the model (t is age). Here, we tell R to predict length at 
#each unique age in our original data, and calculate some bootstrapped confidence
#intervals
boot_preds_male <- data.frame(predict(boot_fit_male, vbO, t = sort(unique(metamorphs_growth_male$AGE))))
names(boot_preds_male) <- c("AGE", "fit", "lwr", "upr")

growth_preds_male <- merge(metamorphs_growth_male, boot_preds_male, by = "AGE")

#plot
growth_curve_plot_male = 
  ggplot(growth_preds_male, aes(x = AGE, y = SVL)) +
  geom_jitter(width = 0.1, alpha = 0.15, size = 5) +
  geom_line(aes(y = fit), linewidth = 2) +
  geom_ribbon(
    aes(x = AGE, ymin = lwr, ymax = upr, color = NULL), alpha = 0.3) +
  scale_y_continuous(limits = c(59,120)) +
  xlab("Age (years)") +
  ylab("Snout-Vent Length (mm)") +
  theme_bw(40)

growth_curve_plot_male

#ggsave(growth_curve_plot_male, filename = "Outputs/growth_curve_plot_male.png",  width = 14, height = 12, dpi = 100)


#calculate residuals 
growth_preds_male = growth_preds_male[,c("DATEOFYEAR", "Year", "INDIV_ID", "INDIV_ID_Year",
                                             "MORPH", "POND", "COHORT", "SEX", "WEIGHT", "SVL",
                                             "fit", "lwr", "upr")]

names(growth_preds_male)[names(growth_preds_male) == 'fit'] <- 'SVL_fit'
names(growth_preds_male)[names(growth_preds_male) == 'lwr'] <- 'SVL_fit_lwr'
names(growth_preds_male)[names(growth_preds_male) == 'upr'] <- 'SVL_fit_upr'

growth_preds_male$SVL_residual = growth_preds_male$SVL - growth_preds_male$SVL_fit 
hist(growth_preds_male$SVL_residual)

growth_preds_male = droplevels(growth_preds_male)
nlevels(growth_preds_male$INDIV_ID)

num_captures = growth_preds_male %>% group_by(INDIV_ID) %>% summarise(N = n())
range(num_captures$N)
mean(num_captures$N)
sd(num_captures$N)

#create df for SVL residual data for metamorphs that were captured at least twice. In other words, eliminate
#metas that were captured in one pond in one year only. These individuals don't have the opportunity to learn
#about pond resources, and they're more likely to end up in ponds based on landscape factors (e.g., biggest pond)
growth_preds_male_multiple_obs = growth_preds_male %>% 
  add_count(INDIV_ID, name = "Num_Captures") %>% 
  filter(Num_Captures > 1)

#add in mean SVL residual data for each meta over their lifetime
growth_preds_male_avg = growth_preds_male_multiple_obs %>% group_by(INDIV_ID) %>% summarise(SVL_residual_mean = mean(SVL_residual))
hist(growth_preds_male_avg$SVL_residual_mean)
growth_preds_male = inner_join(growth_preds_male, growth_preds_male_avg, by = "INDIV_ID")  
growth_preds_male = droplevels(growth_preds_male)

metamorphs_SVL_immigration_factors_male = merge(growth_preds_male, meta_mvmt_factors,
                                                  by = "POND", all.x = T)

#model
m1_SVL_fs_rank_M = lmer(SVL_residual ~ Mean_Historical_FS_Rank + (1|Year) + (1|INDIV_ID), 
                        data = metamorphs_SVL_immigration_factors_male)

m2_SVL_landscape_M = lmer(SVL_residual ~ scale(Size) + scale(Dist_Avg) + (1|Year) + (1|INDIV_ID), 
                          data = metamorphs_SVL_immigration_factors_male)

m3_SVL_hydroperiod_M = lmer(SVL_residual ~ Hydroperiod + (1|Year) + (1|INDIV_ID), 
                            data = metamorphs_SVL_immigration_factors_male)


AIC(m1_SVL_fs_rank_M, m3_SVL_hydroperiod_M, m2_SVL_landscape_M)
car::Anova(m1_SVL_fs_rank_M)
car::Anova(m2_SVL_landscape_M)
car::Anova(m3_SVL_hydroperiod_M)


m1_SVL_fs_rank_predict_M <- ggpredict(m1_SVL_fs_rank_M, terms = "Mean_Historical_FS_Rank")
plot(m1_SVL_fs_rank_predict_M)
plot(m1_SVL_fs_rank_predict_M, rawdata = T)

m2_SVL_landscape_predict_M <- ggpredict(m2_SVL_landscape_M, terms = "Size")
plot(m2_SVL_landscape_predict_M)
plot(m2_SVL_landscape_predict_M, rawdata = T)

m3_SVL_hydroperiod_predict_M <- ggpredict(m3_SVL_hydroperiod_M, terms = "Hydroperiod")
plot(m3_SVL_hydroperiod_predict_M)
plot(m3_SVL_hydroperiod_predict_M, rawdata = T)



#### Delta Condition ~ FS ####
delta_condition = metamorphs %>% group_by(INDIV_ID, Year, POND) %>% filter(n() > 1)

delta_condition$Condition = delta_condition$WEIGHT / delta_condition$SVL

delta_condition = delta_condition %>%
  group_by(INDIV_ID, Year, POND) %>%
  arrange(DATEOFYEAR, .by_group = TRUE) %>%
  mutate(Condition_Change = Condition - lag(Condition, default = first(Condition))) %>%
  mutate(Date_Diff = DATEOFYEAR - lag(DATEOFYEAR, default = first(DATEOFYEAR)))

delta_condition = delta_condition[,c("POND", "INDIV_ID", "Condition_Change", "Date_Diff")]

delta_condition = subset(delta_condition, Condition_Change != 0 & POND != "1A" & !is.na(POND))
delta_condition = droplevels(delta_condition)
delta_condition = subset(delta_condition, Condition_Change < 0.6 & Condition_Change > -0.5) #outliers?

delta_condition$Date_Diff = as.numeric(delta_condition$Date_Diff)
delta_condition$Condition_Change_Rate = delta_condition$Condition_Change / delta_condition$Date_Diff
delta_condition = subset(delta_condition, !is.infinite(Condition_Change_Rate))

delta_condition = subset(delta_condition, POND == "10" | POND == "11" | POND == "15" |
                                  POND == "51" | POND == "52" | POND == "6" | POND == "8" |  
                                  POND == "13")

delta_condition$POND <- factor(delta_condition$POND,levels = c("6", "8", "10", "11", "51", "52", "13", "15"))
delta_condition = merge(delta_condition, unique(metamorphs_SVL_immigration_factors[,c("POND", "Hydroperiod")]),
                        by = "POND", all.x = T)

m1_delta_condition = lm(Condition_Change_Rate ~ Hydroperiod, data = delta_condition)
#plot(m1_delta_condition)
car::Anova(m1_delta_condition)

m1_delta_condition_predict <- ggpredict(m1_delta_condition, "Hydroperiod")
plot(m1_delta_condition_predict)
plot(m1_delta_condition_predict, add.data = TRUE)


##smalls
quantile(growth_preds$SVL_residual, probs = 0.5)
delta_svl_smallest = subset(growth_preds, SVL_residual < 0.422745)
delta_svl_smallest = droplevels(delta_svl_smallest)

delta_condition_smalls = semi_join(delta_condition, delta_svl_smallest, by = "INDIV_ID")  

m2_delta_condition = lm(Condition_Change_Rate ~ POND, data = delta_condition_smalls)
#plot(m1_delta_condition)
car::Anova(m2_delta_condition)

m2_delta_condition_predict <- ggpredict(m2_delta_condition, "POND")
plot(m2_delta_condition_predict)
plot(m2_delta_condition_predict, add.data = TRUE)


#M vs. F number in ponds
pond15 = subset(metamorphs, POND == "15")
summary(pond15$SEX)
#Female   Male 
# 22       6   =  79% females

pond52 = subset(metamorphs, POND == "52")
summary(pond52$SEX)
#Female   Male 
# 41       16  =  72% females


#all I've found is evidence for competitive exclusion in small males...


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
  scale_x_continuous(limits = c(0,12), breaks = seq(0,12, by = 3), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradientn(colours = paletteer_dynamic("cartography::wine.pal", 20), na.value = "grey50") +
  theme_bw(45) + 
  labs(y = "Pond", x = "Sampling Occasion", fill = "Percent Females") +
  theme(axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        legend.key.size = unit(2, 'cm'))

sexratio_heatmap_plot_2021

#ggsave(sexratio_heatmap_plot_2021, filename = "Outputs/sexratio_heatmap_plot_2021.png",  width = 20, height = 12, dpi = 100)



