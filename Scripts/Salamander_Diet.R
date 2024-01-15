
######### Salamander Gut Flush Data #########

rm(list = ls())

#library(remotes)
#remotes::install_github("Ajfrick/ajfhelpR", force = T)
library(ajfhelpR)
library(ggplot2)
library(tidyr)
library(purrr)
library(dplyr)
library(lubridate)
library(data.table)
library(RColorBrewer)
library(scales)
library(car)
library(lme4)
library(lmerTest)
library(visreg)
library(sjPlot)
library(lattice)
library(vegan)
library(permute)
library(statnet.common)
library(network)
library(sna)
library(bipartite)
library(indicspecies)
library(betareg)
library(glmmTMB)
library(drc)

meta_diet_count <- read.delim("Data/Salamander_Gut_Counts.txt", stringsAsFactors = FALSE, sep = "\t")
meta_diet_biomass <- read.csv(file = "Data/Salamander_Gut_Biomass.csv", stringsAsFactors = FALSE)
diet_unknown <- read.csv(file = "Data/Salamander_Diet_Taxa_Unknown.csv", stringsAsFactors = FALSE)

meta_diet_count <-  meta_diet_count[!(is.na(meta_diet_count$Date_Lab) | meta_diet_count$Date_Lab==""), ] #remove blank rows

#set aside biomass data for caddis bodies and cases from 2022
caddis_body_case = subset(meta_diet_biomass, Diet_Taxa == "Caddis Body" | Diet_Taxa == "Caddis Case")
meta_diet_biomass = subset(meta_diet_biomass, Diet_Taxa != "Caddis Body" & Diet_Taxa != "Caddis Case")

#add leading zeros to Sample_ID so they are all two digits
meta_diet_count$Sample_ID = sprintf("%02d", meta_diet_count$Sample_ID)
meta_diet_biomass$Sample_ID = as.numeric(meta_diet_biomass$Sample_ID)
meta_diet_biomass$Sample_ID = sprintf("%02d", meta_diet_biomass$Sample_ID)
diet_unknown$Sample_ID = sprintf("%02d", diet_unknown$Sample_ID)

#make Sample_ID_Year category since there are repeated Sample_ID numbers across years (whoops)
meta_diet_count$Sample_ID_Year = interaction(meta_diet_count$Sample_ID, meta_diet_count$Year, sep="_")
meta_diet_count = meta_diet_count %>% relocate(Sample_ID_Year, .after = Sample_ID)

meta_diet_biomass$Sample_ID_Year = interaction(meta_diet_biomass$Sample_ID, meta_diet_biomass$Year, sep="_")
meta_diet_biomass = meta_diet_biomass %>% relocate(Sample_ID_Year, .after = Sample_ID)

diet_unknown$Sample_ID_Year = interaction(diet_unknown$Sample_ID, diet_unknown$Year, sep="_")
diet_unknown = diet_unknown %>% relocate(Sample_ID_Year, .after = Sample_ID)

#add in salamander ID and sample data variables to biomass df
meta_diet_biomass = merge(meta_diet_biomass, unique(meta_diet_count[,c("Sample_ID_Year", "Date_Sample", "Salamander_ID")]), by = "Sample_ID_Year", all.x = T)
meta_diet_biomass = meta_diet_biomass %>% relocate(c(Salamander_ID, Date_Sample), .after = Sample_ID)
meta_diet_biomass$Salamander_ID = as.factor(meta_diet_biomass$Salamander_ID)
meta_diet_biomass$Date_Sample = mdy(meta_diet_biomass$Date_Sample)


#change Organic_Debris category to Fairy_Shrimp. Based on notes, this category is mostly fairy shrimp. I
#started gut sample processing with counts (followed by biomass), and so FS parts that were not easy to
#count were included in this Organic_Debris category.
#meta_diet_count$Diet_Taxa = dplyr::if_else(meta_diet_count$Diet_Taxa == "Organic_Debris", "Fairy_Shrimp", meta_diet_count$Diet_Taxa)
#meta_diet_biomass$Diet_Taxa = dplyr::if_else(meta_diet_biomass$Diet_Taxa == "Debirs", "Fairy_Shrimp", meta_diet_biomass$Diet_Taxa)

#on second thought, remove FS debris bc I was counting by heads. Not going to get all FS biomass
#that was in the stomach anyways.
meta_diet_count = subset(meta_diet_count, Diet_Taxa != "Organic_Debris" & Diet_Taxa != "Debris" & Diet_Taxa != "Fairy_Shrimp_eggs")
meta_diet_biomass = subset(meta_diet_biomass, Diet_Taxa != "Organic_Debris" & Diet_Taxa != "Debirs" & Diet_Taxa != "Fairy_Shrimp_eggs")

meta_diet_count$Diet_Taxa = as.factor(meta_diet_count$Diet_Taxa)
meta_diet_count$Sample_ID = as.factor(meta_diet_count$Sample_ID)
meta_diet_count$Sample_ID_Year = as.factor(meta_diet_count$Sample_ID_Year)
meta_diet_biomass$Diet_Taxa = as.factor(meta_diet_biomass$Diet_Taxa)
meta_diet_biomass$Sample_ID = as.factor(meta_diet_biomass$Sample_ID)
meta_diet_biomass$Sample_ID_Year = as.factor(meta_diet_biomass$Sample_ID_Year)

##merge in IDs for initially unknown diet taxa to main dataframes
diet_unknown$X <- NULL
diet_unknown$Notes <- NULL

meta_diet_count = merge(meta_diet_count, diet_unknown, by = c("Sample_ID", "Sample_ID_Year", "Year", "Diet_Taxa"), all = T)
meta_diet_count = meta_diet_count %>% relocate(Diet_Taxa_ID, .after = Diet_Taxa)
meta_diet_count = meta_diet_count %>% mutate(Diet_Taxa_ID = coalesce(Diet_Taxa_ID,Diet_Taxa))
meta_diet_count$Diet_Taxa_ID = as.factor(meta_diet_count$Diet_Taxa_ID)

#make diet categories
levels(factor(meta_diet_count$Diet_Taxa_ID))

meta_diet_count$Diet_Category = with(meta_diet_count, #10 categories
                                ifelse(Diet_Taxa_ID %in% c("Asynarchus", "Caddisfly"), "Caddisfly", 
                                ifelse(Diet_Taxa_ID %in% c("Beetle", "Dytiscus", "Hydroporinae"), "Beetle", 
                                ifelse(Diet_Taxa_ID %in% c("Boatman", "Hemipteran", "Hemiptera"), "Hemiptera",
                                ifelse(Diet_Taxa_ID %in% c("Chaoborus", "Chironomidae", "Chironomidae ", "Mosquito"), "Midge",
                                ifelse(Diet_Taxa_ID %in% c("Copepod", "Daphnia", "Copepod_Daphnia"), "Zooplankton",
                                ifelse(Diet_Taxa_ID %in% c("Dragonfly"), "Dragonfly",
                                ifelse(Diet_Taxa_ID %in% c("Fairy_Shrimp"), "Fairy Shrimp",
                                ifelse(Diet_Taxa_ID %in% c("Hymenoptera", "Ant"), "Hymenoptera",
                                ifelse(Diet_Taxa_ID %in% c("Moth", "Moths"), "Moth",  "Other"))))))))))

meta_diet_count = meta_diet_count %>% relocate(Diet_Category, .after = Diet_Taxa_ID)
meta_diet_count$Diet_Category = as.factor(meta_diet_count$Diet_Category)

#Other category includes counts of: "Crambidae", "Fragment", "Grasshopper",
#"Miscellaneous", "Mite, "Needles", "Pine_Needle", "Spider", "Unknown", "Wood_chunks"

##repeat for meta_diet_biomass dataframe
meta_diet_biomass = merge(meta_diet_biomass, diet_unknown, by = c("Sample_ID", "Sample_ID_Year", "Year", "Diet_Taxa"), all = T)
meta_diet_biomass = meta_diet_biomass %>% relocate(Diet_Taxa_ID, .after = Diet_Taxa)
meta_diet_biomass = meta_diet_biomass %>% mutate(Diet_Taxa_ID = coalesce(Diet_Taxa_ID,Diet_Taxa))
meta_diet_biomass$Diet_Taxa_ID = as.factor(meta_diet_biomass$Diet_Taxa_ID)

#make diet categories
levels(factor(meta_diet_biomass$Diet_Taxa_ID))

meta_diet_biomass$Diet_Category = with(meta_diet_biomass, #12 categories (11 in meta samples; 1 Grasshopper in a paedo)
                                     ifelse(Diet_Taxa_ID %in% c("Asynarchus", "Caddisfly"), "Trichoptera", 
                                     ifelse(Diet_Taxa_ID %in% c("Beetle", "Dytiscus", "Hydroporinae"), "Coleoptera", 
                                     ifelse(Diet_Taxa_ID %in% c("Boatman", "Hemipteran", "Hemiptera"), "Hemiptera",
                                     ifelse(Diet_Taxa_ID %in% c("Chaoborus", "Chironomidae", "Chironomidae ", "Mosquito"), "Diptera",
                                     ifelse(Diet_Taxa_ID %in% c("Copepod", "Daphnia", "Copepod_Daphnia"), "Copepoda & Cladocera",
                                     ifelse(Diet_Taxa_ID %in% c("Dragonfly"), "Odonata",
                                     ifelse(Diet_Taxa_ID %in% c("Fairy_Shrimp"), "Anostraca",
                                     ifelse(Diet_Taxa_ID %in% c("Hymenoptera", "Ant"), "Hymenoptera",
                                     ifelse(Diet_Taxa_ID %in% c("Moth", "Moths", "Crambidae"), "Lepidoptera",  
                                     ifelse(Diet_Taxa_ID %in% c("Grasshopper"), "Orthoptera",
                                     ifelse(Diet_Taxa_ID %in% c("Mite"), "Acari",
                                     ifelse(Diet_Taxa_ID %in% c("Spider"), "Araneae","DROP")))))))))))))

#DROP:"Fragment", "Miscellaneous", "Unknown"

meta_diet_biomass = meta_diet_biomass %>% relocate(Diet_Category, .after = Diet_Taxa_ID)
meta_diet_biomass$Diet_Category = as.factor(meta_diet_biomass$Diet_Category)

meta_diet_biomass = subset(meta_diet_biomass, Diet_Category != "DROP")

#separate out metas and paedos
paedos = subset(meta_diet_biomass, Sample_ID_Year == "57_2021" | Sample_ID_Year == "58_2021" | Sample_ID_Year == "59_2021" | Sample_ID_Year == "60_2021" | Sample_ID_Year == "61_2021" | Sample_ID_Year == "62_2021") 
paedos$Sample_ID_Year = factor(paedos$Sample_ID_Year)

meta_diet_count = subset(meta_diet_count, Sample_ID_Year != "57_2021" & Sample_ID_Year != "58_2021" & Sample_ID_Year != "59_2021" & Sample_ID_Year != "60_2021" & Sample_ID_Year != "61_2021" & Sample_ID_Year != "62_2021")
meta_diet_count$Sample_ID_Year = factor(meta_diet_count$Sample_ID_Year)

meta_diet_biomass = subset(meta_diet_biomass, Sample_ID_Year != "57_2021" & Sample_ID_Year != "58_2021" & Sample_ID_Year != "59_2021" & Sample_ID_Year != "60_2021" & Sample_ID_Year != "61_2021" & Sample_ID_Year != "62_2021")
meta_diet_biomass$Sample_ID_Year = factor(meta_diet_biomass$Sample_ID_Year)

#calculate net weight for diet taxa
meta_diet_biomass$Boat_Net = meta_diet_biomass$Boat_Full - meta_diet_biomass$Boat_Empty
paedos$Boat_Net = paedos$Boat_Full - paedos$Boat_Empty

meta_diet_biomass = subset(meta_diet_biomass, Boat_Net > 0) #two obs with negative net weights -- Note for one appears to be associated with wrong observation


#look at Notes for any action items before beginning analysis -- no action items
meta_diet_count$Notes
meta_diet_biomass$Notes

meta_diet_count$Notes <- NULL
meta_diet_biomass$Notes <- NULL

#convert weights in g to mg
meta_diet_biomass$Boat_Net = ifelse(meta_diet_biomass$Units_Full == "g", meta_diet_biomass$Boat_Net*1000, meta_diet_biomass$Boat_Net)
paedos$Boat_Net = ifelse(paedos$Units_Full == "g", paedos$Boat_Net*1000, paedos$Boat_Net)

meta_diet_biomass = meta_diet_biomass[ ,c("Sample_ID_Year", "Salamander_ID", "Date_Sample", "Diet_Category", "Boat_Net")]

#add explicit zero observations and fill in missing values 
meta_diet_biomass = droplevels(meta_diet_biomass)

meta_diet_biomass_full = 
  meta_diet_biomass %>%
  tidyr::complete(Sample_ID_Year, Diet_Category, fill = list(Boat_Net = 0)) %>% 
  dplyr::group_by(Sample_ID_Year, Diet_Category) %>%
  mutate(Boat_Net = sum(Boat_Net))

meta_diet_biomass_full = merge(meta_diet_biomass_full, unique(meta_diet_biomass[,c("Sample_ID_Year", "Salamander_ID", "Date_Sample")]), by = "Sample_ID_Year", all.x = T)
meta_diet_biomass_full = meta_diet_biomass_full[ , -which(names(meta_diet_biomass_full) %in% c("Salamander_ID.x", "Date_Sample.x"))]
names(meta_diet_biomass_full)[names(meta_diet_biomass_full) == 'Salamander_ID.y'] <- 'Salamander_ID'
names(meta_diet_biomass_full)[names(meta_diet_biomass_full) == 'Date_Sample.y'] <- 'Date_Sample'
meta_diet_biomass_full = meta_diet_biomass_full %>% relocate(Boat_Net, .after = Date_Sample)


#meta_diet_biomass <- meta_diet_biomass[!(meta_diet_biomass$Boat_Net == 0 & meta_diet_biomass$Diet_Category != "Fairy Shrimp"),] 
#keep only FS zeros

meta_diet_biomass_full$Year = substr(meta_diet_biomass_full$Sample_ID_Year, start = 4, stop = 7)
meta_diet_biomass_full$Sample_ID = substr(meta_diet_biomass_full$Sample_ID_Year, start = 1, stop = 2)
meta_diet_biomass_full = unique(meta_diet_biomass_full) 
meta_diet_biomass_full$Sample_ID = as.numeric(as.character(meta_diet_biomass_full$Sample_ID))

#exclude weight of asynarchus cases from sample (samples 31 & 33 definitely. sample 32 as well?)
#sample 31_2021 = 5 cases
#sample 33_2021 = 6 cases
caddis_cases = c(22.2, 45.6, 92.5, 125.6) #didn't include largest case (184.7 mg). Mean value would have given negative gut content weights for 31 and 33
mean(caddis_cases) #71.475
sd(caddis_cases) #46.43672

#1 case = 71.475 mg (sample 32_2021)
#5 cases = 357.375 mg (sample 31_2021)
#6 cases = 428.85 mg (sample 33_2021)

meta_diet_biomass_full$Boat_Net <- ifelse(meta_diet_biomass_full$Sample_ID_Year == "31_2021" 
                                          & meta_diet_biomass_full$Diet_Category == "Trichoptera", 
                                          (409.841 - 357.375), meta_diet_biomass_full$Boat_Net)

meta_diet_biomass_full$Boat_Net <- ifelse(meta_diet_biomass_full$Sample_ID_Year == "32_2021" 
                                          & meta_diet_biomass_full$Diet_Category == "Trichoptera", 
                                          (112.138 - 71.475), meta_diet_biomass_full$Boat_Net)

meta_diet_biomass_full$Boat_Net <- ifelse(meta_diet_biomass_full$Sample_ID_Year == "33_2021" 
                                          & meta_diet_biomass_full$Diet_Category == "Trichoptera", 
                                          (502.272 - 428.85), meta_diet_biomass_full$Boat_Net)


#need to have the Sample_ID numbers be continuous across years and not overlap
meta_diet_biomass_2021 = subset(meta_diet_biomass_full, Year == 2021)

meta_diet_biomass_2022 = subset(meta_diet_biomass_full, Year == 2022)
meta_diet_biomass_2022 = meta_diet_biomass_2022 %>% arrange(Sample_ID)
meta_diet_biomass_2022 = meta_diet_biomass_2022 %>% group_by(Sample_ID) %>% mutate(Sample_ID_True = cur_group_id())
meta_diet_biomass_2022$Sample_ID_True = meta_diet_biomass_2022$Sample_ID_True + 56
meta_diet_biomass_2022$Sample_ID <- NULL
names(meta_diet_biomass_2022)[names(meta_diet_biomass_2022) == 'Sample_ID_True'] <- 'Sample_ID'
meta_diet_biomass_2022 = meta_diet_biomass_2022 %>% relocate(Sample_ID, .before = Year)

meta_diet_biomass_full = rbind(meta_diet_biomass_2021, meta_diet_biomass_2022)
meta_diet_biomass_full = meta_diet_biomass_full %>% arrange(Sample_ID)

##plot

#all diet categories
#metas

display.brewer.all(colorblindFriendly = TRUE) #second set of colors = good for categorical datasets

levels(droplevels(meta_diet_biomass)$Diet_Category) #11 bins
nb.cols <- 11
mycolors <- colorRampPalette(brewer.pal(11, "Set2"))(nb.cols)
mycolors
show_col(mycolors)

meta_diet_biomass_plot = 
  ggplot(meta_diet_biomass_full, aes(fill=Diet_Category, y=Boat_Net, x=Sample_ID)) + 
  scale_fill_manual(values = mycolors) +
  geom_bar(position="stack", stat="identity") +
  ylim(0,600) +
  scale_x_continuous(breaks = 1:100) +
  ylab("Diet Taxa Biomass (mg)") +
  theme_bw(55) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5),
        axis.ticks.x = element_blank())

print(meta_diet_biomass_plot)

ggsave(meta_diet_biomass_plot, filename = "Outputs/meta_diet_biomass_plot.png",  width = 24, height = 12, dpi = 100)


#paedos

ggplot(paedos, aes(fill=Diet_Category, y=Boat_Net, x=Sample_ID)) + 
  geom_bar(position="stack", stat="identity") +
  ylim(0,150) +
  ylab("Diet Taxa Biomass (mg)") +
  theme_bw(35) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


#just FS (metas)
meta_diet_biomass_FS = subset(meta_diet_biomass_full, Diet_Category == "Anostraca")

meta_diet_biomass_plot_FS = 
  ggplot(meta_diet_biomass_FS, aes(y=Boat_Net, x=Sample_ID)) + 
  geom_bar(position="stack", stat="identity") +
  ylim(0,600) +
  ylab("Fairy Shrimp Biomass (mg)") +
  theme_bw(35) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

print(meta_diet_biomass_plot_FS)

ggsave(meta_diet_biomass_plot_FS, filename = "Outputs/meta_diet_biomass_plot_FS.png",  width = 18, height = 12, dpi = 100)


##histogram of fairy shrimp count in meta stomachs
meta_diet_count_fs = subset(meta_diet_count, Diet_Taxa_ID == "Fairy_Shrimp" & !is.na(Count))
meta_diet_count_fs$Count = as.numeric(meta_diet_count_fs$Count)
hist(meta_diet_count_fs$Count, col = 'skyblue3', breaks = 20, xlim = c(0,400), ylim = c(0, 20),
     xlab = "Number of Fairy Shrimp in Metamorph Stomach", ylab = "Frequency", main = "")


####### Gut FS ~ Pond FS #######

###pond-level dataframe
fs_pond <- read.csv(file = "Data/fs_pond.csv", stringsAsFactors = FALSE) #may have changed!
fs_pond$X <- NULL

fs_pond$lower_biomass = fs_pond$log_FS_Biomass_Pond_m3 - fs_pond$SE_Biomass*1.96
fs_pond$upper_biomass = fs_pond$log_FS_Biomass_Pond_m3 + fs_pond$SE_Biomass*1.96

meta_diet_biomass_FS = merge(meta_diet_biomass_FS, unique(meta_diet_count[,c("Sample_ID_Year", "Pond")]), by = "Sample_ID_Year", all.x = T)

fs_pond$Date = ymd(fs_pond$Date)
fs_pond$Pond = as.factor(fs_pond$Pond)
meta_diet_biomass_FS$Pond = as.factor(meta_diet_biomass_FS$Pond)
meta_diet_biomass_FS$Year = as.factor(meta_diet_biomass_FS$Year)

names(meta_diet_biomass_FS)[names(meta_diet_biomass_FS) == 'Date_Sample'] <- 'Date_Sample_GF'
names(fs_pond)[names(fs_pond) == 'Date'] <- 'Date_Sample_FS'

#find nearest dates of FS pond sampling to Meta gut flush date, for each pond
meta_diet_biomass_FS = subset(meta_diet_biomass_FS, Sample_ID != 94)
meta_diet_biomass_FS[,c("Pond")] <- factor(meta_diet_biomass_FS[,c("Pond")], levels=c("6", "8", "10", "11", "13", "15", "51", "52", "55"))

temp <- list()

for (i in 1:nrow(meta_diet_biomass_FS)) {       
              fs_pond_subset <- subset(fs_pond, Pond==(meta_diet_biomass_FS[i, 8]))
              temp[[i]] <- ajfhelpR::date_near(fs_pond_subset$Date_Sample_FS, meta_diet_biomass_FS[i, 4], sidepref = 'l')
}

temp

dates = temp %>% purrr::reduce(c) #make list into vector
rm(temp)

meta_diet_biomass_FS$Date_Sample_FS <- dates
meta_diet_biomass_FS$Date_Sample_FS - meta_diet_biomass_FS$Date_Sample_GF #time difference btwn GF and FS sampling

meta_diet_biomass_FS = merge(meta_diet_biomass_FS, fs_pond[,c("Date_Sample_FS", "Pond", "FS_Biomass_Pond_m3", "log_FS_Biomass_Pond_m3",
                                                              "SE_Biomass", "lower_biomass", "upper_biomass")], 
                                                              by = c("Date_Sample_FS", "Pond"), all.x = T)

#need to add a small amount to densities = 0 (can't log transform)
meta_diet_biomass_FS$Boat_Net = meta_diet_biomass_FS$Boat_Net + 1
meta_diet_biomass_FS$log_Boat_Net = log(meta_diet_biomass_FS$Boat_Net)
meta_diet_biomass_FS$Boat_Net = meta_diet_biomass_FS$Boat_Net - 1

meta_diet_biomass_FS = meta_diet_biomass_FS %>% relocate(log_Boat_Net, .after = Boat_Net)

#merge in sex data
salamander = read.csv("Data/Salamander_Master_cleaned.csv", stringsAsFactors = FALSE) #convert to .txt file if issues with INDIV_ID and scientific notation
salamander$X <- NULL
names(salamander)[names(salamander) == 'INDIV_ID'] <- 'Salamander_ID'
names(salamander)[names(salamander) == 'DATEOFYEAR'] <- 'Date_Sample'
salamander$Date_Sample = mdy(salamander$Date_Sample)
salamander$Salamander_ID = as.factor(salamander$Salamander_ID)

meta_diet_biomass_FS = merge(meta_diet_biomass_FS, unique(salamander[,c("Salamander_ID", "SEX")]))
meta_diet_biomass_FS = meta_diet_biomass_FS %>% relocate(SEX, .after = Salamander_ID)
meta_diet_biomass_FS$SEX = car::recode(meta_diet_biomass_FS$SEX,"'Female?'='Female'")
meta_diet_biomass_FS$SEX = as.factor(meta_diet_biomass_FS$SEX)

#code for sex-specific graph. switch out to "Male" for male plots.
meta_diet_biomass_FS = subset(meta_diet_biomass_FS, SEX == "Female")

#average FS in gut over all individuals per pond
meta_diet_biomass_FS = meta_diet_biomass_FS %>% group_by(Date_Sample_GF, Pond) %>% mutate(log_Boat_Net_Avg = mean(log_Boat_Net))
meta_diet_biomass_FS = meta_diet_biomass_FS %>% group_by(Date_Sample_GF, Pond) %>% mutate(log_Boat_Net_SE = sd(log_Boat_Net)/sqrt(n()))
meta_diet_biomass_FS$FS_total_lower = meta_diet_biomass_FS$log_Boat_Net_Avg - meta_diet_biomass_FS$log_Boat_Net_SE*1.96
meta_diet_biomass_FS$FS_total_upper = meta_diet_biomass_FS$log_Boat_Net_Avg + meta_diet_biomass_FS$log_Boat_Net_SE*1.96

#bound 95% CI
meta_diet_biomass_FS$FS_total_lower[meta_diet_biomass_FS$FS_total_lower < 0] <- 0

meta_diet_biomass_FS_avg_total = meta_diet_biomass_FS[,c("Date_Sample_GF", "Year", "Pond", "FS_Biomass_Pond_m3", 
                                                    "log_FS_Biomass_Pond_m3", "SE_Biomass", "lower_biomass", 
                                                    "upper_biomass", "log_Boat_Net_Avg", 
                                                    "log_Boat_Net_SE", "FS_total_lower", "FS_total_upper")]

meta_diet_biomass_FS_avg_total = unique(meta_diet_biomass_FS_avg_total)

#remove one point (which represents two metas) that had FS in their gut despite no FS present
#during sampling. Probably a result of sampling for FS a day after gut flush. So the metas ate
#the last remaining FS or the temperature in the pond reached the CT max for FS when we sampled. 
#Metas could have moved from another pond where they ate FS, but seems unlikely given that both 
#were caught in pond 13 previously
meta_diet_biomass_FS_avg_total = subset(meta_diet_biomass_FS_avg_total, !(Year == "2021" & Pond == "13"))


##build logistic curve before plotting
#great resource for making more custom sigmoid curves than standard logisitic regression:
#http://www.darrenkoppel.com/2020/09/04/dose-response-modelling-and-model-selection-in-r/

getMeanFunctions() #use LL.4 -- log-logistic

mod1 <- drm(log_Boat_Net_Avg ~ FS_Biomass_Pond_m3, data = meta_diet_biomass_FS_avg_total, 
            fct=LL.4(fixed=c(-1, 0, 5.4124484, NA),
                     names = c("Slope", "Lower Limit", "Upper Limit", "Midpoint")), type = "continuous")

#female: fixed=c(-1, 0, 5.766789, NA)
#male: fixed=c(-1, 0, 5.4124484, NA)

summary(mod1)
plot(mod1)


#create new data
newdata = seq(min(meta_diet_biomass_FS_avg_total$FS_Biomass_Pond_m3), max(meta_diet_biomass_FS_avg_total$FS_Biomass_Pond_m3) + 500000, length.out = 100000)
newdata = as.data.frame(newdata)
colnames(newdata) = c("FS_Biomass_Pond_m3")
fit1_predict = predict(mod1, type="response", se.fit = F, newdata=newdata)
fit1_predict = as.data.frame(fit1_predict)
fit1_predict = cbind(fit1_predict, newdata)
colnames(fit1_predict) = c("log_Boat_Net_Avg", "FS_Biomass_Pond_m3")
plot(log(fit1_predict$FS_Biomass_Pond_m3 + 1), fit1_predict$log_Boat_Net_Avg)

fs_stomach_v_pond_plot =
  ggplot(meta_diet_biomass_FS_avg_total, aes(x=log_FS_Biomass_Pond_m3, y=log_Boat_Net_Avg)) + 
  geom_errorbarh(aes(xmin=lower_biomass, xmax=upper_biomass), color = "black", height = 0.15, size = 1.5) +
  geom_errorbar(aes(ymin=FS_total_lower, ymax=FS_total_upper), color = "black", width = 0.15, size = 1.5) +
  geom_point(aes(x=log_FS_Biomass_Pond_m3, y=log_Boat_Net_Avg, color = Pond), size=10, alpha = 0.9) +
  geom_line(data = fit1_predict, aes(x = log(FS_Biomass_Pond_m3 + 1), y = log_Boat_Net_Avg), size = 3) +
  ylab("Fairy Shrimp Biomass in\n Stomach (mg) ± 95% CI") +
  xlab(expression(paste("Fairy Shrimp Biomass in Pond", " ", "(mg/", m^{3}, ")", " ", "± 95% CI"))) +
  coord_cartesian(xlim = c(log(1),log(300000)), ylim = c(-0.2,7)) +
  scale_x_continuous(breaks=c(log(1),log(10),log(100),log(1000),log(10000),log(100000)), labels=c("0", "10", "100", "1,000", "10,000", "100,000")) +
  scale_y_continuous(breaks=c(log(1),log(10),log(100),log(1000)), labels=c("0", "10", "100", "1,000")) +
  scale_color_manual(values = c("#66C2A5", "#A6D854", "#8DA0CB", "#E78AC3", "#FFD92F", "#FC8D62", "#E5C494", "#B3B3B3")) +
  theme_bw(42)

print(fs_stomach_v_pond_plot)

ggsave(fs_stomach_v_pond_plot, filename = "Outputs/fs_stomach_v_pond.png",  width = 18, height = 12, dpi = 100)

#what's going on with the points in the upper left - interpretations?

#stats
#mixed model
fs_stomach_v_pond = lmer(log_Boat_Net ~ log_FS_Biomass_Pond_m3 + (1|Pond), data = meta_diet_biomass_FS)

plot(fs_stomach_v_pond, type=c("p","smooth"), col.line=1) #residuals look really weird
qqmath(fs_stomach_v_pond)

anova(fs_stomach_v_pond, type = 2, ddf = "Kenward-Roger")

visreg(fs_stomach_v_pond, xvar = "log_FS_Biomass_Pond_m3", type = "conditional",
       whitespace = 0.4, points.par = list(cex = 1.1, col = "red"),
       ylab="log(Fairy Shrimp Biomass in Gut)", xlab = "log(Fairy Shrimp Biomass in Pond)")

tab_model(fs_stomach_v_pond, show.reflvl = TRUE, show.re.var = TRUE)
#https://cran.r-project.org/web/packages/sjPlot/vignettes/tab_mixed.html
#marginal R2 is the proportion of the total variance explained by the fixed effects
#conditional R2 is the proportion of the variance explained by both fixed and random effects

meta_diet_biomass_full_total = meta_diet_biomass_full %>% group_by(Sample_ID_Year) %>% mutate(Total_Diet_Biomass = sum(Boat_Net))


####### Percent Diet FS ~ Pond FS  #######
meta_diet_biomass_full = meta_diet_biomass_full %>% group_by(Sample_ID_Year) %>% mutate(Percent_Diet = Boat_Net/sum(Boat_Net))
meta_diet_biomass_full = subset(meta_diet_biomass_full, Sample_ID != 94)

meta_diet_biomass_FS_Percent_Diet = subset(meta_diet_biomass_full, Diet_Category == "Anostraca")

meta_diet_biomass_FS = merge(x = meta_diet_biomass_FS, y = meta_diet_biomass_FS_Percent_Diet[,c("Sample_ID_Year", "Percent_Diet")])
names(meta_diet_biomass_FS)[names(meta_diet_biomass_FS) == 'Percent_Diet'] <- 'Percent_Diet_FS'
meta_diet_biomass_FS = meta_diet_biomass_FS %>% relocate(Percent_Diet_FS, .after = log_Boat_Net)

#write.csv(meta_diet_biomass_FS, "Data/Meta_Percent_FS_Diet.csv")

#code for sex-specific graph. switch out to "Male" for male plots.
meta_diet_biomass_FS = subset(meta_diet_biomass_FS, SEX == "Female") 

#average FS in gut over all individuals per pond
meta_diet_biomass_FS = meta_diet_biomass_FS %>% group_by(Date_Sample_GF, Pond) %>% mutate(Percent_Diet_FS_Avg = mean(Percent_Diet_FS))
meta_diet_biomass_FS = meta_diet_biomass_FS %>% group_by(Date_Sample_GF, Pond) %>% mutate(Percent_Diet_FS_SE = sd(Percent_Diet_FS)/sqrt(n()))
meta_diet_biomass_FS$Percent_Diet_FS_lower = meta_diet_biomass_FS$Percent_Diet_FS_Avg - meta_diet_biomass_FS$Percent_Diet_FS_SE*1.96
meta_diet_biomass_FS$Percent_Diet_FS_upper = meta_diet_biomass_FS$Percent_Diet_FS_Avg + meta_diet_biomass_FS$Percent_Diet_FS_SE*1.96

#bound 95% CI (0,1)
meta_diet_biomass_FS$Percent_Diet_FS_upper[meta_diet_biomass_FS$Percent_Diet_FS_upper > 1] <- 1
meta_diet_biomass_FS$Percent_Diet_FS_lower[meta_diet_biomass_FS$Percent_Diet_FS_lower < 0] <- 0


meta_diet_biomass_FS_avg_percent = meta_diet_biomass_FS[,c("Date_Sample_GF", "Year", "Pond", "FS_Biomass_Pond_m3", "log_FS_Biomass_Pond_m3",
                                                   "SE_Biomass", "lower_biomass", "upper_biomass",
                                                   "Percent_Diet_FS_Avg", "Percent_Diet_FS_SE", "Percent_Diet_FS_lower",
                                                   "Percent_Diet_FS_upper")]

meta_diet_biomass_FS_avg_percent = unique(meta_diet_biomass_FS_avg_percent)
  
#meta_diet_biomass_FS_tran = merge(meta_diet_biomass_FS_tran, meta_diet_biomass_FS[,c("Sample_ID_Year", "Percent_Diet_FS")], by = "Sample_ID_Year")


#remove one point (which represents two metas) that had FS in their gut despite no FS present
#during sampling. Probably a result of sampling for FS a day after gut flush. So the metas ate
#the last remaining FS or the temperature in the pond reached the CT max for FS when we sampled. 
#Metas could have moved from another pond where they ate FS, but seems unlikely given that both 
#were caught in pond 13 previously
meta_diet_biomass_FS_avg_percent = subset(meta_diet_biomass_FS_avg_percent, !(Year == "2021" & Pond == "13"))


##build logistic curve before plotting

mod2 <- drm(Percent_Diet_FS_Avg ~ FS_Biomass_Pond_m3, data = meta_diet_biomass_FS_avg_percent, 
            fct=LL.4(fixed=c(-0.7, 0, 1, NA),
                     names = c("Slope", "Lower Limit", "Upper Limit", "Midpoint")), type = "continuous")

#female: fixed=c(-1.5, 0, 1, NA)
#male: fixed=c(-0.7, 0, 1, NA)

summary(mod2)
plot(mod2)


#create new data with 100 points with min and max bounded by observed values 
newdata = seq(min(meta_diet_biomass_FS_avg_percent$FS_Biomass_Pond_m3), max(meta_diet_biomass_FS_avg_percent$FS_Biomass_Pond_m3) + 500000, length.out = 500000)
newdata = as.data.frame(newdata)
colnames(newdata) = c("FS_Biomass_Pond_m3")
fit2_predict = predict(mod2, type="response", se.fit = F, newdata=newdata)
fit2_predict = as.data.frame(fit2_predict)
fit2_predict = cbind(fit2_predict, newdata)
colnames(fit2_predict) = c("Percent_Diet_FS_Avg", "FS_Biomass_Pond_m3")
plot(log(fit2_predict$FS_Biomass_Pond_m3 + 1), fit2_predict$Percent_Diet_FS_Avg)


percent_fs_diet_plot =
  ggplot(meta_diet_biomass_FS_avg_percent, aes(x=log_FS_Biomass_Pond_m3, y=Percent_Diet_FS_Avg)) + 
  #geom_vline(xintercept = log(64), linetype="dashed", color = "black", size=3) + #these are to show microcosm experiment densities
  #geom_vline(xintercept = log(637), linetype="dashed", color = "black", size=3) +
  #geom_vline(xintercept = log(1599), linetype="dashed", color = "black", size=3) +
  #geom_vline(xintercept = log(3187), linetype="dashed", color = "black", size=3) +
  geom_errorbarh(aes(xmin=lower_biomass, xmax=upper_biomass), color = "black", height = 0.02, size = 1.5) +
  geom_errorbar(aes(ymin=Percent_Diet_FS_lower, ymax=Percent_Diet_FS_upper), color = "black", width = 0.15, size = 1.5) +
  geom_point(aes(x=log_FS_Biomass_Pond_m3, y=Percent_Diet_FS_Avg, color = Pond), size=10, alpha = 0.9) +
  geom_line(data = fit2_predict, aes(x = log(FS_Biomass_Pond_m3 + 1), y = Percent_Diet_FS_Avg), size = 3) +
  coord_cartesian(xlim = c(log(1),log(300000)), ylim=c(0,1)) +
  ylab("Percent Contribution of\n Fairy Shrimp to Diet ± 95% CI") +
  xlab(expression(paste("Fairy Shrimp Biomass in Pond", " ", "(mg/", m^{3}, ")", " ", "± 95% CI"))) +
  scale_x_continuous(breaks=c(log(1),log(10),log(100),log(1000),log(10000), log(100000)), labels=c("0", "10", "100", "1,000", "10,000", "100,000")) +
  scale_color_manual(values = c("#66C2A5", "#A6D854", "#8DA0CB", "#E78AC3", "#FFD92F", "#FC8D62", "#E5C494", "#B3B3B3")) +
  theme_bw(42) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)))

print(percent_fs_diet_plot)

ggsave(percent_fs_diet_plot, filename = "Outputs/percent_fs_diet_plot.png",  width = 18, height = 12, dpi = 100)


####### Community Diet Analyses #######

####### (i) Bipartite Plot #######


names(meta_diet_biomass_full)[names(meta_diet_biomass_full) == 'Boat_Net'] <- 'Biomass_mg'

#merge in metadata for gut samples and get data into correct format
meta_diet_biomass_full = merge(meta_diet_biomass_full, unique(salamander[,c("Salamander_ID", "Date_Sample", "Pond", "COHORT", "SEX", "Hydroperiod")]), 
           by = c("Salamander_ID", "Date_Sample"), all.x = T)

names(meta_diet_biomass_full)[names(meta_diet_biomass_full) == 'COHORT'] <- 'Cohort'
names(meta_diet_biomass_full)[names(meta_diet_biomass_full) == 'SEX'] <- 'Sex'

meta_diet_biomass_full$Date_Sample_Block <- ifelse(meta_diet_biomass_full[,c("Sample_ID")] < 36, "Early", 
                                            ifelse(meta_diet_biomass_full[,c("Sample_ID")] > 35 & meta_diet_biomass_full[,c("Sample_ID")] < 57, "Late", 
                                            ifelse(meta_diet_biomass_full[,c("Sample_ID")] > 56 & meta_diet_biomass_full[,c("Sample_ID")] < 80, "Early", 
                                            ifelse(meta_diet_biomass_full[,c("Sample_ID")] > 79, "Late", "ERROR"))))

meta_diet_biomass_full$Pond = as.factor(meta_diet_biomass_full$Pond)
meta_diet_biomass_full$Year = as.factor(meta_diet_biomass_full$Year)
meta_diet_biomass_full$Sex = as.factor(meta_diet_biomass_full$Sex)
meta_diet_biomass_full$Hydroperiod = as.factor(meta_diet_biomass_full$Hydroperiod)
meta_diet_biomass_full$Date_Sample_Block = as.factor(meta_diet_biomass_full$Date_Sample_Block)

meta_diet_biomass_full = meta_diet_biomass_full[,c("Sample_ID", "Date_Sample", "Date_Sample_Block", 
                                                   "Year", "Pond", "Hydroperiod", "Salamander_ID",
                                                   "Cohort", "Sex", "Diet_Category", "Biomass_mg")]

meta_diet_biomass_full = meta_diet_biomass_full %>% arrange(Sample_ID, Diet_Category)

#drop two samples (sample_ID 90 & 94) that were not included in the salamander_master_cleaned df
#bc they were from pond 56 (90), which was no FS data, and didn't have a SVL measurement (90)
meta_diet_biomass_full = subset(meta_diet_biomass_full, !is.na(Pond))

#in Sex change "Female?" to "Female" 
meta_diet_biomass_full$Sex = car::recode(meta_diet_biomass_full$Sex,"'Female?'='Female'")

#categorize Cohort data (including NAs)
meta_diet_biomass_full$Cohort[is.na(meta_diet_biomass_full$Cohort)] <- "Unknown"
meta_diet_biomass_full$Cohort = as.factor(meta_diet_biomass_full$Cohort)
meta_diet_biomass_full$Cohort_Class = with(meta_diet_biomass_full,
                                      ifelse(Cohort %in% c("1996", "1998", "1999", "2000", "2001"), "Oldest", 
                                      ifelse(Cohort %in% c("2007", "2008", "2009", "2010", "2011"), "Middle", 
                                      ifelse(Cohort %in% c("2012", "2013", "2014", "2015", "2020"), "Youngest",
                                      ifelse(Cohort %in% c("Unknown"), "Unknown",
                                      "ERROR")))))

meta_diet_biomass_full$Cohort_Class = as.factor(meta_diet_biomass_full$Cohort_Class)
meta_diet_biomass_full[,c("Cohort_Class")] <- factor(meta_diet_biomass_full[,c("Cohort_Class")], levels=c("Youngest", "Middle", "Oldest", "Unknown"))
meta_diet_biomass_full = meta_diet_biomass_full %>% relocate(Cohort_Class, .after = Cohort)


#any NA's in df?
apply(meta_diet_biomass_full, 2, function(x) any(is.na(x) | is.infinite(x))) #no

#don't need to transform biomass data bc NMDS uses species ranks (not abundance) and their 
#variation among sites.

#spread data into wide format
meta_diet_biomass_ord = meta_diet_biomass_full %>% spread(Diet_Category, Biomass_mg)


#make column with time period variable for bipartite graph
meta_diet_biomass_ord = meta_diet_biomass_ord %>% arrange(Date_Sample, Pond)

time_period <- c(rep("Mid July 2021 Pond 6", 6), rep("Mid July 2021 Pond 8", 6), rep("Mid July 2021 Pond 10", 4), 
                 rep("Mid July 2021 Pond 11", 6), rep("Mid July 2021 Pond 13", 2), rep("Mid July 2021 Pond 15", 1),
                 rep("Mid July 2021 Pond 51", 1), rep("Mid July 2021 Pond 52", 9), rep("Late July 2021 Pond 6", 1),
                 rep("Late July 2021 Pond 8", 5), rep("Late July 2021 Pond 10", 1), rep("Late July 2021 Pond 11", 5),
                 rep("Late July 2021 Pond 15", 2), rep("Late July 2021 Pond 52", 7), rep("Mid July 2022 Pond 8", 21),
                 rep("Mid July 2022 Pond 15", 2),  rep("Late July 2022 Pond 8", 4), rep("Late July 2022 Pond 10", 1),
                 rep("Late July 2022 Pond 11", 1), rep("Late July 2022 Pond 13", 3), rep("Late July 2022 Pond 15", 4),
                 rep("Late July 2022 Pond 52", 3))

meta_diet_biomass_ord$Time_Period = time_period
meta_diet_biomass_ord$Time_Period = as.factor(meta_diet_biomass_ord$Time_Period)
meta_diet_biomass_ord = meta_diet_biomass_ord %>% relocate(Time_Period, .after = Date_Sample)
meta_diet_biomass_ord$Elevation = with(meta_diet_biomass_ord,
                                     ifelse(Pond %in% c("6", "8", "10", "11", "13", "15"), "Lower Cut", 
                                            "Upper Cut"))
meta_diet_biomass_ord$Elevation = as.factor(meta_diet_biomass_ord$Elevation)
meta_diet_biomass_ord = meta_diet_biomass_ord %>% relocate(Elevation, .after = Pond)


meta_diet_biomass_gather = meta_diet_biomass_ord %>% gather(key = Diet_Category, value = Biomass_mg, Acari:Trichoptera)
meta_diet_biomass_gather = meta_diet_biomass_gather %>% group_by(Time_Period, Diet_Category) %>% summarise(Total_Biomass_mg = sum(Biomass_mg))
meta_diet_biomass_gather = meta_diet_biomass_gather %>% group_by(Time_Period) %>% mutate(Biomass_Percent = Total_Biomass_mg / sum(Total_Biomass_mg))
meta_diet_biomass_gather$Total_Biomass_mg <- NULL

meta_diet_biomass_gather = meta_diet_biomass_gather %>% spread(Diet_Category, Biomass_Percent)

meta_diet_biomass_gather$Time_Period <- factor(meta_diet_biomass_gather$Time_Period, levels = c("Mid July 2021 Pond 6", "Mid July 2021 Pond 8", "Mid July 2021 Pond 10", 
                                                    "Mid July 2021 Pond 11", "Mid July 2021 Pond 13", "Mid July 2021 Pond 15",
                                                    "Mid July 2021 Pond 51", "Mid July 2021 Pond 52", "Late July 2021 Pond 6",
                                                    "Late July 2021 Pond 8", "Late July 2021 Pond 10", "Late July 2021 Pond 11",
                                                    "Late July 2021 Pond 15", "Late July 2021 Pond 52", "Mid July 2022 Pond 8",
                                                    "Mid July 2022 Pond 15", "Late July 2022 Pond 8", "Late July 2022 Pond 10",
                                                    "Late July 2022 Pond 11", "Late July 2022 Pond 13", "Late July 2022 Pond 15",
                                                    "Late July 2022 Pond 52"))

meta_diet_biomass_gather = meta_diet_biomass_gather %>% arrange(Time_Period)
meta_diet_biomass_gather = as.data.frame(meta_diet_biomass_gather)
rownames(meta_diet_biomass_gather) <- meta_diet_biomass_gather[,c("Time_Period")]
meta_diet_biomass_gather$Time_Period <- NULL

meta_diet_biomass_gather = t(meta_diet_biomass_gather)
meta_diet_biomass_gather = as.data.frame(meta_diet_biomass_gather)

plotweb(meta_diet_biomass_gather, method = "normal", labsize = 0.5, 
        col.low = c("#042333ff", 
        "#13306dff", "#253582ff", "#403891ff", "#6b4596ff", "#90548bff", "#b8627dff", "#de7065ff",
        "#f68f46ff", "#f9b641ff", "#efe350ff", "#ccebc5"),
        col.high=c(rep("grey5",8),rep("grey35",6),rep("grey65",2),rep("grey95",6)),  
        bor.col.interaction ="black", bor.col.high="black", bor.col.low="black", 
        text.rot=90, text.high.col=c("black"), text.low.col="black")

#make matrices
meta_diet_comm = meta_diet_biomass_ord[,13:ncol(meta_diet_biomass_ord)]
meta_diet_env = meta_diet_biomass_ord[,c(4:8,11:12)] #environmental variables that will add to plot later

#turn abundance data frame into a matrix
matrix_diet = as.matrix(meta_diet_comm)

#plot bray-curtis & jaccard nmds for date_sample_block and sex

####### (ii) Bray-Curits Sex#####
#https://jkzorz.github.io/2019/07/02/Indicator-species-analysis.html
#https://jkzorz.github.io/2019/06/11/ANOSIM-test.html

#good nmds explanation: https://jonlefcheck.net/2012/10/24/nmds-tutorial-in-r/

#nmds is iterative and will give slightly different results each time it is run on the same dataset

meta_nmds_bc_sex = metaMDS(matrix_diet, distance = "bray", binary = FALSE)
meta_nmds_bc_sex

plot(meta_nmds_bc_sex)

#add environmental information
en = envfit(meta_nmds_bc_sex, meta_diet_env, permutations = 999, na.rm = TRUE)

en

#extract NMDS scores (x and y coordinates), so can make a better plot in ggplot2
data.scores.bc.sex = as.data.frame(scores(meta_nmds_bc_sex)$sites)
species.scores.bc.sex <- as.data.frame(scores(meta_nmds_bc_sex)$species)

#add columns to data frame 
data.scores.bc.sex$Year = meta_diet_biomass_ord$Year
data.scores.bc.sex$Date_Sample_Block = meta_diet_biomass_ord$Date_Sample_Block
data.scores.bc.sex$Hydroperiod = meta_diet_biomass_ord$Hydroperiod
data.scores.bc.sex$Elevation = meta_diet_biomass_ord$Elevation
data.scores.bc.sex$Sex = meta_diet_biomass_ord$Sex
data.scores.bc.sex$Cohort_Class = meta_diet_biomass_ord$Cohort_Class
data.scores.bc.sex$Salamander_ID = meta_diet_biomass_ord$Salamander_ID
#pond won't work in the modeling part below -- too many levels leads to overfit model


head(data.scores.bc.sex)

#extract information from the envfit result
en_coord_cat = as.data.frame(scores(en, "factors"))


#plot (switch out grouping variable to generate other plots in manuscript [year, time period, bench, hydroperiod, and cohort])
NMDS = data.frame(NMDS1 = meta_nmds_bc_sex$points[,1], NMDS2=meta_nmds_bc_sex$points[,2],group=meta_diet_env$Sex)
NMDS.mean=aggregate(NMDS[,1:2],list(group=NMDS$group),mean)

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

plot.new()

ord<-ordiellipse(meta_nmds_bc_sex, meta_diet_env$Sex, display = "sites", 
                 kind = "sd", conf = .95, label = T)

df_ell <- data.frame()
for(g in levels(NMDS$group)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$group==g,],
                                                   veganCovEllipse(ord[[g]]$cov,ord[[g]]$center)))
                                ,group=g))
}

bc_sex_plot = 
  ggplot(data = data.scores.bc.sex, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores.bc.sex, aes(colour = Sex), size = 10, alpha = 0.5) + 
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,colour=group), size=3, linetype=1)+
  geom_text(data = species.scores.bc.sex, aes(x = NMDS1, y = NMDS2), colour = "darkblue", 
            fontface = "bold", label = row.names(species.scores.bc.sex), size = 8) + 
  #annotate("text",x=NMDS.mean$NMDS1,y=NMDS.mean$NMDS2,label=NMDS.mean$group)+
  scale_colour_manual(values = c("orange", "steelblue"))  + 
  #scale_colour_manual(values = c("orange", "steelblue", "#87CA6A", "#DF92B6"))  + 
  #scale_colour_manual(values = c("grey65"))  + 
  #scale_x_continuous(limits = c(-7,6))+
  #scale_y_continuous(limits = c(-2,2))+
  theme(
        axis.title = element_text(size = 40, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30", size = 3), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 30, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 25, colour = "grey30"),
        axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0))) + 
  labs(colour = "Sex")

bc_sex_plot

ggsave(bc_sex_plot, filename = "Outputs/bc_sex_plot.png",  width = 16, height = 12, dpi = 100)

#stats
anosim_bc_sex = anosim(matrix_diet, meta_diet_biomass_ord$Sex, distance = "bray", permutations = 10000)
anosim_bc_sex

Sex = meta_diet_biomass_ord$Sex

inv_bc = multipatt(meta_diet_comm, Sex, func = "r.g", control = how(nperm=10000))

summary(inv_bc, alpha = 1)

inv_bc$sign


####### (ii) Pie Chart#####

#overall contribution of each category to total diet biomass across all metas

total_diet_biomass = meta_diet_biomass_full %>% group_by(Diet_Category) %>% summarise(Total_Biomass = sum(Biomass_mg))
total_diet_biomass = total_diet_biomass %>% arrange(Diet_Category)
total_diet_biomass = as.data.frame(total_diet_biomass)
rownames(total_diet_biomass) <- total_diet_biomass[,1]

png("pie_biomass_plot.png", bg = "transparent")

pie(total_diet_biomass$Total_Biomass, init.angle=90, clockwise = TRUE, labels = "",
  col=c("#042333ff", "#13306dff", "#253582ff", "#403891ff", "#6b4596ff", "#90548bff", "#b8627dff", 
        "#de7065ff", "#f68f46ff", "#f9b641ff", "#efe350ff", "#ccebc5"))

dev.off()

#calculate mean and range for a) total stomach sample biomasses and b) observed sample richness
#a)
meta_diet_biomass_full$Sample_ID = as.factor(meta_diet_biomass_full$Sample_ID)
total_diet_biomass_samples = meta_diet_biomass_full %>% group_by(Sample_ID) %>% summarise(Total_Biomass = sum(Biomass_mg))

mean(total_diet_biomass_samples$Total_Biomass)
sd(total_diet_biomass_samples$Total_Biomass)
range(total_diet_biomass_samples$Total_Biomass)

#b) 
diet_richness = subset(meta_diet_biomass_full, Biomass_mg > 0)
diet_richness = diet_richness %>% group_by(Sample_ID) %>% summarise(Richness = n())

mean(diet_richness$Richness)
sd(diet_richness$Richness)
range(diet_richness$Richness)


##### Repeat Metas ####
#need to see if diets of repeat metas are correlated

salamander_recap_diet = data.scores.bc.sex %>% group_by(Salamander_ID) %>% filter(n()>1)

#two groups with three observations. For each, drop the one 2021 obs and keep the two 2022 obs
salamander_recap_diet = subset(salamander_recap_diet, !(Salamander_ID == "900118000002813" & Year == "2021") &
                                 !(Salamander_ID == "900118000002907" & Year == "2021"))

#NMDS x-values
salamander_recap_diet_x = salamander_recap_diet[,c("NMDS1", "Salamander_ID")]

salamander_recap_diet_x1 = salamander_recap_diet_x %>% group_by(Salamander_ID) %>% filter(row_number()==1)
salamander_recap_diet_x2 = salamander_recap_diet_x %>% group_by(Salamander_ID) %>% filter(row_number()==n())

salamander_recap_diet_x = merge(salamander_recap_diet_x1, salamander_recap_diet_x2, by = "Salamander_ID")


#NMDS y-values
salamander_recap_diet_y = salamander_recap_diet[,c("NMDS2", "Salamander_ID")]

salamander_recap_diet_y1 = salamander_recap_diet_y %>% group_by(Salamander_ID) %>% filter(row_number()==1)
salamander_recap_diet_y2 = salamander_recap_diet_y %>% group_by(Salamander_ID) %>% filter(row_number()==n())

salamander_recap_diet_y = merge(salamander_recap_diet_y1, salamander_recap_diet_y2, by = "Salamander_ID")

#plot
meta_recap_corr = 
  ggplot(data = salamander_recap_diet_y, aes(x= NMDS2.x, y = NMDS2.y)) +
  geom_point(color = "#FC766AFF", size = 7, alpha = 0.8) +
  geom_point(data = salamander_recap_diet_x, aes(x= NMDS1.x, y = NMDS1.y), color = "#5B84B1FF", size = 7, alpha = 0.8) +
  geom_abline(slope=1, intercept=0, lwd = 2) +
  scale_y_continuous(limits = c(-3, 3), breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
  scale_x_continuous(limits = c(-3, 3), breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
  labs(x = "First Capture", y = "Second Capture") +
  theme_bw(42) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0)))

meta_recap_corr

ggsave(meta_recap_corr, filename = "Outputs/meta_recap_corr.png",  width = 12, height = 12, dpi = 100)
  
cor(salamander_recap_diet_x$NMDS1.x, salamander_recap_diet_x$NMDS1.y)
cor(salamander_recap_diet_y$NMDS2.x, salamander_recap_diet_y$NMDS2.y)


#### Meta Diet & Mvmt ####

#need to assess if metas were previously captured in same pond on date before gut flush

all_captures = salamander[,c("Salamander_ID", "Date_Sample", "Pond")]
all_captures$Pond = as.factor(all_captures$Pond)

gut_flush_captures = meta_diet_biomass_FS[,c("Salamander_ID", "Date_Sample_GF", "Pond")]
names(gut_flush_captures)[names(gut_flush_captures) == 'Date_Sample_GF'] <- 'Date_Sample'

gf_metas = merge(all_captures, gut_flush_captures, by = "Salamander_ID", all.y = T)
gf_metas[,c("Date_Sample.y", "Pond.y")] <- NULL
colnames(gf_metas) = c("Salamander_ID", "Date_Sample", "Pond")
