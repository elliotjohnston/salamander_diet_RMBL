
##### Metamorph Foraging Experiment ######

rm(list = ls())
require(lubridate)
require(reshape2)
require(dplyr)
require(ggplot2)
require(scales)
require(stringr)
require(lme4)
require(emmeans)


##import and format data
meta_diet_choice <- read.csv(file = "Data/Foraging_Experiment.csv", stringsAsFactors = TRUE)
meta_fs_rate <- read.csv(file = "Data/Metamorph_Foraging_Trials.csv", stringsAsFactors = FALSE)
salamander <- read.csv(file = "Data/Salamander_Master_cleaned.csv", stringsAsFactors = TRUE)

meta_diet_choice = subset(meta_diet_choice, Date != "XX") #did this so Salamander_ID wouldn't import as scientific notation
meta_diet_choice = droplevels(meta_diet_choice)
meta_diet_choice$FS = as.numeric(as.character(meta_diet_choice$FS))
meta_diet_choice$Mosquito = as.numeric(as.character(meta_diet_choice$Mosquito))
meta_diet_choice$Copepod = as.numeric(as.character(meta_diet_choice$Copepod))
meta_diet_choice$Asynarchus = as.numeric(as.character(meta_diet_choice$Asynarchus))

meta_fs_rate$Treatment = as.factor(meta_fs_rate$Treatment)
meta_fs_rate$Replicate = as.factor(meta_fs_rate$Replicate)
meta_fs_rate$Consumption_Time = as_datetime(ms(meta_fs_rate$Consumption_Time))

salamander$X <- NULL
names(salamander)[names(salamander) == "INDIV_ID"] <- "Salamander_ID"
salamander = subset(salamander, DATEOFYEAR == "7/11/22")
salamander$Condition = (salamander$WEIGHT / salamander$SVL) / 1000

#### Meta Diet Choice Plots ####
#merge in meta demographic information
meta_diet_choice = merge(x = meta_diet_choice, y = salamander[,c("Salamander_ID", "COHORT", "SEX", "SVL", "WEIGHT", "Condition")], by = "Salamander_ID", all.x = T)
meta_diet_choice = meta_diet_choice %>% relocate("Salamander_ID", .after = "Notes")

#drop Time and Notes columns
meta_diet_choice = meta_diet_choice[ , -which(names(meta_diet_choice) %in% c("Time", "Notes"))] 

#melt df
meta_diet_choice_melt = meta_diet_choice[ , -which(names(meta_diet_choice) %in% c("Salamander_ID", "COHORT", "SEX", "SVL", "WEIGHT", "Condition"))] 
meta_diet_choice_melt = melt(meta_diet_choice_melt, value.name = "Count")
names(meta_diet_choice_melt)[names(meta_diet_choice_melt) == 'variable'] <- 'Diet_Taxon'

#summarize data (mean, SE)
meta_diet_means = meta_diet_choice_melt %>% group_by(Date, Treatment, Diet_Taxon) %>% dplyr::summarize(Mean_Count = mean(Count))
meta_diet_SE = meta_diet_choice_melt %>% group_by(Date, Treatment, Diet_Taxon) %>% dplyr::summarize(SE_Count = sd(Count)/sqrt(n()))

meta_diet_means = merge(meta_diet_means, meta_diet_SE, by = c("Date", "Treatment", "Diet_Taxon"))
meta_diet_means$lower = meta_diet_means$Mean_Count - meta_diet_means$SE_Count
meta_diet_means$upper = meta_diet_means$Mean_Count + meta_diet_means$SE_Count

rm(meta_diet_SE)

meta_diet_means$Diet_Taxon <- factor(meta_diet_means$Diet_Taxon, levels = c("FS", "Mosquito", "Copepod", "Asynarchus"))


#plot - PO
meta_diet_means_PO = subset(meta_diet_means, Treatment == "PO")

hline_PO <- data.frame(Diet_Taxon=c("FS", "Mosquito", "Copepod", "Asynarchus"), Mean_Count=c(25, 10, 10, 10))
hline_PO$Diet_Taxon <- factor(hline_PO$Diet_Taxon, levels = c("FS", "Mosquito", "Copepod", "Asynarchus"))


diet_choice_PO = 
ggplot(meta_diet_means_PO, aes(x=Diet_Taxon, y=Mean_Count, shape = Date)) + 
geom_point(data=hline_PO, aes(x=Diet_Taxon, y=Mean_Count), shape=95, size=90, color = "blue") +
geom_errorbar(aes(ymin=lower, ymax=upper), width=0.75, size=6, color = "black", position = position_dodge(1)) +
geom_point(size=20, stroke = 2, position = position_dodge(1)) +
scale_y_continuous(breaks=c(0,5,10,15,20,25), limits=c(0,25)) +
scale_x_discrete(labels=c("FS", "MOSQ", "COPE", "CADDIS")) +
theme(legend.position = "none",
      axis.title = element_blank(),
      axis.text = element_text(size = 50),
      axis.ticks = element_line(size = 2),
      axis.ticks.length = unit(.25, "cm"),
      axis.text.x=element_text(margin=margin(t=10)),
      axis.text.y=element_text(margin=margin(r=10)))

diet_choice_PO

ggsave(diet_choice_PO, filename = "Outputs/meta_diet_choice_PO.png",  width = 18, height = 12, dpi = 100)


#plot - NFS
hline_NFS <- data.frame(Diet_Taxon=c("FS", "Mosquito", "Copepod", "Asynarchus"), Mean_Count=c(0, 10, 10, 10))
hline_NFS$Diet_Taxon <- factor(hline_NFS$Diet_Taxon, levels = c("FS", "Mosquito", "Copepod", "Asynarchus"))

meta_diet_means_NFS = subset(meta_diet_means, Treatment == "N FS")

diet_choice_NFS = 
  ggplot(meta_diet_means_NFS, aes(x=Diet_Taxon, y=Mean_Count, shape = Date)) + 
  geom_point(data=hline_NFS, aes(x=Diet_Taxon, y=Mean_Count), shape=95, size=90, color = "blue") +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.75, size=6, color = "black", position = position_dodge(1)) +
  geom_point(size=20, stroke = 2, position = position_dodge(1)) +
  scale_y_continuous(breaks=c(0,5,10,15,20,25), limits=c(0,25)) +
  scale_x_discrete(labels=c("FS", "MOSQ", "COPE", "CADDIS")) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_text(size = 50),
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x=element_text(margin=margin(t=10)),
        axis.text.y=element_text(margin=margin(r=10)))

diet_choice_NFS

ggsave(diet_choice_NFS, filename = "Outputs/meta_diet_choice_NFS.png",  width = 18, height = 12, dpi = 100)


#plot - LFS
hline_LFS <- data.frame(Diet_Taxon=c("FS", "Mosquito", "Copepod", "Asynarchus"), Mean_Count=c(10, 10, 10, 10))
hline_LFS$Diet_Taxon <- factor(hline_LFS$Diet_Taxon, levels = c("FS", "Mosquito", "Copepod", "Asynarchus"))

meta_diet_means_LFS = subset(meta_diet_means, Treatment == "L FS")

diet_choice_LFS = 
  ggplot(meta_diet_means_LFS, aes(x=Diet_Taxon, y=Mean_Count, shape = Date)) + 
  geom_point(data=hline_LFS, aes(x=Diet_Taxon, y=Mean_Count), shape=95, size=90, color = "blue") +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.75, size=6, color = "black", position = position_dodge(1)) +
  geom_point(size=20, stroke = 2, position = position_dodge(1)) +
  scale_y_continuous(breaks=c(0,5,10,15,20,25), limits=c(0,25)) +
  scale_x_discrete(labels=c("FS", "MOSQ", "COPE", "CADDIS")) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_text(size = 50),
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x=element_text(margin=margin(t=10)),
        axis.text.y=element_text(margin=margin(r=10)))

diet_choice_LFS

ggsave(diet_choice_LFS, filename = "Outputs/meta_diet_choice_LFS.png",  width = 18, height = 12, dpi = 100)


#plot - MFS
hline_MFS <- data.frame(Diet_Taxon=c("FS", "Mosquito", "Copepod", "Asynarchus"), Mean_Count=c(25, 10, 10, 10))
hline_MFS$Diet_Taxon <- factor(hline_MFS$Diet_Taxon, levels = c("FS", "Mosquito", "Copepod", "Asynarchus"))

meta_diet_means_MFS = subset(meta_diet_means, Treatment == "M FS")

diet_choice_MFS = 
  ggplot(meta_diet_means_MFS, aes(x=Diet_Taxon, y=Mean_Count, shape = Date)) + 
  geom_point(data=hline_MFS, aes(x=Diet_Taxon, y=Mean_Count), shape=95, size=90, color = "blue") +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.75, size=6, color = "black", position = position_dodge(1)) +
  geom_point(size=20, stroke = 2, position = position_dodge(1)) +
  scale_y_continuous(breaks=c(0,5,10,15,20,25), limits=c(0,25)) +
  scale_x_discrete(labels=c("FS", "MOSQ", "COPE", "CADDIS")) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_text(size = 50),
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x=element_text(margin=margin(t=10)),
        axis.text.y=element_text(margin=margin(r=10)))

diet_choice_MFS

ggsave(diet_choice_MFS, filename = "Outputs/meta_diet_choice_MFS.png",  width = 18, height = 12, dpi = 100)


#plot - HFS
hline_HFS <- data.frame(Diet_Taxon=c("FS", "Mosquito", "Copepod", "Asynarchus"), Mean_Count=c(30, 10, 10, 10)) #30 is a dummy value for 50
hline_HFS$Diet_Taxon <- factor(hline_HFS$Diet_Taxon, levels = c("FS", "Mosquito", "Copepod", "Asynarchus"))

meta_diet_means_HFS = subset(meta_diet_means, Treatment == "H FS")

diet_choice_HFS = 
  ggplot(meta_diet_means_HFS, aes(x=Diet_Taxon, y=Mean_Count, shape = Date)) + 
  geom_point(data=hline_HFS, aes(x=Diet_Taxon, y=Mean_Count), shape=95, size=90, color = "blue") +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.75, size=6, color = "black", position = position_dodge(1)) +
  geom_point(size=20, stroke = 2, position = position_dodge(1)) +
  scale_y_continuous(breaks=c(0,5,10,15,20,25,30), limits=c(0,30), labels = c("0", "5", "10", "15", "20", "25", "50")) +
  scale_x_discrete(labels=c("FS", "MOSQ", "COPE", "CADDIS")) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_text(size = 50),
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x=element_text(margin=margin(t=10)),
        axis.text.y=element_text(margin=margin(r=10)))

diet_choice_HFS

ggsave(diet_choice_HFS, filename = "Outputs/meta_diet_choice_HFS.png",  width = 18, height = 12, dpi = 100)


##### Diet Choice - Alt. Viz. ####
#divide plots by prey taxa rather than FS treatment

#do I need to recalculate SEs for this alternative approach...??

#mosq
meta_diet_mosq = subset(meta_diet_means, Diet_Taxon == "Mosquito")
meta_diet_mosq$Treatment <- factor(meta_diet_mosq$Treatment, levels = c("PO", "N FS", "L FS", "M FS", "H FS"))
hline_mosq <- data.frame(Treatment=c("PO", "N FS", "L FS", "M FS", "H FS"), Mean_Count=c(10, 10, 10, 10, 10))
hline_mosq$Treatment <- factor(hline_mosq$Treatment, levels = c("PO", "N FS", "L FS", "M FS", "H FS"))

meta_diet_mosq_plot = 
  ggplot(meta_diet_mosq, aes(x=Treatment, y=Mean_Count, shape = Date)) + 
  geom_point(data=hline_mosq, aes(x=Treatment, y=Mean_Count), shape=95, size=90, color = "blue") +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.75, size=6, color = "black", position = position_dodge(1)) +
  geom_point(size=20, stroke = 2, position = position_dodge(1)) +
  scale_y_continuous(breaks=c(0,2,4,6,8,10), limits=c(0,10), labels = c("0", "2", "4", "6", "8", "10")) +
  scale_x_discrete(labels=c("Control", "N FS", "L FS", "M FS", "H FS")) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_text(size = 70),
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x=element_text(margin=margin(t=10)),
        axis.text.y=element_text(margin=margin(r=10)))

meta_diet_mosq_plot

ggsave(meta_diet_mosq_plot, filename = "Outputs/meta_diet_choice_mosq.png",  width = 18, height = 12, dpi = 100)


#cope
meta_diet_cope = subset(meta_diet_means, Diet_Taxon == "Copepod")
meta_diet_cope$Treatment <- factor(meta_diet_cope$Treatment, levels = c("PO", "N FS", "L FS", "M FS", "H FS"))
hline_cope <- data.frame(Treatment=c("PO", "N FS", "L FS", "M FS", "H FS"), Mean_Count=c(10, 10, 10, 10, 10))
hline_cope$Treatment <- factor(hline_cope$Treatment, levels = c("PO", "N FS", "L FS", "M FS", "H FS"))

meta_diet_cope_plot = 
  ggplot(meta_diet_cope, aes(x=Treatment, y=Mean_Count, shape = Date)) + 
  geom_point(data=hline_cope, aes(x=Treatment, y=Mean_Count), shape=95, size=90, color = "blue") +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.75, size=6, color = "black", position = position_dodge(1)) +
  geom_point(size=20, stroke = 2, position = position_dodge(1)) +
  scale_y_continuous(breaks=c(0,2,4,6,8,10), limits=c(0,10), labels = c("0", "2", "4", "6", "8", "10")) +
  scale_x_discrete(labels=c("Control", "N FS", "L FS", "M FS", "H FS")) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_text(size = 70),
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x=element_text(margin=margin(t=10)),
        axis.text.y=element_text(margin=margin(r=10)))

meta_diet_cope_plot

ggsave(meta_diet_cope_plot, filename = "Outputs/meta_diet_choice_cope.png",  width = 18, height = 12, dpi = 100)


#caddis
meta_diet_caddis = subset(meta_diet_means, Diet_Taxon == "Asynarchus")
meta_diet_caddis$Treatment <- factor(meta_diet_caddis$Treatment, levels = c("PO", "N FS", "L FS", "M FS", "H FS"))
hline_caddis <- data.frame(Treatment=c("PO", "N FS", "L FS", "M FS", "H FS"), Mean_Count=c(10, 10, 10, 10, 10))
hline_caddis$Treatment <- factor(hline_caddis$Treatment, levels = c("PO", "N FS", "L FS", "M FS", "H FS"))

meta_diet_caddis_plot = 
  ggplot(meta_diet_caddis, aes(x=Treatment, y=Mean_Count, shape = Date)) + 
  geom_point(data=hline_caddis, aes(x=Treatment, y=Mean_Count), shape=95, size=90, color = "blue") +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.75, size=6, color = "black", position = position_dodge(1)) +
  geom_point(size=20, stroke = 2, position = position_dodge(1)) +
  scale_y_continuous(breaks=c(0,2,4,6,8,10), limits=c(0,10), labels = c("0", "2", "4", "6", "8", "10")) +
  scale_x_discrete(labels=c("Control", "N FS", "L FS", "M FS", "H FS")) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_text(size = 70),
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x=element_text(margin=margin(t=10)),
        axis.text.y=element_text(margin=margin(r=10)))

meta_diet_caddis_plot

ggsave(meta_diet_caddis_plot, filename = "Outputs/meta_diet_choice_caddis.png",  width = 18, height = 12, dpi = 100)


#fs
meta_diet_fs = subset(meta_diet_means, Diet_Taxon == "FS")
meta_diet_fs$Treatment <- factor(meta_diet_fs$Treatment, levels = c("PO", "N FS", "L FS", "M FS", "H FS"))
hline_fs <- data.frame(Treatment=c("PO", "N FS", "L FS", "M FS", "H FS"), Mean_Count=c(25, 0, 10, 25, 50))
hline_fs$Treatment <- factor(hline_fs$Treatment, levels = c("PO", "N FS", "L FS", "M FS", "H FS"))

meta_diet_fs_plot = 
  ggplot(meta_diet_fs, aes(x=Treatment, y=Mean_Count, shape = Date)) + 
  geom_point(data=hline_fs, aes(x=Treatment, y=Mean_Count), shape=95, size=90, color = "blue") +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.75, size=6, color = "black", position = position_dodge(1)) +
  geom_point(size=20, stroke = 2, position = position_dodge(1)) +
  #scale_y_continuous(breaks=c(0,2,4,6,8,10), limits=c(0,10), labels = c("0", "2", "4", "6", "8", "10")) +
  scale_x_discrete(labels=c("Control", "N FS", "L FS", "M FS", "H FS")) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_text(size = 70),
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x=element_text(margin=margin(t=10)),
        axis.text.y=element_text(margin=margin(r=10)))

meta_diet_fs_plot

ggsave(meta_diet_fs_plot, filename = "Outputs/meta_diet_choice_fs.png",  width = 18, height = 12, dpi = 100)




#### Meta FS Foraging Rate ####

xlimits <- as.POSIXct(strptime(c("1970-01-01 00:00:00", "1970-01-01 00:20:00"), 
                            format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))

xbreaks <- as.POSIXct(strptime(c("1970-01-01 00:00:00", "1970-01-01 00:05:00", 
                                 "1970-01-01 00:10:00","1970-01-01 00:15:00", 
                                 "1970-01-01 00:20:00"), 
                               format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))

#plot - N FS
meta_fs_rate_NFS = subset(meta_fs_rate, Treatment == "N FS")

meta_rate_NFS =
  ggplot(meta_fs_rate_NFS, aes(x=Consumption_Time, y=FS_Eaten)) + 
  geom_line(aes(group = Replicate), lwd = 2) +
  geom_point(aes(group = Replicate), size = 10, alpha = 0.8) +
  scale_y_continuous(breaks=c(0,2,4,6,8,10), limits=c(0,10)) +
  scale_x_datetime(breaks = xbreaks,
                   labels = c("0", "5", "10", "15", "20"),
                   limits = xlimits) +
  theme_bw(65) +
  theme(axis.title = element_blank())

meta_rate_NFS

ggsave(meta_rate_NFS, filename = "Outputs/meta_fs_rate_NFS.png",  width = 12, height = 12, dpi = 100)


#plot - L FS
meta_fs_rate_LFS = subset(meta_fs_rate, Treatment == "L FS")

meta_rate_LFS =
  ggplot(meta_fs_rate_LFS, aes(x=Consumption_Time, y=FS_Eaten)) + 
  geom_line(aes(group = Replicate), lwd = 2) +
  geom_point(aes(group = Replicate), size = 10, alpha = 0.8) +
  scale_y_continuous(breaks=c(0,2,4,6,8,10), limits=c(0,10)) +
  scale_x_datetime(breaks = xbreaks,
                   labels = c("0", "5", "10", "15", "20"),
                   limits = xlimits) +
  theme_bw(65) +
  theme(axis.title = element_blank())


meta_rate_LFS

ggsave(meta_rate_LFS, filename = "Outputs/meta_fs_rate_LFS.png",  width = 12, height = 12, dpi = 100)


#plot - M FS
meta_fs_rate_MFS = subset(meta_fs_rate, Treatment == "M FS")

meta_rate_MFS =
  ggplot(meta_fs_rate_MFS, aes(x=Consumption_Time, y=FS_Eaten)) + 
  geom_line(aes(group = Replicate), lwd = 2) +
  geom_point(aes(group = Replicate), size = 10, alpha = 0.8) +
  scale_y_continuous(breaks=c(0,2,4,6,8,10), limits=c(0,10)) +
  scale_x_datetime(breaks = xbreaks,
                   labels = c("0", "5", "10", "15", "20"),
                   limits = xlimits) +
  theme_bw(65) +
  theme(axis.title = element_blank())


meta_rate_MFS

ggsave(meta_rate_MFS, filename = "Outputs/meta_fs_rate_MFS.png",  width = 12, height = 12, dpi = 100)


#plot - H FS
meta_fs_rate_HFS = subset(meta_fs_rate, Treatment == "H FS")

meta_rate_HFS =
  ggplot(meta_fs_rate_HFS, aes(x=Consumption_Time, y=FS_Eaten)) + 
  geom_line(aes(group = Replicate), lwd = 2) +
  geom_point(aes(group = Replicate), size = 10, alpha = 0.8) +
  scale_y_continuous(breaks=c(0,2,4,6,8,10), limits=c(0,10)) +
  scale_x_datetime(breaks = xbreaks,
                   labels = c("0", "5", "10", "15", "20"),
                   limits = xlimits) +
  theme_bw(65) +
  theme(axis.title = element_blank())


meta_rate_HFS

ggsave(meta_rate_HFS, filename = "Outputs/meta_fs_rate_HFS.png",  width = 12, height = 12, dpi = 100)

##stats
meta_fs_rate_max = meta_fs_rate %>% group_by(Treatment, Replicate) %>% slice(which.max(FS_Eaten))

m1_meta_fs_rate = lm(FS_Eaten ~ Treatment, data = meta_fs_rate_max)
car::Anova(m1_meta_fs_rate)

emmeans(m1_meta_fs_rate, ~ Treatment)


#### Manly-Chesson Index ####

#format
manly_chesson = subset(meta_diet_choice_melt, Treatment != "PO" & !is.na(Count) & Date == "7/12/22")
manly_chesson$Date <- NULL
names(manly_chesson)[names(manly_chesson) == "Count"] <- "Count_End"

manly_chesson$Treatment <- factor(manly_chesson$Treatment, levels = c("N FS", "L FS", "M FS", "H FS"))
manly_chesson$Diet_Taxon <- factor(manly_chesson$Diet_Taxon, levels = c("FS", "Mosquito", "Copepod", "Asynarchus"))
manly_chesson$Replicate <- factor(manly_chesson$Replicate, levels = c("1", "2", "3", "4", "5"))

manly_chesson = arrange(manly_chesson, Diet_Taxon, Treatment, Replicate)

Count_Start <- c(rep("10", 5), rep("25", 5), rep("50", 5), rep("10", 60))
manly_chesson$Count_Start = Count_Start
manly_chesson$Count_Start = as.numeric(as.character(manly_chesson$Count_Start))
manly_chesson = manly_chesson %>% relocate("Count_Start", .before = "Count_End")
manly_chesson = arrange(manly_chesson, Treatment, Replicate, Diet_Taxon)
manly_chesson$Num_Consumed = manly_chesson$Count_Start - manly_chesson$Count_End

#calculate index

manly_chesson = manly_chesson %>% group_by(Treatment, Replicate) %>% mutate(Prop_Avail = Count_Start / sum(Count_Start))
manly_chesson = manly_chesson %>% group_by(Treatment, Replicate) %>% mutate(Prop_Used = Num_Consumed / sum(Num_Consumed))
manly_chesson = manly_chesson %>% group_by(Treatment, Replicate) %>% mutate(Used_Avail_Ratio = Prop_Used / Prop_Avail)
manly_chesson = manly_chesson %>% group_by(Treatment, Replicate) %>% mutate(Used_Avail_Ratio_Total = sum(Used_Avail_Ratio))
manly_chesson = manly_chesson %>% group_by(Treatment, Replicate) %>% mutate(Alpha = Used_Avail_Ratio/Used_Avail_Ratio_Total)
manly_chesson = droplevels(manly_chesson)

#average index value for each prey taxa in each treatment across replicates
manly_chesson_means = manly_chesson %>% group_by(Treatment, Diet_Taxon) %>% dplyr::summarize(Mean_Alpha = mean(Alpha))
manly_chesson_SE = manly_chesson %>% group_by(Treatment, Diet_Taxon) %>% dplyr::summarize(SE_Alpha = sd(Alpha)/sqrt(n()))

manly_chesson_means = merge(manly_chesson_means, manly_chesson_SE, by = c("Treatment", "Diet_Taxon"))
manly_chesson_means$CI_lower = manly_chesson_means$Mean_Alpha - manly_chesson_means$SE_Alpha*1.96
manly_chesson_means$CI_lower[manly_chesson_means$CI_lower < 0] <- 0
manly_chesson_means$CI_upper = manly_chesson_means$Mean_Alpha + manly_chesson_means$SE_Alpha*1.96

rm(manly_chesson_SE)

#no selection threshold equals 1/n. 1/3 for N FS treatment and 1/4 for other three treatments
mc_lfs = subset(manly_chesson_means, Treatment == "L FS")

mc_lfs_plot = 
  ggplot(mc_lfs, aes(x = Diet_Taxon, y = Mean_Alpha)) +
  geom_hline(yintercept=0.25, linetype="dashed", color = "black", size = 2.5) +
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), width=.2, size = 3) +
  geom_point(size = 12) +
  scale_y_continuous(limits = c(0, 0.8), breaks = c(0.0, 0.2, 0.4, 0.6, 0.8)) +
  scale_x_discrete(labels=c("FS", "MQ", "CP", "AS")) +
  theme_bw(65) +
  theme(axis.title = element_blank())

mc_lfs_plot

ggsave(mc_lfs_plot, filename = "Outputs/mc_lfs.png",  width = 12, height = 12, dpi = 100)

mc_mfs = subset(manly_chesson_means, Treatment == "M FS")

mc_mfs_plot = 
  ggplot(mc_mfs, aes(x = Diet_Taxon, y = Mean_Alpha)) +
  geom_hline(yintercept=0.25, linetype="dashed", color = "black", size = 2.5) +
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), width=.2, size = 3) +
  geom_point(size = 12) +
  scale_y_continuous(limits = c(0, 0.8), breaks = c(0.0, 0.2, 0.4, 0.6, 0.8)) +
  scale_x_discrete(labels=c("FS", "MQ", "CP", "AS")) +
  theme_bw(65) +
  theme(axis.title = element_blank())

mc_mfs_plot

ggsave(mc_mfs_plot, filename = "Outputs/mc_mfs.png",  width = 12, height = 12, dpi = 100)


mc_hfs = subset(manly_chesson_means, Treatment == "H FS")

mc_hfs_plot = 
  ggplot(mc_hfs, aes(x = Diet_Taxon, y = Mean_Alpha)) +
  geom_hline(yintercept=0.25, linetype="dashed", color = "black", size = 2.5) +
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), width=.2, size = 3) +
  geom_point(size = 12) +
  scale_y_continuous(limits = c(0, 0.8), breaks = c(0.0, 0.2, 0.4, 0.6, 0.8)) +
  scale_x_discrete(labels=c("FS", "MQ", "CP", "AS")) +
  theme_bw(65) +
  theme(axis.title = element_blank())

mc_hfs_plot

ggsave(mc_hfs_plot, filename = "Outputs/mc_hfs.png",  width = 12, height = 12, dpi = 100)


mc_nfs = subset(manly_chesson_means, Treatment == "N FS")

mc_nfs_plot = 
  ggplot(mc_nfs, aes(x = Diet_Taxon, y = Mean_Alpha)) +
  geom_hline(yintercept=0.33, linetype="dashed", color = "black", size = 2.5) +
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), width=.2, size = 3) +
  geom_point(size = 12) +
  scale_y_continuous(limits = c(0, 0.8), breaks = c(0.0, 0.2, 0.4, 0.6, 0.8)) +
  scale_x_discrete(labels=c("MQ", "CP", "AS")) +
  theme_bw(65) +
  theme(axis.title = element_blank())

mc_nfs_plot

ggsave(mc_nfs_plot, filename = "Outputs/mc_nfs.png",  width = 12, height = 12, dpi = 100)

##stats
manly_chesson$Box = str_c(manly_chesson$Treatment, "_", manly_chesson$Replicate)
manly_chesson$Box = as.factor(manly_chesson$Box)
manly_chesson = manly_chesson %>% relocate("Box", .after = "Replicate")


m1_mc = lmer(Alpha ~ Treatment*Diet_Taxon + (1|Box), data = manly_chesson) 
car::Anova(m1_mc)

emmeans_trt = emmeans(m1_mc, ~ Treatment)
pairs(emmeans_trt)

emmeans_prey = emmeans(m1_mc, ~ Diet_Taxon)
pairs(emmeans_prey)

emmeans_interaction = emmeans(m1_mc, ~ Treatment | Diet_Taxon)
pairs(emmeans_interaction)


manly_chesson_FS = subset(manly_chesson, Treatment != "N FS")
m1_mc_v2 = lmer(Alpha ~ Treatment*Diet_Taxon + (1|Box), data = manly_chesson_FS) 
car::Anova(m1_mc_v2)

emmeans_trt2 = emmeans(m1_mc_v2, ~ Treatment)
pairs(emmeans_trt2)

emmeans_prey2 = emmeans(m1_mc_v2, ~ Diet_Taxon)
pairs(emmeans_prey2)

emmeans_interaction2 = emmeans(m1_mc_v2, ~ Treatment | Diet_Taxon)
pairs(emmeans_interaction2)




