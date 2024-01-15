############################################### 
######### FS Length-Weight Regression #########
############################################### 

rm(list = ls())

library(ggplot2)

##import and format data
fs_lw <- read.csv(file = "Data/FS_LW_Regression.csv", stringsAsFactors = FALSE)

#change blank cells in Sex column to UNK
fs_lw$Sex <- sub("^$", "UNK", fs_lw$Sex) 

#remove rows that have a string in the Notes column -- these are all measurements that may be biased for various reasons
#no tail (5), short tail (1), Hard to tell where abdomen ends (6), full tin weight less than empty (1)
fs_lw = fs_lw[!grepl("a", fs_lw$Notes),] #choosing 'a' because it is contained in all of the strings 

#one outlier (obs 60). Identified as outlier by residuals vs leverage plot.
fs_lw = subset(fs_lw, Length_cm > 0.20)
#133 obs


#plot
fs_lw$Length_cm = fs_lw$Length_cm*10
names(fs_lw)[names(fs_lw) == 'Length_cm'] <- 'Length_mm'

lw_reg_curve_output = 
  ggplot(fs_lw, aes(x = log(Length_mm), y = log(Boat_Net))) +
  geom_point(aes(color = Sex), size = 6) +
  stat_smooth(method = "lm", formula = y ~ x, size = 3, color = "black") +
  ylab("log(Dry Weight [mg])") +
  xlab("log(Length [mm])") +
  theme_bw(42) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))

lw_reg_curve_output

ggsave(lw_reg_curve_output, filename = "Outputs/lw_reg_curve.png",  width = 14, height = 12, dpi = 100)

#stats
lw_model = lm(log(Boat_Net) ~ log(Length_mm), data = fs_lw)
#plot(lw_model)

summary(lw_model) 
#Multiple R-squared: 0.91, p-value < 0.001

#Formula: log(y) = 2.22706*log(x) - 3.69561

#check that this formula yields the correct curve
curve(2.22706*x - 3.69561, from=2, to=13, n=300, xlab="xvalue", ylab="yvalue")



