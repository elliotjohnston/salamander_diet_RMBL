#### Whiteman Graphs #####

library(ggplot2)

#Meta calories perm vs nonperm

pond_type = c("Permanent", "Nonpermanent")
data1_est = c(178, 525)
data1_high = c(231, 580)
data1_low = c(133, 475)

data1 = as.data.frame(cbind(pond_type, data1_est, data1_high, data1_low))
data1$pond_type = as.factor(data1$pond_type)
data1$pond_type = factor(data1$pond_type, levels = c("Permanent", "Nonpermanent"))
data1$data1_est = as.numeric(data1$data1_est)
data1$data1_high = as.numeric(data1$data1_high)
data1$data1_low = as.numeric(data1$data1_low)

whiteman_plot1 = 
  ggplot(data = data1, aes(x = pond_type, y = data1_est)) +
  theme_bw(60) +
  geom_point(size = 15) +
  geom_errorbar(ymin = data1_low, ymax = data1_high, width = 0.2, linewidth = 3) +
  scale_y_continuous(limits = c(0, 600), breaks = c(0, 100, 200, 300, 400, 500, 600)) +
  labs(x = "Pond Type", y = "Calories in Metamorph Stomach") +
  theme(axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 45, margin = margin(t = 0, r = 30, b = 0, l = 0)))

whiteman_plot1

#ggsave(whiteman_plot1, filename = "Outputs/whiteman_plot1.png",  width = 14, height = 12, dpi = 100)


###metamorph diet composition in nonpermanent ponds#####

prey_taxa = c("FS", "ZOO", "BEN", "TER")
data2_est = c(89, 2.1, 5.7, 1.8)
data2_high = c(97, 2.9, 8.3, 2.8)
data2_low = c(81, 1.3, 3.1, 0.8)

data2 = as.data.frame(cbind(prey_taxa, data2_est, data2_high, data2_low))
data2$prey_taxa = as.factor(data2$prey_taxa)
data2$prey_taxa = factor(data2$prey_taxa, levels = c("FS", "ZOO", "BEN", "TER"))
data2$data2_est = as.numeric(data2$data2_est)
data2$data2_high = as.numeric(data2$data2_high)
data2$data2_low = as.numeric(data2$data2_low)

whiteman_plot2 = 
  ggplot(data = data2, aes(x = prey_taxa, y = data2_est)) +
  theme_bw(60) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_errorbar(ymin = data2_low, ymax = data2_high, width = 0.2, linewidth = 3) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  labs(x = "Prey Taxa", y = "Calories in Metamorph Stomach") +
  theme(axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 45, margin = margin(t = 0, r = 30, b = 0, l = 0)))

whiteman_plot2

#ggsave(whiteman_plot2, filename = "Outputs/whiteman_plot2.png",  width = 14, height = 12, dpi = 100)


###metamorph diet composition in nonpermanent ponds#####

prey_type = c("MOSQ", "MDG", "COPE", "FS")
data3_est = c(0.07, -0.64, -0.3, 0.83)
data3_high = c(0.2, -0.53, -0.24, 0.87)
data3_low = c(-0.06, -0.75, -0.36, 0.79)

data3 = as.data.frame(cbind(prey_type, data3_est, data3_high, data3_low))
data3$prey_type = as.factor(data3$prey_type)
data3$prey_type = factor(data3$prey_type, levels = c("MOSQ", "MDG", "COPE", "FS"))
data3$data3_est = as.numeric(data3$data3_est)
data3$data3_high = as.numeric(data3$data3_high)
data3$data3_low = as.numeric(data3$data3_low)

whiteman_plot3 = 
  ggplot(data = data3, aes(x = prey_type, y = data3_est)) +
  theme_bw(60) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_hline(yintercept=0, size=2) +
  geom_errorbar(ymin = data3_low, ymax = data3_high, width = 0.2, linewidth = 3) +
  scale_y_continuous(limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1)) +
  labs(x = "Prey Type", y = "Metamorph Electivity") +
  theme(axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 45, margin = margin(t = 0, r = 30, b = 0, l = 0)))

whiteman_plot3

#ggsave(whiteman_plot3, filename = "Outputs/whiteman_plot3.png",  width = 14, height = 12, dpi = 100)


