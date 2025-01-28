# R version 4.3.2 (2023-10-31 ucrt)


library(ggeffects) #  ggeffects_1.5.0 
library(ggplot2)   #  ggplot2_3.5.1
library(dplyr)     #  dplyr_1.1.4
library(glmmTMB)   #  glmmTMB_1.1.8
library(lsmeans)   #  lsmeans_2.30-0
library(patchwork) #  patchwork_1.2.0
library(forcats)   #  forcats_1.0.0
library(qpcR)      #  qpcR_1.4-1
library(ordinal)   #  ordinal_2023.12-4.1

rm(list=ls())
set.seed(12)

# Set working directory to a folder that contains two additional folders: 
# a "Data" folder (with the bleaching data csv file) and an empty "Outputs" folder
# to store figures generated here.

path = getwd()

# Load bleaching data
dat = read.csv("Data/bleaching_data.csv", sep = ",")


dat$Bleaching = as.factor(dat$Bleaching)


################################################################################
#                                  Plot raw data                               #
################################################################################

# Colour scheme for bleaching scores
colours_sev = c("#FFFFCC", "#F0E442","#E69F00", "#D55E00" , "#D52E00","#990000")
names(colours_sev) = c("1", "2", "3", "4", "5", "6")

# Depth labels
depth_names <- c(
  `D` = "Deep",
  `S` = "Shallow"
)


Taxa_plot = ggplot(data = dat, aes(x = Taxa, fill = factor(Bleaching)))+
  geom_bar(position = "fill", stat = "count", col = "black") +
  scale_fill_manual(values = colours_sev, name = "Bleaching category",
                    breaks = c("1", "2", "3", "4", "5", "6"),
                    labels = c("no bleaching",
                               "minor partial bleaching",
                               "major partial bleaching",
                               "fully bleached",
                               "partial mortality",
                               "whole colony mortality")) +
  geom_hline(yintercept = c(0.25, 0.5, 0.75), col = "grey", linetype = "dashed") +
  theme_classic() +
  scale_x_discrete(limits = c("Abra", "Acor", "Atab", "Adgt", "Gspp", "Lobo_sp",
                              "Menc", "Pdam", "Plat", "Pmas", "Pver", "Shys", "Spis"),
                   labels = c("Acropora staghorn", "Acropora corymbose",
                              "Acropora tabular", "Acropora digitate",
                              "Goniastrea sp.", "Lobophyllia sp.",
                              "Montipora encrusting", "Pocillopora damicornis",
                              "Platygyra sp.", "Porites massive", 
                              "Pocillopora verrucosa", "Seriatopora hystrix",
                              "Stylophora pistillata")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 1)) +
  xlab("") +
  ylab("Proportion of colonies")
Taxa_plot



Inshore_plot =
  ggplot(data = dat[dat$Cross_shelf == "inshore", ], aes(x = Exposure, fill = factor(Bleaching)))+
  geom_bar(position = "fill", stat = "count", col = "black") +
  scale_fill_manual(values = colours_sev, name = "Bleaching category",
                    breaks = c("1", "2", "3", "4", "5", "6"),
                    labels = c("no bleaching",
                               "minor partial bleaching",
                               "major partial bleaching",
                               "fully bleached",
                               "partial mortality",
                               "whole colony mortality")) +
  geom_hline(yintercept = c(0.25, 0.5, 0.75), col = "grey", linetype = "dashed") +
  theme_classic() +
  facet_wrap(~ Zone, labeller = as_labeller(depth_names))+
  scale_x_discrete(limits = c("FR", "FL", "BA", "LA"),
                   labels = c("Front", "Flank", "Back", "Lagoon")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 1),
        strip.background =element_rect(fill="white")) +
  xlab("") +
  ylab("Proportion of colonies")
Inshore_plot  

Offshore_plot =
  ggplot(data = dat[dat$Cross_shelf == "offshore", ], aes(x = Exposure, fill = factor(Bleaching)))+
  geom_bar(position = "fill", stat = "count", col = "black") +
  scale_fill_manual(values = colours_sev, name = "Bleaching category",
                    breaks = c("1", "2", "3", "4", "5", "6"),
                    labels = c("no bleaching",
                               "minor partial bleaching",
                               "major partial bleaching",
                               "fully bleached",
                               "partial mortality",
                               "whole colony mortality")) +
  geom_hline(yintercept = c(0.25, 0.5, 0.75), col = "grey", linetype = "dashed") +
  theme_classic() +
  facet_wrap(~ Zone, labeller = as_labeller(depth_names))+
  scale_x_discrete(limits = c("FR", "FL", "BA", "LA"),
                   labels = c("Front", "Flank", "Back", "Lagoon")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 1),
        strip.background = element_rect(fill="white")) +
  xlab("") +
  ylab("Proportion of colonies")


Offshore_plot  

Figure2 = (Taxa_plot + ggtitle("(a) Taxa") + theme(plot.title = element_text(face="bold"))) /
              ((Inshore_plot  + theme(legend.position = "none", plot.title = element_text(face="bold")) +  ggtitle("(b) Inshore")) + 
                 (Offshore_plot + theme(legend.position = "none", plot.title = element_text(face="bold")) + ylab("") + ggtitle("(c) Offshore")))

Figure2
# save plots of raw data
#ggsave(paste(path, "/Outputs/Figure2.png", sep = ""), Figure2, width = 25, height = 20, units = "cm", dpi = 300)



################################################################################
#                    Model selection - shallow sites only                      #
################################################################################

dat$Bleaching = factor(dat$Bleaching, levels = c("1", "2", "3","4", "5", "6"))

# Depth sites are only available in the Palms and Davies Reef
d_Zone = dat[dat$Reef %in% c("PA", "OCDA"), ]

# Select taxa that are abundant at both depths
d_Zone = d_Zone[d_Zone$Taxa %in% c("Abra", "Acor", "Lobo_sp", "Menc", "Pmas",
                                   "Plat"), ]


# Modifty to "run_models = TRUE" if it is the first time running the script
# and the models have not been saved in the Outputs folder
run_models = FALSE

if (run_models) {
  
  # TABLE S2 - COMPARISON 1

  m0 = clmm(Bleaching ~ 
              log_area +
              Taxa +
              log_area : Taxa +
              ubedmean +
              Cross_shelf +
              ubedmean : Cross_shelf +
              fd_rugosity_plane +
              slope +
              rugosity_plane +
              DHW +
              (1|Site) + (1|Plot),
            data =  dat[dat$Zone == "S", ])
  
  saveRDS(m0, paste(path,"/Outputs/m0.rds", sep = ""))
  
 
  
  # no DHW
  m0a = clmm(Bleaching ~ 
               log_area +
               Taxa +
               log_area : Taxa +
               ubedmean +
               Cross_shelf +
               ubedmean : Cross_shelf +
               fd_rugosity_plane +
               slope +
               rugosity_plane +
               (1|Site) + (1|Plot),
             data =  dat[dat$Zone == "S", ])
  saveRDS(m0a, paste(path,"/Outputs/m0a.rds", sep = ""))
  
  # no rugosity
  m0b = clmm(Bleaching ~ 
               log_area +
               Taxa +
               log_area : Taxa +
               ubedmean +
               Cross_shelf +
               ubedmean : Cross_shelf +
               fd_rugosity_plane +
               slope +
               DHW +
               (1|Site) + (1|Plot),
             data =  dat[dat$Zone == "S", ])
  saveRDS(m0b, paste(path,"/Outputs/m0b.rds", sep = ""))
  
  
  # Compare models
  akaike.weights(c(AIC(m0), AIC(m0a), AIC(m0b)))
  
  
  # TABLE S2 - COMPARISON 2
  
  # no DHW
  m0b1 = clmm(Bleaching ~ 
               log_area +
               Taxa +
               log_area : Taxa +
               ubedmean +
               Cross_shelf +
               ubedmean : Cross_shelf +
               fd_rugosity_plane +
               slope +
               (1|Site) + (1|Plot),
             data =  dat[dat$Zone == "S", ])
  saveRDS(m0b1, paste(path,"/Outputs/m0b1.rds", sep = ""))

  akaike.weights(c(AIC(m0b), AIC(m0b1)))
  
 
  
  
  m = m0b
  
} else {

  m <- readRDS( paste(path,"/Outputs/m0b.rds", sep = ""))

  
}


################################################################################
#                       Setting up plot schemes and labels                     #
################################################################################


# Taxa labels for plotting
taxa_labels =  c(
  `Abra` = "Acropora staghorn",
  `Acor` = "Acropora corymbose",
  `Lobo_sp` = "Lobophyllia sp.",
  `Menc` = "Montipora encrusting",
  `Plat` = "Platygyra sp.",
  `Pmas` = "Porites massive",
  `Adgt` = "Acropora digitate",
  `Atab` = "Acropora tabular",
  `Gspp` = "Goniastrea sp.",
  `Pdam` = "Pocillopora damicornis",
  `Pver` = "Pocillopora verrucosa",
  `Shys` = "Seriatopora hystrix",
  `Spis` = "Stylophora pistillata"
)


################################################################################
#                Plot model predictions - shallow sites only                   #
################################################################################



# Predictions for fractal dimension

# Modifty to "run_preds = TRUE" if it is the first time running the script
# and the predictions have not been saved in the Outputs folder
run_preds = FALSE

if (run_preds) {
  
  # Fractal dimension
  fd_df = predict_response(m,
                           terms = list("fd_rugosity_plane" =  seq(min(dat$fd_rugosity, na.rm = TRUE), max(dat$fd_rugosity_plane, na.rm = TRUE), length.out = 100)), margin = "marginalmeans")
  
  write.table(fd_df, paste(path, "/Outputs/fd_df.csv", sep = ""), sep = ",", row.names = FALSE)
  
  
  # Colony area by Taxa
  area_df = predict_response(m, terms = list("log_area" =  seq(min(dat$log_area, na.rm = TRUE), 
                                                               max(dat$log_area, na.rm = TRUE), length.out = 100),
                                             "Taxa" = unique(dat$Taxa)),
                             margin = "marginalmeans")
  
  
  area_df$Taxa = area_df$group
  write.table(area_df, paste(path,"/Outputs/area_df.csv", sep = ""), sep = ",", row.names = FALSE)
  
  # Cross shelf by water velocity
  ubed_df = predict_response(m, terms = list("ubedmean" =  seq(min(dat$ubedmean, na.rm = TRUE), 
                                                              max(dat$ubedmean, na.rm = TRUE), length.out = 100),
                                            "Cross_shelf" = unique(dat$Cross_shelf)),
                             margin = "marginalmeans")
  write.table(ubed_df, paste(path,"/Outputs/ubed_df.csv", sep = ""), sep = ",", row.names = FALSE)
  
  dhw_df = predict_response(m, terms = list("DHW" = seq(min(dat$DHW, na.rm = TRUE),
                                                        max(dat$DHW, na.rm = TRUE), l = 100)), margin = "marginalmeans")
  write.table(dhw_df, paste(path,"/Outputs/dhw_df.csv", sep = ""), sep = ",", row.names = FALSE)
  
  slope_df = predict_response(m, terms = list("slope" = seq(min(dat$slope, na.rm = TRUE),
                                                        max(dat$slope, na.rm = TRUE), l = 100)), margin = "marginalmeans")
  write.table(slope_df, paste(path,"/Outputs/slope_df.csv", sep = ""), sep = ",", row.names = FALSE)
  
} else{
  fd_df = read.csv(paste(path, "/Outputs/fd_df.csv", sep = ""), sep = ",")
  fd_df$response.level = as.factor(fd_df$response.level)
  area_df = read.csv(paste(path, "/Outputs/area_df.csv", sep = ""), sep = ",")
  area_df$response.level = as.factor(area_df$response.level)
  ubed_df = read.csv(paste(path, "/Outputs/ubed_df.csv", sep = ""), sep = ",")
  ubed_df$response.level = as.factor(ubed_df$response.level)
  dhw_df = read.csv(paste(path, "/Outputs/dhw_df.csv", sep = ""), sep = ",")
  dhw_df$response.level = as.factor(dhw_df$response.level)
  slope_df = read.csv(paste(path, "/Outputs/slope_df.csv", sep = ""), sep = ",")
  slope_df$response.level = as.factor(slope_df$response.level)
}





fd_plot = ggplot() +
  geom_bar(data = fd_df, aes(x = x, y = predicted, fill = response.level), position = "stack", 
           stat = "identity", width = 0.005) +
  scale_fill_manual(values = colours_sev, limits = c("1", "2", "3", "4", "5", "6"),
                    labels = c("no bleaching",
                               "minor partial bleaching",
                               "major partial bleaching",
                               "fully bleached",
                               "partial mortality",
                               "whole colony mortality"),
                    name = "") +
  geom_hline(yintercept = c(0.25, 0.5, 0.75), col = "grey", linetype = "dashed") +
  xlab("Fractal dimension") +
  ylab("Bleaching probability") +
  theme_classic() +
  theme(text = element_text(size = 6), panel.grid.major.y = element_line(linetype = "dashed")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Bleaching == "1", ], aes(x = fd_rugosity_plane, y = -0.05, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Bleaching == "2", ], aes(x = fd_rugosity_plane, y = -0.075, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Bleaching == "3", ], aes(x = fd_rugosity_plane, y = -0.1, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Bleaching == "4", ], aes(x = fd_rugosity_plane, y = -0.125, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Bleaching == "5", ], aes(x = fd_rugosity_plane, y = -0.15, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Bleaching == "6", ], aes(x = fd_rugosity_plane, y = -0.175, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  scale_x_continuous(limits = c(min(dat[dat$Zone == "S" , ]$fd_rugosity_plane, na.rm = TRUE), 
                                max(dat[dat$Zone == "S", ]$fd_rugosity_plane, na.rm = TRUE))) 

fd_plot




# Predictions for area
# Exclude colony sizes not observed in the data set
for (taxa in unique(dat$Taxa)) {
  area_df = area_df[! (area_df$Taxa == taxa & 
                         area_df$x < min(dat[dat$Zone == "S" & dat$Taxa == taxa, ]$log_area, na.rm = TRUE)), ]
  area_df = area_df[! (area_df$Taxa == taxa & 
                         area_df$x > max(dat[dat$Zone == "S" & dat$Taxa == taxa, ]$log_area, na.rm = TRUE)), ]
}


dat$Taxa = as.factor(dat$Taxa)



Figure4 = ggplot() +
  geom_bar(data = area_df, aes(x = x, y = predicted, fill = response.level), position = "stack", 
           stat = "identity", width = 0.05) +
  scale_fill_manual(values = colours_sev, limits = c("1", "2", "3", "4", "5", "6"),
                    labels = c("no bleaching",
                               "minor partial bleaching",
                               "major partial bleaching",
                               "fully bleached",
                               "partial mortality",
                               "whole colony mortality"),
                    name = "") +
  geom_hline(yintercept = c(0.25, 0.5, 0.75), col = "grey", linetype = "dashed") +
  facet_wrap(~ fct_relevel(Taxa,"Shys", "Spis", "Pdam", "Abra", "Gspp", "Plat", "Menc",
                           "Adgt", "Lobo_sp", "Acor", "Pver", "Atab", "Pmas"), 
             labeller = as_labeller(taxa_labels), scales = "free_x", ncol = 5)+
  theme_classic() +
  ylab("Bleaching probability") + 
  xlab(expression(paste("Colony size (", m^2, ")"))) +
  theme(text = element_text(size = 12),
        strip.background =element_rect(fill="white")) +
  scale_x_continuous(breaks = seq(-3, 1, by = 1),
                     labels = c(expression(10^-3),
                                expression(10^-2),
                                expression(10^-1),
                                1,
                                10), limits = c(-3.51, 1.2)) +
  theme(text = element_text(size = 10, family = "Calibri"), 
        panel.grid.major.y = element_line(linetype = "dashed")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Bleaching == "1", ], aes(x = log_area, y = -0.05, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Bleaching == "2", ], aes(x = log_area, y = -0.075, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Bleaching == "3", ], aes(x = log_area, y = -0.1, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Bleaching == "4", ], aes(x = log_area, y = -0.125, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Bleaching == "5", ], aes(x = log_area, y = -0.15, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Bleaching == "6", ], aes(x = log_area, y = -0.175, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7)


Figure4 


# Save area plot 
#ggsave(paste(path, "/Outputs/Figure4.png", sep = ""), Figure4, width = 25, height = 15, units = "cm", dpi = 300)



# Predictions for water velocity
ubed_df$Cross_shelf = ubed_df$group

ubed_plot_in = ggplot() +
  geom_bar(data = ubed_df[ubed_df$Cross_shelf == "inshore", ], aes(x = x, y = predicted, fill = response.level), position = "stack", 
           stat = "identity", width = 0.004) +
  scale_fill_manual(values = colours_sev, limits = c("1", "2", "3", "4", "5", "6"),
                    labels = c("no bleaching",
                               "minor partial bleaching",
                               "major partial bleaching",
                               "fully bleached",
                               "partial mortality",
                               "whole colony mortality"),
                    name = "") +
  geom_hline(yintercept = c(0.25, 0.5, 0.75), col = "grey", linetype = "dashed") +
  xlab(expression(paste("Water velocity at bed (m/s", ")"))) +
  ylab("Bleaching probability") +
  theme_classic() +
  theme(text = element_text(size = 6), panel.grid.major.y = element_line(linetype = "dashed")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Cross_shelf == "inshore" & dat$Bleaching == "1", ], aes(x = ubedmean, y = -0.05, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Cross_shelf == "inshore" & dat$Bleaching == "2", ], aes(x = ubedmean, y = -0.075, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Cross_shelf == "inshore" & dat$Bleaching == "3", ], aes(x = ubedmean, y = -0.1, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Cross_shelf == "inshore" & dat$Bleaching == "4", ], aes(x = ubedmean, y = -0.125, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Cross_shelf == "inshore" & dat$Bleaching == "5", ], aes(x = ubedmean, y = -0.15, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Cross_shelf == "inshore" & dat$Bleaching == "6", ], aes(x = ubedmean, y = -0.175, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  scale_x_continuous(limits = c(min(dat[dat$Zone == "S" & dat$Cross_shelf == "inshore", ]$ubedmean, na.rm = TRUE), 
                                max(dat[dat$Zone == "S"& dat$Cross_shelf == "inshore", ]$ubedmean, na.rm = TRUE))) 

ubed_plot_in


ubed_plot_of = ggplot() +
  geom_bar(data = ubed_df[ubed_df$Cross_shelf == "offshore", ], aes(x = x, y = predicted, fill = response.level), position = "stack", 
           stat = "identity", width = 0.005) +
  scale_fill_manual(values = colours_sev, limits = c("1", "2", "3", "4", "5", "6"),
                    labels = c("no bleaching",
                               "minor partial bleaching",
                               "major partial bleaching",
                               "fully bleached",
                               "partial mortality",
                               "whole colony mortality"),
                    name = "") +
  geom_hline(yintercept = c(0.25, 0.5, 0.75), col = "grey", linetype = "dashed") +
  xlab(expression(paste("Water velocity at bed (m/s", ")"))) +
  ylab("Bleaching probability") +
  theme_classic() +
  theme(text = element_text(size = 6), panel.grid.major.y = element_line(linetype = "dashed")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Cross_shelf == "offshore" & dat$Bleaching == "1", ], aes(x = ubedmean, y = -0.05, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Cross_shelf == "offshore" & dat$Bleaching == "2", ], aes(x = ubedmean, y = -0.075, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Cross_shelf == "offshore" & dat$Bleaching == "3", ], aes(x = ubedmean, y = -0.1, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Cross_shelf == "offshore" & dat$Bleaching == "4", ], aes(x = ubedmean, y = -0.125, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Cross_shelf == "offshore" & dat$Bleaching == "5", ], aes(x = ubedmean, y = -0.15, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Cross_shelf == "offshore" & dat$Bleaching == "6", ], aes(x = ubedmean, y = -0.175, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  scale_x_continuous(limits = c(min(dat[dat$Zone == "S" & dat$Cross_shelf == "offshore", ]$ubedmean, na.rm = TRUE), 
                                max(dat[dat$Zone == "S"& dat$Cross_shelf == "offshore", ]$ubedmean, na.rm = TRUE))) 

ubed_plot_of


# Predictions for DHW

dhw_plot = ggplot() +
  geom_bar(data = dhw_df, aes(x = x, y = predicted, fill = response.level), position = "stack", 
           stat = "identity", width = 0.015) +
  scale_fill_manual(values = colours_sev, limits = c("1", "2", "3", "4", "5", "6"),
                    labels = c("no bleaching",
                               "minor partial bleaching",
                               "major partial bleaching",
                               "fully bleached",
                               "partial mortality",
                               "whole colony mortality"),
                    name = "") +
  geom_hline(yintercept = c(0.25, 0.5, 0.75), col = "grey", linetype = "dashed") +
  xlab("Degree Heating Weeks (°C-weeks)") +
  ylab("Bleaching probability") +
  theme_classic() +
  theme(text = element_text(size = 6), panel.grid.major.y = element_line(linetype = "dashed")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Bleaching == "1", ], aes(x = DHW, y = -0.05, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Bleaching == "2", ], aes(x = DHW, y = -0.075, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Bleaching == "3", ], aes(x = DHW, y = -0.1, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Bleaching == "4", ], aes(x = DHW, y = -0.125, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Bleaching == "5", ], aes(x = DHW, y = -0.15, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Bleaching == "6", ], aes(x = DHW, y = -0.175, fill = Bleaching), pch = 21,
                height = 0.005, width = 0 , alpha = 0.7) +
  scale_x_continuous(limits = c(min(dat[dat$Zone == "S", ]$DHW, na.rm = TRUE), 
                                max(dat[dat$Zone == "S", ]$DHW, na.rm = TRUE))) 
dhw_plot




# Predictions for slope

slope_plot = ggplot() +
  geom_bar(data = slope_df, aes(x = x, y = predicted, fill = response.level), position = "stack", 
           stat = "identity", width = 1.5) +
  scale_fill_manual(values = colours_sev, limits = c("1", "2", "3", "4", "5", "6"),
                    labels = c("no bleaching",
                               "minor partial bleaching",
                               "major partial bleaching",
                               "fully bleached",
                               "partial mortality",
                               "whole colony mortality"),
                    name = "") +
  geom_hline(yintercept = c(0.25, 0.5, 0.75), col = "grey", linetype = "dashed") +
  xlab("Slope (°)") +
  ylab("Bleaching probability") +
  theme_classic() +
  theme(text = element_text(size = 6), panel.grid.major.y = element_line(linetype = "dashed")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Bleaching == "1", ], aes(x = slope, y = -0.05, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Bleaching == "2", ], aes(x = slope, y = -0.075, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Bleaching == "3", ], aes(x = slope, y = -0.1, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Bleaching == "4", ], aes(x = slope, y = -0.125, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Bleaching == "5", ], aes(x = slope, y = -0.15, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  geom_jitter(data = dat[dat$Zone == "S" & dat$Bleaching == "6", ], aes(x = slope, y = -0.175, fill = Bleaching), pch = 21,
              height = 0.005, width = 0, alpha = 0.7) +
  scale_x_continuous(limits = c(min(dat[dat$Zone == "S", ]$slope, na.rm = TRUE), 
                                max(dat[dat$Zone == "S", ]$slope, na.rm = TRUE))) +
  ggtitle("c")
slope_plot



Figure3 = (fd_plot +ggtitle("(a)") + theme(text = element_text(family = "Calibri"), legend.position = "none", plot.title = element_text(face="bold")) ) + 
             (dhw_plot +ggtitle("(b)") +  ylab("") + theme(text = element_text(family = "Calibri"),  legend.position = "none", plot.title = element_text(face="bold")) ) +
             (slope_plot +ggtitle("(c)") +  ylab("") + theme(text = element_text(family = "Calibri"), plot.title = element_text(face="bold")) ) +
  (ubed_plot_in + ggtitle("(d) Inshore")  + theme(text = element_text(family = "Calibri"), legend.position = "none", plot.title = element_text(face="bold"))) + 
     (ubed_plot_of + ggtitle("(e) Offshore") + ylab("")+ theme(text = element_text(family = "Calibri"), legend.position = "none", plot.title = element_text(face="bold"))) +plot_layout(ncol = 3) 
Figure3

# Save Figure 3
#ggsave(paste(path,"/Outputs/Figure3.png",sep = ""), Figure3, width = 15, height = 10, units = "cm", dpi = 300)




################################################################################
#                     Model selection - deep vs. shallow                      #
################################################################################


# Depth sites are only available in the Palms and Davies Reef
d_Zone = dat[dat$Reef %in% c("PA", "OCDA"), ]

# Select taxa that are abundant at both depths
d_Zone = d_Zone[d_Zone$Taxa %in% c("Abra", "Acor", "Lobo_sp", "Menc", "Pmas",
                                   "Plat"), ]

# Modifty to "run_models = TRUE" if it is the first time running the script
# and the models have not been saved in the Outputs folder
run_models = FALSE

if (run_models) {
  
  # Initial model, includes all interactions 
  m1 = clmm(Bleaching ~ 
               log_area +
               Taxa +
               log_area : Taxa +
               ubedmean +
               Cross_shelf +
               ubedmean : Cross_shelf +
               fd_rugosity_plane +
               DHW +
               rugosity_plane +
               slope +
               Zone +
               Zone : Taxa +
               Zone : Cross_shelf +
               (1|Site) + (1|Plot),
             data =  d_Zone)
  saveRDS(m1, paste(path, "/Outputs/m1.rds", sep = ""))
  
  
  # TABLE S3 - COMPARISON 1
  # Check if the interactions are important

  # no Zone: cross_shelf interaction
  m1a = clmm(Bleaching ~ 
              log_area +
              Taxa +
              log_area : Taxa +
              ubedmean +
              Cross_shelf +
              ubedmean : Cross_shelf +
              fd_rugosity_plane +
              DHW +
              rugosity_plane +
              slope +
              Zone +
              Zone : Taxa +
              (1|Site) + (1|Plot),
            data =  d_Zone)
  saveRDS(m1a, paste(path, "/Outputs/m1a.rds", sep = ""))
  
  
  akaike.weights(c(AIC(m1), AIC(m1a)))

 
  # TABLE S3 - COMPARISON 2
  # Model selection on the rest of the fixed effects (i.e., no interactions)
  
  # no FD
  m1a1 = clmm(Bleaching ~ 
               log_area +
               Taxa +
               log_area : Taxa +
               ubedmean +
               Cross_shelf +
               ubedmean : Cross_shelf +
               DHW +
               rugosity_plane +
               slope +
               Zone +
               Zone : Taxa +
               (1|Site) + (1|Plot),
             data =  d_Zone)
  saveRDS(m1a1, paste(path, "/Outputs/m1a1.rds", sep = ""))
  
  
  # no DHW
  m1a2 = clmm(Bleaching ~ 
               log_area +
               Taxa +
               log_area : Taxa +
               ubedmean +
               Cross_shelf +
               ubedmean : Cross_shelf +
               fd_rugosity_plane +
               rugosity_plane +
               slope +
               Zone +
               Zone : Taxa +
               (1|Site) + (1|Plot),
             data =  d_Zone)
  saveRDS(m1a2, paste(path, "/Outputs/m1a2.rds", sep = ""))
  
  # no rugosity
  m1a3 = clmm(Bleaching ~ 
               log_area +
               Taxa +
               log_area : Taxa +
               ubedmean +
               Cross_shelf +
               ubedmean : Cross_shelf +
               fd_rugosity_plane +
               DHW +
               slope +
               Zone +
               Zone : Taxa +
               (1|Site) + (1|Plot),
             data =  d_Zone)
  saveRDS(m1a3, paste(path, "/Outputs/m1a3.rds", sep = ""))

  
  # no slope
  m1a4 = clmm(Bleaching ~ 
               log_area +
               Taxa +
               log_area : Taxa +
               ubedmean +
               Cross_shelf +
               ubedmean : Cross_shelf +
               fd_rugosity_plane +
               DHW +
               rugosity_plane +
               Zone +
               Zone : Taxa +
               (1|Site) + (1|Plot),
             data =  d_Zone)
  saveRDS(m1a4, paste(path, "/Outputs/m1a4.rds", sep = ""))
  
  akaike.weights(c(AIC(m1a), AIC(m1a1), AIC(m1a2), AIC(m1a3), AIC(m1a4)))
  
  
  # TABLE S3 - COMPARISON 3
  # Model selection on additional fixed effects
  
  # no FD
  m1a3a = clmm(Bleaching ~ 
                log_area +
                Taxa +
                log_area : Taxa +
                ubedmean +
                Cross_shelf +
                ubedmean : Cross_shelf +
                DHW +
                slope +
                Zone +
                Zone : Taxa +
                (1|Site) + (1|Plot),
              data =  d_Zone)
  saveRDS(m1a3a, paste(path, "/Outputs/m1a3a.rds", sep = ""))
  
  
  # no DHW
  m1a3b = clmm(Bleaching ~ 
                log_area +
                Taxa +
                log_area : Taxa +
                ubedmean +
                Cross_shelf +
                ubedmean : Cross_shelf +
                fd_rugosity_plane +
                slope +
                Zone +
                Zone : Taxa +
                (1|Site) + (1|Plot),
              data =  d_Zone)
  saveRDS(m1a3b, paste(path, "/Outputs/m1a3b.rds", sep = ""))
  
  
  # no slope
  m1a3c = clmm(Bleaching ~ 
                log_area +
                Taxa +
                log_area : Taxa +
                ubedmean +
                Cross_shelf +
                ubedmean : Cross_shelf +
                fd_rugosity_plane +
                DHW +
                Zone +
                Zone : Taxa +
                (1|Site) + (1|Plot),
              data =  d_Zone)
  saveRDS(m1a3c, paste(path, "/Outputs/m1a3c.rds", sep = ""))
  
  akaike.weights(c(AIC(m1a3), AIC(m1a3a), AIC(m1a3b), AIC(m1a3c)))
  
  
 
  # TABLE S3 - COMPARISON 4
  # Model selection on additional fixed effects
  
  # no FD
  m1a3c1 = clmm(Bleaching ~ 
                 log_area +
                 Taxa +
                 log_area : Taxa +
                 ubedmean +
                 Cross_shelf +
                 ubedmean : Cross_shelf +
                 DHW +
                 Zone +
                 Zone : Taxa +
                 (1|Site) + (1|Plot),
               data =  d_Zone)
  saveRDS(m1a3c1, paste(path, "/Outputs/m1a3c1.rds", sep = ""))
  
  # no DHW
  m1a3c2 = clmm(Bleaching ~ 
                 log_area +
                 Taxa +
                 log_area : Taxa +
                 ubedmean +
                 Cross_shelf +
                 ubedmean : Cross_shelf +
                 fd_rugosity_plane +
                 Zone +
                 Zone : Taxa +
                 (1|Site) + (1|Plot),
               data =  d_Zone)
  saveRDS(m1a3c2, paste(path, "/Outputs/m1a3c2.rds", sep = ""))
  
  
  akaike.weights(c(AIC(m1a3c), AIC(m1a3c1), AIC(m1a3c2)))
  
  m1 = m1a3c

} else{
  m1 = readRDS(paste(path, "/Outputs/m1a3c.rds", sep = ""))
}




################################################################################
#                  Plot model predictions - deep vs shallow                    #
################################################################################

# Modifty to "run_preds = TRUE" if it is the first time running the script
# and the predictions have not been saved in the Outputs folder
run_preds = FALSE

if (run_preds) {
  zone_taxa = predict_response(m1, terms = list(
    "Taxa" = unique(d_Zone$Taxa),
    "Zone" = c("S", "D")
  ),
  margin = "marginalmeans")
  
  zone_taxa$Taxa = zone_taxa$x
  zone_taxa$Zone = zone_taxa$group
  
  
  write.table(zone_taxa, paste(path,"/Outputs/zone_taxa.csv", sep = ""), sep = ",", row.names = FALSE)
  
  
} else{
  zone_taxa = read.csv(paste(path,"/Outputs/zone_taxa.csv", sep = ""), sep = ",")
 
}


zone_taxa$response.level = as.factor(zone_taxa$response.level)


Figure5 = ggplot() +
  geom_bar(data = zone_taxa, aes(x = Zone, y = predicted, fill = response.level), position = "stack", 
           stat = "identity", col = "black") +
  scale_fill_manual(values = colours_sev, limits = c("1", "2", "3", "4", "5", "6"),
                    labels = c("no bleaching",
                               "minor partial bleaching",
                               "major partial bleaching",
                               "fully bleached",
                               "partial mortality",
                               "whole colony mortality"),
                    name = "") +
  geom_hline(yintercept = c(0.25, 0.5, 0.75), col = "grey", linetype = "dashed") +
  facet_wrap(~ Taxa, labeller =  as_labeller(taxa_labels), scale = "free_x" )+
  theme_classic() +
  ylab("Bleaching probability") + 
  theme(text = element_text(size = 12),
        strip.background =element_rect(fill="white")) +
  scale_x_discrete(limits = c("S" , "D"), labels = c("shallow", "deep")) +
  xlab("") +
 
  theme(text = element_text(size = 12, family = "Calibri"), panel.grid.major.y = element_line(linetype = "dashed")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  geom_jitter(data = d_Zone[ d_Zone$Bleaching == "1", ], aes(x = Zone, y = -0.05, fill = Bleaching), pch = 21,
              height = 0.005, width = 0.4, alpha = 0.7) +
  geom_jitter(data = d_Zone[d_Zone$Bleaching == "2", ], aes(x = Zone,  y = -0.075, fill = Bleaching), pch = 21,
              height = 0.005, width = 0.4, alpha = 0.7) +
  geom_jitter(data = d_Zone[ d_Zone$Bleaching == "3", ], aes(x = Zone,  y = -0.1, fill = Bleaching), pch = 21,
              height = 0.005, width = 0.4, alpha = 0.7) +
  geom_jitter(data = d_Zone[ d_Zone$Bleaching == "4", ], aes(x = Zone,  y = -0.125, fill = Bleaching), pch = 21,
              height = 0.005, width = 0.4, alpha = 0.7) +
  geom_jitter(data = d_Zone[ d_Zone$Bleaching == "5", ], aes(x = Zone,  y = -0.15, fill = Bleaching), pch = 21,
              height = 0.005, width = 0.4, alpha = 0.7) +
  geom_jitter(data = d_Zone[ d_Zone$Bleaching == "6", ], aes(x = Zone,  y = -0.175, fill = Bleaching), pch = 21,
              height = 0.005, width = 0.4, alpha = 0.7)


Figure5


# save depth plot
#ggsave(paste(path,"/Outputs/Figure5.png", sep = ""), Figure5, width = 17, height = 15, units = "cm", dpi = 300)







################################################################################
#                               Save model outputs                             #
################################################################################

#write.table(summary(m)$coefficients, paste(path, "/Outputs/m_shallow.csv", sep = ""), sep = ",", row.names = TRUE )

#write.table(summary(m1)$coefficients, paste(path, "/Outputs/m_shallow_deep.csv", sep = ""), sep = ",", row.names = TRUE )


