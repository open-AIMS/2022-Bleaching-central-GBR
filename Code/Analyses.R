# R version 4.3.2 (2023-10-31 ucrt)

library(ggeffects) #  ggeffects_1.5.0 
library(ggplot2)   #  ggplot2_3.5.1
library(dplyr)     #  dplyr_1.1.4
library(glmmTMB)   #  glmmTMB_1.1.8
library(lsmeans)   #  lsmeans_2.30-0
library(patchwork) #  patchwork_1.2.0
library(forcats)   #  forcats_1.0.0
library(qpcR)      #  qpcR_1.4-1


set.seed(12)

# Load bleaching data
dat = read.csv("Data/bleaching_data.csv", sep = ",")

################################################################################
#                                  Plot raw data                               #
################################################################################

# Colour scheme for bleaching scores
colours_taxa = c("#FFFFFC", "#F0E442","#E69F00", "#D55E00" , "#D52E00","#990000")
names(colours_taxa) = c("1", "2", "3", "4", "5", "6")

# Depth labels
depth_names <- c(
  `D` = "Deep",
  `S` = "Shallow"
)


Taxa_plot = ggplot(data = dat, aes(x = Taxa, fill = factor(Bleaching)))+
  geom_bar(position = "fill", stat = "count", col = "black") +
  scale_fill_manual(values = colours_taxa, name = "Bleaching category",
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
  scale_fill_manual(values = colours_taxa, name = "Bleaching category",
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
  scale_fill_manual(values = colours_taxa, name = "Bleaching category",
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
Offshore_plot  

Figure2 = (Taxa_plot + ggtitle("a- Taxa"))/ ((Inshore_plot  + theme(legend.position = "none") +  ggtitle("b- Inshore")) + (Offshore_plot + theme(legend.position = "none") + ylab("") +
                                                                                                                             ggtitle("c- Offshore")))
# save plots of raw data
#ggsave("Outputs/Figure2.png", Figure2, width = 25, height = 20, units = "cm", dpi = 300)



################################################################################
#                    Model selection - shallow sites only                      #
################################################################################


# Reef is estimated to be ~0 and gives problems with convergence, so cannot
# include as a random effect in the model, site is also ~0

# TABLE S2 - COMPARISON 1
# Check whether to include FD or rugosity
m0 = glmmTMB(Bleaching_res ~ 
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
               (1|Site/Plot),
             data =  dat[dat$Zone == "S", ],
             family = "binomial",
             start = list(theta = c(1,  1)))

# no FD
m0a = glmmTMB(Bleaching_res ~ 
                log_area +
                Taxa +
                log_area : Taxa +
                ubedmean +
                Cross_shelf +
                ubedmean : Cross_shelf +
                slope +
                rugosity_plane +
                DHW +
                (1|Site/Plot),
              data =  dat[dat$Zone == "S", ],
              family = "binomial",
              start = list(theta = c(1,  1)))

# no slope
m0b = glmmTMB(Bleaching_res ~ 
                log_area +
                Taxa +
                log_area : Taxa +
                ubedmean +
                Cross_shelf +
                ubedmean : Cross_shelf +
                fd_rugosity_plane +
                rugosity_plane +
                DHW +
                (1|Site/Plot),
              data =  dat[dat$Zone == "S", ],
              family = "binomial",
              start = list(theta = c(1,  1)))

# no DHW
m0c = glmmTMB(Bleaching_res ~ 
                log_area +
                Taxa +
                log_area : Taxa +
                ubedmean +
                Cross_shelf +
                ubedmean : Cross_shelf +
                fd_rugosity_plane +
                slope +
                rugosity_plane +
                (1|Site/Plot),
              data =  dat[dat$Zone == "S", ],
              family = "binomial",
              start = list(theta = c(1,  1)))

# no rugosity
m0d = glmmTMB(Bleaching_res ~ 
                log_area +
                Taxa +
                log_area : Taxa +
                ubedmean +
                Cross_shelf +
                ubedmean : Cross_shelf +
                fd_rugosity_plane +
                slope +
                DHW +
                (1|Site/Plot),
              data =  dat[dat$Zone == "S", ],
              family = "binomial",
              start = list(theta = c(1,  1)))

# Compare models
akaike.weights(c(AIC(m0), AIC(m0a), AIC(m0b), AIC(m0c), AIC(m0d)))


# TABLE S2 - COMPARISON 2

# no FD
m0d1 = glmmTMB(Bleaching_res ~ 
                 log_area +
                 Taxa +
                 log_area : Taxa +
                 ubedmean +
                 Cross_shelf +
                 ubedmean : Cross_shelf +
                 slope +
                 DHW +
                 (1|Site/Plot),
               data =  dat[dat$Zone == "S", ],
               family = "binomial",
               start = list(theta = c(1,  1)))

# no slope
m0d2 = glmmTMB(Bleaching_res ~ 
                 log_area +
                 Taxa +
                 log_area : Taxa +
                 ubedmean +
                 Cross_shelf +
                 ubedmean : Cross_shelf +
                 fd_rugosity_plane +
                 DHW +
                 (1|Site/Plot),
               data =  dat[dat$Zone == "S", ],
               family = "binomial",
               start = list(theta = c(1,  1)))

# no DHW
m0d3 = glmmTMB(Bleaching_res ~ 
                 log_area +
                 Taxa +
                 log_area : Taxa +
                 ubedmean +
                 Cross_shelf +
                 ubedmean : Cross_shelf +
                 fd_rugosity_plane +
                 slope +
                 (1|Site/Plot),
               data =  dat[dat$Zone == "S", ],
               family = "binomial",
               start = list(theta = c(1,  1)))

# Compare models
akaike.weights(c(AIC(m0d), AIC(m0d1), AIC(m0d2), AIC(m0d3)))

# TABLE S2 - COMPARISON 3
# no FD
m0d2a = glmmTMB(Bleaching_res ~ 
                  log_area +
                  Taxa +
                  log_area : Taxa +
                  ubedmean +
                  Cross_shelf +
                  ubedmean : Cross_shelf +
                  DHW +
                  (1|Site/Plot),
                data =  dat[dat$Zone == "S", ],
                family = "binomial",
                start = list(theta = c(1,  1)))

# no DHW
m0d2b = glmmTMB(Bleaching_res ~ 
                  log_area +
                  Taxa +
                  log_area : Taxa +
                  ubedmean +
                  Cross_shelf +
                  ubedmean : Cross_shelf +
                  fd_rugosity_plane +
                  (1|Site/Plot),
                data =  dat[dat$Zone == "S", ],
                family = "binomial",
                start = list(theta = c(1,  1)))

# Compare models
akaike.weights(c( AIC(m0d2), AIC(m0d2a), AIC(m0d2b)))


# best-fit model is m0d2. Includes log_area, Taxa, log_area : Taxa,
# ubedmean, Cross_shelf, ubedmean : Cross_shelf, fd_rugosity_plane, and DHW

m0 = m0d2


################################################################################
#         Model fitting - other response variables, shallow sites only         #
################################################################################

# Use this same structure to fit additional models with different response variables.

# Minor bleaching
m1 = glmmTMB(Bleaching_minor ~ 
               log_area +
               Taxa +
               log_area : Taxa +
               ubedmean +
               Cross_shelf +
               ubedmean : Cross_shelf +
               fd_rugosity_plane +
               DHW +
               (1|Site/Plot),
             data =  dat[dat$Zone == "S", ],
             family = "binomial",
             start = list(theta = c(1,  1)))

# Major bleaching
m2 = glmmTMB(Bleaching_major ~ 
               log_area +
               Taxa +
               log_area : Taxa +
               ubedmean +
               Cross_shelf +
               ubedmean : Cross_shelf +
               fd_rugosity_plane +
               DHW +
               (1|Site/Plot),
             data =  dat[dat$Zone == "S", ],
             family = "binomial",
             start = list(theta = c(1,  1)))

# Severe bleaching
m3 = glmmTMB(Bleaching_fully ~ 
               log_area +
               Taxa +
               log_area : Taxa +
               ubedmean +
               Cross_shelf +
               ubedmean : Cross_shelf +
               fd_rugosity_plane +
               DHW +
               (1|Site/Plot),
             data =  dat[dat$Zone == "S", ],
             family = "binomial",
             start = list(theta = c(1,  1)))
 
################################################################################
#                       Setting up plot schemes and labels                     #
################################################################################

# Set colour scheme for bleaching categories
colours_sev = c("black", "#F0E442",  "#E69F00", "#D55E00")
names(colours_sev) = c("all", "minor", "major", "fully")


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
fd_df = predict_response(m0, terms = "fd_rugosity_plane [all]", margin = "marginalmeans")
fd_df$model = "all"
fd_df1 = predict_response(m1, terms = "fd_rugosity_plane [all]", margin = "marginalmeans")
fd_df1$model = "minor"
fd_df2 = predict_response(m2, terms = "fd_rugosity_plane [all]", margin = "marginalmeans")
fd_df2$model = "major"
fd_df3 = predict_response(m3, terms = "fd_rugosity_plane [all]", margin = "marginalmeans")
fd_df3$model = "fully"
fd_df = rbind(fd_df, fd_df1, fd_df2, fd_df3)


fd_plot = ggplot() +
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_minor == 1, ], 
              aes(x = fd_rugosity_plane, y = Bleaching_minor - 0.055), alpha = 0.1, size = 1, width = 0,
              height = 0.025, col = "#F0E442") +
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_minor == 0, ], 
              aes(x = fd_rugosity_plane, y = Bleaching_minor + 0.055), alpha = 0.1, size = 1, width = 0,
              height = 0.025, col = "#F0E442") +
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_major == 1, ], 
              aes(x = fd_rugosity_plane, y = Bleaching_major ), alpha = 0.1, size = 1, width = 0,
              height = 0.025, col =  "#E69F00") +
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_major == 0, ], 
              aes(x = fd_rugosity_plane, y = Bleaching_major ), alpha = 0.1, size = 1, width = 0,
              height = 0.025, col =  "#E69F00") +
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_fully== 1, ], 
              aes(x = fd_rugosity_plane, y = Bleaching_fully + 0.055), alpha = 0.1, size = 1, width = 0,
              height = 0.025, col =  "#D55E00" ) +
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_fully == 0, ], 
              aes(x = fd_rugosity_plane, y = Bleaching_fully - 0.055), alpha = 0.1, size = 1, width = 0,
              height = 0.025, col =  "#D55E00" ) +
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_res == 1, ], 
              aes(x = fd_rugosity_plane, y = Bleaching_res + 0.055*2), alpha = 0.1, size = 1, width = 0,
              height = 0.025, col =  "black" ) +
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_res == 0, ], 
              aes(x = fd_rugosity_plane, y = Bleaching_res - 0.055*2), alpha = 0.1, size = 1, width = 0,
              height = 0.025, col =  "black" ) +
  geom_ribbon(data = fd_df, aes(x = x, ymin = conf.low, ymax = conf.high, fill = model, group = model), 
              alpha = 0.2) +
  geom_line(data = fd_df, aes(x = x, y = predicted, col = model , group = model), size = 1.05) + 
  scale_colour_manual(values = colours_sev, limits = c("all", "minor", "major", "fully"),
                      labels = c("any bleaching", "minor bleaching", "major bleaching",
                                 "severe bleaching"),
                      name = "") +
  scale_fill_manual(values = colours_sev, limits = c("all", "minor", "major", "fully"),
                    labels = c("any bleaching", "minor bleaching", "major bleaching",
                               "severe bleaching"),
                    name = "") +
  xlab("Fractal dimension") +
  ylab("Bleaching probability") +
  theme_classic() +
  theme(text = element_text(size = 10), panel.grid.major.y = element_line(linetype = "dashed")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25))
fd_plot


# Predictions for area
area_df = predict_response(m0, terms = c("log_area [all]", "Taxa [Adgt, Pmas]"), margin = "marginalmeans")

area_df = rbind(area_df,
                predict_response(m0, terms = c("log_area [all]", "Taxa [Pver, Acor]"), margin = "marginalmeans") )
area_df = rbind(area_df,
                predict_response(m0, terms = c("log_area [all]", "Taxa [Lobo_sp, Menc]"), margin = "marginalmeans") )
area_df = rbind(area_df,
                predict_response(m0, terms = c("log_area [all]", "Taxa [Spis, Abra]"), margin = "marginalmeans") )
area_df = rbind(area_df,
                predict_response(m0, terms = c("log_area [all]", "Taxa [Gspp, Atab]"), margin = "marginalmeans") )
area_df = rbind(area_df,
                predict_response(m0, terms = c("log_area [all]", "Taxa [Pdam, Shys]"), margin = "marginalmeans") )
area_df = rbind(area_df,
                predict_response(m0, terms = c("log_area [all]", "Taxa [Plat]"), margin = "marginalmeans") )

area_df$model = "all"

area_df1 = predict_response(m1, terms = c("log_area [all]", "Taxa [Adgt, Pmas]"), margin = "marginalmeans")

area_df1 = rbind(area_df1,
                 predict_response(m1, terms = c("log_area [all]", "Taxa [Pver, Acor]"), margin = "marginalmeans") )
area_df1 = rbind(area_df1,
                 predict_response(m1, terms = c("log_area [all]", "Taxa [Lobo_sp, Menc]"), margin = "marginalmeans") )
area_df1 = rbind(area_df1,
                 predict_response(m1, terms = c("log_area [all]", "Taxa [Spis, Abra]"), margin = "marginalmeans") )
area_df1 = rbind(area_df1,
                 predict_response(m1, terms = c("log_area [all]", "Taxa [Gspp, Atab]"), margin = "marginalmeans") )
area_df1 = rbind(area_df1,
                 predict_response(m1, terms = c("log_area [all]", "Taxa [Pdam, Shys]"), margin = "marginalmeans") )
area_df1 = rbind(area_df1,
                 predict_response(m1, terms = c("log_area [all]", "Taxa [Plat]"), margin = "marginalmeans") )

area_df1$model = "minor"


area_df2 = predict_response(m2, terms = c("log_area [all]", "Taxa [Adgt, Pmas]"), margin = "marginalmeans")

area_df2 = rbind(area_df2,
                 predict_response(m2, terms = c("log_area [all]", "Taxa [Pver, Acor]"), margin = "marginalmeans") )
area_df2 = rbind(area_df2,
                 predict_response(m2, terms = c("log_area [all]", "Taxa [Lobo_sp, Menc]"), margin = "marginalmeans") )
area_df2 = rbind(area_df2,
                 predict_response(m2, terms = c("log_area [all]", "Taxa [Spis, Abra]"), margin = "marginalmeans") )
area_df2 = rbind(area_df2,
                 predict_response(m2, terms = c("log_area [all]", "Taxa [Gspp, Atab]"), margin = "marginalmeans") )
area_df2 = rbind(area_df2,
                 predict_response(m2, terms = c("log_area [all]", "Taxa [Pdam, Shys]"), margin = "marginalmeans") )
area_df2 = rbind(area_df2,
                 predict_response(m2, terms = c("log_area [all]", "Taxa [Plat]"), margin = "marginalmeans") )

area_df2$model = "major"

area_df3 = predict_response(m3, terms = c("log_area [all]", "Taxa [Adgt, Pmas]"), margin = "marginalmeans")

area_df3 = rbind(area_df3,
                 predict_response(m3, terms = c("log_area [all]", "Taxa [Pver, Acor]"), margin = "marginalmeans") )
area_df3 = rbind(area_df3,
                 predict_response(m3, terms = c("log_area [all]", "Taxa [Lobo_sp, Menc]"), margin = "marginalmeans") )
area_df3 = rbind(area_df3,
                 predict_response(m3, terms = c("log_area [all]", "Taxa [Spis, Abra]"), margin = "marginalmeans") )
area_df3 = rbind(area_df3,
                 predict_response(m3, terms = c("log_area [all]", "Taxa [Gspp, Atab]"), margin = "marginalmeans") )
area_df3 = rbind(area_df3,
                 predict_response(m3, terms = c("log_area [all]", "Taxa [Pdam, Shys]"), margin = "marginalmeans") )
area_df3 = rbind(area_df3,
                 predict_response(m3, terms = c("log_area [all]", "Taxa [Plat]"), margin = "marginalmeans") )

area_df3$model = "fully"

area_df = rbind(area_df, area_df1, area_df2, area_df3)
area_df$Taxa = area_df$group

for (taxa in unique(dat$Taxa)) {
  area_df = area_df[! (area_df$Taxa == taxa & 
                         area_df$x < min(dat[dat$Zone == "S" & dat$Taxa == taxa, ]$log_area, na.rm = TRUE)), ]
  area_df = area_df[! (area_df$Taxa == taxa & 
                         area_df$x > max(dat[dat$Zone == "S" & dat$Taxa == taxa, ]$log_area, na.rm = TRUE)), ]
}


dat$Taxa = as.factor(dat$Taxa)

dat_blmin = dat[ , c("log_area", "Zone", "slope", "ubedmean", "Taxa", "Bleaching_minor",
                     "fd_rugosity_plane")]
dat_blmin = dat_blmin[complete.cases(dat_blmin), ]

dat_blmaj = dat[ , c("log_area","Zone", "slope", "ubedmean", "Taxa", "Bleaching_major",
                      "fd_rugosity_plane")]
dat_blmaj = dat_blmaj[complete.cases(dat_blmaj), ]

dat_blful = dat[ , c("log_area", "Zone","slope", "ubedmean", "Taxa", "Bleaching_fully",
                     "fd_rugosity_plane")]
dat_blful = dat_blful[complete.cases(dat_blful), ]

Figure4 = ggplot() +
  
  geom_jitter(data = dat_blmin[dat_blmin$Zone == "S"  & dat_blmin$Bleaching_minor == 1, ], 
              aes(x = log_area, y = Bleaching_minor + 0.055), na.rm = TRUE, alpha = 0.1, size = 1, width = 0,
              height = 0.01, col = "#F0E442") +
  geom_jitter(data = dat_blmin[dat_blmin$Zone == "S"  & dat_blmin$Bleaching_minor == 0, ], 
              aes(x = log_area, y = Bleaching_minor - 0.055), alpha = 0.1, size = 1, width = 0,
              height = 0.01, col = "#F0E442") +
  geom_jitter(data = dat_blmaj[dat_blmaj$Zone == "S"  & dat_blmaj$Bleaching_major == 1, ], 
              aes(x = log_area, y = Bleaching_major + 0.11), alpha = 0.1, size = 1, width = 0,
              height = 0.01, col =  "#E69F00") +
  geom_jitter(data = dat_blmaj[dat_blmaj$Zone == "S"  & dat_blmaj$Bleaching_major == 0, ], 
              aes(x = log_area, y = Bleaching_major  -0.11), alpha = 0.1, size = 1, width = 0,
              height = 0.01, col =  "#E69F00") +
  geom_jitter(data = dat_blful[dat_blful$Zone == "S"  & dat_blful$Bleaching_fully== 1, ], 
              aes(x = log_area, y = Bleaching_fully + 0.165), alpha = 0.1, size = 1, width = 0,
              height = 0.01, col =  "#D55E00" ) +
  geom_jitter(data = dat_blful[dat_blful$Zone == "S"  & dat_blful$Bleaching_fully == 0, ], 
              aes(x = log_area, y = Bleaching_fully - 0.165), alpha = 0.1, size = 1, width = 0,
              height = 0.01, col =  "#D55E00" ) +
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_res== 1, ], 
              aes(x = log_area, y = Bleaching_res+ 0.165+0.055), alpha = 0.1, size = 1, width = 0,
              height = 0.01, col =  "black" ) +
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_res== 0, ], 
              aes(x = log_area, y = Bleaching_res - 0.165-0.055), alpha = 0.1, size = 1, width = 0,
              height = 0.01, col =  "black" ) +
  geom_line(data = area_df, aes(x = x, y = predicted, col = model, group = model), size = 1.05) +
  #geom_ribbon(data = area_df, aes(x = x, ymin = conf.low, ymax = conf.high, fill = model , group = model), alpha = 0.5) +
  facet_wrap(~ fct_relevel(Taxa,"Shys", "Spis", "Pdam", "Abra", "Gspp", "Plat", "Menc",
                           "Adgt", "Lobo_sp", "Acor", "Pver", "Atab", "Pmas"), 
             labeller = as_labeller(taxa_labels), scales = "fixed", ncol = 5)+
  theme_classic() +
  ylab("Bleaching probability") + 
  scale_colour_manual(values = colours_sev, limits = c("all", "minor", "major", "fully"),
                      labels = c("any bleaching", "minor bleaching", "major bleaching",
                                 "severe bleaching"),
                      name = "") +
  scale_fill_manual(values = colours_sev, limits = c("all", "minor", "major", "fully"),
                    labels = c("any bleaching", "minor bleaching", "major bleaching",
                               "severe bleaching"),
                    name = "") +
  xlab(expression(paste("Colony size (", m^2, ")"))) +
  theme(text = element_text(size = 12),
        strip.background =element_rect(fill="white")) +
  scale_x_continuous(breaks = seq(-3, 1, by = 1),
                     labels = c(expression(10^-3),
                                expression(10^-2),
                                expression(10^-1),
                                1,
                                10)) +
  theme(text = element_text(size = 12, family = "Calibri"), panel.grid.major.y = element_line(linetype = "dashed")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25))

Figure4 

# Save area plot with/without ribbons
#ggsave("Outputs/FigureS2.png", Figure4, width = 25, height = 15, units = "cm", dpi = 300)
#ggsave("Outputs/Figure4.png", Figure4, width = 25, height = 15, units = "cm", dpi = 300)



# Predictions for water velocity
ubed_df = predict_response(m0, terms = c("ubedmean [all]", "Cross_shelf"), margin = "marginalmeans")
ubed_df$model = "all"
ubed_df1 = predict_response(m1, terms = c("ubedmean [all]", "Cross_shelf"), margin = "marginalmeans")
ubed_df1$model = "minor"
ubed_df2 = predict_response(m2, terms = c("ubedmean [all]", "Cross_shelf"), margin = "marginalmeans")
ubed_df2$model = "major"
ubed_df3 = predict_response(m3, terms = c("ubedmean [all]", "Cross_shelf"), margin = "marginalmeans")
ubed_df3$model = "fully"

ubed_df = rbind(ubed_df, ubed_df1, ubed_df2, ubed_df3)
ubed_df$Cross_shelf = ubed_df$group

# (inshore)
ubed_plot_in = ggplot() +
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_minor == 1 & dat$Cross_shelf == "inshore", ], 
              aes(x = ubedmean, y = Bleaching_minor - 0.055), alpha = 0.1, size = 1, width = 0,
              height = 0.025, col = "#F0E442") +
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_minor == 0  & dat$Cross_shelf == "inshore", ], 
              aes(x = ubedmean, y = Bleaching_minor + 0.055), alpha = 0.1, size = 1, width = 0,
              height = 0.025, col = "#F0E442") +
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_major == 1  & dat$Cross_shelf == "inshore", ], 
              aes(x = ubedmean, y = Bleaching_major ), alpha = 0.1, size = 1, width = 0,
              height = 0.025, col =  "#E69F00") +
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_major == 0  & dat$Cross_shelf == "inshore", ], 
              aes(x = ubedmean, y = Bleaching_major ), alpha = 0.1, size = 1, width = 0,
              height = 0.025, col =  "#E69F00") +
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_fully== 1  & dat$Cross_shelf == "inshore", ], 
              aes(x = ubedmean, y = Bleaching_fully + 0.055), alpha = 0.1, size = 1, width = 0,
              height = 0.025, col =  "#D55E00" ) +
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_fully == 0 & dat$Cross_shelf == "inshore", ], 
              aes(x = ubedmean, y = Bleaching_fully - 0.055), alpha = 0.1, size = 1, width = 0,
              height = 0.025, col =  "#D55E00" ) +
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_res == 1  & dat$Cross_shelf == "inshore", ], 
              aes(x = ubedmean, y = Bleaching_res + 0.055 *2), alpha = 0.1, size = 1, width = 0,
              height = 0.025, col =  "black" ) +
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_res == 0 & dat$Cross_shelf == "inshore", ], 
              aes(x = ubedmean, y = Bleaching_res - 0.055 * 2), alpha = 0.1, size = 1, width = 0,
              height = 0.025, col =  "black" ) +
  geom_ribbon(data = ubed_df[ubed_df$Cross_shelf == "inshore", ], aes(x = x, ymin = conf.low, ymax = conf.high,
                                                                      fill = model, group = model),
              alpha = 0.2) +
  geom_line(data = ubed_df[ubed_df$Cross_shelf == "inshore", ], aes(x = x, y = predicted, col = model, group = model),  size = 1.05) +
  theme_classic() +
  ylab("Bleaching probability") + 
  scale_colour_manual(values = colours_sev, limits = c("all", "minor", "major", "fully"),
                      labels = c("any bleaching", "minor bleaching", "major bleaching",
                                 "severe bleaching"),
                      name = "") +
  scale_fill_manual(values = colours_sev, limits = c("all", "minor", "major", "fully"),
                    labels = c("any bleaching", "minor bleaching", "major bleaching",
                               "severe bleaching"),
                    name = "") +
  xlab(expression(paste("Water velocity at bed (m",s^-1, ")"))) +
  theme(text = element_text(size = 10), panel.grid.major.y = element_line(linetype = "dashed")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25))
ubed_plot_in



# (offshore)
ubed_plot_of = ggplot() +
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_minor == 1 & dat$Cross_shelf == "offshore", ], 
              aes(x = ubedmean, y = Bleaching_minor - 0.055), alpha = 0.1, size = 1, width = 0,
              height = 0.025, col = "#F0E442") +
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_minor == 0  & dat$Cross_shelf == "offshore", ], 
              aes(x = ubedmean, y = Bleaching_minor + 0.055), alpha = 0.1, size = 1, width = 0,
              height = 0.025, col = "#F0E442") +
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_major == 1  & dat$Cross_shelf == "offshore", ], 
              aes(x = ubedmean, y = Bleaching_major ), alpha = 0.1, size = 1, width = 0,
              height = 0.025, col =  "#E69F00") +
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_major == 0  & dat$Cross_shelf == "offshore", ], 
              aes(x = ubedmean, y = Bleaching_major ), alpha = 0.1, size = 1, width = 0,
              height = 0.025, col =  "#E69F00") +
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_fully== 1  & dat$Cross_shelf == "offshore", ], 
              aes(x = ubedmean, y = Bleaching_fully + 0.055), alpha = 0.1, size = 1, width = 0,
              height = 0.025, col =  "#D55E00" ) +
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_fully == 0 & dat$Cross_shelf == "offshore", ], 
              aes(x = ubedmean, y = Bleaching_fully - 0.055), alpha = 0.1, size = 1, width = 0,
              height = 0.025, col =  "#D55E00" ) +
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_res == 1  & dat$Cross_shelf == "offshore", ], 
              aes(x = ubedmean, y = Bleaching_res + 0.055 *2), alpha = 0.1, size = 1, width = 0,
              height = 0.025, col =  "black" ) +
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_res == 0 & dat$Cross_shelf == "offshore", ], 
              aes(x = ubedmean, y = Bleaching_res - 0.055 * 2), alpha = 0.1, size = 1, width = 0,
              height = 0.025, col =  "black" ) +
  geom_ribbon(data = ubed_df[ubed_df$Cross_shelf == "offshore", ], aes(x = x, ymin = conf.low, ymax = conf.high,
                                                                       fill = model, group = model),
              alpha = 0.2) +
  geom_line(data = ubed_df[ubed_df$Cross_shelf == "offshore", ], aes(x = x, y = predicted, col = model, group = model),  size = 1.05) +
  theme_classic() +
  ylab("Bleaching severity") + 
  scale_colour_manual(values = colours_sev, limits = c("all", "minor", "major", "fully"),
                      labels = c("any bleaching", "minor bleaching", "major bleaching",
                                 "severe bleaching"),
                      name = "") +
  scale_fill_manual(values = colours_sev, limits = c("all", "minor", "major", "fully"),
                    labels = c("any bleaching", "minor bleaching", "major bleaching",
                               "severe bleaching"),
                    name = "") +
  theme_classic() +
  ylab("Bleaching probability") + 
  xlab(expression(paste("Water velocity at bed (m",s^-1, ")"))) +
  theme(text = element_text(size = 10), panel.grid.major.y = element_line(linetype = "dashed")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25))
ubed_plot_of

# Predictions for DHW
dhw_df = predict_response(m0, terms = "DHW [all]", margin = "marginalmeans")
dhw_df$model = "all"
dhw_df1 = predict_response(m1, terms = "DHW [all]", margin = "marginalmeans")
dhw_df1$model = "minor"
dhw_df2 = predict_response(m2, terms = "DHW [all]", margin = "marginalmeans")
dhw_df2$model = "major"
dhw_df3 = predict_response(m3, terms = "DHW [all]", margin = "marginalmeans")
dhw_df3$model = "fully"

dhw_df = rbind(dhw_df, dhw_df1, dhw_df2, dhw_df3)


dhw_plot = ggplot() +
  
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_minor == 1, ], 
              aes(x = DHW, y = Bleaching_minor - 0.055), alpha = 0.1, size = 1, width = 0,
              height = 0.025, col = "#F0E442") +
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_minor == 0, ], 
              aes(x = DHW, y = Bleaching_minor + 0.055), alpha = 0.1, size = 1, width = 0,
              height = 0.025, col = "#F0E442") +
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_major == 1, ], 
              aes(x = DHW, y = Bleaching_major ), alpha = 0.1, size = 1, width = 0,
              height = 0.025, col =  "#E69F00") +
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_major == 0, ], 
              aes(x = DHW, y = Bleaching_major ), alpha = 0.1, size = 1, width = 0,
              height = 0.025, col =  "#E69F00") +
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_fully== 1, ], 
              aes(x = DHW, y = Bleaching_fully + 0.055), alpha = 0.1, size = 1, width = 0,
              height = 0.025, col =  "#D55E00" ) +
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_fully == 0, ], 
              aes(x = DHW, y = Bleaching_fully - 0.055), alpha = 0.1, size = 1, width = 0,
              height = 0.025, col =  "#D55E00" ) +
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_res== 1, ], 
              aes(x = DHW, y = Bleaching_res + 0.055*2), alpha = 0.1, size = 1, width = 0,
              height = 0.025, col =  "black" ) +
  geom_jitter(data = dat[dat$Zone == "S"  & dat$Bleaching_res == 0, ], 
              aes(x = DHW, y = Bleaching_res - 0.055*2), alpha = 0.1, size = 1, width = 0,
              height = 0.025, col =  "black" ) +
  geom_ribbon(data = dhw_df, aes(x = x, ymin = conf.low, ymax = conf.high, fill = model, group = model), 
              alpha = 0.2) +
  geom_line(data = dhw_df, aes(x = x, y = predicted, col = model , group = model), size = 1.05) + 
  scale_colour_manual(values = colours_sev, limits = c("all", "minor", "major", "fully"),
                      labels = c("any bleaching", "minor bleaching", "major bleaching",
                                 "severe bleaching"),
                      name = "") +
  scale_fill_manual(values = colours_sev, limits = c("all", "minor", "major", "fully"),
                    labels = c("any bleaching", "minor bleaching", "major bleaching",
                               "severe bleaching"),
                    name = "") +
  xlab("Degree Heating Weeks (°C-weeks)") +
  ylab("Bleaching probability") +
  theme_classic() +
  theme(text = element_text(size = 10), panel.grid.major.y = element_line(linetype = "dashed")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25))
dhw_plot




Figure3 = ((fd_plot +ggtitle("a") + theme(text = element_text(family = "Calibri"), legend.position = "none") ) + 
             (dhw_plot +ggtitle("b") +  ylab("") + theme(text = element_text(family = "Calibri")) )) /
  ((ubed_plot_in + ggtitle("c - inshore")  + theme(text = element_text(family = "Calibri"), legend.position = "none")) + 
     (ubed_plot_of + ggtitle("d - offshore") + ylab("")+ theme(text = element_text(family = "Calibri"), legend.position = "none") ))


# Save Figure 3
#ggsave("Outputs/Figure3.png", Figure3, width = 15, height = 10, units = "cm", dpi = 300)




################################################################################
#                     Model selection - deep vs. shallow                      #
################################################################################


# Depth sites are only available in the Palms and Davies Reef
d_Zone = dat[dat$Reef %in% c("PA", "OCDA"), ]

# Select taxa that are abundant at both depths
d_Zone = d_Zone[d_Zone$Taxa %in% c("Abra", "Acor", "Lobo_sp", "Menc", "Pmas",
                                   "Plat"), ]

# Initial model, includes all interactions
m1z = glmmTMB(Bleaching_res ~ 
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
                Zone : Cross_shelf +
                Zone : Taxa +
                Zone : ubedmean +
                Zone : DHW +
                (1|Site/Plot),
              data =  d_Zone,
              family = "binomial",
              start = list(theta = c(1,  1)))


# TABLE S3 - COMPARISON 1
# Check if the interactions are important
# no Zone : DHW
m2z = glmmTMB(Bleaching_res ~ 
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
                Zone : Cross_shelf +
                Zone : Taxa +
                Zone : ubedmean +
                (1|Site/Plot),
              data =  d_Zone,
              family = "binomial",
              start = list(theta = c(1,  1)))

# no Zone : ubedmean
m3z = glmmTMB(Bleaching_res ~ 
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
                Zone : Cross_shelf +
                Zone : Taxa +
                Zone : DHW +
                (1|Site/Plot),
              data =  d_Zone,
              family = "binomial",
              start = list(theta = c(1,  1)))

# no Zone : Cross_shelf
m4z = glmmTMB(Bleaching_res ~ 
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
                Zone : ubedmean +
                Zone : DHW +
                (1|Site/Plot),
              data =  d_Zone,
              family = "binomial",
              start = list(theta = c(1,  1)))

# Compare models. Remove Zone : ubedmean
akaike.weights(c(AIC(m1z), AIC(m2z), AIC(m3z), AIC(m4z)))

# TABLE S3 - COMPARISON 2
# no Zone : DHW 
m3za = glmmTMB(Bleaching_res ~ 
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
                 Zone : Cross_shelf +
                 Zone : Taxa +
                 (1|Site/Plot),
               data =  d_Zone,
               family = "binomial",
               start = list(theta = c(1,  1)))

# no Zone : Cross_shelf
m3zb = glmmTMB(Bleaching_res ~ 
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
                 Zone : DHW +
                 (1|Site/Plot),
               data =  d_Zone,
               family = "binomial",
               start = list(theta = c(1,  1)))

# Compare models. Remove  Zone : DHW 
akaike.weights(c(AIC(m3z), AIC(m3za), AIC(m3zb)))


# TABLE S3 - COMPARISON 3

# no Zone : Cross_shelf
m3za1 = glmmTMB(Bleaching_res ~ 
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
                  (1|Site/Plot),
                data =  d_Zone,
                family = "binomial",
                start = list(theta = c(1,  1)))

# Compare models. Remove Zone : Cross_shelf
akaike.weights(c(AIC(m3za), AIC(m3za1)))


# TABLE S3 - COMPARISON 4
# Model selection on the rest of the fixed effects (i.e., no interactions)

# no slope
m3za1a = glmmTMB(Bleaching_res ~ 
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
                   (1|Site/Plot),
                 data =  d_Zone,
                 family = "binomial",
                 start = list(theta = c(1,  1)))

# no rugosity
m3za1b = glmmTMB(Bleaching_res ~ 
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
                   (1|Site/Plot),
                 data =  d_Zone,
                 family = "binomial",
                 start = list(theta = c(1,  1)))

# no DHW
m3za1c = glmmTMB(Bleaching_res ~ 
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
                   (1|Site/Plot),
                 data =  d_Zone,
                 family = "binomial",
                 start = list(theta = c(1,  1)))

# Compare models.  Remove slope.
akaike.weights(c( AIC(m3za1), AIC(m3za1a), AIC(m3za1b), AIC(m3za1c)))


# TABLE S3 - COMPARISON 5
# no rugosity
m3za1a1 = glmmTMB(Bleaching_res ~ 
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
                    (1|Site/Plot),
                  data =  d_Zone,
                  family = "binomial",
                  start = list(theta = c(1,  1)))

# no DHW
m3za1a2 = glmmTMB(Bleaching_res ~ 
                    log_area +
                    Taxa +
                    log_area : Taxa +
                    ubedmean +
                    Cross_shelf +
                    ubedmean : Cross_shelf +
                    fd_rugosity_plane +
                    rugosity_plane +
                    Zone +
                    Zone : Taxa +
                    (1|Site/Plot),
                  data =  d_Zone,
                  family = "binomial",
                  start = list(theta = c(1,  1)))

# Compare models. Remove rugosity
akaike.weights(c( AIC(m3za1a),AIC(m3za1a1), AIC(m3za1a2)))


# TABLE S3 - COMPARISON 6

# no DHW
m3za1a1a = glmmTMB(Bleaching_res ~ 
                     log_area +
                     Taxa +
                     log_area : Taxa +
                     ubedmean +
                     Cross_shelf +
                     ubedmean : Cross_shelf +
                     fd_rugosity_plane +
                     Zone +
                     Zone : Taxa +
                     (1|Site/Plot),
                   data =  d_Zone,
                   family = "binomial",
                   start = list(theta = c(1,  1)))

# no FD
m3za1a1b = glmmTMB(Bleaching_res ~ 
                     log_area +
                     Taxa +
                     log_area : Taxa +
                     ubedmean +
                     Cross_shelf +
                     ubedmean : Cross_shelf +
                     DHW +
                     Zone +
                     Zone : Taxa +
                     (1|Site/Plot),
                   data =  d_Zone,
                   family = "binomial",
                   start = list(theta = c(1,  1)))

# Compare models. Remove DHW
akaike.weights(c(AIC(m3za1a1), AIC(m3za1a1a), AIC(m3za1a1b)))


# TABLE S3 - COMPARISON 7
# no FD
m3za1a1a1 = glmmTMB(Bleaching_res ~ 
                      log_area +
                      Taxa +
                      log_area : Taxa +
                      ubedmean +
                      Cross_shelf +
                      ubedmean : Cross_shelf +
                      Zone +
                      Zone : Taxa +
                      (1|Site/Plot),
                    data =  d_Zone,
                    family = "binomial",
                    start = list(theta = c(1,  1)))

# Compare models. Retain FD.
akaike.weights(c( AIC(m3za1a1a), AIC(m3za1a1a1)))


################################################################################
#           Model fitting - other response variables, deep vs. shallow         #
################################################################################


# Best fit model includes log_area, taxa, ubedmean, cross_shelf, FD, depth, 
# an interaction between taxa and log_area, and an interaction between taxa and
# depth
m0z = m3za1a1a 

# Use this same structure to fit additional models with different response variables.

m1z = glmmTMB(Bleaching_minor ~ 
                log_area +
                Taxa +
                log_area : Taxa +
                ubedmean +
                Cross_shelf +
                ubedmean : Cross_shelf +
                fd_rugosity_plane +
                Zone +
                Zone : Taxa +
                (1|Site/Plot),
              data =  d_Zone,
              family = "binomial",
              start = list(theta = c(1,  1)))

m2z = glmmTMB(Bleaching_major ~ 
                log_area +
                Taxa +
                log_area : Taxa +
                ubedmean +
                Cross_shelf +
                ubedmean : Cross_shelf +
                fd_rugosity_plane +
                Zone +
                Zone : Taxa +
                (1|Site/Plot),
              data =  d_Zone,
              family = "binomial",
              start = list(theta = c(1,  1)))

m3z = glmmTMB(Bleaching_fully ~ 
                log_area +
                Taxa +
                log_area : Taxa +
                ubedmean +
                Cross_shelf +
                ubedmean : Cross_shelf +
                fd_rugosity_plane +
                Zone +
                Zone : Taxa +
                (1|Site/Plot),
              data =  d_Zone,
              family = "binomial",
              start = list(theta = c(1,  1)))



################################################################################
#                  Plot model predictions - deep vs shallow                    #
################################################################################

zone_df = predict_response(m0z, terms = c("Zone", "Taxa"), margin = "marginalmeans")
zone_df$model = "all"
zone_df1 = predict_response(m1z, terms = c("Zone", "Taxa"), margin = "marginalmeans")
zone_df1$model = "minor"
zone_df2 = predict_response(m2z, terms = c("Zone", "Taxa"), margin = "marginalmeans")
zone_df2$model = "major"
zone_df3 = predict_response(m3z, terms = c("Zone", "Taxa"), margin = "marginalmeans")
zone_df3$model = "fully"

zone_df = rbind(zone_df, zone_df1, zone_df2, zone_df3)
zone_df = zone_df[!is.na(zone_df$group), ]
zone_df$Taxa = zone_df$group


d_Zone$Zone_x = ifelse(d_Zone$Zone == "S", 0, 1)
d_Zone0 = d_Zone[complete.cases(d_Zone$Bleaching_minor), ]
d_Zone1 = d_Zone[complete.cases(d_Zone$Bleaching_major), ]

zone_df$Zone_x = ifelse(zone_df$x == "S", 0, 1 )
zone_df$Zone_x = ifelse(zone_df$model == "all", zone_df$Zone_x - 0.15,
                        ifelse(zone_df$model == "minor", zone_df$Zone_x - 0.05,
                               ifelse(zone_df$model == "major", zone_df$Zone_x + 0.05,    
                                      zone_df$Zone_x + 0.15)))

Figure5 <- ggplot() +
  geom_jitter(data = d_Zone0[d_Zone0$Bleaching_minor == 1, ], 
              aes(x = Zone_x, y = Bleaching_minor - 0.055), alpha = 0.1, size = 1, width = 0.25,
              height = 0.025, col = "#F0E442") +
  geom_jitter(data =  d_Zone0[d_Zone0$Bleaching_minor == 0, ], 
              aes(x = Zone_x, y = Bleaching_minor + 0.055), alpha = 0.1, size = 1, width = 0.25,
              height = 0.025, col = "#F0E442") +
  geom_jitter(data = d_Zone1[d_Zone1$Bleaching_major == 1, ], 
              aes(x = Zone_x, y = Bleaching_major ), alpha = 0.1, size = 1, width = 0.25,
              height = 0.025, col =  "#E69F00") +
  geom_jitter(data =  d_Zone1[d_Zone1$Bleaching_major == 0, ], 
              aes(x = Zone_x, y = Bleaching_major ), alpha = 0.1, size = 1, width = 0.25,
              height = 0.025, col =  "#E69F00") +
  geom_jitter(data = d_Zone[d_Zone$Bleaching_fully == 1, ], 
              aes(x = Zone_x, y = Bleaching_fully + 0.055), alpha = 0.1, size = 1, width = 0.25,
              height = 0.025, col =  "#D55E00" ) +
  geom_jitter(data = d_Zone[d_Zone$Bleaching_fully == 0, ], 
              aes(x = Zone_x, y = Bleaching_fully - 0.055), alpha = 0.1, size = 1, width = 0.25,
              height = 0.025, col =  "#D55E00" ) +
  geom_jitter(data = d_Zone[d_Zone$Bleaching_res == 1, ], 
              aes(x = Zone_x, y = Bleaching_res + 0.055*2), alpha = 0.1, size = 1, width = 0.25,
              height = 0.025, col =  "black" ) +
  geom_jitter(data = d_Zone[d_Zone$Bleaching_res == 0, ], 
              aes(x = Zone_x, y = Bleaching_res - 0.055*2), alpha = 0.1, size = 1, width = 0.25,
              height = 0.025, col =  "black" ) +
  geom_errorbar(data = zone_df, aes(x = Zone_x, ymin = conf.low, ymax = conf.high, 
                                    group = paste(model, Taxa)), width = 0.1, linewidth = 0.8) +
  geom_point(data = zone_df, aes(x = Zone_x, y = predicted, fill = model, group = paste(model, Taxa)),
             size = 3, pch = 21) + 
  scale_colour_manual(values = colours_sev, limits = c("all", "minor", "major", "fully"),
                      labels = c("any bleaching", "minor bleaching", "major bleaching",
                                 "severe bleaching"),
                      name = "") +
  scale_fill_manual(values = colours_sev, limits = c("all", "minor", "major", "fully"),
                    labels = c("any bleaching", "minor bleaching", "major bleaching",
                               "severe bleaching"),
                    name = "") +
  xlab("Depth") +
  ylab("Bleaching severity") +
  theme_classic() +
  theme(text = element_text(size = 10), panel.grid.major.y = element_line(linetype = "dashed")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  scale_x_continuous(breaks = c(0,1), labels = c("shallow", "deep")) +
  facet_wrap(~ Taxa, labeller =  as_labeller(taxa_labels), ) 
Figure5

# save depth plot
#ggsave("Outputs/Figure5.png", Figure5, width = 17, height = 15, units = "cm", dpi = 300)







################################################################################
#                               Save model outputs                             #
################################################################################

#write.table(summary(m0)$coefficients$cond, "Outputs/m0.csv", sep = ",", row.names = TRUE )
#write.table(summary(m1)$coefficients$cond, "Outputs/m1.csv", sep = ",", row.names = TRUE )
#write.table(summary(m2)$coefficients$cond, "Outputs/m2.csv", sep = ",", row.names = TRUE )
#write.table(summary(m3)$coefficients$cond, "Outputs/m3.csv", sep = ",", row.names = TRUE )

#write.table(summary(m0z)$coefficients$cond, "Outputs/m0z.csv", sep = ",", row.names = TRUE )
#write.table(summary(m1z)$coefficients$cond, "Outputs/m1z.csv", sep = ",", row.names = TRUE )
#write.table(summary(m2z)$coefficients$cond, "Outputs/m2z.csv", sep = ",", row.names = TRUE )
#write.table(summary(m3z)$coefficients$cond, "Outputs/m3z.csv", sep = ",", row.names = TRUE )

