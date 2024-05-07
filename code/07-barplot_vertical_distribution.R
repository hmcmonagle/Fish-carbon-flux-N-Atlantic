# Stacked bar plot of fish abundance by depth with 
# each bar representing the taxonomic breakdown at each depth

# Nov 14, 2023, HM
# Thanks to Tor Mowatt-Larssen for help with this plotting code!

# load functions and data )
source("code/01-calc_C_flux.R")

# (or use csv file fishdata_Cook_summary.csv commented out below)
# fishdata_Cook_summary <- read.csv("data/MOCNESS_data/MOCNESS 1 Cook/fishdata_Cook_summary.csv", header = TRUE)

# use data for just the 3 most abundant families, which make up over 90% of the catch
fishdata <- fishdata_3_families

# Use Cook data for this plot because we have n = 3 night and n = 3 day tows, 
# and the depth intervals have higher resolution (8 depth intervals instead of just 4)

# filter just Cook data
fishdata_Cook <- fishdata %>% filter(Ship == "Cook")

# add column for depth midpoint
fishdata_Cook$DepthMid <- (fishdata_Cook$DepthEnd+fishdata_Cook$DepthStart) /2

fishdata_Cook$biomass_g_m3 <- fishdata_Cook$weight.to.use / fishdata_Cook$VolFiltM3 

# make table that provides sum of catch biomass by family
# and depth midpoint and day/night
fishdata_Cook_summary <- fishdata_Cook %>% group_by(DepthMid, DepthStart, DepthEnd, DepthInterval, DayNight, Family) %>% 
  summarise(biomass_g_m3 = sum(biomass_g_m3))

# write.csv(fishdata_Cook_summary, "fishdata_Cook_summary.csv")

# make Day catch values negative so that they plot separately from night catch
dn_plot <- fishdata_Cook_summary %>%
  mutate(biomass_g_m3 = ifelse(DayNight == "Night",
                                biomass_g_m3,
                                -1*biomass_g_m3))

# plot stacked bar plot showing density by depth and family
# (standardized by volume filtered, g/m^3)
# Depth interval sampled is reflected by the plot so no need to multipply
# by depth interval to get g/m^2)
diel_plot <- ggplot(dn_plot, aes(x = DepthMid, y = biomass_g_m3, fill = Family)) +
  geom_bar(stat = "identity", width = dn_plot$DepthInterval) +
  theme_bw(base_size = 26) +
  coord_flip() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=-30) +
  scale_x_reverse(position = "bottom", breaks=c(0,50,100,150,200,300,400,500,750,1000), labels = 
                    c(0, 50, 100, 150, 200, 300, 400, 500, 750, 1000), limits = c(1000, -80)) + 
  scale_y_continuous(position = "left", breaks = c(-0.02, -0.01, 0, 0.01, 0.02),
                     labels = c(0.02, 0.01, 0, 0.01, 0.02), limits = c(-0.02, 0.025)) + 
  ylab(expression(paste("Fish biomass (g m"^"-3", ")"))) +
  scale_fill_manual(values = c("#7570B3", "#D95F02", "#1B9E77")) +
  xlab("Depth (m)") +
  theme(legend.position = "top", panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(), legend.title=element_blank()) + 
  
  # add annotations for Day vs Night labels
  annotate("text", x = -80, y = -0.011, label = "Day", size = 7.5) +
  annotate("text", x = -80, y = 0.013, label = "Night", size = 7.5) 

diel_plot

ggsave("graphics/submission_graphics/vertical_distribution.jpg", diel_plot, width = 11, height = 9, units = "in", dpi = 500)

## use Cook data in manuscript (depth distributions are higher resolution because there are more nets)

# make same plot but with Sarmiento data
# filter just Sarmiento data
fishdata_Sarmiento <- fishdata %>% filter(Ship == "Sarmiento")

# filter to just include tows 2, 4, 5, 6, 7, 8 and 9
fishdata_Sarmiento <- fishdata_Sarmiento %>% filter(Tow_number %in% c(2, 4, 5, 6, 7, 8, 9))

# add column for depth midpoint
fishdata_Sarmiento$DepthMid <- (fishdata_Sarmiento$DepthEnd+fishdata_Sarmiento$DepthStart) /2

fishdata_Sarmiento$biomass_g_m3 <- fishdata_Sarmiento$weight.to.use / fishdata_Sarmiento$VolFiltM3

# make table that provides sum of catch biomass by family
# and depth midpoint and day/night
fishdata_Sarmiento_summary <- fishdata_Sarmiento %>% group_by(DepthMid, DepthStart, DepthEnd, DepthInterval, DayNight, Family) %>% 
  summarise(biomass_g_m3 = sum(biomass_g_m3))

# write.csv(fishdata_Sarmiento_summary, "fishdata_Sarmiento_summary.csv")

# make Day catch values negative so that they plot separately from night catch
dn_plot_Sarmiento <- fishdata_Sarmiento_summary %>%
  mutate(biomass_g_m3 = ifelse(DayNight == "Night",
                                biomass_g_m3,
                                -1*biomass_g_m3))

# plot stacked bar plot showing density by depth and family
# (standardized by volume filtered, g/m^3)
# Depth interval sampled is reflected by the plot so no need to multipply
# by depth interval to get g/m^2)
diel_plot_Sarmiento <- ggplot(dn_plot_Sarmiento, aes(x = DepthMid, y = biomass_g_m3, fill = Family)) +
  geom_bar(stat = "identity", width = dn_plot_Sarmiento$DepthInterval) +
  theme_bw(base_size = 26) +
  coord_flip() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=-30) +
  scale_x_reverse(position = "bottom", breaks=c(0,100,300,500,1000), labels = 
                    c(0, 100, 300, 500, 1000), limits = c(1000, -80)) + 
  scale_y_continuous(position = "left", "Biomass (g/m^3)", breaks = c(-0.015, 0, 0.015),
                     labels = c(0.015, 0, 0.015), limits = c(-0.015, 0.015)) + 
  scale_fill_manual(values = c("#7570B3", "#D95F02", "#1B9E77")) + 
  xlab("Depth (m)") +
  theme(legend.position = "top", panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(), legend.title=element_blank()) + 
  
  # add annotations for Day vs Night labels
  annotate("text", x = -80, y = -0.01, label = "Day", size = 7.5) +
  annotate("text", x = -80, y = 0.01, label = "Night", size = 7.5)

diel_plot_Sarmiento

