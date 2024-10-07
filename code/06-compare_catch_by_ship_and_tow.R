
# load data ####
# 1) Load functions and fishdata filtered for 3 most abundant families
source("code/01-calc_c_flux.R")

# the plots below will show results for all fish identified to family level. 
# note that carbon flux results are based only on 3 families, so if desired, 
# use line below to filter results for just these 3 families

length(fishdata$Ship) # all fish with weights available
length(fishdata_3_families$Ship) # just 3 families

# use fishdata_3_families for plotting results of just Myctophidae, 
# Gonostomatidae and Sternoptychidae (which make up ~ 90% of individuals in the
# catch across both ships)

# note that results are slightly different if plotting all fish families. 
# Plotting all fish families is useful for looking at catch by depth, ship
# and tow, but not as useful for fish C flux calculations because mass
# is needed and we used a length-weight regression for finding weights, which
# required a sufficient number of individuals to make that regression. DVM
# behavior must also be known for C flux calculations, which is why we just use
# these more abundant and well-studied fish families in the C flux analysis. 

# include tows 1, 2, 4, 5, 6, 7, 8 and 9 in plotting Sarmiento tows. Tow
# 3 only fished to 500 m (note that net 1 cannot be used for plotting distribution
# by depth because nets fished just from 1000-500 m and then from 500 m to surface)
# so tows 1, 5 and 7 were the Sarmiento day tows that worked the best
# from Cook, tows 10, 11, 47, 48, 89 and 90 worked well (3 paired day/night tows)

# filter by tow number
fishdata <- fishdata %>% filter(Tow_number %in% c(1, 2, 4, 5, 6, 7, 8, 9, 
                                                  10, 11, 47, 48, 89, 90))

# filter by tow number
fishdata_3_families <- fishdata_3_families %>% filter(Tow_number %in% c(1, 2, 4, 5, 6, 7, 8, 9, 
                                                  10, 11, 47, 48, 89, 90))


# for use later: filter by night tows only
fishdata_night <- fishdata %>% filter(DayNight == "Night")

# 2) Add a column for net type, for plotting
fishdata <- fishdata %>% mutate(Net_type = case_when(
  Ship == "Sarmiento"  ~
    "10 sq. meter MOCNESS",
  Ship == "Cook"  ~
    "1 sq. meter MOCNESS"))

fishdata_night <- fishdata_night %>% mutate(Net_type = case_when(
  Ship == "Sarmiento"  ~
    "10 sq. meter MOCNESS",
  Ship == "Cook"  ~
    "1 sq. meter MOCNESS"))

fishdata_3_families <- fishdata_3_families %>% mutate(Net_type = case_when(
  Ship == "Sarmiento"  ~
    "10 sq. meter MOCNESS",
  Ship == "Cook"  ~
    "1 sq. meter MOCNESS"))


# plot catch by tow and ship ####

# plot Cook and Sarmiento tows together (night and day tows separately) in g/m^2
# (grams of fish wet weight per meter squared of sea surface)

catch_by_tow_all_m2_day_night <- ggplot(fishdata, aes(x = Tow_number, y = weight.to.use/VolFiltM3*DepthInterval, fill = Net_type)) +
  geom_bar(stat = "identity") +
  labs(x = "Tow number", y = expression(paste("Fish biomass (g m"^"-2", ")"))) + theme(legend.title=element_blank(), legend.text = element_text(size = 26), legend.position = "right", strip.background.x = element_blank(), strip.text.x = element_text(size = 26), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 30), panel.spacing = unit(1.5, 'lines'), axis.title.y = element_text(vjust=2), axis.title.x = element_text(vjust=-2), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "inches")) + 
  facet_grid(~DayNight, scales = "free", space = "free") + 
  scale_fill_manual(values = c("#7570B3", "#66A61E")) + scale_color_manual(values = c("#7570B3", "#66A61E")) +
  scale_y_continuous(breaks = c(0, 1, 2, 3), limits = c(0, 3)) 

catch_by_tow_all_m2_day_night


# now make same plot as above but with just the 3 fish families we're focused on
# use the figure below in paper ** 

catch_by_tow_3families_m2_day_night <- ggplot(fishdata_3_families, aes(x = Tow_number, y = weight.to.use/VolFiltM3*DepthInterval, fill = Net_type)) +
  geom_bar(stat = "identity") +
  labs(x = "Tow number", y = expression(paste("Fish biomass (g m"^"-2", ")"))) + 
  theme(legend.title=element_blank(), legend.text = element_text(size = 20), 
        legend.position = "right", strip.background.x = element_blank(), 
        strip.text.x = element_text(size = 26), panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), text = element_text(size = 24), 
        panel.spacing = unit(1.5, 'lines'), axis.title.y = element_text(vjust=2), 
        axis.title.x = element_text(vjust=-2), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             "inches")) + 
  facet_grid(~DayNight, scales = "free", space = "free") + 
  scale_fill_manual(values = c("#7570B3", "#66A61E")) + scale_color_manual(values = c("#7570B3", "#66A61E")) +
  scale_y_continuous(breaks = c(0, 1, 2, 3), limits = c(0, 3)) 

catch_by_tow_3families_m2_day_night

# repeat but now manually change tow numbers to days that correspond to date of tow
# for use in manuscript

# could add a column to fishdata_3_families with date of tow where day is the
# date in May when the tow began 
# Note: have to leave a space (" ") before the Sarmiento tows to get this to plot correctly 
fishdata_3_families_date_in_May <- fishdata_3_families %>% mutate(Tow_start_date_in_May = case_when(
  Tow_number == 1  ~
    " 06",
  Tow_number == 2  ~
    " 06",
  Tow_number == 4  ~
    " 13",
  Tow_number == 5  ~
    " 13",
  Tow_number == 6  ~
    " 14",
  Tow_number == 7  ~
    " 17",
  Tow_number == 8  ~
    " 17",
  Tow_number == 9  ~
    " 18",
  Tow_number == 10  ~
    "06",
  Tow_number == 11  ~
    "07", # at 2am on May 7, so this pairs with day tow on May 6
  Tow_number == 47  ~
    "17",
  Tow_number == 48  ~
   "17",
  Tow_number == 89  ~
    "26",
  Tow_number == 90  ~
    "26"))

# or, for manuscript, add a column with date corresponding to 
# day 1, 2, 3, 4 or 5 of a paired tow where May 6 = day 1, May 13 = day 2, 
# May 14 = Day 3, May 17 = Day 4, May 18 = Day 5, May 26 = Day 6
#fishdata_3_families_date_labels <- fishdata_3_families %>% mutate(Tow_start_date_in_May = case_when(
#  Tow_number == 1  ~
#    "1",
#  Tow_number == 2  ~
#    "1",
#  Tow_number == 4  ~
#    "2",
#  Tow_number == 5  ~
#    "2",
#  Tow_number == 6  ~
#    "3",
#  Tow_number == 7  ~
#    "4",
#  Tow_number == 8  ~
#    "4",
#  Tow_number == 9  ~
#    "5",
#  Tow_number == 10  ~
#    "1",
#  Tow_number == 11  ~
#    "1",
#  Tow_number == 47  ~
#    "4",
#  Tow_number == 48  ~
#    "4",
#  Tow_number == 89  ~
#    "6",
#  Tow_number == 90  ~
#    "6"))

catch_by_tow_3families_m2_day_night_date_labels <- ggplot(
  fishdata_3_families_date_in_May, aes(x = Tow_start_date_in_May, y = weight.to.use/VolFiltM3*DepthInterval, 
                                       fill = Net_type)) +
  geom_bar(stat = "identity") +
  labs(x = "Date in May", y = expression(paste("Fish biomass (g m"^"-2", ")"))) + 
  theme(legend.title=element_blank(), legend.text = element_text(size = 20), 
        legend.position = "right", strip.background.x = element_blank(), 
        strip.text.x = element_text(size = 26), panel.background = element_blank(), panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), text = element_text(size = 24), 
        panel.spacing = unit(1.5, 'lines'), axis.title.y = element_text(vjust=2), 
        axis.title.x = element_text(vjust=-2), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),
                                                                  "inches")) + 
  facet_grid(~DayNight, scales = "free", space = "free") +
  scale_fill_manual(breaks = c("10 sq. meter MOCNESS", "1 sq. meter MOCNESS"), values = c("#7570B3", "#66A61E")) +
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5), limits = c(0, 2.5))  

catch_by_tow_3families_m2_day_night_date_labels

ggsave("graphics/submission_graphics/catch_by_tow_3_families.jpg", catch_by_tow_3families_m2_day_night_date_labels, width = 13, height = 7, units = "in", dpi = 500)


# ^ this version of the plot is used in manuscript ** 

### two-way ANOVA for ship and day/night effects ####

# determine whether DayNight or Ship categories (or their interaction) have a significant effect on catch

# make summary table with columns for Ship, DayNight, and weight.to.use/VolFiltM3*DepthInterval
# and rows for each tow
fishdata_3_families_summary <- fishdata_3_families %>% 
  group_by(Ship, DayNight, Tow_number) %>% 
  summarise(weight.to.use.VolFiltM3.DepthInterval = sum(weight.to.use)/sum(VolFiltM3*DepthInterval))

# run two-way ANOVA
fishdata_3_families_summary_aov <- aov(weight.to.use.VolFiltM3.DepthInterval ~ Ship + DayNight, 
                                       data = fishdata_3_families_summary)
summary(fishdata_3_families_summary_aov)

# check if interaction term is significant
summary(aov(weight.to.use.VolFiltM3.DepthInterval ~ Ship * DayNight + Ship:DayNight, data = fishdata_3_families_summary))

# Yes, interaction term is significant (interaction between net size and time of day). 
# Effects on biomass of net size, time of day (day or night), and the interaction 
# between net size and time of day were all significant (two-way ANOVA, p < 0.05). 


# Cook's smaller MOC-1 unexpectedly caught far more during the day than the MOC-10.
# Catches were more comparable at night for both nets. 

# plot histogram of size distribution ####

# plot histogram with Sarmiento and Cook fish data overlapping in a single plot
# with counts per m squared sampled. Includes all tows with reliable volume and depth
# interval sampled (even though some of these tows sampled different depth intervals), 
# which are tows 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 47, 48, 89 and 90.

hist_Sarmiento_Cook_stand_m2 <- ggplot(data = fishdata, 
                                       aes(x = Std_length_mm, color = Net_type, 
                                           fill = Net_type, weight = 1/VolFiltM3*DepthInterval)) + 
  geom_histogram(binwidth = 2, alpha = 0.5, position = "identity") + theme_bw() + 
  theme(legend.title=element_blank(), legend.text = element_text(size = 26), legend.position = "right", 
        strip.background.x = element_blank(), strip.text.x = element_text(size = 26), 
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), text = element_text(size = 30), panel.spacing = unit(1.5, 'lines'), 
        axis.title.y = element_text(vjust=2), axis.title.x = element_text(vjust=-2), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       "inches")) +
  xlab(bquote("Fish standard length (mm)")) + 
  ylab(expression(Number~of~fish~caught~per~m^{2})) + scale_fill_manual(values = c("#7570B3", "#66A61E")) + scale_color_manual(values = c("#47cffc", "#66A61E")) +
  scale_y_continuous(breaks = c(0, 1, 2, 3), limits = c(0, 3.5)) 

# plot
hist_Sarmiento_Cook_stand_m2

# across all tows, including day and night, Cook caught more AND larger individuals! 
# largely because the MOC-1 caught more during the day than the MOC-10. 

# night catches only, plotted in m^2 (use just night tows because these are more similar
# between large and small)

hist_Sarmiento_Cook_stand_night_m2 <- ggplot(data = fishdata_night, 
                                             aes(x = Std_length_mm, color = Net_type, 
                                                 fill = Net_type, weight = 1/VolFiltM3*DepthInterval)) + 
  geom_histogram(binwidth = 2, alpha = 0.5, position = "identity") + theme_bw() + 
  theme(legend.title=element_blank(), legend.text = element_text(size = 26), strip.background.x = element_blank(), 
        strip.text.x = element_text(size = 26), panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), text = element_text(size = 25), 
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)
  ) +
  xlab(bquote("\n Fish standard length (mm)")) + 
  ylab(bquote(Number~of~fish~caught~at~night~per~m^2)) + scale_fill_manual(values = c("#7570B3", "#66A61E")) + 
  scale_color_manual(values = c("#7570B3", "#66A61E")) + ylim(0,2.0)

# plot
hist_Sarmiento_Cook_stand_night_m2

# save
ggsave("graphics/submission_graphics/histogram_size_distribution.jpg", hist_Sarmiento_Cook_stand_night_m2, width = 9, height = 7, units = "in", dpi = 500)


# add column to fishdata_night for g/m2
fishdata_night$standardized_weight <- fishdata_night$weight.to.use/ fishdata_night$VolFiltM3 * 
  fishdata_night$DepthInterval

# summarize fishdata_night by ship and tow_number, summing across weight.to-use 
# for each tow
fishdata_night_summary <- fishdata_night %>% group_by(Ship, Tow_number) %>% 
  summarize(Tow_weight_tow = sum(standardized_weight, na.rm = TRUE))

# factor variation between min and max catch by tow
max(fishdata_night_summary$Tow_weight_tow) / min(fishdata_night_summary$Tow_weight_tow)
# there is over a factor of 3 difference between min and max tow biomass/m^2

# mean biomass/m^2 across nighttime tows
mean(fishdata_night_summary$Tow_weight_tow)

# Cook still seems to catch larger fish and misses the smaller ones, 
# even if included just night catches which had more similar catch among tows 
# between the two ships (except for a few very large individuals from Sarmiento)

# Data suggest that net avoidance is higher for MOC-10 than for 
# MOC-1, which was not our expectation with using a larger net for comparison. 

# best tows for comparison: 
# Sarmiento tow 2, Cook tow 11 (these were conducted at night in space time/place)
# Sarmiento tow 7, Cook tow 47 (these were conducted during day in space time/place)
# Sarmiento tows 2, 4, 6, 8 and 9 were the only night tows conducted on Sarmiento

# Tow 2 is the only Sarmiento tow that fished at night with Cook where all 
# nets fished and sampling extended to 1000 m. 
# Tow 7 was conducted during day but at least was in the same time/place as
# one of the the Cook's daytime tows so we can directly compare catch 
# (avoidance for both nets may be higher because towed during the day). 
# Note that
# tows 7 and 8 may have an abnormally high catch for the top net 9 because 
# this net was held at the surface for ~10 minutes or more before recovery. 
# Nets 8 and 9 were during the edge experiment
# (so, catch may be distinct from catch rates inside the eddy core)

# plot just the best tows for direct comparison between the two nets
fishdata_comparison_tows <- fishdata %>%
  filter(Tow_number == "2" | Tow_number == "7" | Tow_number == "11" | 
           Tow_number == "47")

# plot in g/m^2
same_time_place_plot_m2 <-  ggplot(fishdata_comparison_tows, aes(x = Tow_number, y = weight.to.use/VolFiltM3*DepthInterval, fill = Ship)) +
  geom_bar(stat = "identity") +
  labs(x = "Tow number", y = expression(Fish~biomass~(g~m^{-2}))) +
  facet_wrap(~DayNight, scales = "free") + 
  theme(legend.title=element_blank(), strip.background.x = element_blank(), strip.text.x = element_text(size = 18), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 18)) + scale_fill_manual(values = c("#33638DFF", "#B8DE29FF")) + scale_color_manual(values = c("#33638DFF", "#B8DE29FF")) 

same_time_place_plot_m2


# Observations: 
# - daytime tows from Sarmiento have lower biomass than nighttime tows. Either
# many fish are swimming deeper than 1000 m during the day only, and this is not
# evident in acoustics (no gas filled swim bladders in those fish), or there is 
# much higher net avoidance during the day. in g/m^3
#
# - Before standardizing by volume filtered, Tow 3 that did not fish below 500 m 
# is quite a small catch, as expected. This difference is much less after
# standardizing by volume filtered
#
# Tows 8 and 9, which were held at the surface (where there was a high density 
# of myctophids, at night) have a much higher biomass than other night tows before and
# after standardizing by volume filtered. However, after multiplying by depth
# strata sampled, the biomass for these two tows becomes more comparable to the 
# other night tows. Tows 8 and 9 were also conducted near the eddy edge, not core
#
# Tow 7 was  accidentally held at surface but was conducted during the day, 
# so being held at surface should not have increased catch biomass much (few
# epipelagics were caught in the nets). Not sure why tow 7 is a little higher 
# than the others in catch biomass in g/m^3 and g/m^2--fish patchiness?


# now that dfs from Cook and Sarmiento data are combined, show both in histogram. 
# Plot Sarmiento data first, then Sarmiento plus Cook data

# filter Sarmiento data only
fishdata_Sarmiento <- fishdata %>%
  filter(Ship == "Sarmiento")

# plot all fish lengths, in separate subplot by lowest taxon possible, for all tows
fish_length_all_Sarm <- ggplot(fishdata_Sarmiento, aes(x = Std_length_mm)) +
  geom_histogram(binwidth = 2) +
  labs(x = "Fish length (mm)", y = "Number of fish caught \n Sarmiento only") + 
  facet_wrap(~Lowest_taxon) 

fish_length_all_Sarm

# next, plot all fish lengths, in separate subplot by fish family, for all tows
# first, add column for family

# plot by family and limit x axis to 125 mm for easier data vis (excludes n = 50 fish)
fish_length_family_Sarm <- ggplot(fishdata_Sarmiento, aes(x = Std_length_mm)) +
  geom_histogram(binwidth=2) +
  labs(x = "Fish length (mm)", y = "Number of fish caught \n Sarmiento only") + 
  facet_wrap(~Family) 

fish_length_family_Sarm


# by lowest taxon for Cook and Sarmiento, standarized to g/m^2
fish_length_all_Sarm_Cook <- ggplot(fishdata, aes(x = Std_length_mm, weight = 1/VolFiltM3*DepthInterval)) +
  geom_histogram(binwidth = 2) +
  labs(x = "Fish length (mm)", y = "Number of fish caught \n Cook and Sarmiento") + 
  facet_wrap(~Lowest_taxon) 

fish_length_all_Sarm_Cook

# by family for Cook and Sarmiento
fish_length_family_Sarm_Cook <- ggplot(fishdata, aes(x = Std_length_mm, weight = 1/VolFiltM3*DepthInterval)) +
  geom_histogram(binwidth=2) +
  labs(x = "Fish length (mm)", y = "Number of fish caught \n Cook and Sarmiento") + 
  facet_wrap(~Family) 

fish_length_family_Sarm_Cook


# consider different calibration factor for MOC-1 ####


# re-plot comparison of catch by tow and ship but considering a calibration factor
# for the Cook that is 50% higher than used previously (as calibration factors
# can vary between about 4 and 6, Peter Wiebe pers. comm. Jan 26, 2024 via email)

# 6/4.136 = multiplier to adjust to max feasible calibration factor

# edit volume filtered for Cook rows only
fishdata_3_families_date_in_May_higher_calibration <- fishdata_3_families_date_in_May
fishdata_3_families_date_in_May_higher_calibration$VolFiltM3[fishdata_3_families_date_in_May_higher_calibration$Ship ==
                                                               "Cook"] <- fishdata_3_families_date_in_May_higher_calibration$VolFiltM3[fishdata_3_families_date_in_May_higher_calibration$Ship == 
                                                                                                                                         "Cook"]*1.45
catch_by_tow_3families_m2_day_night_date_labels_higher_calibration <- ggplot(
  fishdata_3_families_date_in_May_higher_calibration, aes(x = Tow_start_date_in_May, y = weight.to.use/VolFiltM3*DepthInterval, 
                                       fill = Net_type)) +
  geom_bar(stat = "identity") +
  labs(x = "Date in May", y = expression(paste("Fish biomass (g m"^"-2", ")"))) + 
  theme(legend.title=element_blank(), legend.text = element_text(size = 20), 
        legend.position = "right", strip.background.x = element_blank(), 
        strip.text.x = element_text(size = 26), panel.background = element_blank(), panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), text = element_text(size = 24), 
        panel.spacing = unit(1.5, 'lines'), axis.title.y = element_text(vjust=2), 
        axis.title.x = element_text(vjust=-2), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),
                                                                  "inches")) + 
  facet_grid(~DayNight, scales = "free", space = "free") + 
  scale_fill_manual(values = c("#7570B3", "#66A61E")) + scale_color_manual(values = c("#7570B3", "#66A61E")) +
  scale_y_continuous(breaks = c(0, 1, 2, 3), limits = c(0, 2)) 

catch_by_tow_3families_m2_day_night_date_labels_higher_calibration

ggsave("graphics/submission_graphics/catch_by_tow_3_families_max_MOC1_cal.jpg", catch_by_tow_3families_m2_day_night_date_labels_higher_calibration, width = 13, height = 7, units = "in", dpi = 500)


# re-do stats in case they change with the higher calibration factor for the MOC-1

# two-way ANOVA to determine whether DayNight or Ship categories have a significant effect on catch

# make summary table with columns for Ship, DayNight, and weight.to.use/VolFiltM3*DepthInterval
# and rows for each tow
fishdata_3_families_summary <- fishdata_3_families_date_in_May_higher_calibration %>% 
  group_by(Ship, DayNight, Tow_number) %>% 
  summarise(weight.to.use.VolFiltM3.DepthInterval = sum(weight.to.use)/sum(VolFiltM3*DepthInterval))

# run two-way ANOVA
fishdata_3_families_summary_aov <- aov(weight.to.use.VolFiltM3.DepthInterval ~ Ship + DayNight, 
                                       data = fishdata_3_families_summary)
summary(fishdata_3_families_summary_aov)

# check if interaction term is signficant
summary(aov(weight.to.use.VolFiltM3.DepthInterval ~ Ship * DayNight + Ship:DayNight, data = fishdata_3_families_summary))

# Interaction term is significant (interaction between net size and time of day) 


# Still, even with calibration factor correction, Cook's smaller MOC-1 unexpectedly 
# caught significantly more during the day than the MOC-10.
# Catches were more comparable at night for both nets. 

### Minor calculation for Introduction ####

# if Proud et al. 2019 estimate that mesopelagic fish biomass is 2-16 billion metric ton globally, 
# and if Jennings et al 2008, Wilson et al 2009, Bianchi et al 2021, and Gjosaeter and Kawaguchi 1980
# estimate that global fish biomass in total is about 1 billion metric ton, then mesopelagic fishes
# comprise anywhere from roughly 67% (2 Gt mesopelagic fishes / 3 Gt fish total) and 94%
# (16 Gt mesopelagic fishes / 17 Gt fish total) of global fish biomass.
