# Bar plots comparing fish flux with other flux pathways#####

# load packages
library(dplyr)
library(ggplot2)

# source previous code (just need a function for plot colors)
source("code/05-fish_diversity.R")

## Find fish flux ####

### Load data ####

# load MC results for VM and NM fish
# these results are for a flux boundary of 200 m (default) unless
# 500 is written in the object name
VM_MC_results <- readRDS("analysis/VM_MC_results.RDS")
NM_MC_results <- readRDS("analysis/NM_MC_results.RDS")

# add a column to VM_MC_results that adds VM and NM fish values in C_flux_site_nsims column
VM_plus_NM_MC_results <- VM_MC_results
VM_plus_NM_MC_results$C_flux_site_VM_NM <- VM_MC_results$C_flux_site_nsims + NM_MC_results$C_flux_site_nsims

# load data for 500 m 
VM_MC_results_500 <- readRDS("analysis/VM_MC_results_500.RDS")
# sum(is.na(VM_MC_results_500$C_flux_site_nsims)) # checks for NA values

NM_MC_results_500 <- readRDS("analysis/NM_MC_results_500.RDS")


### 200 m ####

# 90% range (to remove the unlikely extremes among Monte Carlo simulations)
fish_flux_200_range <- signif(quantile(VM_plus_NM_MC_results$C_flux_site_nsims, probs=c(.05,.95)))
fish_flux_200_range # ** use this in manuscript text

fish_flux_min_200 <- fish_flux_200_range[1]
fish_flux_max_200 <- fish_flux_200_range[2]

# could use midpoint between 5th and 95th percentile as nominal estimate for fish C flux, 
mean(fish_flux_200_range)

# or mean across all Monte Carlo simulations as nominal estimate for plotting (use this in manuscript**)
fish_flux_mean_200 <- mean(VM_plus_NM_MC_results$C_flux_site_nsims)

# could also use median (50th percentile) fish flux
fish_flux_median_200 <- signif(quantile(VM_plus_NM_MC_results$C_flux_site_nsims, probs=0.50))


### 500 m ####

# add a column to VM_MC_results that adds VM and NM fish values in C_flux_site_nsims column
VM_plus_NM_MC_results_500 <- VM_MC_results_500
VM_plus_NM_MC_results_500$C_flux_site_VM_NM <- VM_MC_results_500$C_flux_site_nsims + NM_MC_results_500$C_flux_site_nsims

# 90% range 
fish_flux_500_range <- signif(quantile(VM_plus_NM_MC_results_500$C_flux_site_nsims, probs=c(.05,.95)))
fish_flux_500_range # ** use this in manuscript text

fish_flux_min_500 <- fish_flux_500_range[1]
fish_flux_max_500 <- fish_flux_500_range[2]

# could use midpoint between 5th and 95th percentile as nominal estimate for fish C flux, 
mean(fish_flux_500_range) 

# or, could use mean across all Monte Carlo simulations as nominal estimate for plotting (use this)
fish_flux_mean_500 <- mean(VM_plus_NM_MC_results_500$C_flux_site_nsims)

# median (50th percentile) fish flux
fish_flux_median_500 <- signif(quantile(VM_plus_NM_MC_results_500$C_flux_site_nsims, probs=0.50))

## Find sediment trap particle flux ####

# Note: Sinking particle flux results from the sediment trap approach are uploaded 
# to the Github repository for this paper (Natl_trap_fluxes_share_V2.xlsx) 
# and original data are available at Estapa, M. L. 2022. EXPORTSNA. SeaWIFS Bio-optical 
# Archive and Storage System (SeaBASS). NASA. doi:10.5067/SeaBASS/EXPORTS/DATA001

### Load data ####

passive_flux <- read.csv("data/Passive_flux_sediment_traps/NAtl_trap_fluxes_share_V2_mean_tab.csv")

# passive flux results from the datasheet above, shared by Meg Estapa and 
# Colleen Durkin, include data for 3 epochs: 1, 2, and 3.
# These time periods are based on conditions, but we'll consider fish flux
# across these three time periods and compare to a passive flux estimate
# averaging across all 3 epochs.


### 200 m ####

# mean flux where depth nominal is equal to 175 or 195, from sediment traps (~200 m)
passive_flux_200_data <- passive_flux %>% filter(Depth_nominal_m == 195 | Depth_nominal_m == 175)

# create new column and convert from mmol C m^-2 d^-1 to 
# mg C m^-2 d^-1
# 1 mmol C = 12 mg C
passive_flux_200_data$C_flux_mgC_m2_d <- passive_flux_200_data$C_flux_mmolC_m2_d * 12
passive_flux_200_data$C_flux_uncertainty_mgC_m2_d <- passive_flux_200_data$C_flux_uncertainty * 12

# group by epoch and calculate mean
passive_flux_200_by_epoch <- passive_flux_200_data %>% group_by(Epoch) %>% 
  summarise(C_flux_mgC_m2_d = mean(C_flux_mgC_m2_d), 
            uncertainty = mean(C_flux_uncertainty_mgC_m2_d))

# plot mean by epoch with error bars to show measurement uncertainty
ggplot(passive_flux_200_by_epoch, aes(x=Epoch, y=C_flux_mgC_m2_d)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(ymin=C_flux_mgC_m2_d-uncertainty, ymax=C_flux_mgC_m2_d+uncertainty), width=.2, position=position_dodge(.9)) +
  labs(title="Passive flux at 200 m", x="Epoch", y="C flux (mg C m^-2 d^-1)")

# note high variation in passive flux among epochs (compared to relatively small measurement error within epochs,
# though epoch 2 only has one )

# find mean across all 3 epoch estimates
mean_passive_flux_200 <- mean(passive_flux_200_by_epoch$C_flux_mgC_m2_d)

# min/max flux by epoch
passive_flux_200_min <- min(passive_flux_200_by_epoch$C_flux_mgC_m2_d)
passive_flux_200_max <- max(passive_flux_200_by_epoch$C_flux_mgC_m2_d)

# range in flux across epochs
passive_flux_200_range <- c(passive_flux_200_min, passive_flux_200_max)

# range of passive flux estimates at 200 m: 43.68-236.94 mg C m^-2 d^-1


### 500 m ####

# mean flux where depth nominal is equal to 500 m, from sediment traps
passive_flux_500_data <- passive_flux %>% filter(Depth_nominal_m == 500)

# create new column and convert from mmol C m^-2 d^-1 to 
# mg C m^-2 d^-1
# 1 mmol C = 12 mg C
passive_flux_500_data$C_flux_mgC_m2_d <- passive_flux_500_data$C_flux_mmolC_m2_d * 12
passive_flux_500_data$C_flux_uncertainty_mgC_m2_d <- passive_flux_500_data$C_flux_uncertainty * 12

# group by epoch and calculate mean
passive_flux_500_by_epoch <- passive_flux_500_data %>% group_by(Epoch) %>% 
  summarise(C_flux_mgC_m2_d = mean(C_flux_mgC_m2_d), 
            uncertainty = mean(C_flux_uncertainty_mgC_m2_d))

# plot mean by epoch with error bars to show measurement uncertainty
ggplot(passive_flux_500_by_epoch, aes(x=Epoch, y=C_flux_mgC_m2_d)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(ymin=C_flux_mgC_m2_d-uncertainty, ymax=C_flux_mgC_m2_d+uncertainty), width=.2, position=position_dodge(.9)) +
  labs(title="Passive flux at 500 m", x="Epoch", y="C flux (mg C m^-2 d^-1)")

# note high variation in passive flux among epochs (compared to relatively small measurement error within epochs,
# though epoch 2 only has one )

# find mean across all 3 epoch estimates
mean_passive_flux_500 <- mean(passive_flux_500_by_epoch$C_flux_mgC_m2_d)

# min/max flux by epoch
passive_flux_500_min <- min(passive_flux_500_by_epoch$C_flux_mgC_m2_d)
passive_flux_500_max <- max(passive_flux_500_by_epoch$C_flux_mgC_m2_d)

# range in flux across epochs
passive_flux_500_range <- c(passive_flux_500_min, passive_flux_500_max)

# range of passive flux estimates at 500 m: 31.56-71.64 mg C m^-2 d^-1


## Find Th-234 particle flux ####

### Load data ####
# from Clevenger et al 2024, Table 4

# ** edit code later to pull from csv file from Clevenger et al. 2024 if available. 
# For now, these values are manually copied from Table 4 in Clevenger et al. 2024 

# make data frame with depth, epoch, carbon flux and uncertainty (SD)
Th234_flux <- data.frame(Depth = c(200, 200, 200, 500, 500, 500), 
                         Flux_mgC_m2_d = c(10.7, 11.8, 14.7, 8.9, 11.9, 16.7) * 12, 
                         Uncertainty = c(2.0, 3.0, 2.4, 5.6, 6.2, 5.1) * 12)

### 200 m ####
Th234_flux_200 <- Th234_flux %>% filter(Depth == 200)


### 500 m ####
Th234_flux_500 <- Th234_flux %>% filter(Depth == 500)

## Find zooplankton flux ####

### Load data ####

# Results from Amy Maas and Deb Steinberg:
# In zoop flux, include respiration, DOC, and POC, but not mortality flux
# to avoid double counting the zoop flux mortality in the fish flux 

zoop_flux <- read.csv("data/Zoop_biomass_and_flux/Zoop_C_flux_Maas_TableSG2.csv")

### 200 m ####

# sum respiration, DOC, POC fluxes (exclude mortality flux to avoid double-counting
# with fish flux) to get total zoop flux (total C is already in csv file)
zoop_flux_200 <- zoop_flux %>% filter(Flux_boundary == "200")

# find mean zoop flux
mean(zoop_flux_200$Total_C)

# min flux 
zoop_flux_200_min <- min(zoop_flux_200$Total_C)

# max zoop flux 
zoop_flux_200_max <- max(zoop_flux_200$Total_C)


### 500 m ####

zoop_flux_500 <- zoop_flux %>% filter(Flux_boundary == "500")

# find mean zoop flux
mean(zoop_flux_500$Total_C)

# min flux 
zoop_flux_500_min <- min(zoop_flux_500$Total_C)

# max zoop flux 
zoop_flux_500_max <- max(zoop_flux_500$Total_C)


## Generate data frame of results and plot ####

# Mean value for passive flux from sediment traps and Th-234 = average of 3 epochs,
# and error bars will range from epoch 1 to epoch 3 estimate
# For fish flux, error bars will range from the 5th to 95th percentile estimate 
# and colored bar will show median (50th percentile)
# For zoop flux, error bar will show range from 3 tows (3 diff biomass estimates)
# (see figure caption for additional details) 

# divide uncertainty values for sediment trap and Th-234 by 2 because will plot
# error bars with min/max values from the mean shown by colored bar (mean + or - 
# half the uncertainty)
# see code for "uncertainty_min or "uncertainty_max" columns below

flux_pathways <- data.frame(
  pathway = c("particles \n traps \n time 1", 
              "particles \n traps \n time 2",
              "particles \n traps \n time 3",
              "particles \n Th-234 \n time 1", 
              "particles \n Th-234 \n time 2", 
              "particles \n Th-234 \n time 3", 
              "zooplankton", "fish"), 
  carbon_flux = c(passive_flux_200_by_epoch$C_flux_mgC_m2_d, Th234_flux_200$Flux_mgC_m2_d,
                  mean(zoop_flux_200$Total_C), fish_flux_mean_200), 
  uncertainty_min = c(passive_flux_200_by_epoch$C_flux_mgC_m2_d - passive_flux_200_by_epoch$uncertainty/2, 
                      Th234_flux_200$Flux_mgC_m2_d - Th234_flux_200$Uncertainty/2, 
                      zoop_flux_200_min, fish_flux_min_200),
  uncertainty_max = c(passive_flux_200_by_epoch$C_flux_mgC_m2_d + passive_flux_200_by_epoch$uncertainty/2, 
                      Th234_flux_200$Flux_mgC_m2_d + Th234_flux_200$Uncertainty/2, 
                      zoop_flux_200_max, fish_flux_max_200)
  )

colors <- c("#7570B3","#1b9e77","#66A61E", "#D95F02")

# Lighten or darken the original color (made this function in 05-fish_diversity.R script)
darker_colors <- darken_color(colors)

lighter_colors <- lighten_color(colors)

darkest_colors <- darken_color(darker_colors)

# bar plot of fish compared to other C flux pathways
bar_by_flux_pathway_200 <- ggplot(
  data = flux_pathways, aes(x = pathway, y = carbon_flux, fill = pathway)) + 
  geom_bar(stat = "identity") + theme_bw() + 
  theme(plot.margin = margin(10, 10, 5, 50), legend.title=element_blank(), legend.position = "none", 
        strip.background.x = element_blank(), 
        strip.text.x = element_text(size = 26), panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), text = element_text(size = 24)) +
  xlab(bquote("  \n  ")) + 
  ylab(expression("Carbon flux 200m \n     (mg m-2 d-1)")) + 
  scale_fill_manual(values=c("#7570B3", lighter_colors[2], "#1b9e77", darkest_colors[2], 
                             lighter_colors[3], "#66A61E", darkest_colors[3], "#D95F02")) + 
  ylim(0, 260) + annotate("text", x = 0.7, y = 250, label = "a)", size = 8)

# add error bars to plot using passive_flux_200_uncertainty for passive flux, 
# zoop_uncertainty_200 for zoop flux, and fish_flux_uncertainty_200 for fish flux
bar_by_flux_pathway_with_error_200 <- bar_by_flux_pathway_200 + 
  geom_errorbar(aes(ymin = uncertainty_min, ymax = uncertainty_max), 
                width = 0.2, position = position_dodge(0.9))

# plot
bar_by_flux_pathway_with_error_200


#### re-create flux pathway comparison but at 500 m ####

flux_pathways_500 <- data.frame(
  pathway = c("particles \n traps \n time 1", 
              "particles \n traps \n time 2",
              "particles \n traps \n time 3",
              "particles \n Th-234 \n time 1", 
              "particles \n Th-234 \n time 2", 
              "particles \n Th-234 \n time 3", 
              "zooplankton", "fish"), 
  carbon_flux = c(passive_flux_500_by_epoch$C_flux_mgC_m2_d, Th234_flux_500$Flux_mgC_m2_d,
                  mean(zoop_flux_500$Total_C), fish_flux_mean_500), 
  uncertainty_min = c(passive_flux_500_by_epoch$C_flux_mgC_m2_d - passive_flux_500_by_epoch$uncertainty/2, 
                      Th234_flux_500$Flux_mgC_m2_d - Th234_flux_500$Uncertainty/2, 
                      zoop_flux_500_min, fish_flux_min_500),
  uncertainty_max = c(passive_flux_500_by_epoch$C_flux_mgC_m2_d + passive_flux_500_by_epoch$uncertainty/2, 
                      Th234_flux_500$Flux_mgC_m2_d + Th234_flux_500$Uncertainty/2, 
                      zoop_flux_500_max, fish_flux_max_500)
)

# use same colors as for 200 m plot above

# bar plot of fish compared to other C flux pathways
bar_by_flux_pathway_500 <- ggplot(
  data = flux_pathways_500, aes(x = pathway, y = carbon_flux, fill = pathway)) + 
  geom_bar(stat = "identity") + theme_bw() + 
  theme(plot.margin = margin(10, 10, 5, 50), legend.title=element_blank(), legend.position = "none", 
        strip.background.x = element_blank(), 
        strip.text.x = element_text(size = 26), panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), text = element_text(size = 24)) +
  xlab(bquote("  \n  ")) + 
  ylab(expression("Carbon flux 500m \n     (mg m-2 d-1)")) + 
  scale_fill_manual(values=c("#7570B3", lighter_colors[2], "#1b9e77", darkest_colors[2], 
                             lighter_colors[3], "#66A61E", darkest_colors[3], "#D95F02")) + 
  ylim(0, 260) + annotate("text", x = 0.7, y = 250, label = "b)", size = 8)

# add error bars to plot using passive_flux_500_uncertainty for passive flux, 
# zoop_uncertainty_500 for zoop flux, and fish_flux_uncertainty_500 for fish flux
bar_by_flux_pathway_with_error_500 <- bar_by_flux_pathway_500 + 
  geom_errorbar(aes(ymin = uncertainty_min, ymax = uncertainty_max), 
                width = 0.2, position = position_dodge(0.9))

# plot
bar_by_flux_pathway_with_error_500



# put bar_by_flux_pathway_200 and bar_by_flux_pathway_500 side by side using ggarange
# install.packages("egg")
#install.packages("ggpubr")
library(ggpubr)

flux_comparison_panels <- ggarrange(bar_by_flux_pathway_with_error_200, bar_by_flux_pathway_with_error_500, 
                                    ncol=1, nrow=2, common.legend = TRUE, legend="none")
flux_comparison_panels


#### Alternative plot ####

# Bar plot above with error bars might be misleading because the error bars do
# not have the same meaning for all bars, and people tend to make assumptions
# about what error bars mean. Instead, try making a similar plot but with labels
# on the y-axis (for a cleaner look, as it's hard to fit the labels along the x-axis)
# and with each category as a set of 3 dots, with the first dot being the minimum
# carbon flux estimate, the middle dot being the nominal, and the 3rd dot being the maximum. 

##### 200 m #####

# use flux_pathways and make separate rows in the plot for each pathway. Then, 
# for each row, plot 3 dots for the minimum, nominal, and maximum carbon flux
# estimates. Connect the 3 dots with a line between the minimum and maximum

# re-word pathways so they don't create returns in y axis labels
flux_pathways$pathway <- gsub("particles \n traps \n time 1", "sinking particles from sediment traps, t1", flux_pathways$pathway)
flux_pathways$pathway <- gsub("particles \n traps \n time 2", "sinking particles from sediment traps, t2", flux_pathways$pathway)
flux_pathways$pathway <- gsub("particles \n traps \n time 3", "sinking particles from sediment traps, t3", flux_pathways$pathway)
flux_pathways$pathway <- gsub("particles \n Th-234 \n time 1", "sinking particles from Th-234, t1", flux_pathways$pathway)
flux_pathways$pathway <- gsub("particles \n Th-234 \n time 2", "sinking particles from Th-234, t2", flux_pathways$pathway)
flux_pathways$pathway <- gsub("particles \n Th-234 \n time 3", "sinking particles from Th-234, t3", flux_pathways$pathway)

# reorder y axis labels so that time 1 comes at top before time 3
flux_pathways$pathway <- factor(flux_pathways$pathway, levels = c("sinking particles from sediment traps, t3",
                                                                  "sinking particles from sediment traps, t2", 
                                                                  "sinking particles from sediment traps, t1", 
                                                                  "sinking particles from Th-234, t3", 
                                                                  "sinking particles from Th-234, t2",
                                                                  "sinking particles from Th-234, t1",
                                                                  "zooplankton", "fish"))


# make a dot plot with 3 dots for each pathway
dot_plot_200m <- ggplot(data = flux_pathways, aes(x = carbon_flux, y = pathway)) + 
  geom_point(aes(color = pathway, size = 2)) + theme_bw() + 
  theme(plot.margin = margin(10, 10, 5, 50), legend.title=element_blank(), legend.position = "none", 
        strip.background.x = element_blank(), 
        strip.text.x = element_text(size = 26), panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), text = element_text(size = 24)) +
  xlab(bquote("Carbon flux 200m (mg m-2 d-1)")) + 
  ylab(expression("  \n  ")) + 
  scale_color_manual(values=c(darker_colors[2],"#1b9e77", lighter_colors[2], 
                              darker_colors[3], "#66A61E", lighter_colors[3], "#D95F02", "#7570B3")) + 
  annotate("text", x = 240, y = 8, label = "a)", size = 8) + 
  geom_point(aes(x = uncertainty_min, y = pathway, color = pathway, size = 2)) + 
  geom_point(aes(x = uncertainty_max, y = pathway, color = pathway, size = 2)) + 
  geom_segment(aes(x = uncertainty_min, xend = uncertainty_max, y = pathway, yend = pathway, color = pathway), size = 1.5) +
  xlim(0, 250)

# print plot  
dot_plot_200m

# for simplified figure in presentation, plot just fish, zoop and sediment trap data
# and remove the min/max dots

# subset to exclude Th-234 data
flux_pathways_simple <- flux_pathways[flux_pathways$pathway %in% c("sinking particles from sediment traps, t1", 
                                                                   "sinking particles from sediment traps, t2", 
                                                                   "sinking particles from sediment traps, t3", 
                                                                   "zooplankton", "fish"),]
# rename "sinking particles from sediment traps" to just "sinking particles"
# rename flux_pathways_simple$pathway[1] to "sinking particles, time period 1"
# rename flux_pathways_simple$pathway[2] to "sinking particles, time period 2"
# rename flux_pathways_simple$pathway[3] to "sinking particles, time period 3"
flux_pathways_simple$pathway <- gsub("sinking particles from sediment traps", "sinking particles", flux_pathways_simple$pathway)

# make order of pathways consistent with previous plot
flux_pathways_simple$pathway <- factor(flux_pathways_simple$pathway, levels = c("sinking particles, t3",
                                                                  "sinking particles, t2", 
                                                                  "sinking particles, t1", 
                                                                  "zooplankton", "fish"))

dot_plot_200m_simple <- ggplot(data = flux_pathways_simple, aes(x = carbon_flux, y = pathway)) + 
  geom_point(aes(color = pathway, size = 2)) + theme_bw() + 
  theme(plot.margin = margin(10, 10, 5, 50), legend.title=element_blank(), legend.position = "none", 
        strip.background.x = element_blank(), 
        strip.text.x = element_text(size = 26), panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), text = element_text(size = 24)) +
  # make the 2 and one in superscript
  xlab(expression("Carbon transport to 200m (mg m"^-2*" d"^-1*")")) +
  ylab(expression("  \n  ")) + 
  scale_color_manual(values=c(darker_colors[2],"#1b9e77", lighter_colors[2], 
                              "#D95F02", "#7570B3")) + 
  geom_segment(aes(x = uncertainty_min, xend = uncertainty_max, y = pathway, yend = pathway, color = pathway), size = 1.5) +
  xlim(0, 250)

dot_plot_200m_simple

##### 500 m ######

# re-word pathways so they don't create returns in y axis labels
flux_pathways_500$pathway <- gsub("particles \n traps \n time 1", "sinking particles from sediment traps, t1", flux_pathways$pathway)
flux_pathways_500$pathway <- gsub("particles \n traps \n time 2", "sinking particles from sediment traps, t2", flux_pathways$pathway)
flux_pathways_500$pathway <- gsub("particles \n traps \n time 3", "sinking particles from sediment traps, t3", flux_pathways$pathway)
flux_pathways_500$pathway <- gsub("particles \n Th-234 \n time 1", "sinking particles from Th-234, t1", flux_pathways$pathway)
flux_pathways_500$pathway <- gsub("particles \n Th-234 \n time 2", "sinking particles from Th-234, t2", flux_pathways$pathway)
flux_pathways_500$pathway <- gsub("particles \n Th-234 \n time 3", "sinking particles from Th-234, t3", flux_pathways$pathway)

# reorder y axis labels so that time 1 comes at top before time 3
flux_pathways_500$pathway <- factor(flux_pathways_500$pathway, levels = c("sinking particles from sediment traps, t3",
                                                                  "sinking particles from sediment traps, t2", 
                                                                  "sinking particles from sediment traps, t1", 
                                                                  "sinking particles from Th-234, t3", 
                                                                  "sinking particles from Th-234, t2",
                                                                  "sinking particles from Th-234, t1",
                                                                  "zooplankton", "fish"))

# make a dot plot with 3 dots for each pathway
dot_plot_500m <- ggplot(data = flux_pathways_500, aes(x = carbon_flux, y = pathway)) + 
  geom_point(aes(color = pathway, size = 2)) + theme_bw() + 
  theme(plot.margin = margin(10, 10, 5, 50), legend.title=element_blank(), legend.position = "none", 
        strip.background.x = element_blank(), 
        strip.text.x = element_text(size = 26), panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), text = element_text(size = 24)) +
  xlab(bquote("Carbon flux 500m (mg m-2 d-1)")) + 
  ylab(expression("  \n  ")) + 
  scale_color_manual(values=c(darker_colors[2],"#1b9e77", lighter_colors[2], 
                              darker_colors[3], "#66A61E", lighter_colors[3], "#D95F02", "#7570B3")) + 
  annotate("text", x = 240, y = 8, label = "b)", size = 8) + 
  geom_point(aes(x = uncertainty_min, y = pathway, color = pathway, size = 2)) + 
  geom_point(aes(x = uncertainty_max, y = pathway, color = pathway, size = 2)) + 
  geom_segment(aes(x = uncertainty_min, xend = uncertainty_max, y = pathway, yend = pathway, color = pathway), size = 1.5) +
  xlim(0, 250)

# print plot  
dot_plot_500m

# put these two dot plot side by side into one figure
flux_comparison_panels_dot_plots <- ggarrange(dot_plot_200m, dot_plot_500m, 
                                    ncol=1, nrow=2, common.legend = TRUE, legend="none")
#print plot
flux_comparison_panels_dot_plots

# save
# ggsave("graphics/submission_graphics/flux_comparison.jpg", flux_comparison_panels_dot_plots, width = 13, height = 9, units = "in", dpi = 500)


## Percent fish C flux compared to total BCP #### 

# fish % contribution to total biological carbon transport (BCP), for results section
# (used sediment trap for passive flux estimation). Found min contribution by
# fish as min fish flux divided by the max flux for other pathways. 
# Found max contribution by fish as max fish flux divided by the min flux 
# for other pathways.

# 200 m min
max_BCP <- max(passive_flux_200_by_epoch$C_flux_mgC_m2_d) + max(zoop_flux_200$Total_C) + 
  fish_flux_min_200

fish_flux_min_200 / max_BCP * 100 # 0.52% 

# 200 m max
min_BCP <- min(passive_flux_200_by_epoch$C_flux_mgC_m2_d) + min(zoop_flux_200$Total_C) + 
  fish_flux_max_200

fish_flux_max_200 / min_BCP * 100 # 18% 

# repeat for deeper flux boundary

# 500 m min
max_BCP <- max(passive_flux_500_by_epoch$C_flux_mgC_m2_d) + max(zoop_flux_500$Total_C) + 
  fish_flux_min_500

fish_flux_min_500 / max_BCP * 100 # 0.43%

# 500 m max

min_BCP <- min(passive_flux_500_by_epoch$C_flux_mgC_m2_d) + min(zoop_flux_500$Total_C) + 
  fish_flux_max_500

fish_flux_max_500 / min_BCP * 100 # 13.4%


# factor difference between min and max estimate of fish C flux at 200 m
(fish_flux_max_200 - fish_flux_min_200)/fish_flux_min_200 # 11 times higher max than min 

####Comparison of mesopelagic fish C flux contribution compared to higher trophic levels ####


# Question during a meeting from Dr. Simon Thorrold: how does mesopelagic fish 
# flux compare to that of higher trophic levels? 

# Durfort et al 2022 finds that currently (post-exploitation), whales in Southern Ocean contribute 
# about 0.6 x 10^5 t C/yr. Pre-exploitation: 4 x 10^5 t C/yr. If Southern Ocean is about 20,000,000 km^2, 

# 0.6 x 10^5 t C/yr / 20,000,000 km^2 / 1,000,000 m^2/km^2 * 1 x 10^9 mg/metric tonne / 365 d/yr
.6e5 / 2e7 / 1e6 * 1e9 /365 # post exploitation
# = 0.008 mg C m^-2 d^-1

4e5 / 2e7 / 1e6 * 1e9 /365 # pre-exploitation of whales


# based on fish flux ranging from 1.6 to 21 (delete 2.9-38) mg C m^-2 d^-1 

# pre-exploitation, as a percentage: 
0.055/1.6 *100  
0.055/21 *100  


# post-exploitation, as a percentage: 
0.0082/1.6 *100 
0.0082/21 * 100 
# whales in Southern Ocean (PRE-exploitation, though now there is 15x less contribution according to Durfort et al 2022)
# may have contributed roughly 0.3-3% of fish flux in NE Atlantic based on our estimate. Post-exploitation, 
# whales in Southern Ocean may be contributing 0.04-0.5% of fish flux in NE Atlantic based on our estimate. 

