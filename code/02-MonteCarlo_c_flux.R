
# Set up ####
# 1) Load functions and fishdata filtered for 3 most abundant families
source("code/01-calc_c_flux.R")

# 2) Load all data and specify simulation parameters
thedata <- load_data()

fishdata <- fishdata_3_families # includes Myctohpidae, Gonostomatidae, and 
# Sternoptychidae, which make up 91% of the catch and have reliable length-weight
# regressions available from cruise data to fill in missing weights

volumedata <- thedata[[2]]

# load L- W regression and convert L to W
lwreg <- readRDS("analysis/length_weight_parameters.RDS")
nominal.pars <- make_nominal_pars()
nsims <- 100L # L makes this value an integer rather than numeric
qpars <- c(0.05, 0.5) # min and max of q (capture efficiency)
calpars <- c(1, 1.5) # min and max of cal (calibration factor correction)

# add both Sarmiento and Cook fish to towlist (C_t_bar is mean catch for each taxon t)
C_t_bar_obs <- calc_C_t_bar(fishdata = fishdata, volumedata = volumedata, towlist = c(2, 4, 6, 8, 9, 11, 48, 90))$c_t_bar

groupnames <- c("Myctophidae", "Gonostomatidae", "Sternoptychidae") # added Sternoptychidae as VM fish 
migrate <- c("VM", "NM", "VM") # added "VM" as the third value in this vector for Sternoptychidae
par_sims <- get_pars(nsims)

# Scenario 1 ####
simpars = T
simC = T
total_C_flux <- sim_C_flux(simpars, simC, nsims, groupnames, migrate, thedata, 
                           fishdata, volumedata, par_sims,nominal.pars,C_t_bar_obs, qpars, calpars)

# put results in tibble
model_sims_both <- tibble(Scenario = "Both",
                     Cflux = total_C_flux[[1]] + total_C_flux[[2]])

# repeat, holding catch as observed levels
# Scenario 2 ####
simpars = T
simC = F
total_C_flux <- sim_C_flux(simpars, simC, nsims, groupnames, migrate, thedata, fishdata, volumedata, par_sims,nominal.pars,C_t_bar_obs, qpars, calpars) 

# save results to a tibble
model_sims_indpars <- tibble(Scenario = "Bioenergetics",
                          Cflux = total_C_flux[[1]] + total_C_flux[[2]])

# Scenario 3 ####
# repeat, holding individual parameters constant
simpars = F
simC = T
total_C_flux <- sim_C_flux(simpars, simC, nsims, groupnames, migrate, thedata, fishdata, volumedata, par_sims,nominal.pars,C_t_bar_obs, qpars, calpars) 
model_sims_catch <-  tibble(Scenario = "Biomass",
                            Cflux = total_C_flux[[1]] + total_C_flux[[2]])

model_sims_all <- rbind(model_sims_both,
                        model_sims_indpars,
                        model_sims_catch)

# Plot Results ####

theme_set(theme_bw())
theme_update(panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank(),
             strip.background = element_blank())

# add a margin at top of plot
violin_plot <- ggplot(data = model_sims_all, aes(x = Scenario, y = Cflux)) +
  geom_violin(adjust = 1.5, fill = "#7570B3", colour = "black") +
  ylab(expression(Carbon~flux~(mgC~m^{-2}~d^{-1}))) +
  xlab("\nSource of Uncertainty") +
  theme(legend.position = "none", text = element_text(size = 40), plot.margin = unit(c(1, 0.5, 0.5, 0.5), 
                                                                                     "inches"))

# plot with bigger margin at top of figure so that top label isn't cut off
ggsave("graphics/submission_graphics/violin_plot.jpg", violin_plot, width = 11, height = 9, units = "in", dpi = 500)

# calculate standard deviation of c flux estimates across 3 scenarios
model_sims_all %>% 
  group_by(Scenario) %>% 
  summarise(sd = sd(Cflux))
