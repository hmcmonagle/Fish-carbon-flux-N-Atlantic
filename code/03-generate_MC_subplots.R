# check run time--see end of this script for more
start.time <- Sys.time()


# Run 01-calc_c_flux script (to load functions and assign object "nl" below)
source("code/01-calc_c_flux.R")

# uses get_pars() function to make a data frame containing nsims parameter
# draws, for bioenergetic and movement parameters only (generated in chapter 1) 
pars_df <- get_pars(nsims = 1000)

# now turn the current row of the data frame into a list.
# transpose, then rename each row (each par) as an item in the list becuase
# the function make_tau() requires a list.
# later turn row "1" into "i" in a loop
pars.draw.i.df <- t(pars_df[1,]) # later, make 1 into an i
pars.draw.i <- setNames(split(pars.draw.i.df, seq(nrow(pars.draw.i.df))), rownames(pars.draw.i.df))


#### Find C flux at study site for each nsims, VM fish ####
# Produce a df that gives the par draws (bioenergetic/movement pars, q, and catch)
# as well as associated C flux at the study site for each of nsims iterations
# (this df will be needed for plotting)

# for VM fish only. define nsims and make empty vectors for catch (C_t_bar_nsims)
# capture efficiency (q_nsims), calibration factor correction (cal_nsims), 
# carbon flux for each length increment (tau_l), and carbon flux (C_flux_site_nsims)
nsims = 1000
pars_df <- get_pars(nsims = nsims)
C_t_bar_nsims <- rep(NA, times = nsims)
q_nsims <- rep(NA, times = nsims)
cal_nsims <- rep(NA, times = nsims)
C_flux_site_nsims <- rep(NA, times = nsims)
tau_l <- rep(NA, times = nl) # nl has already been generated above


for (i in 1:nsims){
  # for Monte Carlo analysis: add new column to MC_results_VM with C flux at study site (step 11 below)
  
  # step 2 above: load data
  thedata <- load_data()
  fishdata <- thedata[[1]]
  volumedata <- thedata[[2]]
  
  # step 3 above: choose taxon
  taxagroup <- "Myctophidae"
  
  # step 4 above: assign parameters (does not include q and catch)
  pars.draw.i.df <- t(pars_df[i,])
  pars.draw.i <- setNames(split(pars.draw.i.df, seq(nrow(pars.draw.i.df))),
                          rownames(pars.draw.i.df))
  
  # step 5 above: filter data by taxon group
  taxadata <- dplyr::filter(thedata[[1]], Family == taxagroup)
  
  # step 6 above: calculate length bins and P(l)
  # P(l) is the probability of fish being in each length bin
  Pl_and_l <- make_Pl(taxadata$Std_length_mm)

  nl <- length(Pl_and_l$L)
  
  # step 7 above: For each L, calculate tau (total carbon flux for that length fish)
  # tau_l <- rep(NA, times = nl) <- removed this line and placed outside of loop
  for (j in 1:nl) {
    tau_l[j] <- make_tau(L = Pl_and_l$L[j], pars = pars.draw.i, taxagroup = taxagroup, VM = TRUE, lwreg)[1]
  }
  
  # step 8 above: Multiply P(L) times tau_l and sum
  ind_C_flux <- sum(tau_l * Pl_and_l$P)
  
  # step 9 above: find average of all tows, but instead of using tows from
  # the data with calc_C_t_bar() like we did for nominal par calculation, 
  # now use make_random_catch() to draw tows from the PDF for catches
  C_t <- make_random_catch(fishdata = fishdata, volumedata = volumedata, 
                           towlist = c(2, 4, 6, 8, 9, 11, 48, 90), group = "Myctophidae")
  C_t_bar <- mean(C_t) # average across all simulated catches
  C_t_bar_nsims[i] <- C_t_bar
  
  # step 10 above: calculate n_g for VM fish 
  # n_g = number of fish of that taxonomic group, based now on a simulated C_t_bar for myctophids
  # generate a random draw for capture efficiency for this single set of pars
  # and a random draw for calibration factor correction
  q <- runif(n = 1, min = 0.05, max = 0.50)
  cal <- runif(n = 1, min = 1, max = 1.5) # based on pers. comm. with Peter that calibration
  # factors can vary from roughly 4 to 6 (and calibration factors we use are closer to 4, 
  # so we'll allow this to vary up to 50% higher than calibration factors we use) 
  q_nsims[i] <- q
  cal_nsims[i] <- cal
  n_g_myct <- calc_n_g(catch = C_t_bar_nsims[i], capture_efficiency = q_nsims[i], cal = cal_nsims[i])
  
  # step 11 above: calculate C flux for VM fish and (new for MC analysis) add to MC_results_VM df
  C_flux_site_nsims[i] <- ind_C_flux * n_g_myct 
}

MC_results_VM <- cbind(pars_df, C_t_bar_nsims, q_nsims, cal_nsims, C_flux_site_nsims)
write.csv(MC_results_VM, "analysis/MC_results.csv", row.names = FALSE)


#### Repeat MC analysis for NM fish ####
# change make_tau argument "VM" from
# TRUE to FALSE and to change taxagroup from Myctophidae to Gonostomatidae

# changing argument VM from TRUE to false will change the H predictor 
# (habitat depth in modified Ikeda 2016 regression) from 1 to 0

# use same parameter values as used for VM fish where possible, 
# so that parameters that should be the same for any given simulation 
# (for example, temperature) are consistent for both VM and NM fish 
# when the two are added together to find total fish C flux for that simulation
nsims = 1000
# pars_df <- get_pars(nsims = nsims) Commented out because we do not want to 
# re-run pars_df, but rather use the same pars as used previously 
# for VM fish MC analysis where possible

C_t_bar_nsims_NM <- rep(NA, times = nsims)
q_nsims <- rep(NA, times = nsims)
cal_nsims <- rep(NA, times = nsims)
C_flux_site_nsims_NM <- rep(NA, times = nsims)
tau_l <- rep(NA, times = nl) # nl has already been generated above

for (i in 1:nsims){
  # for Monte Carlo analysis: add new column to MC_results_NM with C flux at study site (step 11 below)
  
  # step 2 above: load data
  thedata <- load_data()
  fishdata <- thedata[[1]]
  volumedata <- thedata[[2]]
  
  # step 3 above: choose taxon
  taxagroup <- "Gonostomatidae"
  
  # step 4 above: assign parameters (does not include q and catch)
  pars.draw.i.df <- t(pars_df[i,])
  pars.draw.i <- setNames(split(pars.draw.i.df, seq(nrow(pars.draw.i.df))),
                          rownames(pars.draw.i.df))
  
  # step 5 above: filter data by taxon group
  taxadata <- dplyr::filter(thedata[[1]], Family == taxagroup)
  
  # step 6 above: calculate length bins and P(l)
  # P(l) is the probability of fish being in each length bin
  Pl_and_l <- make_Pl(taxadata$Std_length_mm)
  
  nl <- length(Pl_and_l$L)
  
  # step 7 above: For each L, calculate tau (total carbon flux for that length fish)
  # tau_l <- rep(NA, times = nl) <- removed this line and placed outside of loop
  for (j in 1:nl) {
    tau_l[j] <- make_tau(L = Pl_and_l$L[j], pars = pars.draw.i, taxagroup = taxagroup, VM = FALSE, lwreg)[1]
  }
  
  # step 8 above: Multiply P(L) times tau_l and sum
  ind_C_flux <- sum(tau_l * Pl_and_l$P)
  
  # step 9 above: find average of all tows, but instead of using tows from
  # the data with calc_C_t_bar() like we did for nominal par calculation, 
  # now use make_random_catch() to draw tows from the PDF for catches
  C_t <- make_random_catch(fishdata = fishdata, volumedata = volumedata, 
                           towlist = c(2, 4, 6, 8, 9, 11, 48, 90), group = "Gonostomatidae")
  C_t_bar <- mean(C_t) # average across all simulated catches
  C_t_bar_nsims_NM[i] <- C_t_bar
  
  # step 10 above: calculate n_g for VM fish 
  # n_g = number of fish of that taxonomic group, based now on a simulated C_t_bar for myctophids
  # generate a random draw for capture efficiency for this single set of pars
  # allow capture efficiency to vary between 5% and 50%
  q <- runif(n = 1, min = 0.05, max = 0.50)
  q_nsims[i] <- q
  cal <- runif(n = 1, min = 1, max = 1.5)
  cal_nsims[i] <- cal
  n_g_myct <- calc_n_g(catch = C_t_bar_nsims_NM[i], capture_efficiency = q_nsims[i], cal = cal_nsims[i])
  
  # step 11 above: calculate C flux for NM fish and (new for MC analysis) add to MC_results_NM df
  C_flux_site_nsims_NM[i] <- ind_C_flux * n_g_myct 
}

MC_results_NM <- cbind(pars_df, C_t_bar_nsims_NM, q_nsims, cal_nsims, C_flux_site_nsims_NM)
write.csv(MC_results_NM, "analysis/MC_results_NM.csv", row.names = FALSE)

mean(MC_results_VM$C_flux_site_nsims) # 10.7. without taking into account calibration
# factor variation that might result in lower biomass: 13.0
mean(MC_results_NM$C_flux_site_nsims_NM) # 2.9. without taking into account
# calibration factor: 3.6

# save VM fish (vertical migrator) and NM fish (non-migrator)
# Monte Carlo results so that analysis doesn't need to keep being re-run
# (takes several minutes to run)

saveRDS(MC_results_VM, "analysis/VM_MC_results.RDS")
saveRDS(MC_results_NM, "analysis/NM_MC_results.RDS")



#### 14) plot VM fish MC results ####
# with standardized x axis (Z score) and find quantile GAM lines

# write get_upper_lower function to make length.out = 1000, k = 8, and shade inner 80% percentile later
get_lower_upper_1000 <- function(x, y, k = 8, err = 0.10, lowerbd, upperbd, medianbd) {
  dat <- tibble(x = x, y = y)
  predict.df = tibble(x = seq(min(x), max(x), length.out = 1000),
                      y = 0
  )
  fit.lower <- qgam(y~s(x, k=k, bs="ad"), data = dat, qu = lowerbd)
  pred.lower <- predict(fit.lower, newdata = predict.df)
  fit.upper <- qgam(y~s(x, k=k, bs="ad"), data = dat, qu = upperbd)
  pred.upper <- predict(fit.upper, newdata = predict.df)
  fit.median <- qgam(y~s(x, k=k, bs="ad"), data = dat, qu = medianbd)
  pred.median <- predict(fit.median, newdata = predict.df)
  
  return(cbind(x = predict.df$x, lower = pred.lower, upper = pred.upper, median = pred.median))
}

# re-write function that generates quantile GAM data frames that includes get_lower_upper_1000 function for nsims = 1000
# create generic code so we can repeat this later for all parameters
quantiles_df_func_1000 <- function(df_all_pars, single_par){
  one_par_df <- filter(df_all_pars, par_name== single_par)
  qgam_df <- as.data.frame(get_lower_upper_1000(x=one_par_df$par_value, y=one_par_df$C_flux_site_nsims,
                                                k=8, err = 0.10, lowerbd=0.10, upperbd=0.90, medianbd=0.50))
  quantile_df_1000 <- data.frame(one_par_df, median = qgam_df$median, lower = qgam_df$lower, upper = qgam_df$upper, x = qgam_df$x)
}

# plot with standardized x axis (Z score). exclude the final column, which is the C flux estimate
run_MC_VM_df_stand_1000 <- MC_results_VM
number_variables_VM <- length(run_MC_VM_df_stand_1000[1,])
run_MC_VM_df_stand_1000[,1:(number_variables_VM-1)] <- scale(MC_results_VM[,1:(number_variables_VM-1)])

# VM_NM habitat column is NaN because sd is zero and can't divide by zero. That's 
# fine before we're not going to plot this binary (0 or 1) categorical variable anyway. 
# convert wide to long format 
run_MC_VM_df_stand_long_1000 <- 
  run_MC_VM_df_stand_1000 %>%
  pivot_longer(cols = a0:cal_nsims, names_to = "par_name", values_to = "par_value")

# filter out parameters we do not want to plot
run_MC_VM_df_stand_long_1000 <- filter(run_MC_VM_df_stand_long_1000, 
                                       par_name != "H_VM" & par_name != "H_NM" 
                                       & par_name != "G_NM" & par_name != "C_c_nm"
                                       & par_name != "prop_non_detrital_NM_prey" )

# re-order factors for plotting in same order as parameters appear in equations in text
run_MC_VM_df_stand_long_ordered_1000 <- run_MC_VM_df_stand_long_1000

run_MC_VM_df_stand_long_ordered_1000$par_name <- factor(run_MC_VM_df_stand_long_ordered_1000$par_name, 
                                                        levels = c("C_flux_site_nsims", "a0", "a1", "a2", "a3", 
                                                                   "TT_epipelagic", "TT_mesopelagic", "TT_mean", 
                                                                   "Q_ox_kJ_g_O2", "activity_factor_foraging", "activity_factor_migrating", 
                                                                   "activity_factor_resting", "time_resting", "theta", "Dmax", "Dmin", "B", "S", 
                                                                   "RQ", "G_VM", "prop_SDA", "peakTime", "timeTo90", "prop_egestion", 
                                                                   "energy_in_pellet", "C_in_pellet", "proportion_egestion_below_150m", "prop_excretion", 
                                                                   "X_r", "X_e", "C_c_vm", "mu", "prop_mortality_below_150m", "q_nsims", "cal_nsims", "C_t_bar_nsims"))

# create parameter groups
resp_flux <- c("a0", "a1", "a2", "a3", "TT_epipelagic", "TT_mesopelagic", "TT_mean", "Q_ox_kJ_g_O2", "activity_factor_foraging", "activity_factor_migrating", "activity_factor_resting", "time_resting", "theta", "Dmax", "Dmin", "B", "S", "RQ", "G_VM")
SDA_flux <- c("prop_SDA", "peakTime", "timeTo90")
eg_flux <- c("prop_egestion", "energy_in_pellet", "C_in_pellet", "proportion_egestion_below_150m")
ex_flux <- c("prop_excretion", "X_r", "X_e")
mort_flux <- c("C_c_vm", "mu", "prop_mortality_below_150m")
biomass <- c("q_nsims", "cal_nsims", "C_t_bar_nsims")

# add column for color coding, by parameter group
MC_VM_data_1000 <- run_MC_VM_df_stand_long_ordered_1000 %>%
  mutate(group = factor(case_when(par_name %in% resp_flux ~ "respiratory flux and growth", 
                                  par_name %in% SDA_flux ~ "specific dynamic action flux", 
                                  par_name %in% eg_flux ~ "egestion flux", 
                                  par_name %in% ex_flux ~ "excretion flux", 
                                  par_name %in% mort_flux ~ "mortality flux", 
                                  par_name %in% biomass ~ "biomass", TRUE ~ NA_character_)))

# re-order for plotting
MC_VM_data_1000$group <- factor(MC_VM_data_1000$group, levels = c("respiratory flux and growth", "specific dynamic action flux", "egestion flux", "excretion flux", "mortality flux", "biomass"))

A_f_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "activity_factor_foraging")
A_m_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "activity_factor_migrating")
A_r_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "activity_factor_resting")
a0_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "a0")
a1_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "a1")
a2_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "a2")
a3_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "a3")
T_epi_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "TT_epipelagic")
T_meso_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "TT_mesopelagic")
T_mean_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "TT_mean")
t_r_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "time_resting")
Dmax_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "Dmax")
Dmin_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "Dmin")
B_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "B")
S_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "S")
prop_SDA_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "prop_SDA")
prop_F_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "prop_egestion")
prop_X_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "prop_excretion")
G_VM_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "G_VM")
theta_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "theta")
t_peak_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "peakTime")
t_90_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "timeTo90")
X_r_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "X_r")
X_e_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "X_e")
mu_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "mu")
Q_ox_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "Q_ox_kJ_g_O2")
RQ_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "RQ")
E_p_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "energy_in_pellet")
C_p_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "C_in_pellet")
F_e_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "proportion_egestion_below_150m")
C_c_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "C_c_vm")
M_e_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "prop_mortality_below_150m")
q_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "q_nsims")
cal_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "cal_nsims")
catch_df <- quantiles_df_func_1000(df_all_pars = MC_VM_data_1000, single_par = "C_t_bar_nsims")


# rbind into one data frame
all_pars_df_VM_1000 <- rbind(A_f_df, A_m_df, A_r_df, a0_df, a1_df, a2_df, a3_df, T_epi_df, T_meso_df, T_mean_df, t_r_df, Dmax_df, Dmin_df, B_df, S_df, prop_SDA_df, prop_F_df, prop_X_df, G_VM_df, theta_df, t_peak_df, t_90_df, X_r_df, X_e_df, mu_df, Q_ox_df, RQ_df, E_p_df, C_p_df, F_e_df, C_c_df, M_e_df, q_df, cal_df, catch_df)

# make labels for multiplot
my_labeller <- as_labeller(c(activity_factor_foraging="A[f]", activity_factor_migrating="A[m]", activity_factor_resting="A[r]", a0="alpha[0]", a1="alpha[1]", a2="alpha[2]", a3="alpha[3]", B="B", C_c_vm="C[c]", C_in_pellet="C[p]", Dmax="D[max]", Dmin="D[min]", energy_in_pellet="e[p]", proportion_egestion_below_150m="E[F]", G_VM="G[VM]", prop_mortality_below_150m="E[mu]", mu="mu", prop_egestion="phi[F]", prop_SDA="phi[SDA]", prop_excretion="phi[X]", Q_ox_kJ_g_O2="Q[ox]", RQ="R[Q]", S="S", timeTo90="t[90]", TT_epipelagic="T[epi]", TT_mean="T[mean]", TT_mesopelagic="T[meso]", peakTime="t[peak]", time_resting="t[r]", theta="theta", X_e="E[X]", X_r="X[r]", q_nsims = "q", cal_nsims = "cal", C_t_bar_nsims = "catch"),
                           default = label_parsed)

# reduce tick marks by setting them to a given min and max C flux on y axis, and by setting Z score min and max on x axis to -2 and +2. See scale_y_continuous and scale_x_continuous arguments below

# create multiplot
plot_MC_VM_qgam_1000 <- ggplot(data = all_pars_df_VM_1000, 
                               aes(x = par_value, y = C_flux_site_nsims)) + 
  scale_x_continuous(breaks = c(-2, 0, 2), limits = c(-2, 2)) + 
  scale_y_continuous(breaks = c(0, 30), limits = c(0, 30)) + 
  geom_point(alpha = 0.02, size = .5) + 
  geom_line(data = all_pars_df_VM_1000, aes(x = x, y = lower), alpha = 0.5) + 
  geom_line(data = all_pars_df_VM_1000, aes(x = x, y = upper), alpha = 0.5) + 
  geom_line(data = all_pars_df_VM_1000, aes(x = x, y = median), alpha = 0.5) + 
  geom_ribbon(data = all_pars_df_VM_1000, aes(ymin=lower, ymax=upper, x=x, 
                                              color = group, group = group,  fill = group), 
              show.legend=TRUE, alpha = 0.2, lwd=0.5) + theme_bw() + 
  theme(legend.title=element_blank(), legend.position = "top", strip.background.x = element_blank(), 
        strip.text.x = element_text(size = 24), panel.spacing = unit(1, 'lines'), 
        panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        text = element_text(size = 24)) + 
  facet_wrap(vars(par_name), ncol = 8, labeller = my_labeller) + 
  xlab(bquote("\n Variable values for vertically migrating fish (z-score)")) + 
  ylab(expression("Carbon flux (mg carbon per square meter per day)")) +
  scale_color_manual(values = c("#66A61E", "#7570B3", "#D95F02", "#E6AB02", "#E7298A", "#666666")) +
  scale_fill_manual(values = c("#66A61E", "#7570B3", "#D95F02", "#E6AB02", "#E7298A", "#666666")) +
  guides(color = guide_legend(nrow=2,byrow=TRUE))

# plot
plot_MC_VM_qgam_1000

# export plot to resolution and size needed using. used preview to save as 170 mm wide for journal requirements
ggsave("graphics/submission_graphics/MC_VM_plot.jpg", plot = plot_MC_VM_qgam_1000, device = "jpeg", width = 170*2, height = 130*2, units = c("mm"), dpi = 500)



# Next, make similar plot for NM fish

#### 15) plot NM fish MC results ####
# with standardized x axis (Z score) and find quantile GAM lines

# use same functions get_upper_lower and quantiles_df_func_1000 as used for VM fish--
# just edit quantiles_df_func_1000 use use NM fish results for C flux

# re-write function that generates quantile GAM data frames that includes get_lower_upper_1000 function for nsims = 1000
# create generic code so we can repeat this later for all parameters
quantiles_df_func_1000_NM <- function(df_all_pars, single_par){
  one_par_df <- filter(df_all_pars, par_name== single_par)
  qgam_df <- as.data.frame(get_lower_upper_1000(x=one_par_df$par_value, y=one_par_df$C_flux_site_nsims_NM,
                                                k=8, err = 0.10, lowerbd=0.10, upperbd=0.90, medianbd=0.50))
  quantile_df_1000 <- data.frame(one_par_df, median = qgam_df$median, lower = qgam_df$lower, upper = qgam_df$upper, x = qgam_df$x)
}


# plot with standardized x axis (Z score). exclude the final column, which is the C flux estimate
run_MC_NM_df_stand_1000 <- MC_results_NM
number_variables_NM <- length(run_MC_NM_df_stand_1000[1,])
run_MC_NM_df_stand_1000[,1:(number_variables_NM-1)] <- scale(MC_results_NM[,1:(number_variables_NM-1)])

# VM_NM habitat column is NaN because sd is zero and can't divide by zero. That's 
# fine before we're not going to plot this binary (0 or 1) categorical variable anyway. 
# convert wide to long format 
run_MC_NM_df_stand_long_1000 <- 
  run_MC_NM_df_stand_1000 %>%
  pivot_longer(cols = a0:cal_nsims, names_to = "par_name", values_to = "par_value")

# filter out parameters we do not want to plot
run_MC_NM_df_stand_long_1000 <- filter(run_MC_NM_df_stand_long_1000, 
                                       par_name != "H_VM" & par_name != "H_NM" 
                                       & par_name != "G_VM" & par_name != "C_c_vm" 
                                       & par_name != "TT_mean" & par_name != "peakTime"
                                       & par_name != "TT_epipelagic" & par_name != "activity_factor_migrating" 
                                       & par_name != "Dmax" & par_name != "Dmin"
                                       & par_name != "B" & par_name != "S"
                                       & par_name != "timeTo90" & par_name != "prop_mortality_below_150m"
                                       & par_name != "X_e" & par_name != "proportion_egestion_below_150m")

# ^ note: should have removed T_mean from NM fish MC results in chapter 1, McMonagle et al 2024 (Fig. 6)

# re-order factors for plotting in same order as parameters appear in equations in text
run_MC_NM_df_stand_long_ordered_1000 <- run_MC_NM_df_stand_long_1000

run_MC_NM_df_stand_long_ordered_1000$par_name <- factor(run_MC_NM_df_stand_long_ordered_1000$par_name, 
                                                        levels = c("C_flux_site_nsims_NM", "a0", "a1", "a2", "a3", 
                                                                   "TT_mesopelagic", "Q_ox_kJ_g_O2", "activity_factor_foraging", 
                                                                   "activity_factor_resting", "time_resting", "theta", "RQ", "G_NM", 
                                                                   "prop_SDA", "prop_egestion", "energy_in_pellet", "C_in_pellet", 
                                                                   "prop_excretion", "X_r", "C_c_nm", "mu", "prop_non_detrital_NM_prey", 
                                                                   "q_nsims", "cal_nsims", "C_t_bar_nsims_NM"))

# create parameter groups for color coding in plot later 
resp_flux <- c("a0", "a1", "a2", "a3", "TT_mesopelagic", "Q_ox_kJ_g_O2", "activity_factor_foraging", "activity_factor_resting", "time_resting", "theta", "RQ", "G_NM")
SDA_flux <- c("prop_SDA")
eg_flux <- c("prop_egestion", "energy_in_pellet", "C_in_pellet")
ex_flux <- c("prop_excretion", "X_r")
mort_flux <- c("C_c_nm", "mu")
all_fluxes <- c("prop_non_detrital_NM_prey")
biomass <- c("q_nsims", "cal_nsims", "C_t_bar_nsims_NM")

# add column for color coding, by parameter group
MC_NM_data_1000 <- run_MC_NM_df_stand_long_ordered_1000 %>%
  mutate(group = factor(case_when(par_name %in% resp_flux ~ "respiratory flux and growth", 
                                  par_name %in% SDA_flux ~ "specific dynamic action flux", 
                                  par_name %in% eg_flux ~ "egestion flux", 
                                  par_name %in% ex_flux ~ "excretion flux", 
                                  par_name %in% mort_flux ~ "mortality flux", 
                                  par_name %in% all_fluxes ~ "all fluxes",
                                  par_name %in% biomass ~ "biomass")))

# re-order for plotting
MC_NM_data_1000$group <- factor(MC_NM_data_1000$group, levels = c("respiratory flux and growth", "specific dynamic action flux", "egestion flux", "excretion flux", "mortality flux", "all fluxes", "biomass"))

A_f_df <- quantiles_df_func_1000_NM(df_all_pars = MC_NM_data_1000, single_par = "activity_factor_foraging")
A_r_df <- quantiles_df_func_1000_NM(df_all_pars = MC_NM_data_1000, single_par = "activity_factor_resting")
a0_df <- quantiles_df_func_1000_NM(df_all_pars = MC_NM_data_1000, single_par = "a0")
a1_df <- quantiles_df_func_1000_NM(df_all_pars = MC_NM_data_1000, single_par = "a1")
a2_df <- quantiles_df_func_1000_NM(df_all_pars = MC_NM_data_1000, single_par = "a2")
a3_df <- quantiles_df_func_1000_NM(df_all_pars = MC_NM_data_1000, single_par = "a3")
T_meso_df <- quantiles_df_func_1000_NM(df_all_pars = MC_NM_data_1000, single_par = "TT_mesopelagic")
t_r_df <- quantiles_df_func_1000_NM(df_all_pars = MC_NM_data_1000, single_par = "time_resting")
prop_SDA_df <- quantiles_df_func_1000_NM(df_all_pars = MC_NM_data_1000, single_par = "prop_SDA")
prop_F_df <- quantiles_df_func_1000_NM(df_all_pars = MC_NM_data_1000, single_par = "prop_egestion")
prop_X_df <- quantiles_df_func_1000_NM(df_all_pars = MC_NM_data_1000, single_par = "prop_excretion")
G_NM_df <- quantiles_df_func_1000_NM(df_all_pars = MC_NM_data_1000, single_par = "G_NM")
theta_df <- quantiles_df_func_1000_NM(df_all_pars = MC_NM_data_1000, single_par = "theta")
X_r_df <- quantiles_df_func_1000_NM(df_all_pars = MC_NM_data_1000, single_par = "X_r")
mu_df <- quantiles_df_func_1000_NM(df_all_pars = MC_NM_data_1000, single_par = "mu")
Q_ox_df <- quantiles_df_func_1000_NM(df_all_pars = MC_NM_data_1000, single_par = "Q_ox_kJ_g_O2")
RQ_df <- quantiles_df_func_1000_NM(df_all_pars = MC_NM_data_1000, single_par = "RQ")
E_p_df <- quantiles_df_func_1000_NM(df_all_pars = MC_NM_data_1000, single_par = "energy_in_pellet")
C_p_df <- quantiles_df_func_1000_NM(df_all_pars = MC_NM_data_1000, single_par = "C_in_pellet")
C_c_df <- quantiles_df_func_1000_NM(df_all_pars = MC_NM_data_1000, single_par = "C_c_nm")
Z_df <- quantiles_df_func_1000_NM(df_all_pars = MC_NM_data_1000, single_par = "prop_non_detrital_NM_prey")
q_df <- quantiles_df_func_1000_NM(df_all_pars = MC_NM_data_1000, single_par = "q_nsims")
cal_df <- quantiles_df_func_1000_NM(df_all_pars = MC_NM_data_1000, single_par = "cal_nsims")
catch_df <- quantiles_df_func_1000_NM(df_all_pars = MC_NM_data_1000, single_par = "C_t_bar_nsims_NM")

# rbind into one data frame
all_pars_df_NM_1000 <- rbind(A_f_df, A_r_df, a0_df, a1_df, a2_df, a3_df, T_meso_df, t_r_df, prop_SDA_df, prop_F_df, prop_X_df, G_NM_df, theta_df, X_r_df, mu_df, Q_ox_df, RQ_df, E_p_df, C_p_df, C_c_df, Z_df, q_df, cal_df, catch_df)

# make labels for multiplot
my_labeller <- as_labeller(c(activity_factor_foraging="A[f]", activity_factor_resting="A[r]", a0="alpha[0]", a1="alpha[1]", a2="alpha[2]", a3="alpha[3]", C_c_nm="C[c]", 
                             C_in_pellet="C[p]", energy_in_pellet="e[p]", G_NM="G[NM]", mu="mu", prop_egestion="phi[F]", prop_SDA="phi[SDA]", prop_excretion="phi[X]", Q_ox_kJ_g_O2="Q[ox]", RQ="R[Q]", 
                             TT_mesopelagic="T[meso]", time_resting="t[r]", theta="theta", X_r="X[r]", prop_non_detrital_NM_prey = "Z", q_nsims = "q", cal_nsims = "cal", C_t_bar_nsims_NM = "catch"),
                           default = label_parsed)

# reduce tick marks by setting them to a given min and max C flux on y axis, and by setting Z score min and max on x axis to -2 and +2. See scale_y_continuous and scale_x_continuous arguments below

# create multiplot
plot_MC_NM_qgam_1000 <- ggplot(data = all_pars_df_NM_1000, 
                               aes(x = par_value, y = C_flux_site_nsims_NM)) + 
  scale_x_continuous(breaks = c(-2, 0, 2), limits = c(-2, 2)) + 
  scale_y_continuous(breaks = c(0, 10), limits = c(0, 10)) + 
  geom_point(alpha = 0.02, size = .5) + 
  geom_line(data = all_pars_df_NM_1000, aes(x = x, y = lower), alpha = 0.5) + 
  geom_line(data = all_pars_df_NM_1000, aes(x = x, y = upper), alpha = 0.5) + 
  geom_line(data = all_pars_df_NM_1000, aes(x = x, y = median), alpha = 0.5) + 
  geom_ribbon(data = all_pars_df_NM_1000, aes(ymin=lower, ymax=upper, x=x, 
                                              color = group, group = group,  fill = group), 
              show.legend=TRUE, alpha = 0.2, lwd=0.5) + theme_bw() + 
  theme(legend.title=element_blank(), legend.position = "top", strip.background.x = element_blank(), 
        strip.text.x = element_text(size = 24), panel.spacing = unit(1, 'lines'), 
        panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        text = element_text(size = 24)) + 
  facet_wrap(vars(par_name), ncol = 8, labeller = my_labeller) + 
  xlab(bquote("\n Variable values for non-migrating fish (z-score)")) + 
  ylab(expression("Carbon flux (mg carbon per square meter per day)")) +
  scale_color_manual(values = c("#66A61E", "#7570B3", "#D95F02", "#E6AB02", "#E7298A", "#1b9e77", "#666666")) +
  scale_fill_manual(values = c("#66A61E", "#7570B3", "#D95F02", "#E6AB02", "#E7298A", "#1b9e77", "#666666")) +
  guides(color = guide_legend(nrow=3,byrow=TRUE))

# plot
plot_MC_NM_qgam_1000

# export plot to resolution and size needed using. used preview to save as 170 mm wide for journal requirements
ggsave("graphics/submission_graphics/MC_NM_plot.jpg", plot = plot_MC_NM_qgam_1000, device = "jpeg", width = 170*2, height = 112*2, units = c("mm"), dpi = 500)



#### next, run Monte Carlo for VM fish but with B = 500 m ####
# instead of 200 m (B = flux boundary)

# for VM fish only. define nsims and make empty vectors for catch (C_t_bar_nsims)
# capture efficiency (q_nsims), calibration factor correction (cal_nsims), 
# carbon flux for each length increment (tau_l), and carbon flux (C_flux_site_nsims)
nsims = 1000
pars_df_500 <- get_pars(nsims = nsims)
pars_df_500$B <- runif(nsims, min = 499, max = 501) # B = 500 m
C_t_bar_nsims <- rep(NA, times = nsims)
q_nsims <- rep(NA, times = nsims)
cal_nsims <- rep(NA, times = nsims)
C_flux_site_nsims <- rep(NA, times = nsims)
tau_l <- rep(NA, times = nl) # nl has already been generated above

for (i in 1:nsims){
  # for Monte Carlo analysis: add new column to MC_results_VM with C flux at study site (step 11 below)
  
  # step 2 above: load data
  thedata <- load_data()
  fishdata <- thedata[[1]]
  volumedata <- thedata[[2]]
  
  # step 3 above: choose taxon
  taxagroup <- "Myctophidae"
  
  # step 4 above: assign parameters (does not include q and catch)
  pars.draw.i.df <- t(pars_df[i,])
  pars.draw.i <- setNames(split(pars.draw.i.df, seq(nrow(pars.draw.i.df))),
                          rownames(pars.draw.i.df))
  
  # step 5 above: filter data by taxon group
  taxadata <- dplyr::filter(thedata[[1]], Family == taxagroup)
  
  # step 6 above: calculate length bins and P(l)
  # P(l) is the probability of fish being in each length bin
  Pl_and_l <- make_Pl(taxadata$Std_length_mm)
  
  nl <- length(Pl_and_l$L)
  
  # step 7 above: For each L, calculate tau (total carbon flux for that length fish)
  # Note: for a 500 m flux boundary, based on vertical distribution results plotted 
  # day versus night, we can see that the dominant migrating fish (Myctophidae and 
  # Sternoptychidae) generally do not migrate past 500 m. Therefore, only their 
  # egestion and mortality fluxes could be contributing to carbon flux past this 
  # 500 m flux boundary. For this reason, instead of using the first value [1] 
  # generated by make_tau(), we will use only the sum of values [4] and [6], 
  # which are the egestion and mortality fluxes, respectively. See the last few lines
  # of C.flux.function.VM code for more details on how this function works. 
  for (j in 1:nl) {
    tau_l[j] <- make_tau(L = Pl_and_l$L[j], pars = pars.draw.i, taxagroup = taxagroup, VM = TRUE, lwreg)[4] + 
      make_tau(L = Pl_and_l$L[j], pars = pars.draw.i, taxagroup = taxagroup, VM = TRUE, lwreg)[6]
  }
  
  # step 8 above: Multiply P(L) times tau_l and sum
  ind_C_flux <- sum(tau_l * Pl_and_l$P)
  
  # step 9 above: find average of all tows, but instead of using tows from
  # the data with calc_C_t_bar() like we did for nominal par calculation, 
  # now use make_random_catch() to draw tows from the PDF for catches
  C_t <- make_random_catch(fishdata = fishdata, volumedata = volumedata, 
                           towlist = c(2, 4, 6, 8, 9, 11, 48, 90), group = "Myctophidae")
  C_t_bar <- mean(C_t) # average across all simulated catches
  C_t_bar_nsims[i] <- C_t_bar
  
  # step 10 above: calculate n_g for VM fish 
  # n_g = number of fish of that taxonomic group, based now on a simulated C_t_bar for myctophids
  # generate a random draw for capture efficiency for this single set of pars
  q <- runif(n = 1, min = 0.05, max = 0.50)
  q_nsims[i] <- q
  cal <- runif(n = 1, min = 1, max = 1.5)
  cal_nsims[i] <- cal
  n_g_myct <- calc_n_g(catch = C_t_bar_nsims[i], capture_efficiency = q_nsims[i], cal= cal_nsims[i])
  
  # step 11 above: calculate C flux for VM fish and (new for MC analysis) add to MC_results_VM df
  C_flux_site_nsims[i] <- ind_C_flux * n_g_myct 
}

MC_results_VM_500 <- cbind(pars_df_500, C_t_bar_nsims, q_nsims, cal_nsims, C_flux_site_nsims)
write.csv(MC_results_VM_500, "analysis/MC_results_VM_500.csv", row.names = FALSE)


#### Repeat MC analysis for NM fish at 500 m ####

# change make_tau argument "VM" from
# TRUE to FALSE and to change taxagroup from Myctohpidae to Gonostomatidae

# otherwise, use same parameter values as used for VM fish where possible, 
# so that parameters that should be the same for any given simulation 
# (for example, temperature) are consistent for both VM and NM fish 
# when the two are added together to find total fish C flux for that simulation

# use 500 m flux boundary set of pars again for NM fish, as above for VM fish
# Unlike for VM fish, we do not adjust which carbon flux pathways contribute 
# to flux past 500 m because most NM fish (gonostomatidae) occur > 500 m based
# on cruise catch data (see vertical distribution plot from the R/V Cook)
nsims = 1000
# pars_df <- get_pars(nsims = nsims) Commented out because we do not want to
# re-run pars_df

# re-run pars_df, but rather use the same pars as used previously for VM fish MC analysis
C_t_bar_nsims_NM <- rep(NA, times = nsims)
q_nsims <- rep(NA, times = nsims)
cal_nsims <- rep(NA, times = nsims)
C_flux_site_nsims_NM <- rep(NA, times = nsims)
tau_l <- rep(NA, times = nl) # nl has already been generated above

for (i in 1:nsims){
  # for Monte Carlo analysis: add new column to MC_results_NM with C flux at study site (step 11 below)
  
  # step 2 above: load data
  thedata <- load_data()
  fishdata <- thedata[[1]]
  volumedata <- thedata[[2]]
  
  # step 3 above: choose taxon
  taxagroup <- "Gonostomatidae"
  
  # step 4 above: assign parameters (does not include q and catch)
  pars.draw.i.df <- t(pars_df_500[i,])
  pars.draw.i <- setNames(split(pars.draw.i.df, seq(nrow(pars.draw.i.df))),
                          rownames(pars.draw.i.df))
  
  # step 5 above: filter data by taxon group
  taxadata <- dplyr::filter(thedata[[1]], Family == taxagroup)
  
  # step 6 above: calculate length bins and P(l)
  # P(l) is the probability of fish being in each length bin
  Pl_and_l <- make_Pl(taxadata$Std_length_mm)

  nl <- length(Pl_and_l$L)
  
  # step 7 above: For each L, calculate tau (total carbon flux for that length fish)
  # tau_l <- rep(NA, times = nl) <- removed this line and placed outside of loop
  for (j in 1:nl) {
    tau_l[j] <- make_tau(L = Pl_and_l$L[j], pars = pars.draw.i, taxagroup = taxagroup, VM = FALSE, lwreg)[1]
  }
  
  # step 8 above: Multiply P(L) times tau_l and sum
  ind_C_flux <- sum(tau_l * Pl_and_l$P)
  
  # step 9 above: find average of all tows, but instead of using tows from
  # the data with calc_C_t_bar() like we did for nominal par calculation, 
  # now use make_random_catch() to draw tows from the PDF for catches
  C_t <- make_random_catch(fishdata = fishdata, volumedata = volumedata, 
                           towlist = c(2, 4, 6, 8, 9, 11, 48, 90), group = "Gonostomatidae")
  C_t_bar <- mean(C_t) # average across all simulated catches
  C_t_bar_nsims_NM[i] <- C_t_bar
  
  # step 10 above: calculate n_g for VM fish 
  # n_g = number of fish of that taxonomic group, based now on a simulated C_t_bar for myctophids
  # generate a random draw for capture efficiency for this single set of pars
  # allow capture efficiency to vary between 5% and 50%
  q <- runif(n = 1, min = 0.05, max = 0.50)
  cal <- runif(n = 1, min = 1, max = 1.5)
  q_nsims[i] <- q
  cal_nsims[i] <- cal
  n_g_myct <- calc_n_g(catch = C_t_bar_nsims_NM[i], capture_efficiency = q_nsims[i], cal = cal_nsims[i])
  
  # step 11 above: calculate C flux for NM fish and (new for MC analysis) add to MC_results_NM df
  C_flux_site_nsims_NM[i] <- ind_C_flux * n_g_myct 
}

MC_results_NM_500 <- cbind(pars_df_500, C_t_bar_nsims, q_nsims, cal_nsims, C_flux_site_nsims)
write.csv(MC_results_NM_500, "analysis/MC_results_NM_500.csv", row.names = FALSE)

# save to RDS files so MC doesn't need to be re-run 
saveRDS(MC_results_VM_500, "analysis/VM_MC_results_500.RDS")
saveRDS(MC_results_NM_500, "analysis/NM_MC_results_500.RDS")


# check run time (given in min below)
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken # in minutes (about 21 minutes as of 2024-01-09)


# ** Note: the MC results are used in a later script, 09-sequestration-100-years.R
