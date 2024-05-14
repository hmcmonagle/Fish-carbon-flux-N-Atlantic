# First R script to load for North Atlantic EXPORTS manuscript on fish carbon flux
# Co-authors: Helena McMonagle, Joel Llopiz, Amy Maas, Deborah Steinberg, 
# Annette Govindarajan, Ray Hilborn and Timothy Essington

# 1) Load functions
source("code/00-functions.R")

# 2) Load all data
thedata <- load_data()
fishdata <- thedata[[1]] 
volumedata <- thedata[[2]]

# 3) Calculate missing weights

# Make new column for use for individuals without an empirical weight measurement
fishdata$estimated_weight <- NA

# fill in missing total_weight values using the length weight regression below
lwreg <- readRDS("analysis/length_weight_parameters.RDS")

families <- c("Myctophidae", "Gonostomatidae", "Sternoptychidae")

for (i in 1:length(families)) {
  
  # find rows with fish from each family
  family_index <- which(fishdata$Family == families[i])
  
  # loop this through all (abundant) families: Myctophidae, Gonostomatidae and Sternoptychidae
  loga <- lwreg$loga[lwreg$Group == families[i]]
  b <- lwreg$b[lwreg$Group == families[i]]
  fishdata$estimated_weight[family_index] <- exp(loga) * fishdata$Std_length_mm[family_index] ^ b
  
}

# make new column in df that has the "weight.to.use", which is total_weight if empirical measurement exists,
# or estimated_weight if total_weight is missing
fishdata$weight.to.use <- fishdata$total_weight # n = 1553
fishdata$weight.to.use[is.na(fishdata$weight.to.use)] <- fishdata$estimated_weight[is.na(fishdata$weight.to.use)]

# edit fishdata to just include families for which we can calculate an empirical
# length-weight relationship from the individuals for which we have both length and weight.
# These make up 92% of the total fish count (exclude rare or unidentified taxa in the analysis)
fishdata <- fishdata[!is.na(fishdata$weight.to.use),] # n = 1,427

# for rows that have both an estimated_weight and a total_weight, 
# plot scatter plot of estimated_weight against total_weight
# to ensure length-weight regression output is accurate

weight_plot <- ggplot(fishdata, aes(x = estimated_weight, y = total_weight)) +
  geom_point() +
  facet_wrap(~Family) +
  geom_abline(intercept = 0, slope = 1, color = "blue") +
  labs(x = "Estimated weight (mg)", y = "Total weight (mg)") +
  ylim(0, 10)

weight_plot


# edit fishdata to just include Family==Myctophidae, Sternoptychidae, and Gonostomatidae, 
# which make up vast majority of fish catch. n = 1357
fishdata_3_families <- fishdata[fishdata$Family %in% c("Myctophidae", "Sternoptychidae", "Gonostomatidae"),]
length(fishdata_3_families$Ship)

# number of fish from Sarmiento
length(fishdata_3_families$Ship[fishdata_3_families$Ship == "Sarmiento"])
# number of fish from Cook
length(fishdata_3_families$Ship[fishdata_3_families$Ship == "Cook"])

# this is the fishdata version that will be used in the carbon flux analysis, so
# use this version of the catch data in subsequent analyses

# do Kolmogorov-Smirnov Test to determine whether fishdata and fishdata_3_families are 
# significantly different in terms of length distribution (they are not, so conclude OK to use fishdata_3_families)
# fishdata_3_families is a subset of fishdata
ks.test(x = fishdata$Std_length_mm, y = fishdata_3_families$Std_length_mm)

# test this is working
#ks.test(x = fishdata$Std_length_mm[1:10], y = fishdata$Std_length_mm[11:20])
#ks.test(x = 1:10, y = 100:110)
#ks.test(x = 1:10, y = 1:13)


# 4) Assign taxa to analyze (later move this into a loop for all taxa)
# and assign parameters
taxagroup <- "Myctophidae"
pars <- make_nominal_pars()

# 5) filter length data
taxadata <- dplyr::filter(fishdata_3_families, Family == taxagroup)

# 6) calculate length bins and P(l)
# P(l) is the probability of fish being in each length bin
Pl_and_l <- make_Pl(taxadata$Std_length_mm)

nl <- length(Pl_and_l$L)

# 7) For each L, calculate tau (total carbon flux for that length fish)
# Note: we take the first value of make_tau() using [1] because this is the 
# total C flux estimate--values [2] through [7] are the breakdown of C flux via
# respiration, SDA, egestion,excretion, and mortality 
tau_l <- rep(NA, times = nl)
for (i in 1:nl) {
  tau_l[i] <- make_tau(Pl_and_l$L[i], pars = pars, taxagroup = taxagroup, VM = TRUE, lwreg)[1]
}

# 8) Multiply P(L) times tau_l and sum
ind_C_flux <- sum(tau_l * Pl_and_l$P)

# 9) Calculate C_g_bar: avg count per m^2 for all tows
# note: from output of load_data() called thedata, the first value is fishdata
# and the second values is volumedata in the resulting list
# add "$c_t_bar" at end because calc_C_t_bar now finds mean as well as sd of catch across tows
# added Cook night tows 11, 48, 90 to the list of tows to include
C_t_bar <- calc_C_t_bar(fishdata = fishdata_3_families, volumedata = volumedata, towlist = c(2, 4, 6, 8, 9, 11, 48, 90))$c_t_bar

# 10) Calculate n_g for VM fish
# n_g is the number of fish of that taxonomic group, which is based on 
# average catch across all tows (C_t_bar) and capture efficiency (q). 
# If we assume we catch 10% of what was in the volume filtered prior to 
# net avoidance, then q = 0.1. Allow this to vary in sensitivity analysis. 
# Also assume that the calibration factor does not need to be corrected 
# in the nominal estimate, so cal = 1 
n_g_myct <- calc_n_g(catch = C_t_bar$Myctophidae, capture_efficiency = 0.1, cal = 1)

# 11) Calculate C flux for VM fish at study site
C_flux_study_site_VM <- ind_C_flux * n_g_myct # for q = 0.1, avg of all Sarmiento and Cook tows, myctophids, and 
# nominal parameter values, C flux = 5.63 mg C / m^2 --> updated after revised data processing in load_data(): 12.96 mg C / m^2 
# if only using Sarmiento tows, or 19.35 if also using all Cook tows day and night (10, 11, 47, 48, 89, 90), or 24.48 if using 
# only Cook night tows along with Sarmiento tows. This includes just Myctophidae as VM fish. --> updated Jan 4 with new raw data:
# myctophid flux = 24.69, and hatchet fish (below) contribute an additional 2.37 



# do same for Sternoptychidae, and ultimately add to the C flux for VM fish at study site 
# repeat steps 4 through 11 above


# 4) Assign taxa to analyze (later move this into a loop for all taxa)
# and assign parameters
taxagroup <- "Sternoptychidae"
pars <- make_nominal_pars()

# 5) filter length data
taxadata <- dplyr::filter(fishdata_3_families, Family == taxagroup)

# 6) calculate length bins and P(l)
# P(l) is the probability of fish being in each length bin
Pl_and_l <- make_Pl(taxadata$Std_length_mm)

nl <- length(Pl_and_l$L)

# 7) For each L, calculate tau (total carbon flux for that length fish)
# Note: we take the first value of make_tau() using [1] because this is the 
# total C flux estimate--values [2] through [7] are the breakdown of C flux via
# respiration, SDA, egestion,excretion, and mortality 
tau_l <- rep(NA, times = nl)
for (i in 1:nl) {
  tau_l[i] <- make_tau(Pl_and_l$L[i], pars = pars, taxagroup = taxagroup, VM = TRUE, lwreg)[1]
}

# 8) Multiply P(L) times tau_l and sum
ind_C_flux <- sum(tau_l * Pl_and_l$P)

# 9) Calculate C_g_bar: avg count per m^2 for all tows
# note: from output of load_data() called thedata, the first value is fishdata
# and the second values is volumedata in the resulting list
# add "$c_t_bar" at end because calc_C_t_bar now finds mean as well as sd of catch across tows
# added Cook night tows 11, 48, 90 to the list of tows to include
C_t_bar <- calc_C_t_bar(fishdata = fishdata_3_families, volumedata = volumedata, towlist = c(2, 4, 6, 8, 9, 11, 48, 90))$c_t_bar

# 10) Calculate n_g for VM fish
# n_g is the number of fish of that taxonomic group, which is based on 
# average catch across all tows (C_t_bar) and capture efficiency (q). 
# If we assume we catch 10%of what was in the volume filtered prior to 
# net avoidance, then q = 0.1. Allow this to vary in sensitivity analysis. 
n_g_hatchet <- calc_n_g(catch = C_t_bar$Sternoptychidae, capture_efficiency = 0.1, cal = 1)

# 11) Calculate C flux for VM fish at study site
C_flux_study_site_VM_hatchet <- ind_C_flux * n_g_hatchet # for q = 0.1, avg of all Sarmiento and Cook tows, hatchetfish, and 
# nominal parameter values, C flux = 2.37 mg C / m^2 

# 12) for all VM fish, add myctophid and hatchetfish C fluxes together
C_flux_study_site_VM_all <- C_flux_study_site_VM + C_flux_study_site_VM_hatchet 




#### Do same for NM fish #### 

# For NM fish, need to re-do steps but with the NM C flux function and taxagroup
# 7_NM) 

# set taxagroup
taxagroup <- "Gonostomatidae"

# use same nominal pars as for VM fish (make_nominal_pars generates both VM and NM pars)

# filter length data 
taxadata <- dplyr::filter(fishdata_3_families, Family == taxagroup)

# find length bins and P(l)
Pl_and_l_NM <- make_Pl(taxadata$Std_length_mm)

nl <- length(Pl_and_l_NM$L)

# find tau
tau_l_NM <- rep(NA, times = nl)
for (i in 1:nl) {
  tau_l_NM[i] <- make_tau(Pl_and_l_NM$L[i], pars = pars, taxagroup = taxagroup, VM = FALSE, lwreg)[1]
}

# multiply P(L) times tau_l and sum, for NM fish
# remove the first two tau_l_NM values than produce NaN
# (there are hardly any fish in this smallest size class anyway)
ind_C_flux_NM <- sum(tau_l_NM * Pl_and_l_NM$P)

# C_t_bar produces catch by family, so we don't need to re-calculate this for NM

# calculate n_g for NM fish
n_g_bristle <- calc_n_g(catch = C_t_bar$Gonostomatidae, capture_efficiency = 0.1, cal = 1)

# calculate C flux for NM fish at study site
C_flux_study_site_NM <- ind_C_flux_NM * n_g_bristle # for q = 0.1, avg of all Sarmiento tows, bristlemouths, and 
# nominal parameter values, C flux = 6.14 mg C 


# summarize fishdata_3_families by grouping by family and then summing density
fishdata_for_summary <- fishdata_3_families
fishdata_for_summary$density <- fishdata_for_summary$weight.to.use / fishdata_for_summary$VolFiltM3 * fishdata_for_summary$DepthInterval
fish_summary <- fishdata_for_summary %>% group_by(Family) %>% summarize(sum_density = sum(density))

# saved fishdata and volumedata as .csv files
# commented out write.csv after initially exporting the final analysis dataframe 
# write.csv(fishdata, file = "data/fishdata.csv", row.names = FALSE)
# write.csv(volumedata, file = "data/datavolumedata.csv", row.names = FALSE)
