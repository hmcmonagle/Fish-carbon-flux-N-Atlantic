# Functions used in North Atlantic fish carbon flux code
# Main code that sources these functions is "Chapter 2 EXPORTS fish C flux.Rmd"


# load packages
library(dplyr)
library(MASS)
library(TMB)
library(ggplot2)
library(gridExtra)
library(readxl)
library(Matrix)
library(mvtnorm)
library(KernSmooth)
library(tidyr) 
library(qgam)
library(stringr)


#### calculate vertical swim speed and proportion of day resting from acoustics ####

# acoustic data provided by Dr. Mei Sato (Woods Hole Oceanographic Institution) 
# using EK80 data, 18 kHz, from the R/V Cook can be found in the two
# Excel spreadsheets titled "Time_in_mesopelagic_migrators" and "Swim_speed_from_acoustics"
# in the "data" folder in this repository. 

# vertical swim speeds and proportion of day that migrating fish spend at 
# daytime depths of 300 m to 700 m are calculated in these spreadsheets 

# swim speed calculations are based on depths between about 250 and 50 m, 
# but if vertically migrating fish are swimming as deep as 700 m (an 
# echogram suggests they are) in a similar amount of time, then the 
# max swim speed is calculated below. 

# minium swim speed measured from two days of acoustic recordings:
min_swim_speed <- 2 # min in cm/s

# max swim speed measured from about 50 to 250 m (3.7 cm/s), multiplied by 
# 1.6 because some fish appear to be moving as deep as 500 m 
# not just 300 m in the same amount of migration time (and 500/300 = 1.6)
max_swim_speed <- 5.9 # 3.7 * 1.6 = 5.9 cm/s max 

# use the mean swim speed as nominal swim speed
mean_swim_speed <- (min_swim_speed + max_swim_speed)/2 # 4.0 cm/s mean

# to compare to the 30% +/- of the nominal swim speed used in McMonagle et al. 2023
# fish flux paper, calculate % change between mean and min or max
min_swim_speed_change <- (max_swim_speed - mean_swim_speed)/mean_swim_speed
# so, we use 5.255 cm/s as our nominal swim speed and allow this to vary roughly 62% +/-
# (but we'll just use the empirically-based min/max of 2 and 8.5 cm/s)

# proportion of day "resting" (i.e., at daytime depth for migrators), calculated 
# from acoustic data, varies from 
min_prop_day_resting <- 0.44
max_prop_day_resting <- 0.47
mean_prop_day_resting <- (min_prop_day_resting + max_prop_day_resting)/2 # 0.455

# note that this is accurate as of early May, but would vary throughout the year
# based on hours of sunlight at that time of year (assuming vertical migration
# is triggered by light availability). 


#### indexing function ####

find_index <- function(x,y) y <- which(y == x)

#### general nominal parameter values ####

make_nominal_pars <- function() {
  # nominal values for each parameter
  resp_model_center <- get_ikeda_coefs()$resp_model_coeff_center
  Wref <- 1000
  Tref <- 15
  a0 <- as.numeric(resp_model_center)[1]
  a1 <- as.numeric(resp_model_center)[2]
  a2 <- as.numeric(resp_model_center)[3]
  a3  <- as.numeric(resp_model_center)[4]
  H_VM  <- 1 # ** change to 0 for NM fish. also hard coded into make_tau()
  TT_epipelagic  <- 14
  TT_mesopelagic  <- 11
  TT_mean  <- 12.5
  activity_factor_resting  <- .75
  activity_factor_migrating  <- 2.5
  activity_factor_foraging  <- 2.5
  time_resting  <- mean(c(0.44, 0.47)) # found empirically, see information above about acoustic data. 
  # note that this is only applicable for this time of year--it would presumably
  # be lower if averaging over the course of a year (average would be lower due to 
  # less time spent at mesopelagic depth when there are fewer hours of daylight) 
  # See Excel sheet Time_in_mesopelagic_migrators" for more information
  Dmax  <- 500
  Dmin  <- 50
  B <- 200
  S <- mean_swim_speed # given in cm /s, calculated above
  prop_SDA  <- 0.14
  prop_egestion  <- 0.20
  prop_excretion  <- 0.07
  G_VM <- 0.018
  theta <- 0.77
  peakTime <- 3
  timeTo90 <- 12
  mu <- 0.606
  Q_ox_kJ_g_O2  <- 13.4 
  RQ  <- 1.35
  energy_in_pellet  <- 2 
  C_in_pellet <- 0.075 
  proportion_egestion_below_150m <- 0.7
  C_c_vm <- 0.1
  prop_mortality_below_150m <- 0.65
  X_r <- 0.008
  X_e <- 0.7
  
  # NM parameters
  H_NM  <- 0 # ** changed to 0 for NM fish
  G_NM <- 0.007
  C_c_nm <- 0.1 
  prop_non_detrital_NM_prey <- 0.6 # might in future change to range of zoop SIA values from 
  # Connor Shea's PhD work (Brian Popp's lab), from EXPORTS North Atlantic
  
  allpars <- list(Wref = Wref,
                  Tref = Tref,
                  a0=a0, 
                  a1=a1 , 
                  a2=a2, 
                  a3=a3, 
                  H_VM_or_NM=H_VM, 
                  TT_epipelagic=TT_epipelagic, 
                  TT_mesopelagic=TT_mesopelagic, 
                  TT_mean=TT_mean,  
                  activity_factor_resting=activity_factor_resting,
                  activity_factor_foraging=activity_factor_foraging,
                  activity_factor_migrating=activity_factor_migrating, 
                  time_resting=time_resting, 
                  Dmax=Dmax, 
                  Dmin=Dmin, 
                  B=B, 
                  S=S, 
                  prop_SDA=prop_SDA,
                  prop_egestion=prop_egestion, 
                  prop_excretion=prop_excretion, 
                  G_VM=G_VM, 
                  theta=theta,
                  peakTime=peakTime, 
                  timeTo90=timeTo90, 
                  Q_ox_kJ_g_O2=Q_ox_kJ_g_O2, 
                  RQ=RQ,
                  energy_in_pellet=energy_in_pellet, 
                  C_in_pellet=C_in_pellet, 
                  proportion_egestion_below_150m=proportion_egestion_below_150m, 
                  X_r=X_r, 
                  X_e=X_e, 
                  mu=mu, 
                  C_c_vm=C_c_vm, 
                  prop_mortality_below_150m=prop_mortality_below_150m,
                  H_NM = H_NM,
                  G_NM = G_NM,
                  C_c_nm =C_c_nm,
                  prop_non_detrital_NM_prey = prop_non_detrital_NM_prey)
  
  return(allpars)
}


load_data <- function() {
  
  ### load Sarmiento data, larger fish ####
  SG2105_MOC_data <- read.csv("data/MOCNESS_data/MOCNESS 10 Sarmiento/Sarmiento_MOC_fish_data.csv", header = TRUE)
  
  # make new column that uses taxa_name_final that is simply called "taxa_code" 
  # (see notes below on genetic vs morphological species identification)
  SG2105_MOC_data$taxa_code <- SG2105_MOC_data$taxa_code_final
  
  # filter out zooplankton (decapods, ostracods, shrimp, amphipods, pteropods, etc.)
  # and also filter out squid
  fish_data_Sarmiento <- SG2105_MOC_data %>% filter(
    taxa_code!="DECA" & taxa_code!="AMPH" & taxa_code!="RORO" & taxa_code!="ANSP" & taxa_code!="ACPU" & 
      taxa_code!="CRAN" & taxa_code!="SYDE" & taxa_code!="PNSP" & taxa_code!="ATOL" & taxa_code!="PTER" & 
      taxa_code!="CLSP" & taxa_code!="OSTR" & taxa_code!="PASI" & taxa_code!="EVSP")  	
  
  # also filter out families for which there is only n=1 individual. These are extremely
  # rare and thus affect carbon cycling negligibly, and there are not enough individuals
  # to determine reliable, empirical length-weight relationship for regression if there
  # are missing weights (we don't have weights for smaller Sarmiento fish, or Cook fish). 
  # These excluded families were Paralepididae (Arctozenus risso, n = 1, taxa code ARRI), 
  # Bathylagidae (Bathylagus europs, n = 1, taxa code BAEU), Microstomatidae 
  # (Nansenia sp., n = 1, NNSP), Platytroctidae (no species ID, n = 1, PLAT),  
  # Melamphaidae (Poromitra megalops, n = 1, MELA) and Evermannellidae (n = 1, EVSP).
  # Excluded Sigmops elongatus because although they are Gonostomatidae, there was only 
  # n = 1 and they are much larger than Cyclothone sp. so length-weight relationship may not
  # be accurate (nor assumptions about non-migratory behavior for Cyclothone). 
  # Excluded European hake (taxa code MEME)--not considered a small mesopelagic (n = 7).
  # Included Valenciennellus tripunctalatus (now known as tripunctatus) 
  # with their Sternoptychidae family for length-weight
  # relationship; although body shape is more distinct from other hatchetfish, n = 2
  # so not enough to make a separate lenght-weight relationship for X. copei
  # Total families included are Sternoptychidae, Myctophidae, Gonostomatidae and 
  # Alepocephalidae, although Alepocephalidae only have n = 9 individuals--
  # though analysis many ultimately focus on Myctohpids and bristlemouths as
  # they make up the bulk of biomass and are relatively clear migrators vs non-migrators. 
  
  # exclude some taxa from Sarmiento fish datasets
  fish_data_Sarmiento <- fish_data_Sarmiento %>% filter(taxa_code!="ARRI" & taxa_code!="BAEU" & taxa_code!="MEME" & taxa_code!="MELA" & taxa_code!="NNSP" & taxa_code!="PLAT")
  
  # sorted data above by standard length and noticed some issues in first few rows 
  #(e.g. a tiny cyclothone around .5mm (should this be cm?) long weighing more than a
  # fish in row 3 that is 11 mm long --> update: fixed this in .csv file, and 
  # let Lyndsey Lefebvre in Llopiz Lab know
  
  # remove fish that were so badly damaged that there's no standard length
  # (results in a loss of about 3% of the abundance)
  
  # n = 593 before removing NAs below
  fish_data_Sarmiento <- fish_data_Sarmiento %>% filter(Std_length_mm!="NA")
  # where there were NAs it was due to either an estimated weight/length in the 
  # handwritten datasheet (which I included to avoid removing too many NAs in terms
  # of underestimating abundance, though this only happened ~6 times) or the fish were 
  # in pieces upon thawing so length/weight could not be recorded accurately. 
  # individual 657 had an estimated weight but not length, so used length of another 
  # individual of similar size (24.9 mm for a mass of 0.1894, according to handwritten
  # data sheet)
  
  # n = 582 after removing NAs, so this excludes just 1.9% of the catch. 
  
  # remove unnecessary columns (for more details, see readme_SG2105_data.txt)
  # remove column "fork_length"
  fish_data_Sarmiento <- fish_data_Sarmiento %>% dplyr::select(-c(total_length, fork_length))
  
  # exclude fish that are smaller than 10 mm (we'll consider these larvae not 
  # juveniles/adults, and vertical migration behavior may be different than expected
  # based on studies of adults... though this is more relevant for the small Sarmiento
  # fish that were photographed only and measured in ImageJ, as these are all >10 mm 
  fish_data_Sarmiento <- fish_data_Sarmiento %>% filter(Std_length_mm>=10)
  
  # make a new column titled "Lowest_taxon" which uses the result of "lowest_taxa_final", 
  # which contains the finanl decision about taxonomic assignment after considering any 
  # discrepancies between the visual and genetic ID (in which case we went with genetic), 
  # and which contains any lower taxonomic level classification if available after genetic ID. 
  fish_data_Sarmiento$Lowest_taxon <- fish_data_Sarmiento$lowest_taxa_final
  
  # make new column of taxon that is just to family level (as opposed to lowest
  # taxonomic level classified, so that this matches the Cook taxon column)
  fish_data_Sarmiento$Taxon <- as.vector(rep(NA, times = length(fish_data_Sarmiento$Lowest_taxon)))
  
  # edited individual 78 from handwritten datasheet from GONO Gonostomatidae to 
  # CYSP Cyclothone sp, because Julia was only pulling morphologically identified
  # Cyclothone for ETS assays (first 100 fish numbers)
  
  fish_data_Sarmiento <- fish_data_Sarmiento %>% mutate(Taxon = case_when(
    Lowest_taxon == "Benthosema glaciale" |
      Lowest_taxon =="Lampanyctus sp." |
      Lowest_taxon == "Benthosema sp." |
      Lowest_taxon == "Diaphus subtilis" |
      Lowest_taxon == "Lampanyctus macdonaldi" |
      Lowest_taxon == "Lampanyctus sp." |
      Lowest_taxon == "Lampanyctus ater" |
      Lowest_taxon == "Myctophidae" |
      Lowest_taxon == "Myctophum punctatum" |
      Lowest_taxon == "Nannobrachium sp." |
      Lowest_taxon == "Notoscopelus resplendens" |
      Lowest_taxon ==  "Protomyctophum arcticum" |
      Lowest_taxon == "Symbolophorus" | 
      Lowest_taxon == "Symbolophorus veranyi" |
      Lowest_taxon == "Symbolophorus sp." |
      Lowest_taxon == "Myctophidae" ~
      "Myctophidae",
    Lowest_taxon == "Argyropelecus aculeatus" |
      Lowest_taxon =="Argyropelecus hemigymnus" |
      Lowest_taxon == "Argyropelecus lychnus" |
      Lowest_taxon == "Argyropelecus olfersii" |
      Lowest_taxon == "Argyropelecus sp." | 
      Lowest_taxon == "Valenciennellus tripunctulatus" |
      Lowest_taxon == "Sternoptychidae" ~
      "Sternoptychidae", 
    Lowest_taxon == "Cyclothone braueri" |
      Lowest_taxon =="Cyclothone microdon" |
      Lowest_taxon == "Cyclothone pallida" |
      Lowest_taxon == "Cyclothone sp." |
      Lowest_taxon == "Gonostomatidae" ~
      "Gonostomatidae", 
    Lowest_taxon == "Scopelogadus beanii" |
      Lowest_taxon == "Scopelogadus sp." |
      Lowest_taxon == "Melamphaidae" ~
      "Melamphaidae", 
    Lowest_taxon == "Xenodermichthys copei" |
      Lowest_taxon == "Alepocephalidae" ~ 
      "Alepocephalidae"))
  
  # change fish_data_Sarmiento df name to distinguish between small Sarmiento fish
  fish_Sarmiento_large <- fish_data_Sarmiento
  
  # filter out the one Sigmops elongatus (n = 1, and this bristlemouth
  # behaves differently than other Cyclothone bristlemouths that are assumed 
  # to be non-migrators). Total "large" (separated from zooplankton at sea) Sarmiento
  # fish are reduced from 582 to 581 at this point. 
  fish_Sarmiento_large <- fish_Sarmiento_large %>% dplyr::filter(is.na(Lowest_taxon) | Lowest_taxon != "Sigmops elongatus")
  
  # add a column for Net_Tow, which will be used later for merging fish data with volume data
  fish_Sarmiento_large$Tow_net <- paste(fish_Sarmiento_large$Tow_number, "_", fish_Sarmiento_large$Net_number)
  
  # add a column for Ship for large fish df, which the smaller fish df has (and
  # this will help keep Sarmiento data straight from Cook data in any merged df)
  fish_Sarmiento_large$Ship <- rep("Sarmiento", times = length(fish_Sarmiento_large$Tow_number))
  
  
  ### load smaller Sarmiento fish ####
  
  # now, upload "small" Sarmiento fish, which were those that are generally smaller
  # than those that were sorted from zooplankton at sea. These were sorted from 
  # zooplankton later and then photographed and measured in ImageJ.
  fish_Sarmiento_small <- read.csv("data/MOCNESS_data/MOCNESS 10 Sarmiento/SG2105 All Small Fish Measurements.csv", header = TRUE)
  
  # filter out ichthyoplankton
  fish_Sarmiento_small <- fish_Sarmiento_small %>% filter(Taxon!="Larvae" & Taxon!="Leptocephalus larvae")
  
  # delete rows with no standard length
  fish_Sarmiento_small <- fish_Sarmiento_small %>% filter(Std_length_mm!="#VALUE!")

  # make length numeric
  fish_Sarmiento_small$Std_length_mm <- as.numeric(fish_Sarmiento_small$Std_length_mm)
  
  # exclude fish that are smaller than 10 mm (we'll consider these larvae not juveniles/adults)
  fish_Sarmiento_small <- fish_Sarmiento_small %>% filter(Std_length_mm>=10)
  
  # duplicate each fish 4x since these were quarter splits (and we'll assume
  # capture efficiencies did not differ greatly between fish lengths since these
  # were all relatively small fish--only about 11% of these "small fish" were > 30 mm)
  
  # each fish ID is unique, so duplicate each row with a unique fish ID four times
  fish_Sarmiento_small <- fish_Sarmiento_small[rep(c(1:nrow(fish_Sarmiento_small)), 4),]
  
  # adjust length to be 5% larger that what was measured, because these fish were
  # stored in ethanol, which shrinks them by about 5% in length (see Moku et al. 2004)
  fish_Sarmiento_small$Std_length_mm <- fish_Sarmiento_small$Std_length_mm * 1.05
  
  # remove column with earlier estimation of weights from lengths. These weights
  # were based on a study of Mediterranean fish, which Julia Cox had given a shot 
  # in Llopiz Lab after the cruise. Instead, we'll use the length-weight model
  # above, so we can base estimates on other fish caught in the same time and place
  fish_Sarmiento_small <- fish_Sarmiento_small %>% dplyr::select(-Weight.estimate..mg...Mediterranean.)
  
  # add a column in small fish dataframe that has c("Tow_number", "_", "Net_number")
  fish_Sarmiento_small$Tow_net <- paste(fish_Sarmiento_small$Tow_number, "_", fish_Sarmiento_small$Net_number)
  
  # remove standard length in cm column to avoid confusion with length in mm column
  fish_Sarmiento_small <- fish_Sarmiento_small %>%
    dplyr::select(-Standard.Length..cm.)
  
  # combine large and small fish data
  fish_Sarmiento_all <- bind_rows(fish_Sarmiento_small, fish_Sarmiento_large)
  
  
  ### R/V Sarmiento de Gamboa volume filtered data from MOC-10 ####
  
  # Note that on Sarmiento, we called net 0 the one that fished from 0-1000 m, 
  # and on Cook this was called net 1. 
  
  # MOC-10 small fish have been multiplied by 4 because they came from a split
  
  # load Sarmiento data. has extra columns to match metadata format from R/V Cook MOCNESS-1
  MOC_data_Sarmiento <- read.csv("data/MOCNESS_data/MOCNESS 10 Sarmiento/EXPORTSMocNetDataSarmiento.csv", header = TRUE)
  
  # Remove unnecessary columns so that merged dataframe isn't so large
  Sarmiento_vol_filtered <- MOC_data_Sarmiento %>% dplyr::select(-c(CTNID, Ship, Cruise, DateTimeStart, DateTimeLocal, TimeStart, DateTimeEnd, DateTimeLocalEnd, 
                                                                    TimeEnd, StudyName, CruiseTow, Year, CruiseName, CruiseDateStart, CruiseDateEnd, TZZoneName, 
                                                                    CruisePersID, R2R_Event, Epoch, DayEpoch, SolarElevationCorr, Heading, SOG, WaterDepth, WindSpeed, 
                                                                    WindDirection, Weather, CloudCover, Seas, SST, SSS, TowPersID, NetID, NetName, NetLength, 
                                                                    LengthHorizontal, AreaM2, TowType, WireOut, R2R_EventCTD, Cast, MOCDNPair, PreCast, PostCast, 
                                                                    FirstMCN, PayRetrievRate, FracPreserved, FracETS, FracGutFluor, FracGutContents, CruiseNotes, 
                                                                    NetNotes, CTNNotes, TimeLocalEnd, DurationTow, LengthVertical, DepthTarget, DepthMaximum, 
                                                                    TimeNet, Dpth, AngleNet, FlowCounts))
  
  # make tow number a factor in the metadata dataframe
  Sarmiento_vol_filtered$Tow_number <- as.factor(Sarmiento_vol_filtered$Tow_number)
  
  Sarmiento_vol_filtered$Tow_net <- paste(Sarmiento_vol_filtered$Tow_number, "_", Sarmiento_vol_filtered$Net_number)
  
  # remove rows 3 and 5, which have VolFiltM3 values of NA because these nets did not fish (issue with net bars, see MOCNESS event log)
  Sarmiento_vol_filtered <- Sarmiento_vol_filtered[-c(3,5),]
  
  # remove net 0 because this was towed at a different speed from other nets, 
  # and is redundant, (and, net 0 was not saved from Cook data so we won't be 
  # using that net for cross-comparison of biomass data)
  Sarmiento_vol_filtered <- Sarmiento_vol_filtered %>% filter(Net_number!="0")
  
  # calculate depth interval sampled
  Sarmiento_vol_filtered$DepthStart <- as.numeric(Sarmiento_vol_filtered$DepthStart)
  Sarmiento_vol_filtered$DepthEnd <- as.numeric(Sarmiento_vol_filtered$DepthEnd)
  Sarmiento_vol_filtered$DepthInterval <- Sarmiento_vol_filtered$DepthStart - Sarmiento_vol_filtered$DepthEnd
  

  ### Now, merge fish data with the volume data ###
  
  # merge so that df has both fish data and vol filtered data (object name includes fish_vol)
  # by combining Sarmiento fish data with Sarmiento MOCNESS-10 metadata
  fish_vol_Sarmiento_all <- inner_join(fish_Sarmiento_all, Sarmiento_vol_filtered, by = "Tow_net")
  
  # delete extra tow and net column from one df, so that you don't have duplicates
  # after merging (they'd be renamed Tow_number.x and Tow_number.y, for example)
  fish_vol_Sarmiento_all <- dplyr::select(fish_vol_Sarmiento_all, -c("Tow_number.y", "Net_number.y"))
  
  # rename back to original column name
  fish_vol_Sarmiento_all <- fish_vol_Sarmiento_all %>%
    dplyr::rename("Tow_number" = "Tow_number.x", 
                  "Net_number" = "Net_number.x")
  
  # check how many fish are missing empirically measured weights 
  rows_missing_weights_before <- which(is.na(fish_Sarmiento_large$total_weight))
  # we'll fill in these missing weights later using a length-weight regression:
  # see 01-calc_C_flux.R script. We'll need to fill in all the weights for the Cook
  # fish using this regression too, as Cook fish only had lengths measured. 
  
  
  
  ### load R/V Cook fish data ####
  
  # These are data from the other ship that collected fish, but using a MOCNESS-1 instead of a MOCNESS-10
  # These are the (generally larger) fish that were pulled from the net and photographed
  # at sea, rather than the smaller fish that were imaged and measured later using Ecotaxa and Zooscan
  fish_data_Cook_large <- read.csv("data/MOCNESS_data/MOCNESS 1 Cook/Cook MOC fish data.csv", header = TRUE)
  
  # note: edited raw data file EXPORTSMocNetDataCook.csv to make the two Tow 12 
  # nets called nets 9 and 10 (net 10 being 0-50 m), instead of how it was before 
  # where the final net (intended to be used in place of the faulty net 10 in 
  # Tow 11) was called net 2 (because it was one of only 2 nets deployed in that 
  # short tow: net 1 went from 0-50 m, and then net 2 came back up from 50-0 m). 
  # In other words, tow 12 nets were renamed to be part of tow 11, as the 
  # shallowest net 10 from tow 12 replaced the faulty net 10 from tow 11. 
  
  MOC_data_Cook <- read.csv("data/MOCNESS_data/MOCNESS 1 Cook/EXPORTSMocNetDataCook.csv", header = TRUE)
  
  # tows to use from Cook are 10, 11, 47, 48, 89 and 90 because they are 
  # paired in time and place and did not have issues with deployment (tow 1 only
  # went to 100 m and at the end of tow 62, nets came up severely tangled so 
  # catch was not quantified)
  
  # filter fish_data_Cook to only include Tow_number == 10, 11, 47, 48, 89 and 90
  fish_data_Cook_large <- fish_data_Cook_large %>% filter(Tow_number %in% c(10, 11, 47, 48, 89, 90))
  
  # filter MOC_data_Cook to only include Tow_number == 10, 11, 47, 48, 89 and 90
  MOC_data_Cook <- MOC_data_Cook %>% filter(Tow_number %in% c(10, 11, 47, 48, 89, 90))

  # Remove unnecessary rows so that merged dataframe isn't much larger than we need
  # (these unused columns came from the standard MOCNESS metadata sheet from Amy Maas's lab)
  Cook_vol_filtered <- MOC_data_Cook %>% dplyr::select(-c(CTNID, Ship, Cruise, DateTimeStart, DateTimeLocal, TimeStart, DateTimeEnd, DateTimeLocalEnd, 
                                                          TimeEnd, StudyName, CruiseTow, Year, CruiseName, CruiseDateStart, CruiseDateEnd, TZZoneName, 
                                                          CruisePersID, R2R_Event, Epoch, DayEpoch, SolarElevationCorr, Heading, SOG, WaterDepth, WindSpeed, 
                                                          WindDirection, Weather, CloudCover, Seas, SST, SSS, TowPersID, NetID, NetName, NetLength, 
                                                          LengthHorizontal, AreaM2, TowType, WireOut, R2R_EventCTD, Cast, MOCDNPair, PreCast, PostCast, 
                                                          FirstMCN, PayRetrievRate, FracPreserved, FracETS, FracGutFluor, FracGutContents, CruiseNotes, 
                                                          NetNotes, CTNNotes, TimeLocalEnd, DurationTow, LengthVertical, DepthTarget, DepthMaximum, 
                                                          TimeNet, Dpth, AngleNet, FlowCounts))
  
  # make tow number a factor in the metadata dataframe
  Cook_vol_filtered$Tow_number <- as.factor(Cook_vol_filtered$Tow_number)
  
  # there was no data from tow 11 net 10 (aka net 9, 0-50 m, something went wrong
  # with that net), and instead the net 10 to be used for tow 11 is
  # data from a subsequent tow 12 was used to fill this data gap. Fixed this in original data file
  # by deleting the faulty tow 11 net 10 data and replacing with tow 12 net 10 data in
  # csv files to avoid confusion. In summary: Tow 12 net 2 takes the place of the
  # missing tow 11 net 10, and is renamed tow 11 net 10. 
  
  ### load Cook small fish data ####
  
  ## Upload smaller Cook fish that were measured from Zooscan (Ecotaxa) images
  # See readme on Ecotaxa fish measurements for more information
  fish_data_Cook_small <- read.csv("data/MOCNESS_data/MOCNESS 1 Cook/Ecotaxa_small_Cook_fish/Ecotaxa_small_fish_Cook.csv", header = TRUE)

  
  # *** update April 2024: Joe Cope in Deb's lab found some errors in the raw data, so 
  # Ecotaxa_small_fish_Cook.csv has been updated so that the splits information is correct. 
  # Basically, needed to fix some of the split_correct values and also needed to remove 
  # Ecotaxa images and fish legnths that had already been measured from the photographs
  # taken at sea (in Cook MOC fish data.csv). See "Ecotaxa_small_fish_Cook_red_delete_green_add_yellow_edit.xlsx"
  # for rows that were added or deleted from the original "Ecotaxa_small_fish_Cook.csv" file, and see 
  # "Ecotaxa_small_fish_Cook_including_lessthan_10mm.xlsx" for the full set of Ecotaxa images and a note on the
  # far right column about whether this Ecotaxa image was included or not in the final analysis and why. 
  # Given the code below, these edits to the raw data files will be incorporated into fishdata.csv. 
  
  # To see edits between old version of fishdata.csv prior to April 2024 and after April 2024, see 
  # fishdata-highlighted-edits.xlsx". This file fishdata-highlighted-edits will be used for editing 
  # the version of the data uploaded to SeaBASS that Inia already formatted for SeaBASS. Delete/add rows
  # as instructed "fishdata-highlighted-edits.xlsx", where yellow rows need to be deleted from SeaBASS
  # because they represent and unnecessary split replicate (this originated from errors in split column where
  # fish actually came from a whole sample), and orange rows need to be deleted because this was a duplicate between
  # fish lengths originating from photographs taken at sea and lengths taken from Ecotaxa images (ie same fish
  # were measured and added to datasheet twice, so we'll delete the Ecotaxa rows). Green rows represent those 
  # that were added from the new Ecotaxa image upload in April 2024 

  # remove columns that aren't needed from the Ecotaxa fish data frame
  # remove Tow_vol column and use volumes from MOCNESS metadata
  # Delete DepthStart, DepthEnd and DepthInterval because these will be added
  # back in later for all Cook fish using tow and net number alignment
  # with the MOCNESS metadata
  fish_data_Cook_small <- fish_data_Cook_small %>% dplyr::select(-c(fraction, sub_part_correct, Tow_Vol, gray_mode, area_mm2, major_mm, 
                                                                    minor_mm, feret_mm, esd_mm, vol, Taxa, L_D, DW, O2_umol, 
                                                                    CO2, DepthStart, DepthEnd, DepthInterval))
  
  
  # duplicate each fish according to the split it came from. 
  # Where split_correct==1, no correction is needed. Where split_correct==0.5, 
  # we need to duplicate those fish. When split_correct==0.25, multiply each of those rows by 4. 
  # Where split_correct==0.125, multiply by 8. Do this by dividing 1/split_correct below
  fish_data_Cook_small <- fish_data_Cook_small[rep(c(1:nrow(fish_data_Cook_small)), 1/fish_data_Cook_small$split_correct),]
  
  # add a column in small fish dataframe that has c("Tow_number", "_", "Net_number")
  fish_data_Cook_small$Tow_net <- paste(fish_data_Cook_small$Tow_number, "_", fish_data_Cook_small$Net_number)
  
  # adjust length to be 2% larger that what was measured, because these fish were
  # stored in formalin, which shrinks them by about 2% in length (see Moku et al. 2004)
  fish_data_Cook_small$Std_length_mm <- fish_data_Cook_small$Std_length_mm * 1.02
  
  # add columns in large fish dataframe that the small fish dataframe has, for successful bind_rows() later
  fish_data_Cook_large$Image.Name <- NA
  fish_data_Cook_large$split_correct <- 1 # all these are from a whole sample, so split_correct = 1
  
  # add a column in large fish dataframe that has c("Tow_number", "_", "Net_number")
  fish_data_Cook_large$Tow_net <- paste(fish_data_Cook_large$Tow_number, "_", fish_data_Cook_large$Net_number)

  # combine large and small fish data
  fish_Cook_all <- bind_rows(fish_data_Cook_small, fish_data_Cook_large)
  
  Cook_vol_filtered$Tow_net <- paste(Cook_vol_filtered$Tow_number, "_", Cook_vol_filtered$Net_number)
  
  # merge so that each fish has vol filtered
  # by combining Cook fish data with Cook MOCNESS-1 metadata
  fish_vol_Cook_all <- inner_join(fish_Cook_all, Cook_vol_filtered, by = "Tow_net")
  
  # delete extra tow and net column from one df, so that you don't have duplicates
  # after merging (they'd be renamed Tow_number.x and Tow_number.y, for example)
  fish_vol_Cook_all <- dplyr::select(fish_vol_Cook_all, -c(Tow_number.y, Net_number.y, DayNight.y))
  
  # rename back to original column name
  fish_vol_Cook_all <- fish_vol_Cook_all %>%
    dplyr::rename("Tow_number" = "Tow_number.x", 
           "Net_number" = "Net_number.x", 
           "DayNight" = "DayNight.x")
  
  # edit maxVolFiltM3 in fish_vol_Cook_all to VolFiltM3
  fish_vol_Cook_all <- fish_vol_Cook_all %>%
    dplyr::rename("VolFiltM3" = "maxVolFiltM3")
  
  # rename to be split_correct_Cook to be clear in df that this is different from the Sarmiento splits column
  fish_vol_Cook_all <- fish_vol_Cook_all %>%
    dplyr::rename("split_correct_Cook" = "split_correct")
  
  
  # find columns that are missing from Sarmiento df that are present in Cook df 
  Sarmiento_cols_to_add <- setdiff(colnames(fish_vol_Cook_all), colnames(fish_vol_Sarmiento_all))
  # this is the only column that is missing from Sarmiento df that is present in Cook df. 
  
  # delete the columns in Cook_cols_to_edit from Cook data frame
  fish_vol_Cook_all <- dplyr::select(fish_vol_Cook_all, -c(LatDegStart, LatMinStart, LongDegStart, 
                                                           LongMinStart, LatDegEnd, LatMinEnd, LongDegEnd, 
                                                           LongMinEnd))

  # run again: find columns that are missing from Sarmiento df that are present in Cook df 
  Sarmiento_cols_to_add <- setdiff(colnames(fish_vol_Cook_all), colnames(fish_vol_Sarmiento_all))
  
  # find columns that are missing from Cook df that are present in Sarmiento df
  Cook_cols_to_add <- setdiff(colnames(fish_vol_Sarmiento_all), colnames(fish_vol_Cook_all))
  
  # add these columns (with NA values) into Cook df and save as a new data frame
  fish_vol_Cook_all[, Cook_cols_to_add] <- NA
  
  # add these needed Cook columns to Sarmiento df and save
  fish_vol_Sarmiento_all[, Sarmiento_cols_to_add] <- NA
  
  # fix DepthInterval column in Cook df to be a single value (e.g., 250), 
  # rather than a range that can't be used in calculations (e.g., 1000-750)
  fish_vol_Cook_all$DepthInterval <- fish_vol_Cook_all$DepthStart - fish_vol_Cook_all$DepthEnd
  
  # combine Cook and Sarmiento data
  fish_vol_all <- rbind(fish_vol_Sarmiento_all, fish_vol_Cook_all)
  
  # add calculated weights to df to later fill in gaps where we don't have
  # an empirically measured weight (see script 01-calc_C_flux.R**)
  
  # note that any fish from Sarmiento catch with a Note, e.g. "head missing", 
  # do not have an accurate length because of this missing head. 
  # Better to still include that individual for the 
  # sake of the count and biomass estimate than to remove all broken individuals
  # from the analysis. There are about 128 individuals among small Sarmiento
  # fish with a note, of 828 total small Sarmiento fish. So, we'd lose 
  # about 10% of the catch if we removed all those damaged fish. Therefore, we'll 
  # keep these damaged fish in, removing just the "head only" and "head fragment" 
  # rows, and add a rough 25% of measured length to length of fish that are missing 
  # a head (body only, body fragment, headless)
  
  fish_vol_all$Notes <- as.character(fish_vol_all$Notes)
  
  # number of length measurements from Ship==Sarmiento with a note
  sum(!is.na(fish_vol_all$Notes) & fish_vol_all$Ship=="Sarmiento")
  
  # number of length measurements from Ship==Sarmiento with a note that is not "head only", 
  sum(!is.na(fish_vol_all$Notes) & fish_vol_all$Ship=="Sarmiento" & !str_detect(fish_vol_all$Notes, "head only"))
  # So, 128 - 112 = 16 fish with a note that is "head only"
  
  # number of length measurements from Ship==Sarmiento with a note that is not "head fragments"
  sum(!is.na(fish_vol_all$Notes) & fish_vol_all$Ship=="Sarmiento" & !str_detect(fish_vol_all$Notes, "head fragment"))
  # 128 - 124 = there are just 4 fish that are a "head fragment"
  
  # Remove all the "head only" and "head fragment" fish from the catch, or ~1% 
  # of the catch (presumably these head only
  # or head fragment individuals belong with a body that was found, so we don't 
  # want to duplicate those individuals in the dataset). Note that lengths are
  # of course not as accurate when head is missing, but we decided to leave in
  # to avoid underestimating abundance. 
  fish_vol_all <- fish_vol_all[-c(
    which(fish_vol_all$Notes=="head only")),]
  
  # remove the rows with "head fragments" in the Notes column
  fish_vol_all <- fish_vol_all[-c(
    which(fish_vol_all$Notes=="head fragments")),]
  
  # now there are 1,465 individuals total --> update as of April 2024: these numbers
  # will vary slightly after raw data edits from Deb Steinberg and Joe Cope 
  
  # add 25% to length of individual fish that have "body only", "body fragment", 
  # or "headless" written in the Notes column
  
  # check mean length before adjusting for missing heads (a rough correction
  # for imperfect samples)
  mean(fish_vol_all$Std_length_mm)
  
  fish_vol_all$Std_length_mm <- case_when(
    fish_vol_all$Notes == "body only" ~ fish_vol_all$Std_length_mm * 1.25,
    fish_vol_all$Notes == "body fragment" ~ fish_vol_all$Std_length_mm * 1.25,
    fish_vol_all$Notes == "headless" ~ fish_vol_all$Std_length_mm * 1.25,
    TRUE ~ fish_vol_all$Std_length_mm)
  
  mean(fish_vol_all$Std_length_mm) # mean length increased slightly as intended
  
  # make tow number a factor (for plotting as a discrete variable)
  fish_vol_all$Tow_number <- as.factor(fish_vol_all$Tow_number)
  
  # add column to Sarmiento data for family (for consistent taxonomic
  # naming between Cook and Sarmiento fish)
  fish_vol_all <- fish_vol_all %>% mutate(Family = case_when(
    Taxon == "Myctophidae" |
      Taxon == "myctophid" ~
      "Myctophidae",
    Taxon == "Sternoptychidae" |
      Taxon == "hatchetfish" ~
      "Sternoptychidae", 
    Taxon == "Gonostomatidae" |
      Taxon == "bristlemouth" ~
      "Gonostomatidae", 
    Taxon == "Melamphaidae" |
      Taxon == "melamphid" ~
      "Melamphaidae", 
    Taxon == "Alepocephalidae" ~ 
      "Alepocephalidae", 
    Taxon == "Paralepididae" ~
      "Paralepididae"))
  
  # number of fish without an identified family
  number_NA_as_family <- length(which(is.na(fish_vol_all$Family))) # 71 individuals
  
  # number of fish with an identified family
  number_fish_with_family <- length(which(!is.na(fish_vol_all$Family))) # 1394 individuals
  
  # so, only ~5% of fish are missing a family. Exclude these from analysis--
  # we'll focus on myctophids, bristlemouths, and hatchetfishes
  # which make up the vast majority of the catch and have known behavioural patterns. 
  
  # number of individuals from Cook
  number_fish_Cook <- length(which(fish_vol_all$Ship == "Cook")) # 233 individuals originally, 
  # and then 279 after adding in the smaller Cook fish from Ecotaxa (adjusted for splits) 
  
  # number of individuals from Sarmiento
  number_fish_Sarmiento <- length(which(fish_vol_all$Ship == "Sarmiento")) # 1277 individuals
  
  number_fish_total <- length(fish_vol_all$Ship) # updated total fish after adding Ecotaxa fish: 1510
  
  # number of fish either in Family == "Gonostomatidae" or "Family" == "Myctophidae"
  number_fish_Mycto_Gono <- length(which(fish_vol_all$Family == "Gonostomatidae" |
                                           fish_vol_all$Family == "Myctophidae")) # 1183 individuals
  
  # percent of individuals we leave out if we leave out fish that are not 
  # part of Myctophidae or Gonostomatidae families
  number_fish_Mycto_Gono/number_fish_total # ~80% of individuals
  
  # add hatchetfish into analysis
  number_fish_Mycto_Gono_Sterno <- length(which(fish_vol_all$Family == "Gonostomatidae" |
                                                  fish_vol_all$Family == "Myctophidae" |
                                                  fish_vol_all$Family == "Sternoptychidae")) # 1211 individuals
  number_fish_Mycto_Gono_Sterno/number_fish_total # ~90% of individuals
  
  # summarize percent of total number of fish in each family
  fish_vol_all %>% group_by(Family) %>% summarize(n = n()) %>% mutate(percent = n/number_fish_total)
  
  # hatchetfish, bristlemouths and myctophids together make up ~90% of the catch
  
  # now, combine both Sarmiento and Cook volume filtered data
  
  # first, find columns that are missing from volumes Sarmiento df that are present in Cook volumes df 
  Cook_cols_to_edit <- setdiff(colnames(Sarmiento_vol_filtered), colnames(Cook_vol_filtered))
  
  # rename Cook Volume Filtered column name to match that in Sarmiento df
  Cook_vol_filtered <- Cook_vol_filtered %>%
    dplyr::rename("VolFiltM3" = "maxVolFiltM3")
  
  # find columns that are missing from Cook df that are present in Sarmiento df
  Cook_cols_to_remove <- setdiff(colnames(Cook_vol_filtered), colnames(Sarmiento_vol_filtered))
  
  # delete the columns in Cook_cols_to_remove from Cook data frame
  Cook_vol_filtered <- dplyr::select(Cook_vol_filtered, -c(LatDegStart, LatMinStart, LongDegStart, LongMinStart, LatDegEnd, LatMinEnd, LongDegEnd, LongMinEnd))

  # add column for Ship, so we can separate the dataset out by ships (and associated 
  # net type) if we want to (Sarmiento used MOC-10 and Cook used MOC-1)
  Cook_vol_filtered$Ship <- "Cook"
  Sarmiento_vol_filtered$Ship <- "Sarmiento"
  
  # remove net 1 from Cook data because net 1 (equivalent to net 0 on Sarmiento, 
  # sampling from 0-1000 m on the way down) was not saved from Cook data so, we won't be 
  # using that net for cross-comparison of biomass data
  Cook_vol_filtered <- Cook_vol_filtered %>% filter(Net_number!="1")
  
  # edit DepthInterval in Cook rows to be DepthStart-DepthEnd
  Cook_vol_filtered$DepthInterval <- Cook_vol_filtered$DepthStart - Cook_vol_filtered$DepthEnd
  
  vol_filtered_all <- rbind(Sarmiento_vol_filtered, Cook_vol_filtered)
  
  return(list(fish_vol_all, vol_filtered_all))
}

  
#### Function to find C_g_t_bar ####


calc_C_t_bar <- function(fishdata, volumedata, towlist = c(2, 4, 6, 8, 9, 11, 48, 90)){
  
  # specify which families to loop through
  groups.2.use <- unique(fishdata$Family)
  
  c_t <- list() # placeholder to hold estimated c_g_t for each taxon (group, abbreviated g) and tow (abbreviated t)
  
  # loop through all groups
  for (i in 1:length(groups.2.use)) {
    loop_group <- groups.2.use[i]
    c_g_t <- NULL
    
    #loop through tows, calculate c_g_t_hat for each and place in array
    for (tow in towlist) {
      c_g_t_vol <- fishdata %>% filter(Family == loop_group,  Tow_number == tow) %>% 
        group_by(Net_number) %>% 
        summarize(n = n())
      
      # for each net, get volume filtered
      volfilt_t <- volumedata %>%
        filter(Net_number %in% c_g_t_vol$Net_number, Tow_number == tow) %>%
        dplyr::select(Net_number, VolFiltM3)
      
      # for each net, get Depthinterval
      deltaDt <- volumedata %>%
        filter(Net_number %in% c_g_t_vol$Net_number, Tow_number == tow) %>%
        dplyr::select(Net_number, DepthInterval)
      
      # add onto existing array, inserting estimate for new tow.  If there is no existing array, it will begin one
      c_g_t <- c(c_g_t, sum(c_g_t_vol$n / volfilt_t$VolFiltM3 * deltaDt$DepthInterval))
    }
    # save group array to list
    c_t[[i]] <- c_g_t
  }
  names(c_t) <- groups.2.use
  
  # get mean by group (count per m^2)
  c_t_bar <- lapply(X = c_t, FUN = mean)
  c_t_sd <- lapply(X = c_t, FUN = sd)
  
  return(list(c_t_bar = c_t_bar, c_t_sd = c_t_sd))
  
}



#### Function to calculation n_g ####
# find standardized catch (n_g) in number of fish per group, 
# given that capture efficiency is not 100% (i.e., q does not equal 1).
# cal = calibration factor correction, which we allow to vary from 1 to 1.5
# This is because our calibration factors are close to 4 but they could in reality
# be as high as about 6 or 50% higher, so we allow this correction factor
# to vary from 1 to 1.5 and divide the catch by this factor becuase 
# higher calibration factor --> higher volume --> lower catch

calc_n_g <- function(catch, taxon, capture_efficiency, cal){
  n_g <- catch / capture_efficiency / cal
  return(n_g)
}


#### Functions to calculate total weight of each fish where missing ####

# Where total weight wasn't measured empirically, estimate total weight
# of each individual fish using lenghts and the following length-weight model by taxon.
# This is based on the method in Thorson et al. 2017 and is further explained in Appendix B. 

# Function to estimate L - W conversion parameters via phylogenetic hierarchical model ####
calc_lw <- function() {
  all.dat <- readRDS(file = "data/length_weight_reg_data/alldata_taxonomy.RDS")
  
  # collapse species when genus is present
  spcnaIndex <- which(is.na(all.dat$Species))
  gennaIndex <- which(is.na(all.dat$Genera))
  
  spc_only <- setdiff(spcnaIndex, gennaIndex)
  
  for (i in 1:length(spc_only)) all.dat$Species[spc_only[i]] <- paste0(all.dat$Genera[spc_only[i]], " spc")
  for (i in 1:length(gennaIndex)) {
    all.dat$Species[gennaIndex[i]] <- paste0(all.dat$Family[gennaIndex[i]], " spc")
    all.dat$Genera[gennaIndex[i]] <- paste0(all.dat$Family[gennaIndex[i]], " gen")
    
  }
  
  ### Setup TMB data and parameters ####
  #### Create new ParentChild matrix for reduced taxonomic structure ####
  
  taxa.list <- c("Order", "Family", "Genera", "Species")
  Z_ik_main <- dplyr::select(all.dat, all_of(taxa.list))
  Z_ik <- unique(Z_ik_main, MARGIN = 1)
  ParentChild_gz = NULL
  # 1st column: child taxon name
  # 2nd column: parent taxon name
  # 3rd column: parent row-number in ParentChild_gz
  # 4th column: Taxon level
  # Loop through
  for( colI in 1:ncol(Z_ik)){
    Taxa_Names = apply( Z_ik[,1:colI,drop=FALSE], MARGIN=1, FUN=paste, collapse="_")
    Unique_Taxa = unique(Taxa_Names)
    for( uniqueI in 1:length(Unique_Taxa) ){
      Which = which( Taxa_Names == Unique_Taxa[uniqueI] )
      if( colI==1 ){
        ParentChild_gz = rbind( ParentChild_gz, c(Unique_Taxa[uniqueI], NA, NA, colI) )
      }else{
        if( length(unique(Z_ik[Which,colI-1]))>1 ) stop("Taxa has multiple parents")
        ChildName = Unique_Taxa[uniqueI]
        ParentName = paste(rev(rev(strsplit(ChildName,"_")[[1]])[-1]),collapse="_")
        ParentChild_gz = rbind( ParentChild_gz, c(ChildName, ParentName, match(ParentName,ParentChild_gz[,1]), colI) )
      }
    }
  }
  
  # Relabel
  ParentChild_gz = data.frame( ParentChild_gz )
  colnames(ParentChild_gz) = c("ChildName", "ParentName", "ParentRowNumber", "ChildTaxon")
  ParentChild_gz[,'ParentRowNumber'] = as.numeric(as.character(ParentChild_gz[,'ParentRowNumber']))
  ParentChild_gz[,'ChildTaxon'] = as.numeric(as.character(ParentChild_gz[,'ChildTaxon']))
  PC_gz<- as.matrix(ParentChild_gz[, c('ParentRowNumber', 'ChildTaxon')]) - 1
  
  # Identify location for every observation
  Taxa_Names = apply( Z_ik, MARGIN=1, FUN=paste, collapse="_")
  g_i = match( Taxa_Names, ParentChild_gz[,'ChildName'] )
  n_k = ncol(Z_ik)
  n_j = 2 # two traits
  n_g = nrow(ParentChild_gz)
  n_i <- length(g_i)
  
  #### Create index of data to Parent - Child ####
  #Z_ik_dat <- dplyr::select(all.dat, Class, Order, Family, Species)
  Taxa_Names_dat <-  apply( Z_ik_main, MARGIN=1, FUN=paste, collapse="_")
  g_i_dat = match( Taxa_Names_dat, ParentChild_gz[,'ChildName'] )
  g_i_i <- sapply(FUN = find_index, X = g_i_dat, y = g_i)
  
  # Create index of species to Parent  - Child
  spc_in_PC_gz <- which(PC_gz[,2] == max(PC_gz[,2]))
  
  
  
  # Setup TMB ####
  data <- list(PC_gz = PC_gz,
               g_i = g_i - 1,
               logL = log(all.dat$Std_length_mm),
               logW = log(all.dat$total_weight_g),
               taxa_id = g_i_i -1,
               spc_in_PCgz = spc_in_PC_gz -1
               
  )
  
  parameters = list(alpha_j = rep(0,n_j),
                    L_z = rep(1, 3),
                    log_lambda = rep(0, length(unique(PC_gz[,2])) -1),
                    beta_gj = matrix(0, nrow = n_g, ncol = n_j),
                    logsigma = 0
  )
  Random <- c("beta_gj")
  
  model <- "hierarchical_LW"
  
  compile(paste0("code/length_weight_reg_code/TMB/", model, ".cpp"))
  
  dyn.load(dynlib(paste0("code/length_weight_reg_code/TMB/",model)))
  
  ### Run TMB ####
  obj <-
    MakeADFun(
      data = data,
      parameters = parameters,
      DLL = model,
      random = Random,
      silent = TRUE
    )
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  rep = sdreport( obj,
                  getReportCovariance = TRUE, 
                  getJointPrecision=TRUE)
  re <- summary(rep, "random")
  fixef <- summary(rep, "fixed")
  beta_mle <- matrix(re[grep(rownames(re), pattern = "beta"),1], nrow = n_g, ncol = 2, byrow = F)
  beta_se <- matrix(re[grep(rownames(re), pattern = "beta"),2], nrow = n_g, ncol = 2, byrow = F)
  
  ### Save output in a tibble #
  # get group names
  Groups <- gsub(".*_","",ParentChild_gz$ChildName)
  
  # make output dataframe
  output <- tibble(Group = Groups,
                   loga = beta_mle[1:length(Groups)],
                   b = beta_mle[(length(Groups) + 1) : length(beta_mle)]
  )
  
  # save output
  saveRDS(output, "analysis/length_weight_parameters.RDS")
  return(output)
}



# Function to create array of fish lengths and associated probabilities from catch data ####
make_Pl <- function(larray) { 
  # Function accepts array of lengths (x), returns a list with 100 lengths (L), 
  # and estimated proportion of fish in each length bin (P)
  
  # constrain range of lengths 
  range.x  = c(5, 1.25*max(larray)) # we only consider fish > 10 mm, so we make
  # the min 5 here, and we allow the max to be up to 25% higher than the max length of fish.
  # this is done to ensure that the length values can't go negative. 
  
  smooth <- bkde(x = larray, gridsize = 101, range.x = range.x)
  
  deltal <- smooth$x[2] - smooth$x[1]
  # get midpoint lengths
  midlens <- smooth$x[-101] + deltal * 0.5
  # get p(x) for the midpoints
  middens <- approx(x = smooth$x, y = smooth$y, xout = midlens)$y
  # get P(x) for the midpoints (denominator fixes small approximation errors to ensure that probabilities sum to one)
  midprobs <- middens * deltal / sum(middens * deltal)
  output <- list(L = midlens, P = midprobs)
  return(output)
}

# Function to calculate carbon transport, given length and species group ####
make_tau <- function(L, taxagroup, pars, VM = TRUE, lwreg) {
  
  loga <- lwreg$loga[lwreg$Group == taxagroup]
  b <- lwreg$b[lwreg$Group == taxagroup]
  W <- exp(loga) * L ^ b
  list2env(pars, env = environment())
  # load all parameters in pars into function environment
  #parnames <- names(pars)
  #npars <- length(parnames)
  ##for (i in 1:npars) {
  #  eval(parse(text = paste0(parnames[i],"<-",pars[[i]])))
  #}
  
  if (VM) tau <- C.flux.function.VM(W_w_f_VM = W, Wref, Tref, a0, a1, a2, a3, H_VM_or_NM = 1, TT_epipelagic, TT_mesopelagic, TT_mean, activity_factor_resting, activity_factor_migrating, activity_factor_foraging, time_resting, Dmax, Dmin, B, S, prop_SDA, prop_egestion, prop_excretion, G_VM, theta, peakTime, timeTo90, Q_ox_kJ_g_O2, RQ, energy_in_pellet, C_in_pellet,  proportion_egestion_below_150m, X_r, X_e, mu, C_c_vm, prop_mortality_below_150m)
  if (!VM) tau <- C.flux.function.NM(W_w_f_NM = W, Wref, Tref, a0=a0, a1=a1, a2=a2, a3=a3, H_VM_or_NM=0, TT_mesopelagic=TT_mesopelagic, activity_factor_resting=activity_factor_resting, activity_factor_foraging=activity_factor_foraging, time_resting=time_resting, prop_SDA=prop_SDA, prop_egestion=prop_egestion, prop_excretion=prop_excretion, G_NM=G_NM, theta=theta, X_r=X_r, mu=mu, prop_non_detrital_NM_prey=prop_non_detrital_NM_prey, Q_ox_kJ_g_O2=Q_ox_kJ_g_O2, RQ=RQ, energy_in_pellet=energy_in_pellet, C_in_pellet=C_in_pellet, C_c_nm=C_c_nm)
  
  return(tau)
  
}


#### Function to generate random draws for catch ####

make_random_catch <- function(fishdata, volumedata, towlist = c(2, 4, 6, 8, 9, 11, 48, 90), group) {
  
  ntows = length(towlist)
  C_t_bar_sd <- calc_C_t_bar(fishdata = fishdata, volumedata = volumedata)
  
  # get C_t_bar for group
  C_t_bar <- C_t_bar_sd$c_t_bar
  C_t_sd <- C_t_bar_sd$c_t_sd
  
  # get mean for group
  groupnames <- names(C_t_bar)
  
  C_g_bar = C_t_bar[[which(groupnames == group)]]
  C_g_sd = C_t_sd[[which(groupnames == group)]]
  
  # generate 5 random draws
  c_g_t_random <- rnorm(ntows, mean = C_g_bar, sd = C_g_sd)
  
  # if any negative draws, set to zero
  c_g_t_random[c_g_t_random<0] <- 0
  
  return(c_g_t_random)
}



#### Function to generate VM pars for use in Monte Carlo function ####

get_pars <- function(nsims){
  
  # use Ikeda 2016 values for regression of both VM and NM fish 
  # made range the mean +/- 2 SE 
  
  # Variance-covariance matrix, for limiting MC draws for intercept and slope to only appropriate values based on fit to data
  
  # set seed for reproducibility
  set.seed(123)
  ikeda_data_fit <- get_ikeda_coefs()
  resp_model_coeff_center <- ikeda_data_fit$resp_model_coeff_center
  Sigma_resp_model_center <- ikeda_data_fit$Sigma_resp_model_center
  Wref <- 1000
  Tref <- 15
  
  # generate twice as many values as needed for each coefficient (so that, after filtering,
  # there will still be sufficient values for use in the Monte Carlo simulation)
  # use centered regression
  coeff_table_unfiltered_c <- mvrnorm(n = nsims*2, mu = resp_model_coeff_center, Sigma = Sigma_resp_model_center) # _c = centered regression
  
  # coerce into data frame 
  coeff_table_unfiltered_df_c <- as.data.frame(coeff_table_unfiltered_c)
  
  # edit column names
  colnames(coeff_table_unfiltered_df_c) <- c("Intercept", "log.Wref", "TempRef", "Depth.VM")
  
  # pull out just the standard errors
  std_errors <- sqrt(diag(Sigma_resp_model_center))
  
  ### code to add respiration rate coefficient thresholds ###
  
  # We need to filter out a few extreme values that create unrealistic results (to avoid unrealistic results). 
  # To do this, we first centered the regression on reference temp and fish mass values, which removes much of 
  # the covariance that otherwise occurs between the intercept and slope of the regression. Next, we constrain 
  # the variance for each parameter, treating each parameter as univariate (rather than multivariate, 
  # which we were doing before centering the regression due to high covariance). Below this is done by 
  # constraining each predictor within +/- 2 standard errors. 
  
  # filter out any values greater than 2 SE from the mean coefficient value (removes ~5% of values)
  all_intercepts <- coeff_table_unfiltered_df_c$Intercept[
    -which(abs(coeff_table_unfiltered_df_c$Intercept-resp_model_coeff_center[1]) > std_errors[1]*2)]
  all_mass_coeff <- coeff_table_unfiltered_df_c$log.Wref[
    -which(abs(coeff_table_unfiltered_df_c$log.Wref-resp_model_coeff_center[2]) > std_errors[2]*2)]
  all_temp_coeff <- coeff_table_unfiltered_df_c$TempRef[
    -which(abs(coeff_table_unfiltered_df_c$TempRef-resp_model_coeff_center[3]) > std_errors[3]*2)]
  all_depth_coeff <- coeff_table_unfiltered_df_c$Depth.VM[
    -which(abs(coeff_table_unfiltered_df_c$Depth.VM-resp_model_coeff_center[4]) > std_errors[4]*2)]
  # samples 1000 values from this vector of filtered values (all with equal chance of being
  # selected)
  # shouldn't matter much whether replace = TRUE or FALSE, if there are sufficient values remaining to sample when replace = FALSE
  a0.vec <- sample(x = all_intercepts, size = nsims, replace = FALSE)
  a1.vec <- sample(x = all_mass_coeff, size = nsims, replace = FALSE)
  a2.vec <- sample(x = all_temp_coeff, size = nsims, replace = FALSE)
  a3.vec <- sample(x = all_depth_coeff, size = nsims, replace = FALSE)
  
  # check what % of values were filtered out
  100-(length(all_intercepts)/length(coeff_table_unfiltered_df_c$Intercept)*100)
  # approximately 5% filtered out, as we would expect constraining +/- 2 SE
  
  ### end of code to add respiration rate coefficient thresholds ###
  H_VM.vec <- runif(nsims, min = 1, max = 1) # use 1 for VM fish or 0 for NM fish function 
  TT_epipelagic.vec <- runif(n = nsims, min = 14-0.5, max = 14+0.5) # +/- 1 degrees from measured epipelagic temp in EXPORTS location. Adjust based on environment modeled
  TT_mesopelagic.vec <- runif(n = nsims, min = 11-0.5, max = 11+0.5) # +/- 1 degrees from measured mesopelagic temp in EXPORTS location. Adjust based on environment modeled
  TT_mean.vec <- runif(n = nsims, min = 12.5-0.5, max = 12.5+0.5) # +/- 1 degrees from measured mean(mean(epipelagic), mean(mesopelagic)) temp in EXPORTS location. Adjust as needed
  
  activity_factor_resting.vec <- runif(nsims, min = 0.5, max = 1) # assumes fish is at or just 2x above SMR during daytime in deep water (i.e., fish is either 1/2 of routine metabolic rate RMR or at RMR)
  activity_factor_foraging.vec <- runif(nsims, min = 1, max = 4) # Assumes fish metabolizes at least at 2x SMR (if SMR is 1/2 RMR) and as high as 8x SMR during nighttime foraging
  activity_factor_migrating.vec <- runif(nsims, min = 1, max = 4) # Assumes fish metabolizes at least at 2x SMR (if SMR is 1/2 RMR) and as high as 8x SMR during during migration (relevant for VM fish only)
  time_resting.vec <- runif(nsims, min = .44, max = .47) # found empirically from acoustics, see code earlier in script. Note that this could vary at other time of the year. 
  G_VM.vec <- runif(nsims, min = 0.018-0.018*0.5, max = 0.018+0.018*0.5)
  # 0.0018 for VM fish, 0.0007 for NM fish
  
  theta.vec <- runif(nsims, min = 0.54, max = 1)
  
  # First, assume that oxygen consumption during respirometry experiments is due only to M, not also to SDA (food assimilation). Realistically some of the fish from Ikeda 2016 regression were probably still digesting food, so we will also consider a case where the measured O2 updated is due to both M and SDA combined (respiration rate, in units of energy usage, = M + SDA). 
  
  prop_SDA.vec <- runif(nsims, min = 0.14 * 0.8, max = 0.14 * 1.2) # Nominal value from Brett and Groves (they used 0.14), allowed to vary by 20%
  prop_egestion.vec <- runif(nsims, min = 0.20 * 0.8, max = 0.20 * 1.2) # Nominal value from Brett and Groves (0.20), allowed to vary by 20%
  prop_excretion.vec <- runif(nsims, min = 0.07 * 0.8, max = 0.07 * 1.2) # Nominal value from Brett and Groves (0.07), allowed to vary by 20%
  
  # for SDA function
  timeTo90.vec <- runif(nsims, min = 10, max = 14)
  peakTime.vec <- runif(nsims, min = 3- 3*.8, max = 3+ 3*.8)
  
  Q_ox_kJ_g_O2.vec <- runif(nsims, min = 13.4- (13.4*0.1), max = 13.4 + (13.4*0.1)) # From Jobling 1994. Oxycalorific equivalent (energy obtained from respired oxygen) in units of kJ/gO2. Nelson et al. 2016 states that this may vary by 10% at most (sd of about 1).
  
  RQ.vec <- runif(nsims, min = 0.7, max = 2) # Hernandez-Leon et al. 2019 used respiratory quotient (CO2 respired/O2 consumed) of 0.97 (from Omori and Ikeda, 1984). Davison et al. 2013 and Ariza et al. 2015 used respiratory quotient of 0.90 (Brett and Groves 1979). Hudson let this range from 0.7-1.0, given uncertainty associated with mesopelagic fish. According to Brett and Groves 1979, RQ could = 2 for fish that undergo some amount of anaerobic metabolism.
  
  # Dmax range is from EXPORTS acoustics echogram from Mei Sato (and happens to 
  # be the same as the assumptions made for chapter 1... migrators appear to be 
  # swimming as deep as 700 m during the day, or as shallow as about 300 m during the day)
  # We'll extend to 750 m to match our net depth intervals. 
  Dmax.vec <- runif(nsims, min = 300, max = 500)
  
  # range used for chapter 1 of 10-90 meters for migrators at night seems reasonable too. 
  # We'll extend to 100 m to match our net depth intervals 
  Dmin.vec <- runif(nsims, min = 10, max = 100) 
  
  B.vec <- runif(nsims, min = 199, max = 201) # set this to 200 m, but adjust when 
  # finding the flux at 500 m flux boundary (see below)
  
  B_500.vec <- runif(nsims, min = 499, max = 501) # 500 m flux boundary
  
  S.vec <- runif(nsims, min = 2, max = 5.9) # # from acoustic data described earlier in this script 
  # note that in the code we make swim speed in cm/s and then multiply by 0.01 (unit used in references), 
  # but in the equation shown in McMonagle et al. 2023 and in Table B1 (to remove the 0.01 unit conversion factor)
  # we provide swim speed in m/s. 
  
  energy_in_pellet.vec <- runif(nsims, min = 0.6, max = 3.5)
  
  C_in_pellet.vec <- runif(nsims, min = 0.05, max = 0.10)
  
  proportion_egestion_below_150m.vec <- runif(nsims, min = 0.50, max = 0.90) # 50-90% sinks below 150m

  X_r.vec <- runif(nsims, min = 0.008 - 0.8*0.008, max = 0.008 + 0.8*0.008) # see Table 1--these values came from Wilson et al. 2009 and then were converted to mgC /d/ WMfish. Allow to vary by 80%
  
  X_e.vec <- runif(nsims, min = 0.5, max = 0.9)
  
  # For Pauly 1980 mortality regression, simply allow mu to vary by +/- 0.3
  # log_{10}(\mu) = a_0 + a_1 log_{10}(L_{\infty}) + a_2log_{10}(K) + a_3log_{10}(T_{mean}) 
  aa0 <- -0.0066
  aa1 <- -0.2790
  aa2 <- 0.6543
  aa3 <- 0.4634
  Linfinity <- 8.7 # average of the 3 myctophid Linfinity values in Pauly 1980
  K <- 0.2 # arbitrary value between 0.3 and 0.4 for VM fish (Pauly 1980) 
  # and 0.05 and 0.17 for NM fish (Childress et al 1980)
  # use same parameter values for NM fish for now due to limited data
  TT_mean_nominal_value <- 12.5 # set TT_mean
  mu.mean <- 10^ (aa0 + aa1*log10(Linfinity) + aa2*log10(K) + aa3*log10(TT_mean_nominal_value))
  mu.vec <- runif(nsims, min = mu.mean- mu.mean*.5, max = mu.mean+ mu.mean*.5) 
  C_c_vm.vec <- runif(nsims, min = 0.038, max = 0.248) # Ratio of carbon:wet weight of VM fish, from mean of species captured 
  # in this study, and Childress and Nygaard (1973). In Davison et al. 2013--they used 0.49 for carbon:wet weight ratio, but in Childress and Nygaard paper this was about 0.50 for carbon:DRY weight. Mesopelagic fish ranged from about 60% to 90% water, and carbon:ash free dry weight (assumed here to be the same as dry weight, because % ash is small) ranged from 38-62%. To find bounds: min carbon:wet weight = (1-0.90) * 0.38 = 0.038, and max carbon:wet weight = (1-0.60) * 0.62 = 0.248
  prop_mortality_below_150m.vec <- runif(nsims, min = .40, max = .90)
  
  # NM parameters
  H_NM.vec <- runif(nsims, min = 0, max = 0) # use 1 for VM fish or 0 for NM fish function 
  G_NM.vec <- runif(nsims, min = 0.007-0.007*0.5, max = 0.007+0.007*0.5) 
  C_c_nm.vec <- runif(nsims, min = 0.038, max = 0.248) # Ratio of carbon:wet weight of VM fish, from mean of species captured  
  prop_non_detrital_NM_prey.vec <- runif(nsims, min = 0.6-0.6*.5, max = 0.6+0.6*.5) # based on 
  # Table B1 from previous paper (chapter 1)
  
  
  # make a table with parameter values for each simulation
  pars <- as.data.frame(cbind(Wref, Tref, a0.vec, a1.vec, a2.vec, a3.vec, H_VM.vec, TT_epipelagic.vec, TT_mesopelagic.vec, TT_mean.vec, activity_factor_resting.vec, activity_factor_foraging.vec, activity_factor_migrating.vec, time_resting.vec, Dmax.vec, Dmin.vec, B.vec, S.vec, prop_SDA.vec, prop_egestion.vec, prop_excretion.vec, G_VM.vec, theta.vec, peakTime.vec, timeTo90.vec, Q_ox_kJ_g_O2.vec, RQ.vec, energy_in_pellet.vec, C_in_pellet.vec, proportion_egestion_below_150m.vec, X_r.vec, X_e.vec, mu.vec, C_c_vm.vec, prop_mortality_below_150m.vec, H_NM.vec, G_NM.vec, C_c_nm.vec, prop_non_detrital_NM_prey.vec))
  colnames(pars) <- c("Wref", "Tref", "a0", "a1", "a2", "a3", "H_VM", "TT_epipelagic", "TT_mesopelagic", "TT_mean", "activity_factor_resting", "activity_factor_foraging", "activity_factor_migrating", "time_resting", "Dmax", "Dmin", "B", "S", "prop_SDA", "prop_egestion", "prop_excretion", "G_VM", "theta", "peakTime", "timeTo90", "Q_ox_kJ_g_O2", "RQ", "energy_in_pellet", "C_in_pellet",  "proportion_egestion_below_150m", "X_r", "X_e", "mu", "C_c_vm", "prop_mortality_below_150m", "H_NM", "G_NM", "C_c_nm", "prop_non_detrital_NM_prey")
  
  return(pars)
}



  
### Functions to calculate C flux for VM fish ####
# with weight in g as the input parameter to be passed to "W_w_f"

# First, below are function to use inside C.flux.function.VM, 
# for finding SDA flux

# Main function. pSDA is the proportion of exported SDA, called SDAe in text of manuscript. 
calc.pSDA <- function(pars) {
  
  # temporarily turn off warnings that are annoying
  defaultW <- getOption("warn")
  options(warn = -1)
  
  peakTime <- pars$peakTime
  timeTo90 <- pars$timeTo90
  tmax <- pars$tmax
  tau <- pars$tau
  estGamma <- function(peakTime, timeTo90) {
    setOptim <- function(logalpha, peakTime, timeTo90) {
      alpha <- exp(logalpha) + 1 # adding 1 ensures that peak is not at 0
      beta <- (alpha - 1) / peakTime
      estTimeTo90 <-
        qgamma(0.90,
               shape = alpha,
               rate = beta,
               lower.tail = TRUE)
      return((timeTo90  - estTimeTo90) ^ 2)
    }
    
    #### estimate alpha and beta based on peakTime and timeTo90 
    log.alpha.start <- 5
    fit <- optim(
      par = log.alpha.start,
      fn = setOptim,
      method = "Brent",
      lower = -200,
      upper = 200,
      peakTime = peakTime,
      timeTo90 = timeTo90
    )
    alpha <- exp(fit$par) + 1
    beta <- (alpha - 1) / peakTime
    return (c(alpha = alpha, beta = beta))
  }
  
  #### function to return F(tmeal) #
  foft <- function(tmeal, alpha, beta, tmax, tau) {
    if (tmeal == 0)
      return(dunif(tmeal, min = 0, max = tmax))
    if (tmeal > 0)
      return(dunif(tmeal, min = 0, max = tmax) * (
        1 - pgamma(
          tmax + tau - tmeal,
          shape = alpha,
          rate = beta,
          lower.tail = TRUE
        )
      ))
  }
  ### Use functions to calculate integral
  
  gammaPars <- estGamma(peakTime, timeTo90)
  alpha <- gammaPars["alpha"]
  beta <- gammaPars["beta"]
  
  pSDA <-
    tmax / 6 * (foft(tmeal = 0, alpha, beta, tmax, tau) + 
                  4 * foft(tmeal = tmax / 2, alpha, beta, tmax, tau) + 
                  foft(tmeal =tmax, alpha, beta, tmax, tau)
    )
  # turn warnings back to default
  options(warn = defaultW)
  return(pSDA)
}



#=================VM fish flux function==============================
C.flux.function.VM <- function(W_w_f_VM, Wref, Tref, a0, a1, a2, a3, H_VM_or_NM, TT_epipelagic, TT_mesopelagic, TT_mean, activity_factor_resting, activity_factor_migrating, activity_factor_foraging, time_resting, Dmax, Dmin, B, S, prop_SDA, prop_egestion, prop_excretion, G_VM, theta, peakTime, timeTo90, Q_ox_kJ_g_O2, RQ, energy_in_pellet, C_in_pellet,  proportion_egestion_below_150m, X_r, X_e, mu, C_c_vm, prop_mortality_below_150m) {
  
  # Notes on Tref and Wref arguments above:
  # Tref <- 15 # Reference temperature 
  # Wref <- 1000 # Reference fish wet mass 
  # Pass these values in the MC function that comes later. These reference values
  # are used to center the regression to avoid as much covariance between the
  # intercept and slopes. 
  
  # Equations 1 and 2: oxygen utilization
  
  # center regression on reference fish weight (Wref) and reference temperature (Tref), done to avoid strong covariance 
  # between intercept and slope
  R_O_VM_r <- exp(a0 + a1* (log(W_w_f_VM*1000) - log(Wref)) + a2 * (1000/(TT_mesopelagic+273.15) - 1000/(Tref + 273.15)) + a3 * H_VM_or_NM) # ul O2/ind / hour. Raised both sides by e to remove the ln() before respiration rate on the left side of equation
  
  R_O_VM_m <- exp(a0 + a1 * (log(W_w_f_VM*1000) - log(Wref)) + a2 * (1000/(TT_mean+273.15) - 1000 /(Tref + 273.15)) + a3 * H_VM_or_NM)
  
  R_O_VM_f <- exp(a0 + a1 * (log(W_w_f_VM*1000) - log(Wref)) + a2 * (1000/(TT_epipelagic+273.15) - 1000 / (Tref + 273.15)) + a3 * H_VM_or_NM) 
  
  # Equations 3 and 4: energy expenditure due to non-SDA metabolism for each daily activity cycle.
  # Convert from units of O2 to energy, and from per hour to per day, for each part of daily cycle
  
  # find time (proportion of a day) spent in each behavior
  # time resting will be 0.3475 nominally but allowed to vary based on acoustic data.
  # note that this should be increased if finding expected average annual fish flux, 
  # because the mean number of daylight hours will be lower than during the May spring bloom
  # in the Northern Hemisphere where the EXPORTS cruise took place 
  t_m <- 2 * (Dmax - Dmin) / ((S/100)*(60*60*24)) # S= about 5.9 cm/s, or 0.059 m/s)
  t_f <- 1 - (t_m + time_resting) 
  
  # Calculate energy expenditure for each behavior in a daily cycle (resting, migrating, foraging)
  R_VM_r <- 24 * .000001 * 1.429 * R_O_VM_r * Q_ox_kJ_g_O2 * time_resting *
    activity_factor_resting * theta
  R_VM_m <- 24 * .000001 * 1.429 * R_O_VM_m * Q_ox_kJ_g_O2 * t_m *
    activity_factor_migrating * theta
  R_VM_f <- 24 * .000001 * 1.429 * R_O_VM_f * Q_ox_kJ_g_O2 * t_f *
    activity_factor_foraging * theta
  
  # Units are now in kj/day, converted by multiplying resp rate for 1 individual * 24 hr/day * 1 ml / 1000 ul * 1 L / 1000 ml * Qox in kcal / L O2 * 4.184 kJ / kcal)
  
  # Sum all 3 activities to get daily energy expended on non-SDA metabolism per day
  R_VM_daily <- R_VM_r + R_VM_m + R_VM_f
  
  # Equation 5: energy ingested (kJ/d)
  I_VM <- (G_VM + R_VM_daily) / (1 - (prop_SDA + prop_egestion + prop_excretion))
  
  # Equation 6 and 7: C flux from respiration (note; equation numbers differ in manuscript due to moving some material to supplement)
  
  # Calculate proportion of respired carbon that is exported
  R_exported_VM = R_VM_m * (1- ((B-Dmin) / (Dmax - Dmin))) + R_VM_r 
  
  # Convert from units of energy, to units of oxygen, to units of carbon
  C_r_VM = 12 * RQ * 1000 * R_exported_VM / (32 * Q_ox_kJ_g_O2)
  
  # units: kj/day respired and exported * (1 g O2 / 13.4 kJ from Jobling 1994) * (mol O2 / 32 g O2) * respiratory quotient * (1 mol C / 1 mol CO2) * (12 g C / mol C) * 1000 to convert g C to mg C
  
  # The constant 32 converts from grams to moles of oxygen, and there is a 1:1 conversion between moles of carbon dioxide and carbon. The constant 12 converts from moles of carbon to grams of carbon, and 1000 converts from grams of carbon to milligrams (W_w_f_VM is entered earlier as 1 for a 1 g fish)
  
  # Equation 8.5a  C flux from SDA
  # Some parts of these calculations are explained in greater detail in supplementary info. 
  
  # time between when fish starts and ends foraging
  tmax <- 24*t_f
  
  # time after feeding ends that fish crosses flux boundary
  tau <- (B-Dmin) / (Dmax-Dmin) * 24 * t_m * .5
  
  # make object containing parameters for SDA flux calculations
  pars <- list(peakTime = peakTime,
               timeTo90 = timeTo90,
               tmax = tmax,
               tau = tau) 
  
  # find proportion of respired energy lost to SDA that is exported
  SDA_exported <- calc.pSDA(pars)
  
  # convert to C flux per day lost to SDA below flux boundary
  C_SDA_VM <- (12 * RQ * 1000 * prop_SDA * SDA_exported * I_VM) / 
    (32 * Q_ox_kJ_g_O2)
  # 1000 converts from g to mg (no need to multiply by W_w_f because )
  
  # Equation 8.5b: SDA flux for NM fish (NM fish function is a separate function).
  
  
  # Equation 9: egestion flux
  
  ## Assume that 75-100% of carbon released as fecal pellets ends up below flux boundary. Fecal pellets are thought to be released during daytime depth, due to slow rate of gut evacuation (slow compared to zooplankton), and even if released above flux boundary, presumably most would sink below flux boundary before remineralization. Minimum flux here assumes that 25% gets released above the flux boundary and is remineralized before reaching deeper than flux boundary.
  
  # Amount of carbon released below flux boundary also depends on carbon content of total energy lost from egestion. Use published caloric value of fish fecal pellets: 0.6-3.5 kcal/g in Bailey and Robertson 1982, or 2.7 kcal/g mentioned in Brett and Groves 1979. convert kj --> kcal --> g --> POC content
  # from Saba and Steinberg 2012: pellet C content ranged from 15 ug C / pellet to 31 ug C / pellet, and avg pellet size was 314 ug (dry weight). So, range is 15 ug C / 314 ug pellet DW to 31 ug C / 314 ug pellet DW, or 0.05-0.1 g C / g DW pellet.
  
  C_f_VM <- proportion_egestion_below_150m * I_VM* (1/4.184) * (1/energy_in_pellet) * C_in_pellet * proportion_egestion_below_150m * 1000 # egested energy F_VM, in kj/day/ind, * (kcal/ 4.184 kJ) * (g DW /0.6-3.5 kcal) (0.05-0.10 g C / g DW pellet), then 1000 converts from g to mg 
  # In sensitivity analyses, we allow energy_in_pellet to vary from 0.6-3.5 kcal/ g DW. Make C_p (carbon in pellet) vary from 0.05-0.10 g C/g DW. 
  
  # Equation 10 finds egestion flux for NM fish.
  
  # Equation 10.1: excretion flux
  C_X_VM <- X_r * X_e # nominal value of 0.008 is in mgC per day, so no need to multiply this flux by 1000 to convert from g to mg. 
  
  # Equation 10.2 finds excretion flux for NM fish
  
  # Equation 11: mortality flux
  
  # Pauly 1980 mortality rate regression:
  # set mu as a constant based on aa0, aa1, aa2, aa3, K, and Linfinity from lit but 
  # only allow mu to vary in sensitivity analysis 
  # This allows us to determine overall sensitivity of C flux to mu
  # mu <- exp(aa0 - aa1*log(Linfinity) + aa2*log(K) + aa3*log(TT_mean))
  # kJ/d * g C / g wet mass of fish * g wet mass of fish / kJ * 1000 mg / g
  # mg C /d = 1 g fish * g C/g fish * rate in units per day * 1000 mg / g
  C_mu_VM <- prop_mortality_below_150m * C_c_vm * (1 - exp(-mu/365)) * 1000 * W_w_f_VM
  
  # # in paper, we use C_c units of mg C / g fish wet mass, to avoid needing the 1000 constant in the equation, but in the code we just multiply by 1000 and keep this is mg C / mg fish wet mass
  
  # Final equation 
  # Not a numbered equation because easily described in text
  
  # Add C_r, C_f, C_m and C_x together 
  C_vm_tot <- C_r_VM + C_SDA_VM + C_f_VM + C_X_VM + C_mu_VM # mg C / ind / day (so a 1g fish is on avg transporting something like 1 mg C per day)
  
  # Note that this value can be treated as a rough average--it does not account for seasonality. If desired, could multiply by 365 to get an approximate annual carbon flux rate. 
  
  return(c(C_vm_tot, C_r_VM, C_SDA_VM, C_f_VM, C_X_VM, C_mu_VM, I_VM))
}


#### Functions for finding C flux of NM fish ####

#=================NM fish flux function==============================

C.flux.function.NM <- function(W_w_f_NM, Wref, Tref, a0, a1, a2, a3, H_VM_or_NM, TT_mesopelagic, TT_mean, activity_factor_resting, activity_factor_foraging, time_resting, prop_SDA, prop_egestion, prop_excretion, G_NM, theta, X_r, mu, Q_ox_kJ_g_O2, RQ, prop_non_detrital_NM_prey, energy_in_pellet, C_in_pellet, C_c_nm) {
  
  # Equations 1 and 2
  # Find oxygen utilization, which is dependent on mass of individual, temp of water and habitat type). Use Ikeda 2016 (in ul O2/ind/hr) for all fish, but with modifications to exclude mesopelagic fish and to account for NM or VM fish in the habitat depth predictor
  
  # write out in log space then exponentiate. pass Wref and Tref earlier so they work in this function. do same for VM and NM fish 
  
  R_O_NM_r <- exp(a0 + a1* (log(W_w_f_NM*1000) - log(Wref)) + a2 * (1000/(TT_mesopelagic+273.15) - 1000/(Tref + 273.15)) + a3 * H_VM_or_NM) # ul O2/ind / hour. Raised both sides by e to remove the ln() before respiration rate on the left side of equation
  
  R_O_NM_f <- R_O_NM_r # baseline respiration rate (before accounting for difference in activity level) is the same during foraging and resting for NM fish (unlike for VM fish where the temperature used is different during resting, foraging and migrating)
  
  
  # Equation 3 and 4: convert from units of O2 to energy, and from per hour to per day, for each part of a daily cycle
  
  # find time (proportion of a day) spent in each behavior
  # time resting will be 0.5 nominally but allowed to vary in sensitivity analysis
  t_f <- 1 - time_resting
  
  # Calculate energy expenditure for each behavior in a daily cycle (resting, migrating, foraging)
  R_NM_r <- 24 * .000001 * 1.429 * R_O_NM_r * Q_ox_kJ_g_O2 * time_resting *
    activity_factor_resting * theta
  R_NM_f <- 24 * .000001 * 1.429 * R_O_NM_f * Q_ox_kJ_g_O2 * t_f * 
    activity_factor_foraging * theta
  
  # Sum all 2 activities to get daily energy expended on metabolism
  R_NM_daily <- R_NM_r + R_NM_f # total energy expenditure per day
  
  # kj/day, converted by multiplying resp rate for 1 individual * 24 hr/day * 1 ml / 1000 ul * 1 L / 1000 ml * Qox in kcal / L O2 * 4.184 kJ / kcal)
  
  
  # Equation 5: energy ingested (kJ/d)
  I_NM <- (G_NM + R_NM_daily) / (1 - (prop_SDA + prop_egestion + prop_excretion))
  
  
  # Equation 8
  # Find respiratory flux of NM fish
  
  ## amount of carbon released: convert from kj/day to O2 util/day to C released /day for each individual
  
  C_r_NM = R_NM_daily * prop_non_detrital_NM_prey * (1/Q_ox_kJ_g_O2) * (1/32) * RQ * 12 * 1000
  
  # kj/day respired and exported * (1 g O2 / 13.4 kJ from Jobling 1994) * (mol O2 / 32 g O2) * respiratory quotient * (1 mol CO2 / 1 mol O2) (1 mol C / 1 mol CO2) * (12 g C / mol C) * 1000 to convert g C to mg C
  
  # Equation 8.5b
  # Find C flux from SDA
  C_SDA_NM <- 12 * 1000 * RQ * prop_SDA * I_NM * prop_non_detrital_NM_prey * (1/32) * (1/Q_ox_kJ_g_O2)
  
  # Equation 10
  # Find C flux from egestion
  # For NM fish, assume 100% of carbon egested is released > 150 m, but that amount is reduced by multiplying by Z (proportion of zooplankton of non-detrital prey)
  
  C_f_NM <- 1000 * prop_egestion * C_in_pellet * prop_non_detrital_NM_prey * I_NM * (1/4.184) * (1/energy_in_pellet)
  
  # energy lost to ingestion in kj/day/ind * (kcal/ 4.184 kJ) * (g DW /0.6-3.5 kcal) (0.05-0.10 g C / g DW pellet)
  
  
  # Equation 10.2 finds excretion flux for NM fish
  
  C_X_NM <- X_r * prop_non_detrital_NM_prey # nominal value of 0.008 is in mgC per day, so no need to multiply this flux by 1000 to convert from g to mg. 
  
  # Equation 12
  # Find C flux from mortality
  C_mu_NM <- prop_non_detrital_NM_prey * C_c_nm * (1 - exp(-mu/365)) * W_w_f_NM * 1000 
  
  # Assumes 100% of NM fish mortality (predation too) occurs > 150 m. 
  
  # Final equation 
  # Not a numbered equation because easily described in text
  
  # Add C_r, C_SDA, C_f, C_X, and C_mu together
  C_nm_tot <- C_r_NM + C_SDA_NM + C_f_NM + C_X_NM + C_mu_NM # mg C / ind / day (so a 1g fish is on avg transporting something like 2 mg C per day)
  
  # energy ingested (for checking estimate makes sense) = I_NM
  
  # Note that this value can be treated as a rough average--it does not account for seasonality. If desired, could multiply by 365 to get an approximate annual carbon flux rate. 
  
  return(c(C_nm_tot, C_r_NM, C_SDA_NM, C_f_NM, C_X_NM, C_mu_NM, I_NM))
}

#### Function to run Monte Carlo simulation for VM fish ####
# This code was used to run the Monte Carlo in McMonagle et al. 2023 "High uncertainty 
# in estimates of fish-mediated carbon flux into the ocean's twilight zone", but we use
# a different code for this manuscript as found in 02-MonteCarlo_C_flux.R. 
# 
# ## Monte Carlo simulation for VM fish carbon flux
# ## Note that this returns two data frames within a list. The first data frame contains the total C flux and the parameter vectors for each simulation, and a second df contains the total C flux as well as the C flux associated with each carbon pathway (resp, SDA, egestion, etc.) for each simulation.
# # Enter W_w_f in g and Wref in mg
# MC_func_VM <- function(nsims, W_w_f_VM, Tref, Wref) {
#   output <- rep(NA, nsims)
#   W_w_f_VM <- W_w_f_VM # set this in argument of MC_func_VM
#   
#   # use Ikeda 2016 values for regression of both VM and NM fish 
#   # made range the mean +/- 2 SE 
#   
#   # Variance-covariance matrix, for limiting MC draws for intercept and slope to only appropriate values based on fit to data
#   
#   # set seed for reproducibility
#   set.seed(123)
#   
#   # load data 
#   resp_data <- read.csv(here("data/Individual_fish_flux_data",
#                              "Ikeda_2016_data_VM_NM.csv"), header = TRUE)
#   
#   # make the categorical depth category a factor
#   resp_data$Cat_depth <- as.factor(resp_data$Cat_depth)
#   
#   # re-calculate temp data as in Ikeda 2016 (take inverse, multiple by 1000, and convert to K)
#   resp_data$Temp_var <- 1000/(resp_data$Temp_C + 273.15) # original temperature predictor
#   resp_data$Inv_Temp_var <- 1000/(resp_data$Temp_C + 273.15) - 1000/(Tref+ 273.15) # modification for centering at Tref
#   resp_data$logW_ref <- log(resp_data$mgWM) - log(Wref) # modification for centering at Wref
#   # Note: only need to center temp and fish mass because these are the only two continuous predictors (habitat depth is categorical)
#   
#   # fit linear regression (with log-log transformation for respiration rate and fish wet mass)
#   resp_model <- lm(log(R) ~ log(mgWM) + Temp_var + Cat_depth, resp_data)
#   
#   # same model but with centered temp and mass predictors
#   resp_model_center <- lm(log(R) ~ logW_ref + Inv_Temp_var + Cat_depth, resp_data)
#   
#   # regression summaries for old model and new, centered model
#   summary(resp_model)
#   summary(resp_model_center) # use this model (centered, so less covariance between intercept and coefficients of predictors)
#   
#   # variance-covariance plot
#   Sigma_resp_model_old <- vcov(resp_model)
#   Sigma_resp_model_center <- vcov(resp_model_center)
#   
#   # best fit model coefficients 
#   resp_model_coeff <- coefficients(resp_model)
#   resp_model_coeff_center <- coefficients(resp_model_center)
#   
#   # generate twice as many values as needed for each coefficient (so that, after filtering,
#   # there will still be sufficient values for use in the Monte Carlo simulation)
#   # use centered regression
#   coeff_table_unfiltered_c <- mvrnorm(n = nsims*2, mu = resp_model_coeff_center, Sigma = Sigma_resp_model_center) # _c = centered
#   
#   # coerce into data frame 
#   coeff_table_unfiltered_df_c <- as.data.frame(coeff_table_unfiltered_c)
#   
#   # edit column names
#   colnames(coeff_table_unfiltered_df_c) <- c("Intercept", "log.Wref", "TempRef", "Depth.VM")
#   
#   # pull out just the standard errors
#   std_errors <- summary(resp_model_center)$coefficients[, 2]
#   
#   ### code to add respiration rate coefficient thresholds ###
#   
#   # We need to filter out a few extreme values that create unrealistic results (to avoid unrealistic results). To do this, we first centered the regression on reference temp and fish mass values, which removes much of the covariance that otherwise occurs between the intercept and slope of the regression. Next, we constrain the variance for each parameter, treating each parameter as univariate (rather than multivariate, which we were doing before centering the regression due to high covariance). Below this is done by constraining each predictor within +/- 2 standard errors. 
#   
#   # filter out any values greater than 2 SE from the mean coefficient value (removes ~5% of values)
#   all_intercepts <- coeff_table_unfiltered_df_c$Intercept[
#     -which(abs(coeff_table_unfiltered_df_c$Intercept-resp_model_coeff_center[1]) > std_errors[1]*2)]
#   all_mass_coeff <- coeff_table_unfiltered_df_c$log.Wref[
#     -which(abs(coeff_table_unfiltered_df_c$log.Wref-resp_model_coeff_center[2]) > std_errors[2]*2)]
#   all_temp_coeff <- coeff_table_unfiltered_df_c$TempRef[
#     -which(abs(coeff_table_unfiltered_df_c$TempRef-resp_model_coeff_center[3]) > std_errors[3]*2)]
#   all_depth_coeff <- coeff_table_unfiltered_df_c$Depth.VM[
#     -which(abs(coeff_table_unfiltered_df_c$Depth.VM-resp_model_coeff_center[4]) > std_errors[4]*2)]
#   # samples 1000 values from this vector of filtered values (all with equal chance of being
#   # selected)
#   # shouldn't matter much whether replace = TRUE or FALSE, if there are sufficient values remaining to sample when replace = FALSE
#   a0.vec <- sample(x = all_intercepts, size = nsims, replace = FALSE)
#   a1.vec <- sample(x = all_mass_coeff, size = nsims, replace = FALSE)
#   a2.vec <- sample(x = all_temp_coeff, size = nsims, replace = FALSE)
#   a3.vec <- sample(x = all_depth_coeff, size = nsims, replace = FALSE)
#   
#   # check what % of values were filtered out
#   100-(length(all_intercepts)/length(coeff_table_unfiltered_df_c$Intercept)*100)
#   # approximately 5% filtered out, as we would expect constraining +/- 2 SE
#   
#   ### end of code to add respiration rate coefficient thresholds ###
#   
#   
#   H_VM_or_NM.vec <- runif(nsims, min = 1, max = 1) # use 1 for VM fish or 0 for NM fish function 
#   TT_epipelagic.vec <- runif(n = nsims, min = 14-0.5, max = 14+0.5) # +/- 1 degrees from measured epipelagic temp in EXPORTS location. Adjust based on environment modeled
#   TT_mesopelagic.vec <- runif(n = nsims, min = 11-0.5, max = 11+0.5) # +/- 1 degrees from measured mesopelagic temp in EXPORTS location. Adjust based on environment modeled
#   TT_mean.vec <- runif(n = nsims, min = 12.5-0.5, max = 12.5+0.5) # +/- 1 degrees from measured mean(mean(epipelagic), mean(mesopelagic)) temp in EXPORTS location. Adjust as needed
#   
#   activity_factor_resting.vec <- runif(nsims, min = 0.5, max = 1) # assumes fish is at or just 2x above SMR during daytime in deep water (i.e., fish is either 1/2 of routine metabolic rate RMR or at RMR)
#   activity_factor_foraging.vec <- runif(nsims, min = 1, max = 4) # Assumes fish metabolizes at least at 2x SMR (if SMR is 1/2 RMR) and as high as 8x SMR during nighttime foraging
#   activity_factor_migrating.vec <- runif(nsims, min = 1, max = 4) # Assumes fish metabolizes at least at 2x SMR (if SMR is 1/2 RMR) and as high as 8x SMR during during migration (relevant for VM fish only)
#   time_resting.vec <- runif(nsims, min = .33, max = .66) # assumes migrating fish spend 50% of the time inactive (deep) in mesopelagic zone resting. Same for NM fish.
#   
#   G_VM.vec <- runif(nsims, min = 0.018-0.018*0.5, max = 0.018+0.018*0.5)
#   # 0.0018 for VM fish, 0.0007 for NM fish
#   
#   theta.vec <- runif(nsims, min = 0.54, max = 1)
#   
#   # First, assume that oxygen consumption during respirometry experiments is due only to M, not also to SDA (food assimilation). Realistically some of the fish from Ikeda 2016 regression were probably still digesting food, so we will also consider a case where the measured O2 updated is due to both M and SDA combined (respiration rate, in units of energy usage, = M + SDA). 
#   
#   prop_SDA.vec <- runif(nsims, min = 0.14 * 0.8, max = 0.14 * 1.2) # Nominal value from Brett and Groves (they used 0.14), allowed to vary by 20%
#   prop_egestion.vec <- runif(nsims, min = 0.20 * 0.8, max = 0.20 * 1.2) # Nominal value from Brett and Groves (0.20), allowed to vary by 20%
#   prop_excretion.vec <- runif(nsims, min = 0.07 * 0.8, max = 0.07 * 1.2) # Nominal value from Brett and Groves (0.07), allowed to vary by 20%
#   
#   # for SDAe function
#   timeTo90.vec <- runif(nsims, min = 10, max = 14)
#   peakTime.vec <- runif(nsims, min = 3- 3*.8, max = 3+ 3*.8)
#   
#   Q_ox_kJ_g_O2.vec <- runif(nsims, min = 13.4- (13.4*0.1), max = 13.4 + (13.4*0.1)) # From Jobling 1994. Oxycalorific equivalent (energy obtained from respired oxygen) in units of kJ/gO2. Nelson et al. 2016 states that this may vary by 10% at most (sd of about 1).
#   
#   RQ.vec <- runif(nsims, min = 0.7, max = 2) # Hernandez-Leon et al. 2019 used respiratory quotient (CO2 respired/O2 consumed) of 0.97 (from Omori and Ikeda, 1984). Davison et al. 2013 and Ariza et al. 2015 used respiratory quotient of 0.90 (Brett and Groves 1979). Hudson let this range from 0.7-1.0, given uncertainty associated with mesopelagic fish. According to Brett and Groves 1979, RQ could = 2 for fish that undergo some amount of anaerobic metabolism.
#   
#   # Dmax range is from EXPORTS acoustics echogram from Mei Sato (and happens to 
#   # be the same as the assumptions made for chapter 1... migrators appear to be 
#   # swimming as deep as 700 m during the day, or as shallow as about 300 m during the day)
#   # We'll extend to 750 m to match our net depth intervals. 
#   Dmax.vec <- runif(nsims, min = 300, max = 750)
#   
#   # range used for chapter 1 of 10-90 meters for migrators at night seems reasonable too. 
#   # We'll extend to 100 m to match our net depth intervals 
#   Dmin.vec <- runif(nsims, min = 10, max = 100) 
#   
#   B.vec <- runif(nsims, min = 199, max = 201) # set this to 200 m, but adjust when 
#   # finding the flux at 500 m flux boundary 
#   
#   S.vec <- runif(nsims, min = 5- (5*.3), max = 5+ (5*.3)) # see table 1 for refs for swimming speed. note that in the code we make swim speed in cm/s and then multiply by 0.01 (unit used in references), but in the equation (to remove the 0.01 unit conversion factor) we provide swim speed in m/s. 
#   
#   energy_in_pellet.vec <- runif(nsims, min = 0.6, max = 3.5)
#   
#   C_in_pellet.vec <- runif(nsims, min = 0.05, max = 0.10)
#   
#   proportion_egestion_below_150m.vec <- runif(nsims, min = 0.5, max = 0.90) # 50-90% sinks below 150m
#   
#   C_c_vm.vec <- runif(nsims, min = 0.038, max = 0.248) # Ratio of carbon:wet weight of VM fish, from mean of species captured 
#   # in this study, and Childress and Nygaard (1973). In Davison et al. 2013--they used 0.49 for carbon:wet weight ratio, but in Childress and Nygaard paper this was about 0.50 for carbon:DRY weight. Mesopelagic fish ranged from about 60% to 90% water, and carbon:ash free dry weight (assumed here to be the same as dry weight, because % ash is small) ranged from 38-62%. To find bounds: min carbon:wet weight = (1-0.90) * 0.38 = 0.038, and max carbon:wet weight = (1-0.60) * 0.62 = 0.248
#   
#   X_r.vec <- runif(nsims, min = 0.008 - 0.8*0.008, max = 0.008 + 0.8*0.008) # see Table 1--these values came from Wilson et al. 2009 and then were converted to mgC /d/ WMfish. Allow to vary by 80%
#   
#   X_e.vec <- runif(nsims, min = 0.5, max = 0.9)
#   
#   # For Pauly 1980 mortality regression, simply allow mu to vary by +/- 0.3
#   # log_{10}(\mu) = a_0 + a_1 log_{10}(L_{\infty}) + a_2log_{10}(K) + a_3log_{10}(T_{mean}) 
#   aa0 <- -0.0066
#   aa1 <- -0.2790
#   aa2 <- 0.6543
#   aa3 <- 0.4634
#   Linfinity <- 8.7 # average of the 3 myctophid Linfinity values in Pauly 1980
#   K <- 0.2 # arbitrary value between 0.3 and 0.4 for VM fish (Pauly 1980) 
#   # and 0.05 and 0.17 for NM fish (Childress et al 1980)
#   # use same parameter values for NM fish for now due to limited data
#   TT_mean_nominal_value <- 12.5 # set TT_mean
#   mu.mean <- 10^ (aa0 + aa1*log10(Linfinity) + aa2*log10(K) + aa3*log10(TT_mean_nominal_value))
#   mu.vec <- runif(nsims, min = mu.mean- mu.mean*.5, max = mu.mean+ mu.mean*.5) 
#   prop_mortality_below_150m.vec <- runif(nsims, min = .40, max = .90) 
#   
#   # with all parameter vectors created, now calculate C flux for each MC iteration in nsims
#   for(i in 1:nsims) {
#     C_flux_i <- C.flux.function.VM(W_w_f_VM=W_w_f_VM, Wref=Wref, Tref=Tref, a0=a0.vec[i], a1=a1.vec[i], a2=a2.vec[i], a3=a3.vec[i], H_VM_or_NM=H_VM_or_NM.vec[i], TT_epipelagic=TT_epipelagic.vec[i], TT_mesopelagic=TT_mesopelagic.vec[i], TT_mean=TT_mean.vec[i],  activity_factor_resting=activity_factor_resting.vec[i], activity_factor_foraging=activity_factor_foraging.vec[i], activity_factor_migrating=activity_factor_migrating.vec[i], time_resting=time_resting.vec[i], Dmax=Dmax.vec[i], Dmin=Dmin.vec[i], B=B.vec[i], S=S.vec[i], prop_SDA=prop_SDA.vec[i], prop_egestion=prop_egestion.vec[i], prop_excretion=prop_excretion.vec[i], G_VM=G_VM.vec[i], theta=theta.vec[i], peakTime=peakTime.vec[i], timeTo90=timeTo90.vec[i], X_r=X_r.vec[i], X_e=X_e.vec[i], mu = mu.vec[i], Q_ox_kJ_g_O2=Q_ox_kJ_g_O2.vec[i], RQ=RQ.vec[i], energy_in_pellet=energy_in_pellet.vec[i], C_in_pellet=C_in_pellet.vec[i],  proportion_egestion_below_150m=proportion_egestion_below_150m.vec[i], C_c_vm=C_c_vm.vec[i], prop_mortality_below_150m=prop_mortality_below_150m.vec[i])[1] 
#     output[i] = C_flux_i # units of mg
#   }
#   output.mat <- cbind(output, a0.vec, a1.vec, a2.vec, a3.vec, H_VM_or_NM.vec, TT_epipelagic.vec, TT_mesopelagic.vec, TT_mean.vec, activity_factor_resting.vec, activity_factor_foraging.vec, activity_factor_migrating.vec, time_resting.vec, Dmax.vec, Dmin.vec, B.vec, S.vec, prop_SDA.vec, prop_egestion.vec, prop_excretion.vec, G_VM.vec, theta.vec, peakTime.vec, timeTo90.vec, X_r.vec, X_e.vec, mu.vec, Q_ox_kJ_g_O2.vec, RQ.vec, energy_in_pellet.vec, C_in_pellet.vec, proportion_egestion_below_150m.vec, C_c_vm.vec, prop_mortality_below_150m.vec)
#   colnames(output.mat) <- c("C_flux", "a0", "a1", "a2", "a3", "VM_NM", "T_epi", "T_meso", "T_mean", "A_r", "A_f", "A_m", "t_r", "Dmax", "Dmin", "B", "S", "prop_SDA", "prop_F","prop_X", "G_VM", "theta", "t_peak", "t_90", "X_r", "X_e", "mu", "Q_ox", "RQ", "E_p", "C_p", "F_e", "C_c", "M_e")
#   
#   # re-run but make new data frame that includes only the C_flux values by pathway (for boxplot later). Note that from C.flux.function, output is c(total, resp, SDA, egestion, excretion, mortality, ingestion). At the end of the C.flux.function call below, the [1], [2], etc selects which pathway to enter into each column
#   pathways_output <- tibble(total_flux = rep(NA, nsims), resp_flux = rep(NA, nsims), SDA_flux = rep(NA, nsims), eg_flux = rep(NA, nsims), ex_flux = rep(NA, nsims), mort_flux = rep(NA, nsims))
#   for (i in 1:nsims){
#     
#     # total C flux
#     pathways_output[i,1] <- C.flux.function.VM(W_w_f_VM=W_w_f_VM, Wref=Wref, Tref=Tref, a0=a0.vec[i], a1=a1.vec[i], a2=a2.vec[i], a3=a3.vec[i], H_VM_or_NM=H_VM_or_NM.vec[i], TT_epipelagic=TT_epipelagic.vec[i], TT_mesopelagic=TT_mesopelagic.vec[i], TT_mean=TT_mean.vec[i],  activity_factor_resting=activity_factor_resting.vec[i], activity_factor_foraging=activity_factor_foraging.vec[i], activity_factor_migrating=activity_factor_migrating.vec[i], time_resting=time_resting.vec[i], Dmax=Dmax.vec[i], Dmin=Dmin.vec[i], B=B.vec[i], S=S.vec[i], prop_SDA=prop_SDA.vec[i], prop_egestion=prop_egestion.vec[i], prop_excretion=prop_excretion.vec[i], G_VM=G_VM.vec[i], theta=theta.vec[i], peakTime=peakTime.vec[i], timeTo90=timeTo90.vec[i], X_r=X_r.vec[i], X_e=X_e.vec[i], mu = mu.vec[i], Q_ox_kJ_g_O2=Q_ox_kJ_g_O2.vec[i], RQ=RQ.vec[i], energy_in_pellet=energy_in_pellet.vec[i], C_in_pellet=C_in_pellet.vec[i],  proportion_egestion_below_150m=proportion_egestion_below_150m.vec[i], C_c_vm=C_c_vm.vec[i], prop_mortality_below_150m=prop_mortality_below_150m.vec[i])[1]
#     
#     # respiratory C flux
#     pathways_output[i,2] <- C.flux.function.VM(W_w_f_VM=W_w_f_VM, Wref=Wref, Tref=Tref, a0=a0.vec[i], a1=a1.vec[i], a2=a2.vec[i], a3=a3.vec[i], H_VM_or_NM=H_VM_or_NM.vec[i], TT_epipelagic=TT_epipelagic.vec[i], TT_mesopelagic=TT_mesopelagic.vec[i], TT_mean=TT_mean.vec[i],  activity_factor_resting=activity_factor_resting.vec[i], activity_factor_foraging=activity_factor_foraging.vec[i], activity_factor_migrating=activity_factor_migrating.vec[i], time_resting=time_resting.vec[i], Dmax=Dmax.vec[i], Dmin=Dmin.vec[i], B=B.vec[i], S=S.vec[i], prop_SDA=prop_SDA.vec[i], prop_egestion=prop_egestion.vec[i], prop_excretion=prop_excretion.vec[i], G_VM=G_VM.vec[i], theta=theta.vec[i], peakTime=peakTime.vec[i], timeTo90=timeTo90.vec[i], X_r=X_r.vec[i], X_e=X_e.vec[i], mu = mu.vec[i], Q_ox_kJ_g_O2=Q_ox_kJ_g_O2.vec[i], RQ=RQ.vec[i], energy_in_pellet=energy_in_pellet.vec[i], C_in_pellet=C_in_pellet.vec[i],  proportion_egestion_below_150m=proportion_egestion_below_150m.vec[i], C_c_vm=C_c_vm.vec[i], prop_mortality_below_150m=prop_mortality_below_150m.vec[i])[2]
#     
#     # SDA C flux
#     pathways_output[i,3] <- C.flux.function.VM(W_w_f_VM=W_w_f_VM, Wref=Wref, Tref=Tref, a0=a0.vec[i], a1=a1.vec[i], a2=a2.vec[i], a3=a3.vec[i], H_VM_or_NM=H_VM_or_NM.vec[i], TT_epipelagic=TT_epipelagic.vec[i], TT_mesopelagic=TT_mesopelagic.vec[i], TT_mean=TT_mean.vec[i],  activity_factor_resting=activity_factor_resting.vec[i], activity_factor_foraging=activity_factor_foraging.vec[i], activity_factor_migrating=activity_factor_migrating.vec[i], time_resting=time_resting.vec[i], Dmax=Dmax.vec[i], Dmin=Dmin.vec[i], B=B.vec[i], S=S.vec[i], prop_SDA=prop_SDA.vec[i], prop_egestion=prop_egestion.vec[i], prop_excretion=prop_excretion.vec[i], G_VM=G_VM.vec[i], theta=theta.vec[i], peakTime=peakTime.vec[i], timeTo90=timeTo90.vec[i], X_r=X_r.vec[i], X_e=X_e.vec[i], mu = mu.vec[i], Q_ox_kJ_g_O2=Q_ox_kJ_g_O2.vec[i], RQ=RQ.vec[i], energy_in_pellet=energy_in_pellet.vec[i], C_in_pellet=C_in_pellet.vec[i],  proportion_egestion_below_150m=proportion_egestion_below_150m.vec[i], C_c_vm=C_c_vm.vec[i], prop_mortality_below_150m=prop_mortality_below_150m.vec[i])[3]
#     
#     # egestion flux
#     pathways_output[i,4] <- C.flux.function.VM(W_w_f_VM=W_w_f_VM, Wref=Wref, Tref=Tref, a0=a0.vec[i], a1=a1.vec[i], a2=a2.vec[i], a3=a3.vec[i], H_VM_or_NM=H_VM_or_NM.vec[i], TT_epipelagic=TT_epipelagic.vec[i], TT_mesopelagic=TT_mesopelagic.vec[i], TT_mean=TT_mean.vec[i],  activity_factor_resting=activity_factor_resting.vec[i], activity_factor_foraging=activity_factor_foraging.vec[i], activity_factor_migrating=activity_factor_migrating.vec[i], time_resting=time_resting.vec[i], Dmax=Dmax.vec[i], Dmin=Dmin.vec[i], B=B.vec[i], S=S.vec[i], prop_SDA=prop_SDA.vec[i], prop_egestion=prop_egestion.vec[i], prop_excretion=prop_excretion.vec[i], G_VM=G_VM.vec[i], theta=theta.vec[i], peakTime=peakTime.vec[i], timeTo90=timeTo90.vec[i], X_r=X_r.vec[i], X_e=X_e.vec[i], mu = mu.vec[i], Q_ox_kJ_g_O2=Q_ox_kJ_g_O2.vec[i], RQ=RQ.vec[i], energy_in_pellet=energy_in_pellet.vec[i], C_in_pellet=C_in_pellet.vec[i],  proportion_egestion_below_150m=proportion_egestion_below_150m.vec[i], C_c_vm=C_c_vm.vec[i], prop_mortality_below_150m=prop_mortality_below_150m.vec[i])[4]
#     
#     # excretion flux
#     pathways_output[i,5] <- C.flux.function.VM(W_w_f_VM=W_w_f_VM, Wref=Wref, Tref=Tref, a0=a0.vec[i], a1=a1.vec[i], a2=a2.vec[i], a3=a3.vec[i], H_VM_or_NM=H_VM_or_NM.vec[i], TT_epipelagic=TT_epipelagic.vec[i], TT_mesopelagic=TT_mesopelagic.vec[i], TT_mean=TT_mean.vec[i],  activity_factor_resting=activity_factor_resting.vec[i], activity_factor_foraging=activity_factor_foraging.vec[i], activity_factor_migrating=activity_factor_migrating.vec[i], time_resting=time_resting.vec[i], Dmax=Dmax.vec[i], Dmin=Dmin.vec[i], B=B.vec[i], S=S.vec[i], prop_SDA=prop_SDA.vec[i], prop_egestion=prop_egestion.vec[i], prop_excretion=prop_excretion.vec[i], G_VM=G_VM.vec[i], theta=theta.vec[i], peakTime=peakTime.vec[i], timeTo90=timeTo90.vec[i], X_r=X_r.vec[i], X_e=X_e.vec[i], mu = mu.vec[i], Q_ox_kJ_g_O2=Q_ox_kJ_g_O2.vec[i], RQ=RQ.vec[i], energy_in_pellet=energy_in_pellet.vec[i], C_in_pellet=C_in_pellet.vec[i],  proportion_egestion_below_150m=proportion_egestion_below_150m.vec[i], C_c_vm=C_c_vm.vec[i], prop_mortality_below_150m=prop_mortality_below_150m.vec[i])[5]
#     
#     # mortality flux
#     pathways_output[i,6] <- C.flux.function.VM(W_w_f_VM=W_w_f_VM, Wref=Wref, Tref=Tref, a0=a0.vec[i], a1=a1.vec[i], a2=a2.vec[i], a3=a3.vec[i], H_VM_or_NM=H_VM_or_NM.vec[i], TT_epipelagic=TT_epipelagic.vec[i], TT_mesopelagic=TT_mesopelagic.vec[i], TT_mean=TT_mean.vec[i],  activity_factor_resting=activity_factor_resting.vec[i], activity_factor_foraging=activity_factor_foraging.vec[i], activity_factor_migrating=activity_factor_migrating.vec[i], time_resting=time_resting.vec[i], Dmax=Dmax.vec[i], Dmin=Dmin.vec[i], B=B.vec[i], S=S.vec[i], prop_SDA=prop_SDA.vec[i], prop_egestion=prop_egestion.vec[i], prop_excretion=prop_excretion.vec[i], G_VM=G_VM.vec[i], theta=theta.vec[i], peakTime=peakTime.vec[i], timeTo90=timeTo90.vec[i], X_r=X_r.vec[i], X_e=X_e.vec[i], mu = mu.vec[i], Q_ox_kJ_g_O2=Q_ox_kJ_g_O2.vec[i], RQ=RQ.vec[i], energy_in_pellet=energy_in_pellet.vec[i], C_in_pellet=C_in_pellet.vec[i],  proportion_egestion_below_150m=proportion_egestion_below_150m.vec[i], C_c_vm=C_c_vm.vec[i], prop_mortality_below_150m=prop_mortality_below_150m.vec[i])[6]
#   }
#   return(list(output.mat, pathways_output))
# }
# 
# 
# 
# 
# 
# 
# 
# #### Function to run MC simulation for NM fish ####
# 
# 
# # Find NM fish C flux such that NM fish MC function output is generated using the same parameter values (where applicable) as for VM fish. For argument VM_output, input the name of the first object in the list generated in the MC_func_VM function (in the MC_func_VM function, the first object that is returned in the list is the data frame of C flux values per iteration and parameter values selected for each of those iterations)
# 
# # Note: need to create the appropriate VM fish output and pass this to MC_func_NM each time MC_func_NM is used (e.g. if we want nsims=5000 for MC_func_NM, need to pass the MC_func_VM output that also uses nsims = 5000 so that the parameter vectors are the correct length)
# 
# # Monte Carlo function for NM fish. Always use same parameters as in VM fish function, where applicable (not all VM fish parameters appear in NM fish C flux function)
# 
# # Enter W_w_f in g and Wref in mg
# 
# MC_func_NM <- function(nsims, W_w_f_NM, Wref, Tref, output_VM){
#   output_NM <- rep(NA, nsims)
#   W_w_f_NM <- W_w_f_NM # set this in argument of MC_func_NM
#   
#   # generate parameter vector that are unique to VM fish (no equivalent for VM fish)
#   G_NM.vec <- runif(nsims, min = 0.007-0.007*0.5, max = 0.007+0.007*0.5) 
#   prop_non_detrital_NM_prey.vec <- runif(nsims, min = 0.7-0.7*.2, max = 0.7+0.7*.2) 
#   H_VM_or_NM.vec <- runif(nsims, min = 0, max = 0) # use 1 for VM fish or 0 for NM fish function 
#   
#   # generate parameter vectors that are the same as those generated for MC_func_VM
#   a0.vec <- output_VM$a0
#   a1.vec <- output_VM$a1
#   a2.vec <- output_VM$a2
#   a3.vec <- output_VM$a3
#   TT_mesopelagic.vec <- output_VM$T_meso
#   TT_mean.vec <- output_VM$T_mean 
#   activity_factor_resting.vec <- output_VM$A_r
#   activity_factor_foraging.vec <-output_VM$A_f
#   time_resting.vec <- output_VM$t_r
#   theta.vec <- output_VM$theta
#   prop_SDA.vec <- output_VM$prop_SDA
#   prop_egestion.vec <- output_VM$prop_F
#   prop_excretion.vec <- output_VM$prop_X
#   Q_ox_kJ_g_O2.vec <- output_VM$Q_ox
#   RQ.vec <- output_VM$RQ
#   energy_in_pellet.vec <- output_VM$E_p
#   C_in_pellet.vec <- output_VM$C_p
#   C_c_nm.vec <- output_VM$C_c
#   X_r.vec <- output_VM$X_r
#   mu.vec <- output_VM$mu
#   
#   # with all parameter vectors created, now calculate C flux for each MC iteration in nsimsow fa
#   for(i in 1:nsims) {
#     C_flux_i <- C.flux.function.NM (W_w_f_NM=W_w_f_NM, Wref=Wref, Tref=Tref, a0=a0.vec[i], a1=a1.vec[i], a2=a2.vec[i], a3=a3.vec[i], H_VM_or_NM=H_VM_or_NM.vec[i], TT_mesopelagic=TT_mesopelagic.vec[i], TT_mean=TT_mean.vec[i],  activity_factor_resting=activity_factor_resting.vec[i], activity_factor_foraging=activity_factor_foraging.vec[i], time_resting=time_resting.vec[i], prop_SDA=prop_SDA.vec[i], prop_egestion=prop_egestion.vec[i], prop_excretion=prop_excretion.vec[i], G_NM=G_NM.vec[i], theta=theta.vec[i], X_r=X_r.vec[i], mu = mu.vec[i], Q_ox_kJ_g_O2=Q_ox_kJ_g_O2.vec[i], RQ=RQ.vec[i], prop_non_detrital_NM_prey=prop_non_detrital_NM_prey.vec[i], energy_in_pellet=energy_in_pellet.vec[i], C_in_pellet=C_in_pellet.vec[i], C_c_nm=C_c_nm.vec[i])[1] 
#     output_NM[i] = C_flux_i # units of mg 
#   }
#   output_NM.mat <- cbind(output_NM, a0.vec, a1.vec, a2.vec, a3.vec, H_VM_or_NM.vec, TT_mesopelagic.vec, TT_mean.vec, activity_factor_resting.vec, activity_factor_foraging.vec, time_resting.vec, prop_SDA.vec, prop_egestion.vec, prop_excretion.vec, G_NM.vec, theta.vec, X_r.vec, mu.vec, Q_ox_kJ_g_O2.vec, RQ.vec, prop_non_detrital_NM_prey.vec, energy_in_pellet.vec, C_in_pellet.vec, C_c_nm.vec)
#   colnames(output_NM.mat) <- c("C_flux", "a0", "a1", "a2", "a3", "VM_NM", "T_meso", "T_mean", "A_r", "A_f", "t_r", "prop_SDA", "prop_F","prop_X", "G_NM", "theta", "X_r", "mu", "Q_ox", "RQ", "Z", "E_p", "C_p", "C_c")
#   
#   # re-run but make new data frame that includes only the C_flux values by pathway (for boxplot later). Note that from C.flux.function, output is c(total, resp, SDA, egestion, excretion, mortality, ingestion). At the end of the C.flux.function call below, the [1], [2], etc selects which pathway to enter into each column
#   pathways_output_NM <- tibble(total_flux = rep(NA, nsims), resp_flux = rep(NA, nsims), SDA_flux = rep(NA, nsims), eg_flux = rep(NA, nsims), ex_flux = rep(NA, nsims), mort_flux = rep(NA, nsims))
#   for (i in 1:nsims){
#     
#     # total C flux
#     pathways_output_NM[i,1] <- C.flux.function.NM(W_w_f_NM=W_w_f_NM, Wref=Wref, Tref=Tref, a0=a0.vec[i], a1=a1.vec[i], a2=a2.vec[i], a3=a3.vec[i], H_VM_or_NM=H_VM_or_NM.vec[i], TT_mesopelagic=TT_mesopelagic.vec[i], TT_mean=TT_mean.vec[i],  activity_factor_resting=activity_factor_resting.vec[i], activity_factor_foraging=activity_factor_foraging.vec[i], time_resting=time_resting.vec[i], prop_SDA=prop_SDA.vec[i], prop_egestion=prop_egestion.vec[i], prop_excretion=prop_excretion.vec[i], G_NM=G_NM.vec[i], theta=theta.vec[i], X_r=X_r.vec[i], mu = mu.vec[i], Q_ox_kJ_g_O2=Q_ox_kJ_g_O2.vec[i], RQ=RQ.vec[i], prop_non_detrital_NM_prey=prop_non_detrital_NM_prey.vec[i], energy_in_pellet=energy_in_pellet.vec[i], C_in_pellet=C_in_pellet.vec[i], C_c_nm=C_c_nm.vec[i])[1]
#     
#     # respiratory C flux
#     pathways_output_NM[i,2] <- C.flux.function.NM(W_w_f_NM=W_w_f_NM, Wref=Wref, Tref=Tref, a0=a0.vec[i], a1=a1.vec[i], a2=a2.vec[i], a3=a3.vec[i], H_VM_or_NM=H_VM_or_NM.vec[i], TT_mesopelagic=TT_mesopelagic.vec[i], TT_mean=TT_mean.vec[i],  activity_factor_resting=activity_factor_resting.vec[i], activity_factor_foraging=activity_factor_foraging.vec[i], time_resting=time_resting.vec[i], prop_SDA=prop_SDA.vec[i], prop_egestion=prop_egestion.vec[i], prop_excretion=prop_excretion.vec[i], G_NM=G_NM.vec[i], theta=theta.vec[i], X_r=X_r.vec[i], mu = mu.vec[i], Q_ox_kJ_g_O2=Q_ox_kJ_g_O2.vec[i], RQ=RQ.vec[i], prop_non_detrital_NM_prey=prop_non_detrital_NM_prey.vec[i], energy_in_pellet=energy_in_pellet.vec[i], C_in_pellet=C_in_pellet.vec[i], C_c_nm=C_c_nm.vec[i])[2] 
#     
#     # SDA C flux
#     pathways_output_NM[i,3] <- C.flux.function.NM(W_w_f_NM=W_w_f_NM, Wref=Wref, Tref=Tref, a0=a0.vec[i], a1=a1.vec[i], a2=a2.vec[i], a3=a3.vec[i], H_VM_or_NM=H_VM_or_NM.vec[i], TT_mesopelagic=TT_mesopelagic.vec[i], TT_mean=TT_mean.vec[i],  activity_factor_resting=activity_factor_resting.vec[i], activity_factor_foraging=activity_factor_foraging.vec[i], time_resting=time_resting.vec[i], prop_SDA=prop_SDA.vec[i], prop_egestion=prop_egestion.vec[i], prop_excretion=prop_excretion.vec[i], G_NM=G_NM.vec[i], theta=theta.vec[i], X_r=X_r.vec[i], mu = mu.vec[i], Q_ox_kJ_g_O2=Q_ox_kJ_g_O2.vec[i], RQ=RQ.vec[i], prop_non_detrital_NM_prey=prop_non_detrital_NM_prey.vec[i], energy_in_pellet=energy_in_pellet.vec[i], C_in_pellet=C_in_pellet.vec[i], C_c_nm=C_c_nm.vec[i])[3]
#     
#     # egestion flux
#     pathways_output_NM[i,4] <- C.flux.function.NM(W_w_f_NM=W_w_f_NM, Wref=Wref, Tref=Tref, a0=a0.vec[i], a1=a1.vec[i], a2=a2.vec[i], a3=a3.vec[i], H_VM_or_NM=H_VM_or_NM.vec[i], TT_mesopelagic=TT_mesopelagic.vec[i], TT_mean=TT_mean.vec[i],  activity_factor_resting=activity_factor_resting.vec[i], activity_factor_foraging=activity_factor_foraging.vec[i], time_resting=time_resting.vec[i], prop_SDA=prop_SDA.vec[i], prop_egestion=prop_egestion.vec[i], prop_excretion=prop_excretion.vec[i], G_NM=G_NM.vec[i], theta=theta.vec[i], X_r=X_r.vec[i], mu = mu.vec[i], Q_ox_kJ_g_O2=Q_ox_kJ_g_O2.vec[i], RQ=RQ.vec[i], prop_non_detrital_NM_prey=prop_non_detrital_NM_prey.vec[i], energy_in_pellet=energy_in_pellet.vec[i], C_in_pellet=C_in_pellet.vec[i], C_c_nm=C_c_nm.vec[i])[4] 
#     
#     # excretion flux
#     pathways_output_NM[i,5] <- C.flux.function.NM(W_w_f_NM=W_w_f_NM, Wref=Wref, Tref=Tref, a0=a0.vec[i], a1=a1.vec[i], a2=a2.vec[i], a3=a3.vec[i], H_VM_or_NM=H_VM_or_NM.vec[i], TT_mesopelagic=TT_mesopelagic.vec[i], TT_mean=TT_mean.vec[i],  activity_factor_resting=activity_factor_resting.vec[i], activity_factor_foraging=activity_factor_foraging.vec[i], time_resting=time_resting.vec[i], prop_SDA=prop_SDA.vec[i], prop_egestion=prop_egestion.vec[i], prop_excretion=prop_excretion.vec[i], G_NM=G_NM.vec[i], theta=theta.vec[i], X_r=X_r.vec[i], mu = mu.vec[i], Q_ox_kJ_g_O2=Q_ox_kJ_g_O2.vec[i], RQ=RQ.vec[i], prop_non_detrital_NM_prey=prop_non_detrital_NM_prey.vec[i], energy_in_pellet=energy_in_pellet.vec[i], C_in_pellet=C_in_pellet.vec[i], C_c_nm=C_c_nm.vec[i])[5]
#     
#     # mortality flux
#     pathways_output_NM[i,6] <- C.flux.function.NM(W_w_f_NM=W_w_f_NM, Wref=Wref, Tref=Tref, a0=a0.vec[i], a1=a1.vec[i], a2=a2.vec[i], a3=a3.vec[i], H_VM_or_NM=H_VM_or_NM.vec[i], TT_mesopelagic=TT_mesopelagic.vec[i], TT_mean=TT_mean.vec[i],  activity_factor_resting=activity_factor_resting.vec[i], activity_factor_foraging=activity_factor_foraging.vec[i], time_resting=time_resting.vec[i], prop_SDA=prop_SDA.vec[i], prop_egestion=prop_egestion.vec[i], prop_excretion=prop_excretion.vec[i], G_NM=G_NM.vec[i], theta=theta.vec[i], X_r=X_r.vec[i], mu = mu.vec[i], Q_ox_kJ_g_O2=Q_ox_kJ_g_O2.vec[i], RQ=RQ.vec[i], prop_non_detrital_NM_prey=prop_non_detrital_NM_prey.vec[i], energy_in_pellet=energy_in_pellet.vec[i], C_in_pellet=C_in_pellet.vec[i], C_c_nm=C_c_nm.vec[i])[6]
#   }
#   return(list(output_NM.mat, pathways_output_NM))
# }
# 

#### function to get Ikeda 2016 respiration rate coefficients ####
get_ikeda_coefs <- function() {
  resp_data <- read.csv("data/Individual_fish_flux_data/Ikeda_2016_data_VM_NM.csv", header = TRUE)
  
  # make the categorical depth category a factor
  resp_data$Cat_depth <- as.factor(resp_data$Cat_depth)
  
  # re-calculate temp data as in Ikeda 2016 (take inverse, multiple by 1000, and convert to K)
  Tref <- 15 # Reference temperature (same as used in C flux function)
  Wref <- 1000 # Reference fish wet mass (same as used in C flux function)
  resp_data$Temp_var_T_ref <- 1000/(resp_data$Temp_C + 273.15) - 1000/(Tref+ 273.15) # modification for centering at Tref
  resp_data$logW_ref <- log(resp_data$mgWM) - log(Wref) # modification for centering at Wref
  # Note: only need to center temp and fish mass because these are the only two continuous predictors (habitat depth is categorical)
  
  # fit linear regression (with log-log transformation for respiration rate and fish wet mass)
  resp_model_center <- lm(log(R) ~ logW_ref + Temp_var_T_ref + Cat_depth, resp_data)
  
  # variance-covariance matrix
  Sigma_resp_model_center <- vcov(resp_model_center)
  
  # best fit model coefficients 
  resp_model_coeff_center <- coefficients(resp_model_center)
  se <- resp_model_coeff_center
  output <- list(resp_model_coeff_center = resp_model_coeff_center,
                 Sigma_resp_model_center = Sigma_resp_model_center)
  return(output)
}

#### function to generate Monte Carlo sensitivity analysis results ####

sim_C_flux <- function(simpars, simC, nsims, groupnames, migrate, thedata, fishdata, volumedata,par_sims,nominal.pars,C_t_bar_obs, qpars, calpars) {
  lwreg <- readRDS("analysis/length_weight_parameters.RDS")
  total_C_flux <- list()
  allgroupnames <- names(C_t_bar_obs)
  for (g in 1:length(groupnames)) {
    # 3) Assign taxa to analyze 
    taxagroup <- groupnames[g]
    # 5) filter length data
    taxadata <- dplyr::filter(thedata[[1]], Family == taxagroup)
    # 6) calculate length bins and P(l)
    # P(l) is the probability of fish being in each length bin
    Pl_and_l <- make_Pl(taxadata$Std_length_mm)
    nl <- length(Pl_and_l$L)
    C_flux <- rep(NA, nsims)
    for (sim in 1:nsims) {
      if (simpars) pars <- as.list(par_sims[sim,])
      if (!simpars) pars <- nominal.pars
      # generate random catch
      if (simC)  {  
        C_t_g_sim <- make_random_catch(fishdata = fishdata,
                                       volumedata = volumedata,
                                       group = taxagroup)
        C_t_bar <- mean(C_t_g_sim)
        q <- runif(1,  qpars[1], qpars[2])
        cal <- runif(1,  calpars[1], calpars[2])
      }
      if (!simC) {
        C_t_bar <- C_t_bar_obs[[which(allgroupnames == taxagroup)]]
        q <- mean(qpars)
        cal <- mean(calpars)
      }
      if(migrate[g] == "VM") tau_l <- sapply(X = Pl_and_l$L, FUN = make_tau, pars = pars, taxagroup = taxagroup, VM = TRUE, lwreg = lwreg)[1,]
      if(migrate[g] == "NM") tau_l <- sapply(X = Pl_and_l$L, FUN = make_tau, pars = pars, taxagroup = taxagroup, VM = FALSE, lwreg = lwreg)[1,]
      # 8) Multiply P(L) times tau_l and sum
      ind_C_flux <- sum(tau_l * Pl_and_l$P)
      
      # 9) Calculate mean catch from simulated data
      n_hat <- C_t_bar / q / cal
      
      # 10) Calculate C flux from group
      C_flux[sim] <- n_hat * ind_C_flux
    }
    total_C_flux[[g]] <- C_flux
  }
  names(total_C_flux) <- groupnames
  return(total_C_flux)
}

