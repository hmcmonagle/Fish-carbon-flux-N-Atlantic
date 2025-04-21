
# For Dr. Mei Sato's work adding rare taxa back into fishdata for acoustic analysis 

# make table of only Cook catch 
fishdata_cook <- fishdata[fishdata$Ship == "Cook",]

# summarize the catch by Tow_net, then Taxon
fishdata_cook_summary <- fishdata_cook %>% group_by(Tow_net, Taxon) %>% summarize(n = n(), total_weight = sum(weight.to.use))


# make table of only Sarmiento catch

# first, need to add back in rare Sarmiento fish taxa that were removed in 
# 00-functions.R (in the "load_data()" function)

### load Sarmiento data ####
SG2105_MOC_data_include_rare <- read.csv("data/MOCNESS_data/MOCNESS 10 Sarmiento/Sarmiento_MOC_fish_data.csv", header = TRUE)

# make new column that uses taxa_name_final that is simply called "taxa_code" 
# (see notes below on genetic vs morphological species identification)
SG2105_MOC_data_include_rare$taxa_code <- SG2105_MOC_data_include_rare$taxa_code_final

# filter out zooplankton (decapods, ostracods, shrimp, amphipods, pteropods, etc.)
# and also filter out squid
fishdata_Sarmiento_include_rare <- SG2105_MOC_data_include_rare %>% filter(
  taxa_code!="DECA" & taxa_code!="AMPH" & taxa_code!="RORO" & taxa_code!="ANSP" & taxa_code!="ACPU" & 
    taxa_code!="CRAN" & taxa_code!="SYDE" & taxa_code!="PNSP" & taxa_code!="ATOL" & taxa_code!="PTER" & 
    taxa_code!="CLSP" & taxa_code!="OSTR" & taxa_code!="PASI" & taxa_code!="EVSP" & taxa_code!="THSP")  	

# remove NA (empty) rows by selecting rows for which Ship=="Sarmiento"
fishdata_Sarmiento_include_rare <- fishdata_Sarmiento_include_rare[fishdata_Sarmiento_include_rare$Ship == "Sarmiento",]

# repeat steps in 00-functions.R for data processing, but now df includes rare fish taxa

# remove fish that were so badly damaged that there's no standard length
# (results in a loss of about 3% of the abundance)
fishdata_Sarmiento_include_rare <- fishdata_Sarmiento_include_rare %>% filter(Std_length_mm!="NA")

# remove unnecessary columns (for more details, see readme_SG2105_data.txt)
# remove column "fork_length"
fishdata_Sarmiento_include_rare <- fishdata_Sarmiento_include_rare %>% dplyr::select(-c(total_length, fork_length))

# exclude fish that are smaller than 10 mm (we'll consider these larvae not 
# juveniles/adults, and vertical migration behavior may be different than expected
# based on studies of adults... though this is more relevant for the small Sarmiento
# fish that were photographed only and measured in ImageJ, as these are all >10 mm 
fishdata_Sarmiento_include_rare <- fishdata_Sarmiento_include_rare %>% filter(Std_length_mm>=10)

# make a new column titled "Lowest_taxon" which uses the result of "lowest_taxa_final", 
# which contains the finanl decision about taxonomic assignment after considering any 
# discrepancies between the visual and genetic ID (in which case we went with genetic), 
# and which contains any lower taxonomic level classification if available after genetic ID. 
fishdata_Sarmiento_include_rare$Lowest_taxon <- fishdata_Sarmiento_include_rare$lowest_taxa_final

# make new column of taxon that is just to family level (as opposed to lowest
# taxonomic level classified, so that this matches the Cook taxon column)
fishdata_Sarmiento_include_rare$Taxon <- as.vector(rep(NA, times = length(fishdata_Sarmiento_include_rare$Lowest_taxon)))


fishdata_Sarmiento_include_rare <- fishdata_Sarmiento_include_rare %>% mutate(Taxon = case_when(
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
    Lowest_taxon == "Nannobrachium atrum" | # added this rare taxon (no longer an accepted name)
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
    Lowest_taxon == "Gonostoma sp." | # added this rare taxon (just had Cyclothone in original analysis)
    Lowest_taxon == "Sigmops elongatus" | # added this rare taxon (just had Cyclothone in original analysis)
    Lowest_taxon == "Gonostomatidae" ~ 
    "Gonostomatidae", 
  Lowest_taxon == "Scopelogadus beanii" |
    Lowest_taxon == "Scopelogadus sp." |
    Lowest_taxon == "Poromitra megalops" |
    Lowest_taxon == "Melamphaidae" ~
    "Melamphaidae", 
  Lowest_taxon == "Xenodermichthys copei" |
    Lowest_taxon == "Alepocephalidae" ~ 
    "Alepocephalidae", 
  Lowest_taxon == "Arctozenus risso" |
    Lowest_taxon == "Paralepididae" ~  # added this (relatively) rare taxon
    "Paralepididae", 
  Lowest_taxon == "Bathylagus euryops" ~ # added this rare taxon
    "Bathylagidae",
  Lowest_taxon == "Nansenia sp." ~ # added this rare taxon
    "Microstomatidae",
  Lowest_taxon == "Merluccius merluccius" ~ # added this rare taxon
    "Merlucciidae",
  Lowest_taxon == "Platytroctidae" ~ # added this rare taxon
    "Platytroctidae"
  ))

# change fish_data_Sarmiento df name to distinguish between small Sarmiento fish
fish_Sarmiento_large_include_rare <- fishdata_Sarmiento_include_rare

# add a column for Net_Tow, which will be used later for merging fish data with volume data
fish_Sarmiento_large_include_rare$Tow_net <- paste(fish_Sarmiento_large_include_rare$Tow_number, "_", fish_Sarmiento_large_include_rare$Net_number)

# next add smaller Sarmiento fish data to the large Sarmiento fish data
# (those that were photographed and measured in ImageJ)

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
fish_Sarmiento_all_include_rare <- bind_rows(fish_Sarmiento_small, fish_Sarmiento_large_include_rare)

# remove fish from Net_number==0
fish_Sarmiento_all_include_rare <- fish_Sarmiento_all_include_rare %>% filter(Net_number!=0)


# make summary tables for sharing with Mei

# summarize the catch by Tow_net, then Taxon
fish_Sarmiento_all_include_rare_summary <- fish_Sarmiento_all_include_rare %>% group_by(Tow_net, Taxon) %>% summarize(n = n())

# re-do Sarmiento summary table but without rare taxa
fishdata_sarmiento <- fishdata[fishdata$Ship == "Sarmiento",]

# summarize the catch by Tow_net, then Taxon
fishdata_sarmiento_summary_without_rare <- fishdata_sarmiento %>% group_by(Tow_net, Taxon) %>% summarize(n = n(), total_weight = sum(weight.to.use))
# sum number of fish before adding back in rare taxa
sum(fishdata_sarmiento_summary_without_rare$n)

# create additional summary df with rare Sarmiento fishes that groups by lowest taxon (instead of family)
fish_Sarmiento_all_include_rare_sp <- fish_Sarmiento_all_include_rare %>% group_by(Tow_net, Lowest_taxon) %>% summarize(n = n())

# write csv to share with Mei that includes rare species and both large and smaller Sarmiento fishes
fishdata_Sarmiento_include_rare <- write.csv(fish_Sarmiento_all_include_rare, "data/MOCNESS_data/MOCNESS 10 Sarmiento/Sarmiento_fishdata_include_rare_spp.csv", row.names = FALSE)

# print number of fish picked from zooplankton at sea and more carefully identified beyond family
# level, versus those that were sorted back in lab from zooplankton and identified to family using photos
# rather than previously-frozen fishes 
# those that were sorted at sea and identified beyond species level have a non-NA "taxa_code" 
# (the convension used in Llopiz lab for more thoroughly processed fishes)
lower_taxonomic_ID_Sarmiento_fish <- fish_Sarmiento_all_include_rare %>% filter(taxa_code!="NA")
# find number of rows in lower_taxonomic_ID_Sarmiento_fish
nrow(lower_taxonomic_ID_Sarmiento_fish) # processed fully in Joel's lab
nrow(fish_Sarmiento_all_include_rare)
# nrow(fish_Sarmiento_all_include_rare) - nrow(lower_taxonomic_ID_Sarmiento_fish) # 884 ID-ed from photo

# difference between df with vs without rare taxa: 
(nrow(fish_Sarmiento_all_include_rare) - sum(fishdata_sarmiento_summary_without_rare$n))/nrow(fish_Sarmiento_all_include_rare)
# about 113 fishes were added back in, or about 10%
