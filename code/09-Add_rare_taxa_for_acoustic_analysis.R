
# For Dr. Mei Sato's work adding rare taxa back into fishdata for acoustic analysis 

# make table of only Cook catch 
fishdata_cook <- fishdata[fishdata$Ship == "Cook",]

# summarize the catch by Tow_net, then Taxon
fishdata_cook_summary <- fishdata_cook %>% group_by(Tow_net, Taxon) %>% summarize(n = n(), total_weight = sum(weight.to.use))


# make table of only Sarmiento catch

# first, need to add back in rare Sarmiento fish taxa that were removed in 
# 00-functions.R (in the "load_data()" function)
### load Sarmiento data, larger fish ####
SG2105_MOC_data_include_rare <- read.csv("data/MOCNESS_data/MOCNESS 10 Sarmiento/Sarmiento_MOC_fish_data.csv", header = TRUE)

# make new column that uses taxa_name_final that is simply called "taxa_code" 
# (see notes below on genetic vs morphological species identification)
SG2105_MOC_data_include_rare$taxa_code <- SG2105_MOC_data_include_rare$taxa_code_final

# filter out zooplankton (decapods, ostracods, shrimp, amphipods, pteropods, etc.)
# and also filter out squid
fish_data_Sarmiento_include_rare <- SG2105_MOC_data_include_rare %>% filter(
  taxa_code!="DECA" & taxa_code!="AMPH" & taxa_code!="RORO" & taxa_code!="ANSP" & taxa_code!="ACPU" & 
    taxa_code!="CRAN" & taxa_code!="SYDE" & taxa_code!="PNSP" & taxa_code!="ATOL" & taxa_code!="PTER" & 
    taxa_code!="CLSP" & taxa_code!="OSTR" & taxa_code!="PASI" & taxa_code!="EVSP")  	

# make new df that includes rare fishes
# Sarmiento fish taxa
fish_data_Sarmiento_include_rare <- fish_data_Sarmiento

# remove NA (empty) rows by selecting rows for which Ship=="Sarmiento"
fishdata_sarmiento_include_rare <- fish_data_Sarmiento_include_rare[fish_data_Sarmiento_include_rare$Ship == "Sarmiento",]



# need to fix summary to include Tow_net column so that code below runs *** (come back to this) 

# summarize the catch by Tow_net, then Taxon
fishdata_sarmiento_summary <- fishdata_sarmiento_include_rare %>% group_by(Tow_net, Taxon) %>% summarize(n = n(), total_weight = sum(weight.to.use))
