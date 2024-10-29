Explanation of data for data used in code to generate fishdata.csv
Written by Helena McMonagle, 2023-12-14
Last revised October 8, 2024

fishdata.csv contains the fish from the R/V Sarmiento's MOCNESS-10 and the 
R/V James Cook's MOCNESS-1 collected during the North Atlantic NASA EXPORTS/
WHOI Ocean Twilight Zone Project cruise in May 2021. 

We calculated a weight (or where possible used an empirical weight measurement) for all individuals. The current dataset also excludes fish smaller than 10 mm, which are generally not adult fish and have lesser known migratory behaviors and life histories compared to adults, and thus are not included in the fish carbon flux analysis that this dataset will be used for. 

Fish collections at sea were done by Joel Llopiz, Kayla Gardner, Julia Cox, Cristina Garcia-Fernandez and Helena McMonagle on the Sarmiento and led by Deborah Steinberg and Amy Maas on the Cook. Data collection in the lab and with ImageJ was done primarily by Lyndsey Lefebvre, Julia Cox and Helena McMonagle.

Column descriptions in fishdata.csv, which compiles fish data from both the RRS James Cook and the R/V Sarmiento de Gamboa, are as follows: 

Ship: Ship name, either the R/V Sarmiento or RRS Cook

Tow_number: the number for the entire haul, including all nets on the MOCNESS

Net_number: the depth-stratified net associated with each tow. Note that  on the MOCNESS-1 datasheets, net 1 was defined as that which remained open from the surface to 1000 m, whereas on the MOCNESS-10 datasheets, this net was called net zero. Catch from net 1 on the Cook was not saved, so we were not able to compare fish counts between the two ships for this net. This net was not as useful in any case because it is not depth-stratified. 

Taxon: taxonomic family, as it was originally written on the raw datasheets. 

Fish.ID: for the smaller of the Sarmiento fish specimens, Julia Cox (Llopiz Lab)
assigned each individual fish an ID from the fish she imaged and then measured 
in ImageJ. These fish were preserved in 90% ethanol prior to being photographed. 

Image.Name: This corresponds to the photos taken of the smaller Sarmiento fish. 

ImageJ.Number: This corresponds to the ImageJ measurements taken of the smaller Sarmiento fish. 

Std_length_mm: Standard length of each fish in millimeters. 

Notes: Notes that Julia Cox took while measuring the smaller Sarmiento fish.

Tow_net: used to merging volume data into the fishdata dataframe. A combination of the tow, then net associated with each fish. 

fish_no: Used to assign individual fish numbers to the larger Sarmiento fish that were flash frozen at sea and later processed in Llopiz Lab. 

packet_no: The packet that these larger fish were frozen in at sea. Information about the contents of these packets is available in the Llopiz Lab datasheets. 

Lowest_Taxon: For the larger Sarmiento fish specimens, a much lower taxonomic level was often recorded from laboratory morphological IDs or genetic barcoding IDs. 

taxa_code: Assigned by Lyndsey Lefebvre (Llopiz Lab) to each fish based on taxon.

total_weight: Empirically measured weight taken by Lyndsey Lefebvre for the larger Sarmiento specimens that were flash frozen at sea and then processed individually in Llopiz Lab. 

Date: the date that each fish was collected from the study site

TimeLocal: Local time during collection (which was also UTC time)

DayNight: time of day during collection

LatitudeStart: Latitude at the start of the tow

LongitudeStart: Longitude at the start of the tow

LatitudeEnd: Latitude at the end of the tow

LongitudeEnd: Longitude at the end of the tow

MeshSize: Mesh size of each MOCNESS net system, in micrometers

VolFilt: can delete this column for future analysis—this is the initial (but erroneous) volume filtered calculation from Cook for each net as determined by a flowmeter connected to the shipboard computer during MOCNESS deployment. 

VolFiltM3: This is the corrected volume filtered by each net, to be used for future analysis 

Temp: Temperature reading associated with each tow and net, though this was not always working on the MOCNESS so we can rely on CTD data from the same cruise

Salinity: Similar to temperature, ideally measured by MOCNESS but had some issues with this sensor. 

FracBiom: This column is no longer relevant for future analysis for the MOCNESS-1 samples (see split_correct_Cook below, which was used for generating the final dataset used for analysis). The smaller fish from the MOCNESS-10 samples that Julia Cox measured with ImageJ were all from quarter splits, and were multiplied by 4 accordingly. The smaller fish from the MOCNESS-1 that Helena McMonagle measured with ImageJ using Ecotaxa images (after formalin preservation) were from various different splits, now indicated in split_correct_Cook. 

split_correct_Cook: this number equals the split that Cook fish came from, if they did not come from a whole sample. If this was the case, these rows were replicated according to this split factor for Cook fish. This is the column used to replicate rows as needed in fishdata.csv to represent estimated fish counts and biomass given that some came from splits. 

Family: taxonomic family of each individual. This was the lowest taxonomic resolution that we were able to obtain for smaller Sarmiento fish and for all Cook fish (which were obtained exclusively from photographs for the larger Cook fish, or from Ecotaxa images for the smaller ones that were not sorted and photographed at sea) 

estimated_weight: weight estimated from a length-weight regression for the lowest taxon possible, where empirical weights were not available

weight.to.use: the empirical weight measured for that individual where available, or the estimated weight from the length-weight regression where an empirical weight measurement is not available. 