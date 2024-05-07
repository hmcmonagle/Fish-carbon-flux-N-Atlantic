# written by Helena McMonagle
# Notes created 14 December, 2023
# Last edited January 4, 2024

** = notes for which files to share with SeaBASS

Steps taken to get to final MOCNESS dataset for SG2105 were as follows:

1) downloaded "OTZ_FISH_master.xlsx" from Llopiz Lab Drive (Shared Drives > Llopiz Lab). Saved as "OTZ_FISH_masters_231214.xlsx". Note that other documents from cruise SG2105 (R/V Sarmiento) can be found at Llopiz Lab > OTZ Llopiz Lab > Cruises > SG2105 - Sarmiento 2021. This was downloaded on Dec 14 2023--any more recent updates made on the Google Drive will not be incorporated into the R code and analysis unless this process is repeated after those edits are made. 

2) Removed latter rows in dataset in data_entry tab that were from other cruises. This left 662 rows from cruise_no = SG2105. 

3) Saved "LINKED_taxa_code" tab in "OTZ_FISH_master.xlsx" as a separate xlsx document titled "OTZ_FISH_master_231214_taxa_code_info.xlsx". This links the taxa_code with other taxonomic information if using as a relational database. 

** May want to share OTZ_FISH_master_231214_taxa_code_info.xlsx with Inia at NASA for SeaBASS database, if taxa code is deemed necessary information


4) Saved "data_entry" tab as a separate CSV file ("OTZ_FISH_master_231214_SG2105.csv") for use in downstream data processing for fish carbon flux project for SG2105 data. 

5) Remaining steps relate to OTZ_FISH_master_231214_SG2105.csv. This CSV file is comma delimited. Deleted "Flag for imaging/fec" column from this CSV file. 

6) Renamed "OTZ_FISH_master_231214_SG2105.csv" to "Sarmiento_MOC_fish_data.csv" to match what this file was called (before data updates from Lyndsey Llevebre) in the analysis code. 


7) Edited this data frame "Sarmiento_MOC_fish_data.csv" to remove unnecessary columns. Removed the following columns: 
- "station_no" (this column was blank anyway for SG2105 and was a carry over from other Llopiz Lab cruises that use a similar spreadsheet template)
- "gear_tow" (this is redundant with tow_no, and is a carry over from cruises where multiple net types were used) 
- "tow_type" (all rows contained MOC and this is known for SG2105 as all fish were collected with this gear from this cruise)
- deleted "flash frozen" column. Some rows (individuals used for respirometry) say "N" here, but actually all were flash frozen after respirometry was finished (just not immediately upon being brought on deck after net recovery) 
- "family", "genus" and "species" columns (these were left blank anyway but are found later in the code as needed)
- "sl_estimated" because all lengths are an estimate to some extent and SL that could not be measured as accurately are still used the same in final analysis
- "gut_visual_ID_complete", "gut_comments", "dna_extracted", "dna_extracted_date", "GenPopVial", "genetic_ID", "Flag for imagine/fec" were deleted as they were almost entirely blank anyway and a carry over from previous cruise templates, except for GenPopVial, which had some "Y" for yes filled out but relate to a completely different project 
- "sia_collected", which was blank for all rows and relates to other labs doing stable isotope analysis
- "fully_processed", also a carry over from other cruises where some fish were processed in various ways that others weren't (various tissue collections, otolith dissections, etc.) 
- gonad_stage, which was used for Lyndsey's reference to compare apparent stage with what the histological staging shows
- gonads_collected, which was just a way for Lyndsey to keep track of what gonads might already be preserved

The column "sex" (in which 1 = male, 2 = female, and 9 = unknown) is left in, but these assignments are uncertain particularly for immature fish. 

The column "gut_dispension" was left in, and this indicates whether the gut contents were collected for meta barcoding by Ann Bucklin's lab (UConn) with gut_dispensation=mb for meta barcoding or visual, morphological ID by Sarah Glancy in Llopiz Lab with gut_dispensation=visual. Not all stomachs have been or will be processed, but many morphological gut content IDs have been done and are available. 

photo_ID column was left in, and Helena added the photo IDs for the respirometry fish. These photos were taken at sea rather than back in Llopiz Lab after the cruise. These photos are too large for free storage on GitHub, so they are stored on Llopiz Lab Google Drive and on Helena McMonagle's work computer hard drive (and ext drive backup) and could be submitted to SeaBASS ** Note that there are generally duplicate photos because photos were also taken of respirometry individuals after flash freezing back in Llopiz lab after the cruise before processing for ETS. For the fish that had ETS done but not also respirometry (individuals 33-104), there is just the frozen fish photo from Llopiz lab. All fish were photographed except for individual 78 which seems to be missing the photo. 

8) Renamed some columns in "Sarmiento_MOC_fish_data.csv" to match what is in code (and was in previously used data frame for analysis)
- tow_no was renamed to Tow_number
- net_no was renamed to Net_number
- leave lowest_taxa_final and lowest_taxa_original because these show where genetic ID was used to update a morphological ID (and could be useful if genetic IDs are later found to be incorrect or become update as taxonomic assignment evolves). 
- renamed standard_length to Std_length_mm to shorten and provide unit


9) Added column for Ship, for which all rows are Sarmiento in "Sarmiento_MOC_fish_data.csv"

Changes made later in code, rather than in raw data frame: 
- In code, made a new column in which lowest_taxa_final is duplicated and named Lowest_taxon. 
- In code, removed total_length and fork_length, as all individuals from MOC-10 have standard length and relatively few have total or fork lengths (and only standard length was measured for smaller MOC-10 fish and for all MOC-1 fish, using ImageJ)
- made a new column for taxa_code_final that is just taxa_code 

** Note that "photo_ID" was left blank for the respirometry fish, but these were photographed too and those photos are still available. I have now filled out these image numbers from photos in the "respirometry specimen photos" folder in Llopiz Lab drive. 


** this "Sarmiento_MOC_fish_data.csv" data sheet is probably the most useful to share with Inia as a supplemental file, because we have removed unnecessary columns for the analysis, some of which were carry overs from other cruise datasheete templates in Llopiz Lab or which related to sample organization that is beyond the scope of what EXPORTS will use the data for (e.g., "FLAG_genetic" is just for use between Lyndsey and Annette to determine which samples to prioritize for genetic ID and was mostly left blank anyway). However, I will also share with her a version of fishdata.csv called "fishdata_EXPORTS.csv" that has the following columns removed, which are metadata that are not central to repeating the fish carbon flux analysis: 

- moved Fish.ID for the photographed fish to "fish_no" for fish number, as the photographed fish that were measured in imageJ did not have an individual fish number but rather this fish ID

- deleted "ImageJ.Number" (used for matching lengths with image J features)

- cruise_no (change "Ship" to cruise number if this needs to be in the data frame)

- packet_no (used for sample storage but not in analysis)

- taxa_code_original

- taxa_updated_date

- visual_ID_conf for visual ID confidence

- moved "Image.Name" for Sarmiento and Cook fish that were measured with photographs and imageJ over to "photo_ID" and deleted "Image.Name"

- ultimately delete "taxa_code_final" (will ultimately want to use Worms ID, but for now leave in for matching to this Worms ID) and just edit from "taxa_code_final" to "taxa_code"

- added Worms aphia ID column and called this "data_provider_category_manual" per SEABASS directions

- renamed "taxa_code_final_method" to "lowest_taxa_method"

- renamed "lowest_taxa_final_" to "lowest_taxa" 

- deleted Tow_net (this was just used for merging between individual fish data and MOCNESS metadata) 

- deleted "sex", "sex_confidence", "gonad_weight", "head_collected", "body_location", 

- removed extra "taxa_code" and "lowest_taxon" column (extraneous)

- VolFilt (this was an old calculation that Amy Maas's lab did for their net, but the final value to use for volume filtered for all of catch is in VolViltM3)

- edited "weight.to.use" to "weight_to_use" to match convention in other column names 


- changed NA to -9999 using find and replace (match case to NA, whole cell only so as not to accidentally change any notes text of species names with "na" in the word)

- removed "FracBiom" as fractions (net splits) have already been incorporated into duplicated individuals as necessary depending on these splits of the original samples (e.g., if sample had a split of "4", this means we only measured a quarter of the catch from that net so we multiply all those individuals by 4 to represent the full net's sample) 

- edited "VolFiltM3" to "volume_filtered_m^3"

- changed "photo_ID" to "/associated archives" associatedMedia

- deleted "taxon" (this is the same as either Family or lowest_taxon columns

- added whatever was in Family to lowest_taxon where there was no genus or species ID available in lowest_taxon 

- edited Gonostoma sp. to Sigmops sp. (more recently accepted taxon name)

- edited "lowest_taxon" to scientificName_manual

- edited old file to fishdata_EXPORTS_noAphiaID.csv and the newer file to fishdata_EXPORTS.csv so that the file name we sent to Inia is simpler 

- added fields that I think we need to add at the top, based on Sosik lab SeaBASS files 

e.g.: 

/begin_header									
/received=20231019									
/identifier_product_doi=10.5067/SeaBASS/EXPORTS/DATA001									
/investigators=Heidi_Sosik	Collin_Roesler								
/affiliations=Woods_Hole_Oceanographic_Institution	Bowdoin_College								
/contact=hsosik@whoi.edu	croesler@bowdoin.edu											
/experiment=EXPORTS									
/cruise=EXPORTSNA									
/r2r_event=dy131-SE-20210504.2125.001	dy131-SE-20210504.2211.001								
/data_file_name=EXPORTS-EXPORTSNA_DY131_discrete_IFCB_plankton_and_particles_20210504225054_R1.sb									
/calibration_files=no_cal_files									
/eventID=D20210504T225054_IFCB125									
			

- changed "lowest_taxa_method" to "scientificNameMethod"

- edited all genus and species terms to remove the space and use an underscore _ instead 

- removed column Notes
- removed column fish_comments
- removed taxa_code

- edited header names so that naming convention is consistent (e.g., DayNight --> day_night and TimeLocal --> time_local)


Consider removing/editing the follow columns from fishdata_EXPORTS.csv:

- processed initials
- change "weight_to_use" to "weight" and remove "estimated_weight" and "total_weight" 

Ask Inia: wasn't sure if I was supposed to use /associated_archives or associatedMedia

