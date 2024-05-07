library(tidyr)
library(dplyr)
library(rMR)




### Read in adjusted file that includes Alphia IDs for each species / genera
all.dat <- read.csv(file = "data/data_for_weight_model_with_aphiaid.csv", header = T)
# Remove three outliers
all.dat <- filter(all.dat, !fish_no %in% c(88, 192 ,612))
all.dat <- select(all.dat, Scientific_name, AphiaID, Std_length_mm, total_weight_g, Lowest_taxon)
all.dat <- rename(all.dat, scientific.name = Scientific_name, lowest_taxon = Lowest_taxon)
all.dat$lowest_taxon <- tolower(all.dat$lowest_taxon)

#install.packages("jsonlite", repos="http://cran.r-project.org")
#install.packages("httr")


#Use the libraries
library(jsonlite) #https://cran.r-project.org/web/packages/jsonlite/
library(httr)

#Fill in the AphiaID you need

unique.taxa <- unique(all.dat$scientific.name)
all.dat$Class <- NA
all.dat$Subclass <- NA
all.dat$Superorder <- NA
all.dat$Order <- NA
all.dat$Suborder <- NA
all.dat$Infraorder <- NA
all.dat$Superfamily <- NA
all.dat$Family <- NA
all.dat$Genera <- NA
all.dat$Species <- NA

lookup.taxa <- function(x, level) {
  xunlist <- unlist(x)
  level.index <- which(xunlist == level)
  if(length(level.index)>0) {
    return(unname(xunlist[level.index + 1]))
  } else {return(NA)
  }
}

for (i in 1:length(unique.taxa)) {
  taxa <- unique.taxa[i]
  taxa.index <- which(all.dat$scientific.name == taxa)
  aphiaID <- all.dat$AphiaID[taxa.index[1]]
  #Build the URL to get the data from
  url <- sprintf("https://www.marinespecies.org/rest/AphiaClassificationByAphiaID/%d", aphiaID);
  #Get the actual data from the URL
  classificationTree <- fromJSON(url)
  all.dat$Class[taxa.index] <- lookup.taxa(classificationTree, "Class")
  all.dat$Subclass[taxa.index] <- lookup.taxa(classificationTree, "Subclass")
  all.dat$Superorder[taxa.index] <- lookup.taxa(classificationTree, "Superorder")
  all.dat$Order[taxa.index] <- lookup.taxa(classificationTree, "Order")
  all.dat$Suborder[taxa.index] <- lookup.taxa(classificationTree, "Suborder")
  all.dat$Infraorder[taxa.index] <- lookup.taxa(classificationTree, "Infraorder")
  all.dat$Superfamily[taxa.index] <- lookup.taxa(classificationTree, "Superfamily")
  all.dat$Family[taxa.index] <- lookup.taxa(classificationTree, "Family")
  if(all.dat$lowest_taxon[taxa.index[1]] %in% c("genus", "species")) all.dat$Genera[taxa.index] <- lookup.taxa(classificationTree, "Genus")
  if(all.dat$lowest_taxon[taxa.index[1]] =="species")  all.dat$Species[taxa.index] <- lookup.taxa(classificationTree, "Species")

}
nrow(all.dat)

saveRDS(all.dat, file = "data/alldata_taxonomy.RDS")
all.dat2 <- readRDS(file = "data/alldata_taxonomy.RDS")
nrow(all.dat2)

