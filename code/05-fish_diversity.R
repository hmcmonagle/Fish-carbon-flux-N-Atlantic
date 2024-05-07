# # Stacked bar plot of fish fish diversity by tow
# using Sarmiento fish data (subsample that was IDed
# morphologically and genetically)

# Nov 29, 2023, HM

# 1) Load functions and packages and fishdata filtered for 3 most abundant families
source("code/01-calc_c_flux.R")
library(scales)
library(glue)
library(ggtext)

# 2) Subset fishdata as needed for diversity indices 

# subset to include only the fish from Ship==Sarmiento 
# that also have a non-NA value for Lowest_taxon. add "div" to indicate 
# that this is the subsect of fishdata used for diversity analyses
fishdata_div <- subset(fishdata, Ship=="Sarmiento" & !is.na(Lowest_taxon))

length(fishdata_div$Ship)

fishdata_div$Lowest_taxon <- as.factor(fishdata_div$Lowest_taxon)

levels(fishdata_div$Lowest_taxon)

# create a subset of data for plot that shows only the taxa that are identified to 
# species level by excluding rows for which Lowest_taxon==Cyclothone sp., 
# Lampanyctus sp., Nannobrachium sp., or Scopelogadus sp. 
fishdata_spp <- subset(fishdata_div, Lowest_taxon!="Cyclothone sp." & 
                         Lowest_taxon!="Lampanyctus sp." & 
                         Lowest_taxon!="Nannobrachium sp." & 
                         Lowest_taxon!="Scopelogadus sp.")


# find number of species by finding unique values of Lowest_taxon
unique_Lowest_taxon <- unique(fishdata_spp$Lowest_taxon)

# add up number of unique values to find species richness
length(unique_Lowest_taxon)

# find Shannon-Wiener diversity index using Lowest_taxon ** 


# 3) Re-order taxa factors for plotting

# change the order of Lowest_taxon factors so that Lowest_taxon
# categories are organized by family, 
# with different tints or shades of the same color for each family
order_all_taxa <- c("Argyropelecus hemigymnus",  "Argyropelecus lychnus", "Argyropelecus olfersii", "Argyropelecus sp.",
                    "Valenciennellus tripunctulatus",
                    "Cyclothone braueri", "Cyclothone microdon", "Cyclothone pallida", "Cyclothone sp.",
                    "Scopelogadus beanii", "Scopelogadus sp.", "Xenodermichthys copei", 
                    "Benthosema glaciale", "Diaphus subtilis", "Lampanyctus macdonaldi", "Lampanyctus ater", 
                    "Lampanyctus sp.", "Myctophum punctatum", "Nannobrachium sp.", "Notoscopelus resplendens", 
                    "Protomyctophum arcticum", "Symbolophorus veranyi", "Myctophidae")

# Manually change order of the legend
fishdata_div$Lowest_taxon <- factor(
  fishdata_div$Lowest_taxon,
  levels = order_all_taxa,  # Specify the desired order
  ordered = TRUE
)

# we will included all species in the first plot, but if desired, the rarer species
# could be excluded to simplify the plot (e.g., by removing taxa
# with only n = 1 individual). This will be done later in step 6) in the code below. 


# 4) Create colors for stacked bar plot

# edit colors to match Dark2 color palette. Make all fish of the same
# family different tints or shades of the same color in Dark2 palette. 

# generate lighter colors from Dark2

# Define a function to lighten a color
lighten_color <- function(color, factor = 0.5) {
  rgb <- col2rgb(color)
  new_rgb <- rgb + (255 - rgb) * factor
  new_hex <- rgb(new_rgb[1, ], new_rgb[2, ], new_rgb[3, ], maxColorValue = 255)
  return(new_hex)
}

# make a vector of the colors to be lightened
colors_to_lighten <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666")

# loop through colors_to_lighten to make a lighter version of each one
lighter_colors <- c()
for (i in 1:length(colors_to_lighten)) {
  lighter_colors[i] <- lighten_color(colors_to_lighten[i])
}

# make an even lighter version of each color
lighter_colors2 <- c()
for (i in 1:length(lighter_colors)) {
  lighter_colors2[i] <- lighten_color(lighter_colors[i])
}

# make an even lighter version of each color
lighter_colors3 <- c()
for (i in 1:length(lighter_colors2)) {
  lighter_colors3[i] <- lighten_color(lighter_colors2[i])
}

# lighter colors 4 and up are too light, so instead make some shades darker than Dark2
# make a vector of the colors to be darkened
colors_to_darken <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666")

# Define a function to darken a color
darken_color <- function(color, factor = 0.2) {
  rgb <- col2rgb(color)
  new_rgb <- rgb * (1 - factor)
  new_hex <- rgb(new_rgb[1, ], new_rgb[2, ], new_rgb[3, ], maxColorValue = 255)
  return(new_hex)
}

# Define a function to darken a color
darken_color <- function(color, factor = 0.2) { # adjust factor to make the jump to light or darker more or less
  rgb <- col2rgb(color)
  new_rgb <- rgb * (1 - factor)
  new_hex <- rgb(new_rgb[1, ], new_rgb[2, ], new_rgb[3, ], maxColorValue = 255)
  return(new_hex)
}

# Darken the original color
darker_color <- darken_color(colors_to_darken)

# now, make a loop that darkens all of the following colors in Dark2:

# loop through colors_to_darken to make a darker version of each one
darker_colors <- c()
for (i in 1:length(colors_to_darken)) {
  darker_colors[i] <- darken_color(colors_to_darken[i])
}

# make an even darker version of each color
darker_colors2 <- c()
for (i in 1:length(darker_colors)) {
  darker_colors2[i] <- darken_color(darker_colors[i])
}

# make an even darker color
darker_colors3 <- c()
for (i in 1:length(darker_colors2)) {
  darker_colors3[i] <- darken_color(darker_colors2[i])
}


# 5) Generate plot with all taxa

# plot fish_most_common_18_plot, making font size bigger
fish_all_taxa_plot <- ggplot(fishdata_div, aes(x = Tow_number, fill = Lowest_taxon)) + 
  geom_bar(position="fill") + 
  scale_y_continuous(labels = scales::percent, breaks = c(0, 0.25, 0.5, 0.75, 1), expand = c(0,0)) +
  labs(x="Tow Number \n", y="\n Percent of fish belonging to each taxon") + 
  theme(legend.title=element_blank(), legend.position = "right", strip.background.x = element_blank(), 
        strip.text.x = element_text(size = 24),
        panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.background = element_blank(),
        panel.grid.minor = element_blank(), plot.margin = margin(t = 10, unit = "pt"),
        axis.line = element_line(colour = "black"), 
        text = element_text(size = 24)) + 
  scale_fill_manual(values = c(lighter_colors[2], colors_to_darken[2], darker_colors[2], darker_colors2[2], darker_colors3[2],
                                            lighter_colors[3], colors_to_darken[3], darker_colors[3], darker_colors2[3],
                                            lighter_colors[4], darker_colors[4], colors_to_darken[6], 
                               lighter_colors3[1], lighter_colors2[1], lighter_colors[1], colors_to_darken[1],
                               darker_colors[1], darker_colors2[1], 
                               lighter_colors2[5], lighter_colors[5], colors_to_darken[5], darker_colors[5], 
                               darker_colors2[5]))

# print plot
fish_all_taxa_plot


# 5) Subset data for another version of the plot that excludes the rarest (n=1) taxa.

# find the rarest among all the Lowest_taxon
rarest_taxon <- fishdata_div %>% 
  group_by(Lowest_taxon) %>% 
  summarize(n = n()) %>% 
  arrange(n) %>% 
  head(10)

rarest_taxon

# subset the dataset fishdata_div to remove the rarest taxa with n = 1
fishdata_without_rarest <- fishdata_div %>% 
  filter(Lowest_taxon != "Cyclothone pallida", Lowest_taxon != "Scopelogadus sp.", 
         Lowest_taxon != "Diaphus subtilis", Lowest_taxon != "Lampanyctus ater",
         Lowest_taxon != "Notoscopelus resplendens")

# change the order of Lowest_taxon factors so that Lowest_taxon
# categories are organized by family, 
# with different tints or shades of the same color for each family

order_without_rarest <- c("Argyropelecus hemigymnus",  "Argyropelecus lychnus", "Argyropelecus olfersii", "Argyropelecus sp.",
                            "Valenciennellus tripunctulatus",
                            "Cyclothone braueri", "Cyclothone microdon", "Cyclothone sp.",
                            "Scopelogadus beanii", "Scopelogadus sp.", "Xenodermichthys copei", 
                            "Benthosema glaciale", "Lampanyctus macdonaldi",
                            "Lampanyctus sp.", "Myctophum punctatum", "Nannobrachium sp.",
                            "Protomyctophum arcticum", "Symbolophorus veranyi", "Myctophidae")


# Manually change order of the legend
fishdata_without_rarest$Lowest_taxon <- factor(
  fishdata_without_rarest$Lowest_taxon,
  levels = order_without_rarest,  # Specify the desired order
  ordered = TRUE
)


# 6) Replot with rarest taxa excluded 
fish_without_rarest_plot <- ggplot(fishdata_without_rarest, aes(x = Tow_number, fill = Lowest_taxon)) + 
  geom_bar(position="fill") + 
  scale_y_continuous(labels = scales::percent, breaks = c(0, 0.25, 0.5, 0.75, 1), expand = c(0,0)) +
  labs(x="Tow Number \n", y="\n Percent of fish belonging to each taxon") + 
  theme(legend.title=element_blank(), legend.position = "right", legend.text = element_text(size = 21),
        strip.background.x = element_blank(), 
        strip.text.x = element_text(size = 24),
        panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.background = element_blank(),
        panel.grid.minor = element_blank(), plot.margin = margin(t = 10, unit = "pt"),
        axis.line = element_line(colour = "black"), 
        text = element_text(size = 24)) + 
  scale_fill_manual(values = c(lighter_colors[2], colors_to_darken[2], darker_colors[2], darker_colors2[2], darker_colors3[2],
                               lighter_colors[3], colors_to_darken[3], darker_colors[3],
                               colors_to_darken[4], colors_to_darken[6], 
                               lighter_colors3[1], lighter_colors2[1], lighter_colors[1], colors_to_darken[1], 
                               darker_colors[1], darker_colors2[1], darker_colors3[1], 
                               colors_to_darken[5]))

# print plot
fish_without_rarest_plot

# plot with x axis as date in May instead of tow number
# first need to add date column to data frame

# 7) repeat but now manually change tow numbers to days that correspond to date of tow
# for use in manuscript

# could add a column to fishdata_3_families with date of tow where day is the
# date in May when the tow began 
# Note: have to leave a space (" ") before the Sarmiento tows to get this to plot correctly 
fishdata_3_families_date_in_May <- fishdata_without_rarest %>% mutate(Tow_start_date_in_May = case_when(
  Tow_number == 1  ~
    "06 \nDay",
  Tow_number == 2  ~
    "06 \nNight", 
  Tow_number == 3  ~
    "07 \nDay",
  Tow_number == 4  ~
    "13 \nNight",
  Tow_number == 5  ~
    "13 \nDay",
  Tow_number == 6  ~
    "14 \nNight",
  Tow_number == 7  ~
    "17 \nDay",
  Tow_number == 8  ~
    "17 \nNight",
  Tow_number == 9  ~
    "18 \nNight"))


# 8) plot
fish_without_rarest_plot_date <- ggplot(fishdata_3_families_date_in_May, aes(x = Tow_start_date_in_May, fill = Lowest_taxon)) + 
  geom_bar(position="fill") + 
  scale_y_continuous(labels = scales::percent, breaks = c(0, 0.25, 0.5, 0.75, 1), expand = c(0,0)) +
  labs(x="\n Date in May and Time of Day \n", y="\n Percent of fish belonging to each taxon") + 
  theme(legend.title=element_blank(), legend.position = "right", legend.text = element_text(size = 21),
        strip.background.x = element_blank(), 
        strip.text.x = element_text(size = 24),
        panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.background = element_blank(),
        panel.grid.minor = element_blank(), plot.margin = margin(t = 10, unit = "pt"),
        axis.line = element_line(colour = "black"), 
        text = element_text(size = 24)) + 
  scale_fill_manual(values = c(lighter_colors[2], colors_to_darken[2], darker_colors[2], darker_colors2[2], darker_colors3[2],
                               lighter_colors[3], colors_to_darken[3], darker_colors[3],
                               colors_to_darken[4], colors_to_darken[6], 
                               lighter_colors3[1], lighter_colors2[1], lighter_colors[1], colors_to_darken[1], 
                               darker_colors[1], darker_colors2[1], darker_colors3[1], 
                               colors_to_darken[5]))

# print plot
fish_without_rarest_plot_date

# make legend in italics
fish_without_rarest_plot_date + theme(legend.text = element_text(face = "italic"))

# remove italic formatting from the Lowest_taxon==Myctophidae and wherever Lowest_taxon ends in "_sp."
# so that it can be formatted correctly in the manuscript


fishdata_3_families_date_in_May_formatted <- fishdata_3_families_date_in_May %>% mutate(
  noItalics = str_ends(Lowest_taxon, "ae|sp."),
  Italics = if_else(noItalics, glue("{Lowest_taxon}"), glue("<i>{Lowest_taxon}<i>")))

class(fishdata_3_families_date_in_May_formatted$Italics)
fishdata_3_families_date_in_May_formatted$Italics <- as.factor(fishdata_3_families_date_in_May_formatted$Italics)

# Manually change order of the legend (need to re-do after fixing italics)
order_without_rarest <- c("<i>Argyropelecus hemigymnus<i>",  "<i>Argyropelecus lychnus<i>", "<i>Argyropelecus olfersii<i>", "Argyropelecus sp.",
                          "<i>Valenciennellus tripunctulatus<i>",
                          "<i>Cyclothone braueri<i>", "<i>Cyclothone microdon<i>", "Cyclothone sp.",
                          "<i>Scopelogadus beanii<i>", "<i>Xenodermichthys copei<i>", 
                          "<i>Benthosema glaciale<i>", "<i>Lampanyctus macdonaldi<i>",
                          "Lampanyctus sp.", "<i>Myctophum punctatum<i>", "Nannobrachium sp.",
                          "<i>Protomyctophum arcticum<i>", "<i>Symbolophorus veranyi<i>", "Myctophidae")

fishdata_3_families_date_in_May_formatted$Italics <- factor(
  fishdata_3_families_date_in_May_formatted$Italics,
  levels = order_without_rarest,  # Specify the desired order
  ordered = TRUE
)

fish_without_rarest_plot_date_italics <- ggplot(fishdata_3_families_date_in_May_formatted, aes(x = Tow_start_date_in_May, fill = Italics)) + 
  geom_bar(position="fill") + 
  scale_y_continuous(labels = scales::percent, breaks = c(0, 0.25, 0.5, 0.75, 1), expand = c(0,0)) +
  labs(x="\n Date in May and Time of Day \n", y="\n Percent of fish belonging to each taxon") + 
  theme(legend.title=element_blank(), legend.position = "right", legend.text = element_markdown(size = 21),
        strip.background.x = element_blank(), 
        strip.text.x = element_text(size = 24),
        panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.background = element_blank(),
        panel.grid.minor = element_blank(), plot.margin = margin(t = 10, unit = "pt"),
        axis.line = element_line(colour = "black"), 
        text = element_text(size = 24)) + 
  scale_fill_manual(values = c(lighter_colors[2], colors_to_darken[2], darker_colors[2], darker_colors2[2], darker_colors3[2],
                               lighter_colors[3], colors_to_darken[3], darker_colors[3],
                               colors_to_darken[4], colors_to_darken[6], 
                               lighter_colors3[1], lighter_colors2[1], lighter_colors[1], colors_to_darken[1], 
                               darker_colors[1], darker_colors2[1], darker_colors3[1], 
                               colors_to_darken[5]))

# print plot
fish_without_rarest_plot_date_italics

# find n= for plot
length(fishdata_3_families_date_in_May_formatted$Taxon)

ggsave("graphics/submission_graphics/fish_diversity.jpg", fish_without_rarest_plot_date_italics, width = 13, height = 9, units = "in", dpi = 500)


