
# Set up ####

# 1) Load functions and fishdata filtered for 3 most abundant families
source("code/01-calc_c_flux.R")

# 2) Modify volumedata so there is just one row for each unique tow (for
# plotting lat/long of each tow)
volumedata_for_plotting <- volumedata[!duplicated(volumedata$Tow_number),]

# filter to include just night tows (plot is too busy otherwise)
volumedata_for_plotting <- volumedata_for_plotting[volumedata_for_plotting$DayNight == "Night",]

# 2) Generate map of study site and tow start/end location

library(sf)
library(rnaturalearth)
library(cowplot)

world <- ne_countries(scale = "medium", returnclass = "sf")

labels_df <- data.frame(
  longitude = c(-22, -6.5, -8, -19, -12, -20.3),
  latitude = c(54, 42.8, 53.6, 51, 49, 47.2),
  label = c("Atlantic\nOcean", "Spain", "Ireland", "Porcupine\nAbyssal\nPlain", 
            "Study site", "Study site with tows")
)

# offset to make study area box
long_lims <- c(-15.3, -14.6)
lat_lims <- c(48.7, 49.2)

# zoomed out to show general study area
main_map <- ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-25, -5), ylim = c(41, 55)) +
  geom_rect(data = volumedata_for_plotting, aes(xmin = long_lims[1], 
                                                xmax = long_lims[2], 
                                                ymin = lat_lims[1], 
                                                ymax = lat_lims[2]),
            fill = NA, color = "black", linewidth = 0.5) +
  theme_bw(base_size=23) + theme(panel.grid = element_blank()) +
  labs(x="\nLongitude", y="Latitude\n") + geom_text(data = labels_df, 
                                                    aes(x = longitude, y = latitude, label = label), size = 7)

# zoomed in to show tow locations (start in circles and end in triangles)
# if desired, set starting lat/long shape to 3 to view tow numbers 

inset_map <- ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = long_lims, ylim = lat_lims) +
  geom_point(data = volumedata_for_plotting, aes(x=LongitudeStart, y=LatitudeStart, color = Ship), size = 4, shape = 1) +
  geom_point(data = volumedata_for_plotting, aes(x=LongitudeEnd, y=LatitudeEnd, color = Ship), size = 4, shape = 2) +
  theme_void() + theme(panel.border = element_rect(fill =NA, colour = "black", linewidth=1.5)) + 
  geom_segment(data = volumedata_for_plotting, aes(x = LongitudeStart, xend = LongitudeEnd, 
                                                   y = LatitudeStart, yend = LatitudeEnd, color = Ship)) +
  scale_color_manual(values = c("Sarmiento" = "#D95F02", "Cook" = "#7570B3")) +
  theme(legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 19))

# combine into one plot with an inset
study_site_map <- ggdraw() + 
  draw_plot(main_map) + 
  draw_plot(inset_map, 
            height = 0.4, 
            x = -0.015, 
            y = 0.16
  )

# save
ggsave("graphics/submission_graphics/study_map.jpg", study_site_map, width = 11, height = 9, units = "in", dpi = 500)

