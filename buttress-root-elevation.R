# ############################################################## ####
#                                                                ####
# Title: A novel use of mangrove buttress root height to infer   ####
#        soil surface elevation change following forest          ####  
#        degradation and clearance                               ####  
#                                                                ####
# Author: Jack W Hill                                            ####
# Collaborators: Jaona Ravelonjatovo                             #### 
#                Ismael Ratefinjanahary                          ####
#                Lisa Benson                                     ####
#                Leah Glass                                      ####
#                                                                ####
# ############################################################## ####

##### 0. set up workspace #####

# packages may need to be installed using 
# install.packages("package-name")
library(broom)
library(cowplot)
library(egg)
library(emmeans)
library(ggspatial)
library(lemon)
library(ragg)
library(raster)
library(readxl)
library(rstatix)
library(sf)
library(terra)

# at end to avoid namespace conflicts 
library(tidyverse)

# set seed for work with random elements
set.seed(20220522)

# feature colours
triL_col <- "#7D31ED" # purple
triM_col <- "#ED7D31" # orange
triR_col <- "#13d863" # green

# function for ggplot custom theme -- general
theme_buttress <- function() {
  theme_classic(base_size=12,base_family="sans") %+replace%
    theme(
      text = element_text(colour = "black"),
      axis.text = element_text(size = rel(0.75)),
      axis.text.x = element_text(margin = margin(2, 0, 3, 0)),
      axis.text.y = element_text(margin = margin(0, 2, 0, 1)),
      axis.ticks = element_line(colour = "black"),
      axis.line = element_line(colour = "black"),
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_blank()
    )
}

# function for ggplot custom theme -- mapping
theme_buttress_map <- function() {
  theme_classic(base_size=12,base_family="sans") %+replace%
    theme(
      axis.text = element_text(size = rel(0.75)),
      axis.title = element_blank(),
      axis.ticks = element_line(colour = "black"),
      axis.line = element_blank(),
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_blank(),
      panel.border = element_rect(fill = NA,
                                  size = 1)
    )
}

#
##### 1. import soil data #####

soil_dat <- read_xlsx("data/bv-soil-loss-data.xlsx") %>% 
  rename(species = Species,
         tree.stump = "Tree/Stump",
         loss = "Loss Measurement (cm)",
         type = "smart class",
         plotID = ClusterPlot) %>% 
  mutate(type = str_remove(type, " mangrove"))

soil_dat$species <- as_factor(soil_dat$species)
soil_dat$type <- as_factor(soil_dat$type)
soil_dat$type <- fct_relevel(soil_dat$type,
                             "intact", "degraded", "deforested")
soil_dat$tree.stump <- as_factor(soil_dat$tree.stump)

# extract level of degradation for each plot
plot_deg <- soil_dat %>% 
  group_by(plotID) %>% 
  distinct(type, .keep_all = TRUE) %>%
  ungroup() %>% 
  select(plotID, type, 
         canopy_perc = `Canopy cover (%)`,
         degrade_perc = `Degradartion level (%)`)

# get mean soil loss for each species at each plot, 
# using stumps only
plot_soil_loss <- soil_dat %>% 
  filter(tree.stump == "Stump") %>% 
  select(plotID,
         species,
         loss) %>% 
  group_by(plotID, species) %>% 
  summarise(loss = mean(loss)) %>% 
  ungroup()

#
##### 2. import tree data #####

tree_dat_raw <- read_xlsx("data/bv-tree-size-data.xlsx") %>% 
  rename(species = Species,
         dbh = Circumference_measured,
         stump_diam_top = Stump_C_Top,
         stump_diam_base = Stump_C_base,
         stump_h = Stump_height,
         plotID = PlotID) %>% 
  # correct tree size measurements --
  # change to numeric and convert circumference to diameter
  mutate_at(vars(stump_diam_top,
                 stump_diam_base,
                 dbh),
            ~as.numeric(.)/pi) %>% 
  mutate(stump_h = as.numeric(stump_h))

# add the level of degradation for each plot
tree_dat <- left_join(tree_dat_raw, plot_deg,
                      by = "plotID")
  
# get mean stump size for each species at each plot
plot_stump_sizes <- tree_dat %>% 
  filter(!is.na(stump_h)) %>% 
  select(plotID,
         species,
         starts_with("stump",
                     ignore.case = FALSE)) %>% 
  group_by(plotID, species) %>% 
  summarise_at(vars(starts_with("stump",
                                ignore.case = FALSE)),
               mean) %>% 
  ungroup()

# create plotting dataset for tree size
tree_plot <- tree_dat %>%
  filter(!is.na(type),
         stump_diam_base > 0) %>% 
  group_by(type, species) %>% 
  summarise(mean_diam = mean(stump_diam_base),
            se_diam = sd(stump_diam_base)/sqrt(n()))

# save tree_plot for access out of R
write_csv(tree_plot,
          "table-output/tree-summary.csv")

#
##### 3. import GPS point data #####

# shapefiles required below are available at www.gadm.org,
# specifically (as at 18/07/2023) at www.gadm.org/download_country.html,
# by selecting "Madagascar" and downloading the "Shapefile" ZIP, which contains
# four shapefiles of different detail (levels). Shapefiles required are 
# indicated below

# read country and province shapefiles
mada_country <- st_read("data/mada-shps/NAME OF LEVEL 0 SHAPEFILE")
mada_provinces <- st_read("data/mada-shps/NAME OF LEVEL 1 SHAPEFILE")

# read Blue Ventures sampling point GPS data
mada_points_raw <- read_csv("data/bv-gps-data.csv")

# remove one error point;
# rename ID col for joining, join with degradation level info
mada_points_full <- mada_points_raw %>% 
  filter(Year_PlotID != "2017_PX_10-04") %>% 
  rename(plotID = PlotID) %>% 
  inner_join(plot_deg,
             by = "plotID")

# retain only the central plot point for each grouped set of plots
mada_points <- mada_points_full %>% 
  filter(str_detect(Year_PlotID, "-01") == TRUE)

crop_mada_provinces <- st_crop(mada_provinces,
                               xmin = 48.41, xmax = 48.585,
                               ymin = -13.48, ymax = -13.6)

# import data to show coring points of Arias-Ortiz et al. 2020 -- 
# this dataset can be downloaded from the "Data Availability" section of 
# Arias-Ortiz et al. 2020. No modifications are required to their raw data. 

    # Arias-Ortiz, A., Masqué, P., Glass, L. et al. Losses of Soil Organic
    # Carbon with Deforestation in Mangroves of Madagascar. Ecosystems 24,
    # 1–19 (2021). https://doi.org/10.1007/s10021-020-00500-z

ao_raw_dat <- read_xlsx("data/NAME OF ARIAS-ORTIZ RAW DATA FILE",
                        sheet = "Soil properties",
                        skip = 1) 

ao_clean_dat <- ao_raw_dat %>%
  slice(-1) %>% 
  select(core_id = `...1`,
         lat = Latitude,
         lon = Longitude) %>% 
  separate(col = core_id,
           into = c("major_id",
                    "depth_id"),
           sep = 6) %>% 
  filter(!is.na(major_id))

ao_core_gps <- ao_clean_dat %>% 
  distinct(major_id, .keep_all = TRUE) %>% 
  select(major_id, lat, lon)

#
##### 4. analyse soil data #####

# model development and assumption testing has been removed

# use square root transformation to normalise data
model1 <- lm(sqrt(loss) ~ type*species*tree.stump,
             data = soil_dat)

summary(model1)
Anova(model1) # publication Table 3
# result: 3-way interaction (sig)
# result: 2-way interaction (type*tree.stump, not sig)
# result: 2-way interaction (type*species, sig, decomposed below)

# simple two-way for trees vs stumps
test1 <- soil_dat %>% 
  group_by(tree.stump) %>%
  anova_test(sqrt(loss) ~ type*species,
             detailed = TRUE)
test1
# result: 2-way interaction (type*species, sig for stumps, decomposed below) 

# group by tree.stump and species to create simple main effects,
# filter for stumps only due to their sig 2-way interaction above
test2 <- soil_dat %>% 
  group_by(tree.stump, species) %>% 
  anova_test(sqrt(loss) ~ type,
             detailed = TRUE) %>% 
  as_tibble() %>% 
  filter(tree.stump == "Stump")
test2
# result: overall main effect of forest type

# pairwise comparisions using emmeans_test()
test3 <- soil_dat %>% 
  group_by(tree.stump, species) %>% 
  emmeans_test(sqrt(loss) ~ type,
               detailed = FALSE) %>% 
  filter(tree.stump == "Stump")
test3
# result: pairwise comparisions for forest type effect

# create marginal means with back-transformation
# including "pairwise" for post-hoc test with corrected p-vals
trans_emmeans <- emmeans(model1, pairwise ~ type*species*tree.stump,
                         type = "response")

# extract and tidy marginal means
mar_means <- tidy(trans_emmeans$emmeans) %>% 
  select(-c(statistic,
            p.value))

# extract and tidy contrasts (post-hoc results)
post_hoc <- tidy(trans_emmeans$contrasts) %>% 
  rename(t.statistic = statistic)

# keep only significant post-hoc pairs
sig_post_hoc <- post_hoc %>% 
  filter(adj.p.value <= 0.05)

# create dataset to plot 3-way interaction
soil_plot <- soil_dat %>% 
  group_by(type, species, tree.stump) %>% 
  summarise(mean_loss = mean(loss),
            se_loss = sd(loss)/sqrt(n())) %>% 
  ungroup()

# save soil_plot for access out of R
write_csv(soil_plot, 
          "table-output/soil-summary.csv")

#
##### 5. analyse tree data #####

# model development and assumption testing has been removed

# filter for stump entries only
tree_dat_stump <- tree_dat %>% 
  filter(!is.na(stump_diam_base))

# test two stump size variables of interest
model2 <- manova(cbind(stump_diam_base, stump_h) ~ type*species,
                 data = tree_dat_stump)
summary(aov(model2))

# single ANOVA for stump diameter
model3 <- lm(stump_diam_base ~ type * species,
              data = tree_dat_stump)
summary(aov(model3))
# interaction not significant

# post-hoc for stump diameter
model3_tukey <- tukey_hsd(model3)
model3_tukey

#
##### 6. visualise soil data #####

figure_3 <- ggplot(data = soil_plot,
                   aes(x = type,
                       y = mean_loss,
                       colour = tree.stump,
                       shape = tree.stump)) + 
  geom_point(aes(size = species),
             position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymax = mean_loss + se_loss,
                    ymin = mean_loss - se_loss),
                width = 0,
                position = position_dodge(width = 0.5),
                show.legend = FALSE) +
  facet_rep_grid(rows = vars(species)) +
  scale_x_discrete(name = "Forest status",
                   labels = c("Intact",
                              "Degraded",
                              "Deforested")) +
  scale_y_continuous(name = "\u0394d\n(difference in root elevation and soil surface, cm)",
                     limits = c(0,33),
                     expand = expansion(mult = c(0, 0.1))) +
  scale_colour_manual(name = "Measurement type",
                      labels = c("Trees",
                                 "Stumps"),
                      values = c(triL_col,
                                 triM_col)) + 
  scale_shape_manual(name = "Measurement type",
                     values = c(20, 18)) +
  scale_size_manual(name = "Measurement type",
                    values = c(3, 3)) +
  guides(colour = guide_legend(title = NULL,
                               override.aes = list(shape = c(16, 18),
                                                   colour = c(triL_col,
                                                              triM_col),
                                                   size = c(2.3,
                                                            2.8))),
         shape = "none",
         size = "none") +
  theme_buttress() %+replace%
  theme(legend.position = c(1, 0),
        legend.justification = c("right", "bottom"),
        legend.text = element_text(size = rel(0.75),
                                   margin = margin(0,0,0,-6)),
        legend.key.height = unit(11, "points"),
        legend.background = element_blank())
figure_3

# setup plot device
agg_png(filename = "figure-output/figure-3.png",
        width = 8, height = 13, 
        units = "cm", res = 1080,
        scaling = 1)

# plot, with facet labels
tag_facet(figure_3,
          open = "", close = "",
          tag_pool = c("B. gymnorrhiza", "C. tagal"),
          hjust = -0.1,
          fontface = 3, size = 3)

# close plot device
dev.off()

#
##### 7. visualise tree data #####

figure_4 <- ggplot(data = tree_plot,
                   aes(x = type,
                       shape = species,
                       colour = species)) +
  geom_errorbar(aes(ymin = mean_diam - se_diam,
                    ymax = mean_diam + se_diam),
                width = 0,
                show.legend = FALSE) +
  geom_point(aes(y = mean_diam,
                 size = species)) +
  scale_colour_manual(name = "Species",
                      labels = c("Bruguiera gymnorrhiza",
                                 "Ceriops tagal"),
                      values = c(triM_col,
                                 triL_col)) + 
  scale_shape_manual(name = "Species",
                     values = c(18, 20)) +
  scale_size_manual(name = "Species",
                    values = c(2,
                               2)) +
  scale_x_discrete(name = "Forest status",
                   labels = c("Intact",
                              "Degraded",
                              "Deforested")) +
  scale_y_continuous(name = "Stump diameter at base (cm)",
                     limits = c(0, 14),
                     expand = expansion(mult = c(0, 0))) +
  guides(colour = guide_legend(title = NULL,
                               override.aes = list(shape = c(18, 16),
                                                   colour = c(triM_col,
                                                              triL_col),
                                                   size = c(2.8,
                                                            2.2))),
         shape = "none",
         size = "none") +
  theme_buttress() %+replace%
  theme(legend.position = c(1.02, -0.01),
        legend.justification = c("right", "bottom"),
        legend.text = element_text(size = rel(0.75),
                                   face = "italic",
                                   margin = margin(0,0,0,-6)),
        legend.key.height = unit(11, "points"),
        legend.background = element_blank())
figure_4

agg_png(filename = "figure-output/figure-4.png",
        width = 8, height = 6.5, 
        units = "cm", res = 900,
        scaling = 1)
figure_4
dev.off()

#
##### 8. visualise GPS point data #####

main_map <- ggplot() +
  geom_sf(data = crop_mada_provinces,
          colour = NA,
          fill = "grey80") +
  geom_point(data = mada_points,
             aes(x = lon,
                 y = lat,
                 colour = type,
                 shape = type,
                 size = type)) +
  geom_point(data = ao_core_gps,
             aes(x = lon,
                 y = lat),
             shape = 4,
             size = 1, 
             colour = "black") +
  geom_rect(aes(xmin = 48.45, xmax = 48.516,
                ymin = -13.4895, ymax = -13.4945),
            fill = "white",
            colour = NA) +
  coord_sf(crs = "+proj=longlat +datum=WGS84",
           expand = FALSE) +
  scale_x_continuous(name = "") +
  scale_y_continuous(name = "") +
  annotation_scale(location = "bl",
                   width_hint = 0.25,
                   text_cex = 0.8) +
  scale_colour_manual(name = "Forest status",
                      labels = c("Intact - stumps <15% and canopy cover >60%",
                                 "Degraded - stumps 15-60% and canopy cover 30-60%",
                                 "Deforested - stumps >60% and canopy cover <30%"),
                      values = c(triR_col,
                                 triL_col,
                                 triM_col)) +
  scale_shape_manual(name = "Forest status",
                      values = c(20,
                                 18,
                                 15)) +
  scale_size_manual(name = "Forest status",
                    values = c(2,
                               2,
                               1.5)) +
  guides(colour = guide_legend(title = NULL,
                               override.aes = list(shape = c(20, 18, 15),
                                                   colour = c(triR_col,
                                                              triL_col,
                                                              triM_col),
                                                   size = c(3.1,
                                                            3,
                                                            2.5))),
         shape = "none",
         size = "none") +
  theme_buttress_map() %+replace%
  theme(axis.text.x = element_text(margin = margin(2.5, 0, 0, 0)),
        legend.position = c(0.0075, 1),
        legend.justification = c("left", "top"),
        legend.text = element_text(size = rel(0.75),
                                   margin = margin(0,0,0,-6)),
        legend.key.height = unit(11, "points"),
        legend.background = element_rect(fill = NA,
                                         colour = NA))

inset_map <- ggplot() +
  geom_sf(data = mada_country,
          colour = NA, 
          fill = "grey80") +
  geom_rect(aes(xmin = 48.228, xmax = 48.730,
                ymin = -13.307, ymax = -13.755),
            colour = "black",
            fill = NA,
            size = 1) +
  theme_nothing() %+replace%
  theme(panel.border = element_rect(fill = NA,
                                    colour = "black",
                                    size = 0.5),
        panel.background = element_rect(fill = "white",
                                        colour = NA))

figure_1 <- ggdraw() +
  draw_plot(main_map) +
  draw_plot(inset_map,
            x = 0.64, y = 0.097,
            width = 0.48, height = 0.48) 

agg_png(filename = "figure-output/figure-1.png",
        width = 16, height = 11, 
        units = "cm", res = 1080,
        scaling = 1)
figure_1
dev.off()

#
##### 9. import spatial data #####

# soil characteristics values from Arias-Ortiz et al. 2020 -- 
# see Table 2, values from "Intact (all)" row.
ao_dbd = 0.58 # dry bulk density
ao_cdens = 0.057 # carbon content 

# import NDVI rasters from Blue Ventures team -- 

# "pre-clearing NDVI" was created from the Landsat image with ID
# "LE07_L2SP_159069_20000613_20200918_02_T1",

# and "post-clearing NDVI" was created from the Landsat image with ID
# "LC08_L2SP_159069_20180709_20200831_02_T1",

# both available at www.earthexplorer.usgs.gov 

raw_pre_ndvi <- raster("data/spatial-analysis/NAME OF PRE-CLEARING NDVI FILE")
raw_post_ndvi <- raster("data/spatial-analysis/NAME OF POST-CLEARING NDVI FILE")

# project rasters to correct coordinate reference system -- often takes some
# time, and may fail on first attempt -- re-run before further troubleshooting
proj_pre_ndvi <- projectRaster(raw_pre_ndvi,
                               crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
proj_post_ndvi <- projectRaster(raw_post_ndvi,
                                crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# extract raster values at Blue Ventures sample points
ndvi_dat_raw <- mada_points %>% 
  mutate(pre_ndvi = raster::extract(proj_pre_ndvi, 
                                    mada_points %>% dplyr::select(lon, lat)),
         post_ndvi = raster::extract(proj_post_ndvi, 
                                     mada_points %>% dplyr::select(lon, lat)),
         change_ndvi = post_ndvi - pre_ndvi)

# combine NDVI pre/post values with deltaD values, keeping only species 
# Bruguiera gymnorrhiza (plot_soil_loss is also stumps only)
ndvi_dat <- plot_soil_loss %>% 
  filter(species == "Bg") %>% 
  left_join(., ndvi_dat_raw,
            by = "plotID")

# test if NDVI predicts deltaD
delta_m1 <- lm(loss ~ change_ndvi, 
               data = ndvi_dat)
Anova(delta_m1) # publication Table 4

# limit raster extent to Blue Ventures supplied limits of study area (Tsimipaika 
# Bay) -- import KML
tbay_limits_raw <- st_read("data/spatial-analysis/study-area-limits.kml")
# drop unusued z dimension
tbay_limits <- st_zm(tbay_limits_raw)
# reproject to same ESPG CRS as NDVI rasters
proj_tbay_limits <- st_transform(tbay_limits, "EPSG:4326")

# crop both rasters -- creates rectangle, downsizes raster;
# mask both rasters -- NAs cells outside the limits, no cells 'removed'
mask_pre_ndvi <- mask(crop(proj_pre_ndvi, proj_tbay_limits),
                      tbay_limits)
mask_post_ndvi <- mask(crop(proj_post_ndvi, proj_tbay_limits),
                       tbay_limits)

change_ndvi <- mask_post_ndvi - mask_pre_ndvi

# also mask by the extent of mangrove ecosystem -- using mangrove extent raster
# from Global Mangrove Watch (v2), available at www.doi.org/10.5281/zenodo.5658808
gmw_extent <- raster("data/spatial-analysis/NAME OF GMW MANGROVE EXTENT FILE")

crop_gmw_extent <- crop(gmw_extent, tbay_limits)

# correct CRS of GMW raster using change_ndvi raster as a template
proj_gmw_extent <- projectRaster(from = crop_gmw_extent,
                                 to = change_ndvi)

# mask change_ndvi with GMW raster (removes water and terrestrial forest area)
change_ndvi <- mask(change_ndvi, proj_gmw_extent)

##### 10. analyse spatial data #####

# calculate mean soil elevation loss for each forest type
loss_by_type <- ndvi_dat %>% 
  group_by(type) %>% 
  summarise(mean_loss = mean(loss),
            se_loss = sd(loss)/sqrt(n()))

# extract mean loss for intact forests 
mean_intact_loss = loss_by_type %>% 
  filter(type == "intact") %>% 
  pull(mean_loss)

# create deltaD raster using model formula
deltad_raster <- ((delta_m1$coefficients[[2]] * change_ndvi) + 
                    delta_m1$coefficients[[1]]) - mean_intact_loss
  
# calculate sum of deltaD across T Bay
deltad_sum <- cellStats(deltad_raster, stat = sum,
                        na.rm = TRUE)

# convert raster to tibble for plotting
deltad_dat <- as_tibble(as.data.frame(deltad_raster, xy = TRUE))

# convert raster::RasterLayer to terra::SpatRaster for area calcs
deltad_spat_raster <- rast(deltad_raster)

# calculate size of each raster cell in m
cell_sizes_m <- cellSize(deltad_spat_raster,
                         mask = TRUE,
                         unit = "m")

# convert cell sizes in m to cm (consistent with AO per cm values)
cell_sizes <- cell_sizes_m * 10000

# multiply deltaD values by cell area = cubic soil loss per cell
cubic_loss <- deltad_spat_raster * cell_sizes

# multiply cubic loss by AO density = gram soil loss per cell
gram_loss <- cubic_loss * ao_dbd

# multiply gram loss by carbon density = gram C loss per cell
c_loss <- gram_loss * ao_cdens

# make a reference column of cell areas for per ha calculation
cell_sizes_dat <- as_tibble(cell_sizes) %>% 
  rename(area_cm = area) %>% 
  # convert cm2 to ha
  mutate(area_ha = area_cm / 100000000)

c_loss_dat <- as_tibble(as.data.frame(c_loss, xy = TRUE)) %>% 
  bind_cols(., cell_sizes_dat) %>% 
  rename(g_loss = layer) %>% 
  mutate(mg_loss = g_loss / 1000000, # convert g to Mg C (Mg C = tonne)
         # scale Mg C per cell area up to a Mg C per ha loss value
         mg_ha_loss = 1 / area_ha * mg_loss)

# summarise C loss across T Bay
c_loss_dat %>% 
  summarise(total_raster_loss = sum(mg_loss),
            mean_cell_loss = mean(mg_loss),
            se_cell_loss = sd(mg_loss)/n(),
            mean_ha_loss = mean(mg_ha_loss),
            se_ha_loss = sd(mg_ha_loss)/n())

# size of study area 
expanse(c_loss,
        unit = "ha")
#
##### 11. visualise spatial data #####

spat_main_map <- ggplot() +
  geom_sf(data = crop_mada_provinces,
          colour = NA,
          fill = "grey80") +
  geom_tile(data = c_loss_dat,
            aes(x = x, y = y,
                fill = mg_ha_loss)) +
  scale_fill_gradient(name = NULL,
                      low = "#4daf4a",
                      high = "#e41a1c",
                      na.value = "white") +
  coord_sf(crs = "+proj=longlat +datum=WGS84",
           expand = FALSE) +
  scale_x_continuous(name = "") +
  scale_y_continuous(name = "") +
  annotation_scale(location = "bl",
                   width_hint = 0.25,
                   text_cex = 0.8) +
  guides(fill = guide_colourbar(title = expression(paste(
    "Soil C loss (Mg ", ha^-1, ")")),
    title.position = "top",
    direction = "horizontal")) +
  theme_buttress_map() %+replace%
  theme(axis.text.x = element_text(margin = margin(2.5, 0, 0, 0)),
        legend.position = c(0.0075, 1),
        legend.justification = c("left", "top"),
        legend.text = element_text(size = rel(0.75),
                                   margin = margin(0,0,0,-6)),
        legend.key.height = unit(11, "points"),
        legend.background = element_rect(fill = NA,
                                         colour = NA))

figure_5 <- ggdraw() +
  draw_plot(spat_main_map) +
  draw_plot(inset_map,
            x = 0.64, y = 0.097,
            width = 0.48, height = 0.48) 

agg_png(filename = "figure-output/figure-5.png",
        width = 16, height = 11, 
        units = "cm", res = 1080,
        scaling = 1)
figure_5
dev.off()

#
