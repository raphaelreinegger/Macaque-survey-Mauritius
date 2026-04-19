# Clear plots
if(!is.null(dev.list())) dev.off()
# Clear console                                                           
cat("\014") 
# Clean workspace
rm(list=ls())
# Set working directory
setwd("C:/path")

#load relevant packages for MaxEnt
library(dismo)
library(sf) #need to install Rtools beforehand on https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html
library(rJava)
#if loading sf fails, also install package 'proxy' (also requires RTools)
library(terra)
#load relevant packages for GLMM
library(spatstat.geom) 
library(ggplot2)
library(glmmTMB)
library(emmeans)
library(dplyr)
library(DHARMa)
library(pROC)
library(cvAUC)
library(modEvA)
library(performance)
library(ggeffects)
library(spaMM)
library(ROI.plugin.glpk)
library(ggpubr)
library(igraph)
library(spThin)

##load macaque occurrence data
macaque_presence <- read.csv("GIS/MaxEnt/Macaque_presence.csv",na.strings=c("", "NA"),header=T,sep=",") #macaque presences
macaque_presence <- macaque_presence[,-c(1:11)] #need only coordinates for maxent
#thin records, keep only those more than 500 m apart (~1 record per sampling plot)
pts <- st_as_sf(macaque_presence, coords = c("x_32740", "y_32740"), crs = 32740)
neighbors <- st_is_within_distance(pts, pts, dist = 500)
keep <- rep(TRUE, nrow(pts))
for (i in seq_len(nrow(pts))) {
  if (keep[i]) {
    idx <- neighbors[[i]]
    idx <- idx[idx > i]
    keep[idx] <- FALSE
  }
}
pts_thinned <- pts[keep, ]
#convert back to dataframe with x and y coordinates
coords <- st_coordinates(pts_thinned)
macaque_presence <- data.frame( x = coords[,1], y = coords[,2])

#################################################################
#################################################################
#################################################################
########                                                 ########
######## Prepare proportion and distance data for MaxEnt ########   
########                                                 ########
#################################################################
#################################################################
#################################################################
#load all shapefiles
layers <- list(
  forest     = vect(st_read("GIS/MaxEnt/Forest_clusters.shp")),
  scrub      = vect(st_read("GIS/MaxEnt/Scrub_clusters.shp")),
  sparse     = vect(st_read("GIS/MaxEnt/Sparse_clusters.shp")),
  barren     = vect(st_read("GIS/MaxEnt/Barren_polygons.shp")),
  waterbody  = vect(st_read("GIS/MaxEnt/Waterbody_polygon.shp")),
  wetmarshy  = vect(st_read("GIS/MaxEnt/Wet-marshy_polygon.shp")),
  river      = vect(st_read("GIS/River_polygon.shp")),
  urban      = vect(st_read("GIS/MaxEnt/Urban_polygon.shp")),
  agri       = vect(st_read("GIS/Agri_polygon.shp"))
)
#create 10 m template raster
all_vect <- Reduce(c, layers)
r_template <- rast( ext(all_vect), resolution = 10, crs = crs(all_vect))
#rasterize (binary: 1 = presence, 0 = absence)
rasters_10m <- lapply(layers, function(x) {
  rasterize(x, r_template, field = 1, background = 0)
})
#create % cover rasters (320 by 320 m = ~10 ha resolution)
fact <- 32  # 320 × 320 m → ~100,000 m2
cover_10ha <- lapply(rasters_10m, function(r) {
  aggregate(r, fact = fact, fun = mean)
})
#create distance rasters again (10 by 10 m resolution)
dist_10m <- lapply(rasters_10m, function(r) {
  r[r == 0] <- NA   #only habitat cells remain
  distance(r)
})
#resample distance rasters from 100 m2 to 25 ha
dist_10ha <- mapply(function(d, ref) {
  resample(d, ref, method = "bilinear")
}, dist_10m, cover_10ha, SIMPLIFY = FALSE)
#load & align mauritius mask
MU_outline <- rast("GIS/MaxEnt/MU_outline.tif")
#Align mask to predictors
MU_outline <- resample(MU_outline, cover_10ha[[1]], method = "near")
#Ensure outside = NA
MU_outline[MU_outline == 0] <- NA
#apply mask
cover_10ha <- lapply(cover_10ha, function(r) mask(r, MU_outline))
dist_10ha  <- lapply(dist_10ha,  function(r) mask(r, MU_outline))
#stack predictors
predictors_stack <- rast(c(cover_10ha, dist_10ha))
names(predictors_stack) <- c(
  paste0(names(cover_10ha), "_cover"),
  paste0(names(dist_10ha), "_dist")
)
#and finally correlation analysis
vals <- values(predictors_stack, na.rm = TRUE)
cor_matrix <- cor(vals, method = "spearman")
print(cor_matrix) #high correlations for everything: forest (-0.95), scrub (-0.73), sparse (-0.91), agri (-0.97),
#river (-0.89), urban (-0.91), waterbody (-0.74), barren (-0.82), wetmarshy (-0.43)

#Export some %cover rasters
selected_layers <- c("forest_cover", "scrub_cover", "sparse_cover")
for (nm in selected_layers) {
  writeRaster(
    predictors_stack[[nm]],
    filename = paste0("GIS/MaxEnt/", nm, ".tif"),
    overwrite = TRUE
  )
}

#now lets check the maxent analysis with % cover rasters
#now stack for analysis
cover_stack <- predictors_stack[[grep("_cover", names(predictors_stack))]]
stack_cover <- raster::stack(as(cover_stack, "Raster"))
plot(stack_cover) 
#lets make a MaxEnt model for macaques
set.seed(123)
me1 <- maxent(stack_cover, macaque_presence, args=c('responsecurves=true'))
me1 #directs you to the browser to show you detailed model output, including response curves
plot(me1) 

#do the same for the distance rasters (don't use from predictors_stack, cause it has 100 m resolution rasters)
dist_stack_10m <- rast(dist_10m)
MU_outline <- resample(MU_outline, dist_stack_10m[[1]], method = "near")
MU_outline[MU_outline == 0] <- NA
dist_stack_10m <- mask(dist_stack_10m, MU_outline)
stack_dist <- stack(as(dist_stack_10m, "Raster"))
plot(stack_dist) #looks good
#lets make a MaxEnt model for macaques
set.seed(123)
me2 <- maxent(stack_dist, macaque_presence, args=c('responsecurves=true'))
me2 #directs you to the browser to show you detailed model output, including response curves
plot(me2) 

#get AUC with confidence intervals for both models
#its easiest and quickest using cvAUC with five-fold cross-validation 
#first start with % cover model
set.seed(123)
bg1 <- randomPoints(stack_cover, 10000, p = macaque_presence)
me1 <- maxent(stack_cover, macaque_presence)
pred_raster1 <- predict(me1, stack_cover)
pres_pred1 <- raster::extract(pred_raster1, macaque_presence)
bg_pred1   <- raster::extract(pred_raster1, bg1)
#now set up cross-validation
preds1  <- c(pres_pred1, bg_pred1)
labels1 <- c(rep(1, length(pres_pred1)),
            rep(0, length(bg_pred1)))
keep1 <- !is.na(preds1) #remove NA predictions
preds1  <- preds1[keep1]
labels1 <- labels1[keep1]
#split data into equal folds
set.seed(123)
folds1 <- dismo::kfold(preds1, k = 5)
#run ci.cvAUC and get confidence intervals
cv1 <- ci.cvAUC(predictions = preds1,labels = labels1,folds = folds1) #0.83 auc with 0.79 - 0.86 lower and upper CI

#now distance model
set.seed(123)
bg2 <- randomPoints(stack_dist, 10000, p = macaque_presence)
me2 <- maxent(stack_dist, macaque_presence)
pred_raster2 <- predict(me2, stack_dist)
pres_pred2 <- raster::extract(pred_raster2, macaque_presence)
bg_pred2   <- raster::extract(pred_raster2, bg2)
#now set up cross-validation
preds2  <- c(pres_pred2, bg_pred2)
labels2 <- c(rep(1, length(pres_pred2)),
             rep(0, length(bg_pred2)))
keep2 <- !is.na(preds2) #remove NA predictions
preds2  <- preds2[keep2]
labels2 <- labels2[keep2]
#split data into equal folds
set.seed(123)
folds2 <- dismo::kfold(preds2, k = 5)
#run ci.cvAUC and get confidence intervals
cv2 <- ci.cvAUC(predictions = preds2,labels = labels2,folds = folds2) #0.85 auc with 0.82 - 0.88 lower and upper CI
#so the distance model performs a little better

#lets make the suitability maps and see how the surfaces compare
#first for % cover
set.seed(123)
me1 <- predict(me1, stack_cover, args="outputformat=raw", sep="", progress="text")
#now export the map
writeRaster(me1, "GIS/MaxEnt/me_cover_prediction.tif",overwrite = TRUE)
#extract the suitability threshold as the 10th percentile training presence
pres_vals1 <- raster::extract(me1, macaque_presence)
quantile(pres_vals1, probs = 0.10, na.rm = TRUE) #0.0001083471 
#now in qgis, enter the following into the raster calculator for me_cover_prediction.tif: "me_cover_prediction@1" >= 0.0001083471
#convert the raster calculation raster to polygon here in R, because polygonize in QGIS has problems
Binary_suitability <- raster('GIS/MaxEnt/Binary_suitability_cover.tif')
bin_terra <- rast(Binary_suitability) 
polygons <- as.polygons(bin_terra)
writeVector(polygons, "GIS/MaxEnt/Binary_suitability_cover.shp", overwrite=TRUE)

#now for distance
set.seed(123)
me2 <- predict(me2, stack_dist, args="outputformat=raw", sep="", progress="text")
#now export the map
writeRaster(me2, "GIS/MaxEnt/me_dist_prediction.tif",overwrite = TRUE)
#extract the suitability threshold as the 10th percentile training presence
pres_vals2 <- raster::extract(me2, macaque_presence)
quantile(pres_vals2, probs = 0.10, na.rm = TRUE) #0.0001006692 
#now in qgis, enter the following into the raster calculator for me_dist_prediction.tif: "me_cover_prediction@1" >= 0.0001006692
#convert the raster calculation raster to polygon here in R, because polygonize in QGIS has problems
Binary_suitability <- raster('GIS/MaxEnt/Binary_suitability_dist.tif')
bin_terra <- rast(Binary_suitability) 
polygons <- as.polygons(bin_terra)
writeVector(polygons, "GIS/MaxEnt/Binary_suitability_dist.shp", overwrite=TRUE)

#the distance suitability has 608 km2 and the cover suitability has 675 km2



#####################################
#####################################
#####                            ####
#####  optional bias correction  ####
#####                            ####
#####################################
#####################################
#optional, add bias correction file based on where the sampling effort is concentrated
#not really necessary though, because sampling effort is quite proportional to forest region size
colnames(macaque_presence) <- c("x", "y")
pts <- st_as_sf(macaque_presence, coords = c("x", "y"), crs = 32740)
#Create bias raster with distance stack as template
r_template <- rast(stack_dist)
#rasterize presence points (counts per cell)
presence_r <- rasterize(vect(pts), r_template, fun = "count", background = 0)
#smooth (≈ 1 km window → 101 cells at 10 m)
bias_r <- focal(presence_r, w = 101, fun = mean, na.rm = TRUE)
#normalize (0–1)
max_val <- global(bias_r, "max", na.rm = TRUE)[1,1]
bias_r <- bias_r / max_val
#mask to study area (important!)
bias_r <- mask(bias_r, r_template)
#Sample background points (bias-weighted)
set.seed(123)
bg <- spatSample(bias_r, size = 10000, method = "weights", na.rm = TRUE, xy = TRUE)
#convert to matrix/data.frame for MaxEnt
bg <- as.data.frame(bg)[, c("x", "y")]
#Fit MaxEnt model
set.seed(123)
me_dist_bias <- maxent(stack_dist, p = macaque_presence, a = bg)
#now make binary suitability map
set.seed(123)
me3 <- predict(me_dist_bias, stack_dist, args="outputformat=raw", sep="", progress="text")
#now export the map
writeRaster(me3, "GIS/MaxEnt/me_dist_prediction_bias.tif",overwrite = TRUE)
#extract the suitability threshold as the 10th percentile training presence
pres_vals3 <- raster::extract(me3, macaque_presence)
quantile(pres_vals3, probs = 0.10, na.rm = TRUE) #0.00007408291 
#now in qgis, enter the following into the raster calculator for me_dist_prediction.tif: "me_cover_prediction@1" >= 0.00007408291
#convert the raster calculation raster to polygon here in R, because polygonize in QGIS has problems
Binary_suitability <- raster('GIS/MaxEnt/Binary_suitability_dist_bias.tif')
bin_terra <- rast(Binary_suitability) 
polygons <- as.polygons(bin_terra)
writeVector(polygons, "GIS/MaxEnt/Binary_suitability_dist_bias.shp", overwrite=TRUE)

#really does not make a difference, now habitat is 619 instead of 608 km2