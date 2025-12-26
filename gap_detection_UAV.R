# gap_detection_UAV.R
# Complete R code for canopy gap detection using DSM change analysis
# Simonetti et al. <paper name and journal> 2026

# Load required packages
library(terra)
library(sf)
library(dplyr)

# ---------------------------
# CONFIGURATION PARAMETERS
# ---------------------------
# Input file paths (update these to your data)
dem_older_path <- "path/to/older_dem.tif"        # e.g., 20210831 DEM
dem_newer_path <- "path/to/newer_dem.tif"        # e.g., 20210928 DEM  
plot_shp_path <- "path/to/plot.shp"              # shapefile or geopackage
snap_raster_path <- "path/to/reference_1m.tif"   # Reference raster for resampling

# Output directory
output_dir <- "output/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Parameters
buffer_m <- 25          # Buffer around plots (meters)
target_res_m <- 1       # Target resolution (meters)
focal_window <- 99      # Focal median window size (cells)
height_threshold_m <- -5 # Gap threshold (newer - older < threshold)
min_gap_area_m2 <- 5    # Minimum gap area (m²)
min_area_perim_ratio <- 0.6  # Minimum area/perimeter ratio

# Derived filenames (FIXED: added missing flat_difference)
intermediate_files <- list(
  clipped_older = file.path(output_dir, "01_older_dem_clipped.tif"),
  resampled_older = file.path(output_dir, "02_older_dem_1m.tif"),
  dem_difference = file.path(output_dir, "03_dem_diff_newer_minus_older.tif"),
  focal_median = file.path(output_dir, "04_diff_median_focal.tif"),
  flat_difference = file.path(output_dir, "05_flat_difference.tif"),  # ← ADDED
  threshold_mask = file.path(output_dir, sprintf("06_mask_lt_%sm.tif", abs(height_threshold_m))),
  raw_polygons = file.path(output_dir, "07_raw_gaps.shp"),
  selected_gaps = file.path(output_dir, "08_selected_gaps.shp")
)

cat("Pipeline started. Output directory:", output_dir, "\n")

# ---------------------------
# 1. LOAD AND PREPARE DATA
# ---------------------------
cat("1. Loading input data...\n")

# Load rasters and shapefile
dem_older <- rast(dem_older_path)
dem_newer <- rast(dem_newer_path)
plot_shp <- st_read(plot_shp_path, quiet = TRUE)
reference_raster <- rast(snap_raster_path)

# Ensure consistent CRS
target_crs <- crs(reference_raster)
if (st_crs(plot_shp) != target_crs) {
  plot_shp <- st_transform(plot_shp, target_crs)
}

cat("Data loaded. CRS:", target_crs, "\n")

# ---------------------------
# 2. CREATE PLOT BUFFER (25m) 
# ---------------------------
cat("2. Creating plot buffer...\n")

plot_union <- st_union(plot_shp) |> st_make_valid()
plot_buffer <- st_buffer(plot_union, dist = buffer_m)  

# Optional: Save buffer for inspection
# st_write(plot_buffer, file.path(output_dir, "plot_buffer_25m.shp"), 
#          delete_layer = TRUE, quiet = TRUE)

# ---------------------------
# 3. CLIP OLDER DEM TO BUFFERED PLOT 
# ---------------------------
cat("3. Clipping older DEM...\n")

plot_vect <- vect(plot_buffer)  

dem_older_clipped <- crop(dem_older, plot_vect) |> mask(plot_vect)
writeRaster(dem_older_clipped, intermediate_files$clipped_older, overwrite = TRUE)

# ---------------------------
# 4. RESAMPLE OLDER DEM TO MATCH REFERENCE
# ---------------------------
cat("4. Resampling older DEM to 1m...\n")

if (crs(dem_older_clipped) != target_crs) {
  dem_older_clipped <- project(dem_older_clipped, target_crs)
}
dem_older_resampled <- resample(dem_older_clipped, reference_raster, method = "bilinear")
writeRaster(dem_older_resampled, intermediate_files$resampled_older, overwrite = TRUE)

# ---------------------------
# 5. CALCULATE DEM DIFFERENCE (image_date_newer minus image_date_older)
# ---------------------------
cat("5. Calculating DEM difference...\n")

dem_newer_resampled <- resample(dem_newer, dem_older_resampled, method = "bilinear")
dem_diff <- dem_newer_resampled - dem_older_resampled

writeRaster(dem_diff, intermediate_files$dem_difference, overwrite = TRUE)
cat("Difference stats:", paste(round(summary(values(dem_diff)), 2), collapse = " | "), "\n")

# ---------------------------
# 6. FOCAL MEDIAN FILTER (99x99)
# ---------------------------
cat("6. Applying focal median filter...\n")

focal_matrix <- matrix(1, nrow = focal_window, ncol = focal_window)
dem_diff_median <- focal(dem_diff, w = focal_matrix, fun = median, na.rm = TRUE)
writeRaster(dem_diff_median, intermediate_files$focal_median, overwrite = TRUE)

# ---------------------------
# 7. FLAT DIFFERENCE (REMOVE SLOW WARPING)
# ---------------------------
cat("7. Correcting vertical warping...\n")

dem_diff_flat <- dem_diff - dem_diff_median
writeRaster(dem_diff_flat, intermediate_files$flat_difference, overwrite = TRUE)  

# ---------------------------
# 8. THRESHOLD AND PATCHES (< -5m)
# ---------------------------
cat("8. Applying height threshold...\n")

gap_mask <- ifel(dem_diff_flat < height_threshold_m, 1, NA)
gap_patches <- patches(gap_mask, directions = 8, filename = intermediate_files$threshold_mask, overwrite = TRUE)

# ---------------------------
# 9. POLYGONIZE GAPS
# ---------------------------
cat("9. Converting patches to polygons...\n")

raw_polygons <- as.polygons(gap_patches, dissolve = TRUE)
raw_polygons_sf <- st_as_sf(raw_polygons) |> 
  filter(!is.na(patches)) |>
  st_make_valid()

st_write(raw_polygons_sf, intermediate_files$raw_polygons, delete_layer = TRUE, quiet = TRUE)

# ---------------------------
# 10. FILTER GAPS BY GEOMETRY - FIXED: MULTILINESTRING
# ---------------------------
cat("10. Filtering gaps by area and shape...\n")

raw_polygons_sf <- st_set_crs(raw_polygons_sf, target_crs)

raw_polygons_sf <- raw_polygons_sf |>
  mutate(
    area_m2 = as.numeric(st_area(.)),
    perimeter_m = as.numeric(st_length(st_cast(., "MULTILINESTRING"))),  
    area_perim_ratio = area_m2 / perimeter_m
  ) |>
  filter(area_m2 >= min_gap_area_m2, area_perim_ratio >= min_area_perim_ratio)

st_write(raw_polygons_sf, intermediate_files$selected_gaps, delete_layer = TRUE, quiet = TRUE)

# ---------------------------
# SUMMARY 
# ---------------------------
cat("\n=== PIPELINE COMPLETED ===\n")
cat("Total raw gaps found:", nrow(raw_polygons_sf |> filter(area_m2 < min_gap_area_m2 | area_perim_ratio < min_area_perim_ratio)), "\n")  
cat("Selected gaps (area >", min_gap_area_m2, "m², ratio >", min_area_perim_ratio, "):", nrow(raw_polygons_sf), "\n")
cat("Total gap area:", round(sum(raw_polygons_sf$area_m2, na.rm = TRUE), 1), "m²\n")
cat("\nKey outputs:\n")
for (name in names(intermediate_files)) {
  cat(sprintf("  %s: %s\n", name, intermediate_files[[name]]))
}
cat("Final selected gaps:", intermediate_files$selected_gaps, "\n")

cat("\nVisualize results:\n")
cat("plot(rast('", intermediate_files$selected_gaps, "'))\n")
