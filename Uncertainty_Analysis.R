
suppressPackageStartupMessages({
  library(raster)
  library(terra)
  library(sf)
  library(tmap)
  library(openxlsx)
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(cowplot)
  library(grid)
  library(gridExtra)
})

# Base directory paths
dir_base    <- "E:/Modelling"
dir_input   <- file.path(dir_base, "output_maps")          
dir_output  <- file.path(dir_base, "output_maps")    
dir_figures <- file.path(dir_base, "Figures")
dir_data    <- file.path(dir_base, "Data")
dir_maps    <- file.path(dir_base, "maps")

if (!dir.exists(dir_output))  dir.create(dir_output, recursive = TRUE)
if (!dir.exists(dir_figures)) dir.create(dir_figures, recursive = TRUE)
if (!dir.exists(dir_data))    dir.create(dir_data, recursive = TRUE)

# Species definitions
species_keys <- c("krempfii", "dalatensis", "cernua", "fenzeliana", "latteri", "henryi")
species_names <- c("Pinus krempfii", "Pinus dalatensis", "Pinus cernua", 
                   "Pinus fenzeliana", "Pinus latteri", "Pinus henryi")
names(species_names) <- species_keys

panel_letters <- c("krempfii" = "(a)", "dalatensis" = "(b)", "cernua" = "(c)",
                   "fenzeliana" = "(d)", "latteri" = "(e)", "henryi" = "(f)")

# GCMs & Climate Scenarios
gcms <- c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")
ssps <- c("ssp126", "ssp585")

# ------------------------------------------------------------------------------
# LOAD GEOGRAPHICAL BASEMAPS & BOUNDARIES
# ------------------------------------------------------------------------------
cat(">> Loading administrative shapefiles and Vietnam raster mask...\n")
vn_mask_path <- file.path(dir_maps, "VN_raster.tif")
vn_raster    <- raster::raster(vn_mask_path)

crop_bbox  <- sf::st_bbox(c(xmin = 98, ymin = 6, xmax = 118, ymax = 26), crs = sf::st_crs(4326))
tight_bbox <- sf::st_bbox(c(xmin = 102, ymin = 6, xmax = 113, ymax = 24), crs = sf::st_crs(4326))

suppressWarnings({
  load_and_prep <- function(shp_path) {
    x <- sf::st_read(shp_path, quiet = TRUE)
    if (is.na(sf::st_crs(x))) {
      sf::st_crs(x) <- 4326
    } else {
      x <- sf::st_transform(x, 4326)
    }
    return(sf::st_crop(x, crop_bbox))
  }
  
  vietnam       <- load_and_prep(file.path(dir_maps, "VN.shp"))
  china         <- load_and_prep(file.path(dir_maps, "china.shp"))
  laos          <- load_and_prep(file.path(dir_maps, "laos.shp"))
  cambodia      <- load_and_prep(file.path(dir_maps, "cambodia.shp"))
  pas           <- load_and_prep(file.path(dir_maps, "VN_pas.shp"))
  vn_geographic <- load_and_prep(file.path(dir_maps, "vn_geographic.shp"))
  
  # For Truongsa & Hoangsa, transform to 4326
  truongsa      <- sf::st_transform(sf::st_read(file.path(dir_maps, "Truongsa_polyline.shp"), quiet = TRUE), 4326)
  hoangsa       <- sf::st_transform(sf::st_read(file.path(dir_maps, "Hoangsa_polyline.shp"), quiet = TRUE), 4326)
})

# Helper function to find existing raster with flexible case-insensitivity
find_raster_file <- function(base_dir, pattern) {
  matched <- list.files(base_dir, pattern = paste0("^", pattern, "$"), ignore.case = TRUE, full.names = TRUE)
  if (length(matched) > 0) return(matched[1])
  return(file.path(base_dir, pattern))
}

# ==============================================================================
# SECTION 1: GCM AGREEMENT MAPS (0 TO 5 GCMS)
# ==============================================================================
cat("\n==============================================================================\n")
cat("SECTION 1: Computing GCM Agreement Rasters and Maps (0 to 5 GCMs)\n")
cat("==============================================================================\n")

agreement_pal <- c("#FFFFFF", "#B3DDF2", "#4A90E2", "#F5D76E", "#F39C12", "#C0392B")
agreement_labels <- c("0 (Unsuitable)", "1 GCM", "2 GCMs", "3 GCMs", "4 GCMs", "5 GCMs (High Consensus)")

for (sp in species_keys) {
  sp_title <- species_names[sp]
  cat(sprintf("\n--- Processing GCM Agreement for %s ---\n", sp_title))
  
  for (ssp in ssps) {
    ssp_label <- ifelse(ssp == "ssp126", "SSP1-2.6", "SSP5-8.5")
    
    # Locate and load 5 GCM binary prediction rasters (from original output_maps)
    gcm_raster_list <- list()
    for (gcm in gcms) {
      gcm_pattern <- paste0("pinus_", sp, "_", gcm, "_", ssp, "\\.tif")
      gcm_file <- find_raster_file(dir_input, gcm_pattern)
      
      if (!file.exists(gcm_file)) {
        warning(sprintf("GCM binary raster not found: %s", gcm_file))
      } else {
        gcm_raster_list[[gcm]] <- raster::raster(gcm_file)
      }
    }
    
    if (length(gcm_raster_list) == length(gcms)) {
      gcm_stack <- raster::stack(gcm_raster_list)
      agreement_r <- raster::calc(gcm_stack, fun = sum, na.rm = TRUE)
      agreement_r <- agreement_r * vn_raster
      
      out_tif_name <- paste0("pinus_", sp, "_", ssp, "_gcm_agreement_26_8.tif")
      out_tif_path <- file.path(dir_output, out_tif_name)
      raster::writeRaster(agreement_r, filename = out_tif_path, format = "GTiff", overwrite = TRUE)
      cat(sprintf("     Saved raster: %s\n", out_tif_name))
      
      # Generate individual map figure
      map_agree <- tm_shape(vietnam, bbox = tight_bbox) +
        tm_polygons(border.col = "black", col = "grey90", alpha = 1, lwd = 0.7, border.alpha = 0.8) +
        tm_shape(agreement_r) +
        tm_raster(palette = agreement_pal, style = "cat", labels = c("0","1","2","3","4","5"),
                  title = "GCM Agreement", legend.show = TRUE) +
        tm_shape(china) + tm_polygons(border.col = "black", alpha = 0, lwd = 0.7, border.alpha = 0.8) +
        tm_shape(laos) + tm_polygons(border.col = "black", alpha = 0, lwd = 0.7, border.alpha = 0.8) +
        tm_shape(cambodia) + tm_polygons(border.col = "black", alpha = 0, lwd = 0.7, border.alpha = 0.8) +
        tm_shape(pas) + tm_polygons(border.col = "darkgreen", alpha = 0, lwd = 0.6, border.alpha = 0.8) +
        tm_shape(vn_geographic) + tm_polygons(border.col = "black", alpha = 0, lwd = 0.7, border.alpha = 0.8) +
        tm_shape(truongsa) + tm_lines(col = "grey20", lwd = 0.7) +
        tm_shape(hoangsa) + tm_lines(col = "grey20", lwd = 0.7) +
        tm_layout(title = paste0(sp_title, " - ", ssp_label),
                  title.position = c("left", "top"), title.size = 0.9, frame = FALSE,
                  legend.position = c("left", "bottom"), legend.text.size = 0.75, legend.title.size = 0.85,
                  legend.bg.color = "white", legend.bg.alpha = 0.6)
      
      out_jpg_path <- file.path(dir_figures, paste0("pinus_", sp, "_", ssp, "_gcm_agreement_26_8.jpg"))
      tmap::tmap_save(map_agree, filename = out_jpg_path, width = 10, height = 15, units = "cm", dpi = 300)
    }
  }
}

# ==============================================================================
# SECTION 2: HABITAT DYNAMICS CLASSIFICATION (4-CLASS DYNAMICS + UNCERTAINTY)
# ==============================================================================
cat("\n==============================================================================\n")
cat("SECTION 2: Classifying Habitat Dynamics (Unsuitable, Stable, Loss, Gain, Uncertain)\n")
cat("==============================================================================\n")

classify_habitat_dynamics <- function(curr, agr) {
  na_mask <- is.na(curr) | is.na(agr)
  res <- ifelse(na_mask, NA,
         ifelse(agr %in% c(1, 2), 4,      # 4 = Uncertain (1 to 2 GCMs: low model consensus / uncertainty)
         ifelse(curr == 1 & agr >= 3, 1,  # 1 = Stable Refugia (>= 3 GCMs majority consensus)
         ifelse(curr == 1 & agr == 0, 2,  # 2 = Loss (0 GCMs: complete loss across all GCMs)
         ifelse(curr == 0 & agr >= 3, 3,  # 3 = Gain (>= 3 GCMs majority consensus)
         0)))))                           # 0 = Unsuitable (curr == 0 & agr == 0)
  return(res)
}

for (sp in species_keys) {
  sp_title <- species_names[sp]
  cat(sprintf("\n--- Classifying Habitat Dynamics for %s ---\n", sp_title))
  
  curr_file <- find_raster_file(dir_input, paste0("pinus_", sp, "_current\\.tif"))
  if (!file.exists(curr_file)) {
    warning(sprintf("Current binary raster not found: %s", curr_file))
    next
  }
  current_r <- raster::raster(curr_file) * vn_raster
  
  for (ssp in ssps) {
    agree_file <- file.path(dir_output, paste0("pinus_", sp, "_", ssp, "_gcm_agreement_26_8.tif"))
    if (!file.exists(agree_file)) next
    
    agreement_r <- raster::raster(agree_file)
    dynamics_r <- raster::overlay(current_r, agreement_r, fun = classify_habitat_dynamics)
    dynamics_r <- dynamics_r * vn_raster
    
    out_dyn_name <- paste0("pinus_", sp, "_", ssp, "_dynamics_26_8.tif")
    out_dyn_path <- file.path(dir_output, out_dyn_name)
    raster::writeRaster(dynamics_r, filename = out_dyn_path, format = "GTiff", overwrite = TRUE)
    cat(sprintf("  -> Saved dynamics raster: %s\n", out_dyn_name))
  }
}

