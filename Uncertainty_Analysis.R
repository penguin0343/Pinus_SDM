# ==============================================================================
# Script: 02_Uncertainty_Analysis_All_Species_26_8.R
# Purpose: Comprehensive Spatial Uncertainty & Habitat Dynamics Analysis across ALL 6 Pinus Species
#          (PLOS ONE Revision - Preserves original model outputs)
# Author: SDM Research Team
# Date: 2026-08-26
#
# Target Species (All 6 Pinus in Vietnam):
#   1. Pinus krempfii
#   2. Pinus dalatensis
#   3. Pinus cernua
#   4. Pinus fenzeliana
#   5. Pinus latteri
#   6. Pinus henryi
#
# Analyses Included:
#   SECTION 1: GCM Agreement / Consensus Maps (0 to 5 GCMs)
#   SECTION 2: Habitat Dynamics Classification (Unsuitable, Stable Refugia, Loss, Gain, Uncertain)
#   SECTION 3: Quantitative Latitudinally Corrected Area Calculations (km2 & % change)
#   SECTION 4: Individual and Multi-Panel High-Resolution Figure Generation (tmap)
# ==============================================================================

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

cat("==============================================================================\n")
cat("Starting Spatial Uncertainty & Habitat Dynamics Analysis for ALL 6 Pinus Species\n")
cat("==============================================================================\n\n")

# Base directory paths
dir_base    <- "E:/OneDrive/2024/DTCSCL_2024/Modelling"
dir_input   <- file.path(dir_base, "output_maps")          # Original outputs (read-only)
dir_output  <- file.path(dir_base, "output_maps_26_8")     # Revision outputs
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

# ==============================================================================
# SECTION 3: QUANTITATIVE LATITUDINALLY CORRECTED AREA CALCULATIONS
# ==============================================================================
cat("\n==============================================================================\n")
cat("SECTION 3: Calculating Latitudinally Corrected Habitat Areas (km2) and Dynamics\n")
cat("==============================================================================\n")

area_summary_list <- list()

for (sp in species_keys) {
  sp_title <- species_names[sp]
  cat(sprintf("  Calculating area statistics for: %s\n", sp_title))
  
  curr_file <- find_raster_file(dir_input, paste0("pinus_", sp, "_current\\.tif"))
  if (!file.exists(curr_file)) next
  
  current_r <- raster::raster(curr_file) * vn_raster
  area_grid_km2 <- raster::area(current_r) * vn_raster
  
  current_km2 <- raster::cellStats(area_grid_km2 * (current_r == 1), stat = "sum", na.rm = TRUE)
  
  # SSP1-2.6
  dyn_126_file <- file.path(dir_output, paste0("pinus_", sp, "_ssp126_dynamics_26_8.tif"))
  if (file.exists(dyn_126_file)) {
    d126 <- raster::raster(dyn_126_file)
    s126_stable <- raster::cellStats(area_grid_km2 * (d126 == 1), stat = "sum", na.rm = TRUE)
    s126_loss   <- raster::cellStats(area_grid_km2 * (d126 == 2), stat = "sum", na.rm = TRUE)
    s126_gain   <- raster::cellStats(area_grid_km2 * (d126 == 3), stat = "sum", na.rm = TRUE)
    s126_uncert <- raster::cellStats(area_grid_km2 * (d126 == 4), stat = "sum", na.rm = TRUE)
    s126_tot_suit <- s126_stable + s126_gain
    s126_pct_chg <- ifelse(current_km2 > 0, ((s126_tot_suit - current_km2) / current_km2) * 100, NA)
  } else {
    s126_stable <- s126_loss <- s126_gain <- s126_uncert <- s126_tot_suit <- s126_pct_chg <- NA
  }
  
  # SSP5-8.5
  dyn_585_file <- file.path(dir_output, paste0("pinus_", sp, "_ssp585_dynamics_26_8.tif"))
  if (file.exists(dyn_585_file)) {
    d585 <- raster::raster(dyn_585_file)
    s585_stable <- raster::cellStats(area_grid_km2 * (d585 == 1), stat = "sum", na.rm = TRUE)
    s585_loss   <- raster::cellStats(area_grid_km2 * (d585 == 2), stat = "sum", na.rm = TRUE)
    s585_gain   <- raster::cellStats(area_grid_km2 * (d585 == 3), stat = "sum", na.rm = TRUE)
    s585_uncert <- raster::cellStats(area_grid_km2 * (d585 == 4), stat = "sum", na.rm = TRUE)
    s585_tot_suit <- s585_stable + s585_gain
    s585_pct_chg <- ifelse(current_km2 > 0, ((s585_tot_suit - current_km2) / current_km2) * 100, NA)
  } else {
    s585_stable <- s585_loss <- s585_gain <- s585_uncert <- s585_tot_suit <- s585_pct_chg <- NA
  }
  
  area_summary_list[[sp]] <- data.frame(
    Species                  = sp_title,
    Current_Area_km2         = round(current_km2, 2),
    SSP126_Stable_Refugia_km2= round(s126_stable, 2),
    SSP126_Loss_km2          = round(s126_loss, 2),
    SSP126_Gain_km2          = round(s126_gain, 2),
    SSP126_Uncertain_km2     = round(s126_uncert, 2),
    SSP126_Total_Suitable_km2= round(s126_tot_suit, 2),
    SSP126_Net_Change_pct    = round(s126_pct_chg, 2),
    SSP585_Stable_Refugia_km2= round(s585_stable, 2),
    SSP585_Loss_km2          = round(s585_loss, 2),
    SSP585_Gain_km2          = round(s585_gain, 2),
    SSP585_Uncertain_km2     = round(s585_uncert, 2),
    SSP585_Total_Suitable_km2= round(s585_tot_suit, 2),
    SSP585_Net_Change_pct    = round(s585_pct_chg, 2),
    stringsAsFactors         = FALSE
  )
}

all_area_df <- do.call(rbind, area_summary_list)

cat("\n========================================================================================================\n")
cat("                             HABITAT DYNAMICS AREA SUMMARY (ALL 6 SPECIES)                              \n")
cat("========================================================================================================\n")
print(all_area_df, row.names = FALSE)
cat("========================================================================================================\n\n")

# Save master area summary workbook
master_area_xlsx <- file.path(dir_data, "All_6_Species_Habitat_Dynamics_Area_26_8.xlsx")
openxlsx::write.xlsx(all_area_df, file = master_area_xlsx, overwrite = TRUE)
cat(sprintf(">> Master Habitat Dynamics Area table saved to:\n   %s\n\n", master_area_xlsx))
# ==============================================================================
# SECTION 4: COMBINED 6-PANEL HABITAT DYNAMICS & GCM AGREEMENT FIGURES
# ==============================================================================
cat("==============================================================================\n")
cat("SECTION 4: Generating Combined 6-Panel Dynamics & GCM Agreement Figures\n")
cat("==============================================================================\n")

dynamics_pal <- c("0" = "grey90", "1" = "forestgreen", "2" = "firebrick2", "3" = "dodgerblue3", "4" = "gold")
dynamics_labels <- c("Unsuitable", "Stable Refugia", "Loss", "Gain", "Uncertain")

# # Panel generator for Habitat Dynamics (legend hidden per panel, PAS removed, title shifted up)
create_dynamics_panel <- function(sp_key, ssp_code) {
  sp_title <- paste0(panel_letters[sp_key], " ", species_names[sp_key])
  dyn_file <- file.path(dir_output, paste0("pinus_", sp_key, "_", ssp_code, "_dynamics_26_8.tif"))
  if (!file.exists(dyn_file)) return(NULL)
  
  dyn_r <- raster::raster(dyn_file)
  
  p <- tm_shape(vietnam, bbox = tight_bbox) +
    tm_polygons(border.col = "black", col = "grey90", alpha = 1, lwd = 0.6, border.alpha = 0.8) +
    tm_shape(dyn_r) +
    tm_raster(palette = dynamics_pal, style = "cat", labels = dynamics_labels,
              alpha = 1, legend.show = FALSE) +
    tm_shape(china) + tm_polygons(border.col = "black", alpha = 0, lwd = 0.6, border.alpha = 0.8) +
    tm_shape(laos) + tm_polygons(border.col = "black", alpha = 0, lwd = 0.6, border.alpha = 0.8) +
    tm_shape(cambodia) + tm_polygons(border.col = "black", alpha = 0, lwd = 0.6, border.alpha = 0.8) +
    tm_shape(vn_geographic) + tm_polygons(border.col = "black", alpha = 0, lwd = 0.6, border.alpha = 0.8) +
    tm_shape(truongsa) + tm_lines(col = "grey20", lwd = 0.6) +
    tm_shape(hoangsa) + tm_lines(col = "grey20", lwd = 0.6) +
    tm_layout(title = sp_title, title.fontface = "italic",
              title.position = c("left", "top"), title.size = 1.05,
              inner.margins = c(0.01, 0.02, 0.07, 0.02),
              frame = FALSE, legend.show = FALSE)
  return(p)
}

# Bottom Dynamics Legend (shared horizontal)
df_dyn <- data.frame(x = 1:5, y = 1, Dynamics = factor(dynamics_labels, levels = dynamics_labels))
p_dyn <- ggplot(df_dyn, aes(x = x, y = y, fill = Dynamics)) +
  geom_tile(color = "black", linewidth = 0.4) +
  scale_fill_manual(values = setNames(as.character(dynamics_pal), dynamics_labels), name = NULL) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 13, face = "bold", color = "black", margin = margin(r = 18, l = 4)),
    legend.key.size = unit(0.65, "cm"),
    legend.margin = margin(t = 5, b = 5)
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))
leg_dyn_grob <- cowplot::get_legend(p_dyn)

# Assemble 6-panel Dynamics for SSP1-2.6
cat(">> Building 6-panel Dynamics figure for SSP1-2.6...\n")
p126_list <- lapply(species_keys, function(k) create_dynamics_panel(k, "ssp126"))
if (all(!sapply(p126_list, is.null))) {
  g_maps126 <- lapply(p126_list, tmap::tmap_grob)
  g_grid126 <- gridExtra::arrangeGrob(grobs = g_maps126, ncol = 2, nrow = 3)
  fig126 <- gridExtra::arrangeGrob(g_grid126, leg_dyn_grob, ncol = 1, heights = c(16, 1))
  fig126_path <- file.path(dir_figures, "All_species_dynamics_ssp126_26_8.jpg")
  jpeg(fig126_path, width = 30, height = 45, units = "cm", res = 600, quality = 95)
  grid::grid.draw(fig126)
  dev.off()
  cat(sprintf("   Saved SSP1-2.6 6-panel dynamics figure: %s\n", fig126_path))
}

# Assemble 6-panel Dynamics for SSP5-8.5
cat(">> Building 6-panel Dynamics figure for SSP5-8.5...\n")
p585_list <- lapply(species_keys, function(k) create_dynamics_panel(k, "ssp585"))
if (all(!sapply(p585_list, is.null))) {
  g_maps585 <- lapply(p585_list, tmap::tmap_grob)
  g_grid585 <- gridExtra::arrangeGrob(grobs = g_maps585, ncol = 2, nrow = 3)
  fig585 <- gridExtra::arrangeGrob(g_grid585, leg_dyn_grob, ncol = 1, heights = c(16, 1))
  fig585_path <- file.path(dir_figures, "All_species_dynamics_ssp585_26_8.jpg")
  jpeg(fig585_path, width = 30, height = 45, units = "cm", res = 600, quality = 95)
  grid::grid.draw(fig585)
  dev.off()
  cat(sprintf("   Saved SSP5-8.5 6-panel dynamics figure: %s\n", fig585_path))
}

# Panel generator for GCM Agreement (legend hidden per panel, PAS removed, title shifted up)
create_agreement_panel <- function(sp_key, ssp_code) {
  sp_title <- paste0(panel_letters[sp_key], " ", species_names[sp_key])
  agr_file <- file.path(dir_output, paste0("pinus_", sp_key, "_", ssp_code, "_gcm_agreement_26_8.tif"))
  if (!file.exists(agr_file)) return(NULL)
  
  agr_r <- raster::raster(agr_file)
  
  p <- tm_shape(vietnam, bbox = tight_bbox) +
    tm_polygons(border.col = "black", col = "grey90", alpha = 1, lwd = 0.6, border.alpha = 0.8) +
    tm_shape(agr_r) +
    tm_raster(palette = agreement_pal, style = "cat", labels = c("0","1","2","3","4","5"),
              alpha = 1, legend.show = FALSE) +
    tm_shape(china) + tm_polygons(border.col = "black", alpha = 0, lwd = 0.6, border.alpha = 0.8) +
    tm_shape(laos) + tm_polygons(border.col = "black", alpha = 0, lwd = 0.6, border.alpha = 0.8) +
    tm_shape(cambodia) + tm_polygons(border.col = "black", alpha = 0, lwd = 0.6, border.alpha = 0.8) +
    tm_shape(vn_geographic) + tm_polygons(border.col = "black", alpha = 0, lwd = 0.6, border.alpha = 0.8) +
    tm_shape(truongsa) + tm_lines(col = "grey20", lwd = 0.6) +
    tm_shape(hoangsa) + tm_lines(col = "grey20", lwd = 0.6) +
    tm_layout(title = sp_title, title.fontface = "italic",
              title.position = c("left", "top"), title.size = 1.05,
              inner.margins = c(0.01, 0.02, 0.07, 0.02),
              frame = FALSE, legend.show = FALSE)
  return(p)
}

# Bottom GCM Agreement Legend (shared horizontal)
df_agr <- data.frame(x = 1:6, y = 1, Agreement = factor(agreement_labels, levels = agreement_labels))
p_agr <- ggplot(df_agr, aes(x = x, y = y, fill = Agreement)) +
  geom_tile(color = "black", linewidth = 0.4) +
  scale_fill_manual(values = setNames(agreement_pal, agreement_labels), name = NULL) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 13, face = "bold", color = "black", margin = margin(r = 16, l = 4)),
    legend.key.size = unit(0.65, "cm"),
    legend.margin = margin(t = 5, b = 5)
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))
leg_agr_grob <- cowplot::get_legend(p_agr)

# Assemble 6-panel GCM Agreement for SSP1-2.6
cat(">> Building 6-panel GCM Agreement figure for SSP1-2.6...\n")
p126_agr_list <- lapply(species_keys, function(k) create_agreement_panel(k, "ssp126"))
if (all(!sapply(p126_agr_list, is.null))) {
  g_agr_maps126 <- lapply(p126_agr_list, tmap::tmap_grob)
  g_agr_grid126 <- gridExtra::arrangeGrob(grobs = g_agr_maps126, ncol = 2, nrow = 3)
  fig126_agr <- gridExtra::arrangeGrob(g_agr_grid126, leg_agr_grob, ncol = 1, heights = c(16, 1))
  fig126_agr_path <- file.path(dir_figures, "All_species_gcm_agreement_ssp126_26_8.jpg")
  jpeg(fig126_agr_path, width = 30, height = 45, units = "cm", res = 600, quality = 95)
  grid::grid.draw(fig126_agr)
  dev.off()
  cat(sprintf("   Saved SSP1-2.6 6-panel GCM agreement figure: %s\n", fig126_agr_path))
}

# Assemble 6-panel GCM Agreement for SSP5-8.5
cat(">> Building 6-panel GCM Agreement figure for SSP5-8.5...\n")
p585_agr_list <- lapply(species_keys, function(k) create_agreement_panel(k, "ssp585"))
if (all(!sapply(p585_agr_list, is.null))) {
  g_agr_maps585 <- lapply(p585_agr_list, tmap::tmap_grob)
  g_agr_grid585 <- gridExtra::arrangeGrob(grobs = g_agr_maps585, ncol = 2, nrow = 3)
  fig585_agr <- gridExtra::arrangeGrob(g_agr_grid585, leg_agr_grob, ncol = 1, heights = c(16, 1))
  fig585_agr_path <- file.path(dir_figures, "All_species_gcm_agreement_ssp585_26_8.jpg")
  jpeg(fig585_agr_path, width = 30, height = 45, units = "cm", res = 600, quality = 95)
  grid::grid.draw(fig585_agr)
  dev.off()
  cat(sprintf("   Saved SSP5-8.5 6-panel GCM agreement figure: %s\n", fig585_agr_path))
}

cat("\n==============================================================================\n")
cat("ALL SPATIAL UNCERTAINTY ANALYSES COMPLETED SUCCESSFULLY FOR ALL 6 SPECIES!\n")
cat("==============================================================================\n")
