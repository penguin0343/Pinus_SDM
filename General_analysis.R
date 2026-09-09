
# ==============================================================================
# Script: Analysis.R
# Purpose: Post-modeling analyses across 6 native Pinus species in Vietnam:
#          1. Spatial thinning of occurrence records
#          2. Range area calculation and comparison (Current vs SSPs)
#          3. Proportion of suitable habitat in Protected Areas (PAs)
#          4. Elevational shift analysis
#          5. Absolute suitable habitat area in Protected Areas (PAs)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(terra)
  library(raster)
  library(spThin)
  library(openxlsx)
})

# ==============================================================================
# SECTION 1: SPATIAL THINNING OF OCCURRENCE RECORDS
# ==============================================================================
# Thinning occurrence records to reduce spatial autocorrelation (5 km distance)

data_dir <- "E:/OneDrive/2024/DTCSCL_2024/Modelling/Data"
pinus_combine_file <- file.path(data_dir, "pinus_combine.xlsx")

if (file.exists(pinus_combine_file)) {
  pinus_combine <- openxlsx::read.xlsx(pinus_combine_file, sheet = 1)
  
  for (sp in unique(pinus_combine$species)) {
    cat(sprintf(">> Thinning occurrence records for: %s\n", sp))
    pinus_sp <- pinus_combine %>% dplyr::filter(species == sp)
    
    # Coordinate column detection
    lat_col <- if ("latitude" %in% names(pinus_sp)) "latitude" else "Latitude"
    long_col <- if ("longitute" %in% names(pinus_sp)) "longitute" else if ("longitude" %in% names(pinus_sp)) "longitude" else "Longitude"
    
    thinned_results <- spThin::thin(
      loc.data = pinus_sp,
      lat.col = lat_col,
      long.col = long_col,
      spec.col = "species",
      thin.par = 5,                         # Minimum distance threshold (5 km)
      reps = 100,                           # 100 iterations to find optimal thinned subset
      locs.thinned.list.return = TRUE,
      write.files = FALSE,
      verbose = FALSE
    )
    
    records_count <- sapply(thinned_results, nrow)
    best_rep_index <- which.max(records_count)
    pinus_sp_thin <- thinned_results[[best_rep_index]]
    
    out_thin_path <- file.path(data_dir, paste0(sp, "_thin.xlsx"))
    openxlsx::write.xlsx(pinus_sp_thin, file = out_thin_path, overwrite = TRUE)
    cat(sprintf("   Saved thinned records (N = %d) to: %s\n", nrow(pinus_sp_thin), out_thin_path))
  }
}

# ==============================================================================
# SECTION 2: SUITABLE HABITAT AREA CALCULATION AND COMPARISON
# ==============================================================================

species_names <- c("Pinus dalatensis", "Pinus krempfii", "Pinus latteri",
                   "Pinus fenzeliana", "Pinus cernua", "Pinus henryi")
clean_names   <- gsub(" ", "_", tolower(species_names))
scenarios     <- c("current", "ssp126", "ssp585")
gcms          <- c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")
base_path     <- "E:/OneDrive/2024/DTCSCL_2024/Modelling/output_maps/"
fig_dir       <- "E:/OneDrive/2024/DTCSCL_2024/Modelling/Figures"

if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# Process suitable habitat area per species and climate scenario
process_species_area <- function(sp_clean, sp_real) {
  purrr::map_df(scenarios, function(scen) {
    suffix <- if (scen == "current") "_current.tif" else paste0("_", scen, "_ensemble.tif")
    path_main <- file.path(base_path, paste0(sp_clean, suffix))
    
    if (!file.exists(path_main)) {
      alt_clean <- paste0("Pinus_", substring(sp_clean, 7))
      path_main <- file.path(base_path, paste0(alt_clean, suffix))
    }
    
    if (!file.exists(path_main)) {
      return(data.frame(species = sp_real, scenario = scen, area = NA_real_, Min_GCM = NA_real_, Max_GCM = NA_real_))
    }
    
    r_main <- terra::rast(path_main)
    area_val <- as.numeric(terra::global(r_main, fun = "sum", na.rm = TRUE))
    
    # Calculate Min/Max across GCMs for future scenarios
    min_gcm <- NA_real_
    max_gcm <- NA_real_
    
    if (scen != "current") {
      gcm_paths <- file.path(base_path, paste0(sp_clean, "_", gcms, "_", scen, ".tif"))
      existing_paths <- gcm_paths[file.exists(gcm_paths)]
      if (length(existing_paths) == 0) {
        alt_clean <- paste0("Pinus_", substring(sp_clean, 7))
        gcm_paths <- file.path(base_path, paste0(alt_clean, "_", gcms, "_", scen, ".tif"))
        existing_paths <- gcm_paths[file.exists(gcm_paths)]
      }
      
      if (length(existing_paths) > 0) {
        gcm_areas <- purrr::map_dbl(existing_paths, function(p) {
          as.numeric(terra::global(terra::rast(p), fun = "sum", na.rm = TRUE))
        })
        min_gcm <- min(gcm_areas, na.rm = TRUE)
        max_gcm <- max(gcm_areas, na.rm = TRUE)
      }
    }
    
    data.frame(
      species  = sp_real,
      scenario = scen,
      area     = area_val,
      Min_GCM  = min_gcm,
      Max_GCM  = max_gcm
    )
  })
}

pine_data <- purrr::map2_dfr(clean_names, species_names, process_species_area) %>%
  dplyr::mutate(scenario = dplyr::recode(scenario,
                                         "current" = "Current",
                                         "ssp126"  = "SSP1-2.6",
                                         "ssp585"  = "SSP5-8.5"))

p_area <- ggplot(pine_data, aes(x = factor(species, levels = species_names), y = area, fill = scenario)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", alpha = 0.85) +
  geom_errorbar(aes(ymin = Min_GCM, ymax = Max_GCM),
                position = position_dodge(width = 0.8),
                width = 0.25,
                na.rm = TRUE) +
  labs(x = "", y = expression(paste("Suitable Habitat Area (", km^2, ")")), fill = "") +
  scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(face = "italic", angle = 30, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = c(0.85, 0.85),
    legend.background = element_rect(fill = "transparent")
  ) +
  scale_fill_manual(values = c("Current" = "#2ecc71", "SSP1-2.6" = "#3498db", "SSP5-8.5" = "#e74c3c"))

ggsave(p_area, filename = file.path(fig_dir, "Pinus_area.jpg"), width = 12, height = 10, units = 'cm', dpi = 600)

# ==============================================================================
# SECTION 3: PROPORTION OF SUITABLE HABITAT IN PROTECTED AREAS (%)
# ==============================================================================

vn_pas_path <- "E:/OneDrive/2024/DTCSCL_2024/Modelling/maps/VN_pas_raster.tif"
if (file.exists(vn_pas_path)) {
  vn_pas <- terra::rast(vn_pas_path)
  
  process_species_pa_prop <- function(sp_clean, sp_real) {
    purrr::map_df(scenarios, function(scen) {
      suffix <- if (scen == "current") "_current.tif" else paste0("_", scen, "_ensemble.tif")
      path_main <- file.path(base_path, paste0(sp_clean, suffix))
      
      if (!file.exists(path_main)) {
        alt_clean <- paste0("Pinus_", substring(sp_clean, 7))
        path_main <- file.path(base_path, paste0(alt_clean, suffix))
      }
      if (!file.exists(path_main)) {
        return(data.frame(species = sp_real, scenario = scen, area = NA_real_, Min_GCM = NA_real_, Max_GCM = NA_real_))
      }
      
      r_main <- terra::rast(path_main)
      vn_pas_sync <- terra::resample(vn_pas, r_main, method = "near")
      
      # Calculate protected proportion for ensemble
      sp_pas <- r_main * vn_pas_sync
      tot_val <- as.numeric(terra::global(r_main, fun = "sum", na.rm = TRUE))
      pa_val  <- as.numeric(terra::global(sp_pas, fun = "sum", na.rm = TRUE))
      area_val <- ifelse(tot_val > 0, (pa_val / tot_val) * 100, NA_real_)
      
      min_gcm <- NA_real_
      max_gcm <- NA_real_
      
      if (scen != "current") {
        gcm_paths <- file.path(base_path, paste0(sp_clean, "_", gcms, "_", scen, ".tif"))
        existing_paths <- gcm_paths[file.exists(gcm_paths)]
        if (length(existing_paths) == 0) {
          alt_clean <- paste0("Pinus_", substring(sp_clean, 7))
          gcm_paths <- file.path(base_path, paste0(alt_clean, "_", gcms, "_", scen, ".tif"))
          existing_paths <- gcm_paths[file.exists(gcm_paths)]
        }
        
        if (length(existing_paths) > 0) {
          gcm_props <- purrr::map_dbl(existing_paths, function(p) {
            r_gcm <- terra::rast(p)
            num <- as.numeric(terra::global(r_gcm * vn_pas_sync, fun = "sum", na.rm = TRUE))
            den <- as.numeric(terra::global(r_gcm, fun = "sum", na.rm = TRUE))
            if (is.na(den) || den == 0) return(NA_real_)
            return((num / den) * 100)
          })
          gcm_props_clean <- gcm_props[!is.na(gcm_props)]
          if (length(gcm_props_clean) > 0) {
            min_gcm <- min(gcm_props_clean)
            max_gcm <- max(gcm_props_clean)
          }
        }
      }
      
      data.frame(
        species  = sp_real,
        scenario = scen,
        area     = area_val,
        Min_GCM  = min_gcm,
        Max_GCM  = max_gcm
      )
    })
  }
  
  pine_pa_prop_data <- purrr::map2_dfr(clean_names, species_names, process_species_pa_prop) %>%
    dplyr::mutate(scenario = dplyr::recode(scenario,
                                           "current" = "Current",
                                           "ssp126"  = "SSP1-2.6",
                                           "ssp585"  = "SSP5-8.5"))
  
  p_pa_prop <- ggplot(pine_pa_prop_data, aes(x = factor(species, levels = species_names), y = area, fill = scenario)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", alpha = 0.85) +
    geom_errorbar(aes(ymin = Min_GCM, ymax = Max_GCM),
                  position = position_dodge(width = 0.8),
                  width = 0.25,
                  na.rm = TRUE) +
    labs(x = "", y = "Protected Area Coverage (%)", fill = "") +
    theme_classic() +
    theme(
      axis.text.x = element_text(face = "italic", angle = 30, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      legend.position = c(0.2, 0.85),
      legend.background = element_rect(fill = "transparent")
    ) +
    scale_fill_manual(values = c("Current" = "#2ecc71", "SSP1-2.6" = "#3498db", "SSP5-8.5" = "#e74c3c"))
  
  ggsave(p_pa_prop, filename = file.path(fig_dir, "Pinus_pas_proportion.jpg"), width = 12, height = 10, units = 'cm', dpi = 600)
}

# ==============================================================================
# SECTION 4: ELEVATION SHIFT ANALYSIS
# ==============================================================================

path_ele <- "E:/OneDrive/2024/DTCSCL_2024/Modelling/Enviromental/current/elevation.tif"

if (file.exists(path_ele)) {
  elevation <- terra::rast(path_ele)
  
  elevation_files <- c(
    "Pinus_dalatensis_current.tif", "Pinus_dalatensis_ssp126_ensemble.tif", "Pinus_dalatensis_ssp585_ensemble.tif",
    "Pinus_krempfii_current.tif",   "Pinus_krempfii_ssp126_ensemble.tif",   "Pinus_krempfii_ssp585_ensemble.tif",
    "Pinus_latteri_current.tif",    "Pinus_latteri_ssp126_ensemble.tif",    "Pinus_latteri_ssp585_ensemble.tif",
    "Pinus_fenzeliana_current.tif", "Pinus_fenzeliana_ssp126_ensemble.tif", "Pinus_fenzeliana_ssp585_ensemble.tif",
    "Pinus_cernua_current.tif",     "Pinus_cernua_ssp126_ensemble.tif",     "Pinus_cernua_ssp585_ensemble.tif",
    "Pinus_henryi_current.tif",     "Pinus_henryi_ssp126_ensemble.tif",     "Pinus_henryi_ssp585_ensemble.tif"
  )
  
  ele_data_list <- list()
  
  for (f in elevation_files) {
    full_f_path <- file.path(base_path, f)
    if (!file.exists(full_f_path)) {
      full_f_path <- file.path(base_path, tolower(f))
    }
    if (!file.exists(full_f_path)) next
    
    cat(sprintf(">> Processing elevation shift for: %s\n", f))
    r_sp <- terra::rast(full_f_path)
    
    if (!terra::compareGeom(r_sp, elevation, stopOnError = FALSE)) {
      r_sp <- terra::resample(r_sp, elevation, method = "near")
    }
    
    r_combined <- c(r_sp, elevation)
    df <- as.data.frame(r_combined, cells = FALSE, na.rm = TRUE)
    colnames(df) <- c("Suitability", "Elevation")
    df_filtered <- df %>% dplyr::filter(Suitability > 0)
    
    clean_f <- stringr::str_remove(f, "\\.tif")
    sp_name_raw <- stringr::str_extract(clean_f, "(?i)Pinus_[a-z]+")
    sp_name_clean <- stringr::str_replace(sp_name_raw, "_", " ")
    sp_name_clean <- paste0("Pinus ", tolower(substring(sp_name_clean, 7)))
    
    scen_name <- stringr::str_remove_all(clean_f, paste0(sp_name_raw, "_|_ensemble|_current"))
    
    ele_data_list[[f]] <- data.frame(
      Species   = sp_name_clean,
      Scenario  = scen_name,
      Elevation = df_filtered$Elevation
    )
  }
  
  if (length(ele_data_list) > 0) {
    final_ele_df <- dplyr::bind_rows(ele_data_list) %>%
      dplyr::mutate(
        Scenario = factor(Scenario, levels = c("current", "ssp126", "ssp585"),
                          labels = c("Current", "SSP1-2.6", "SSP5-8.5")),
        Species  = factor(Species, levels = species_names)
      )
    
    p_ele <- ggplot(final_ele_df, aes(x = Scenario, y = Elevation, fill = Scenario)) +
      geom_boxplot(outlier.size = 0.4, alpha = 0.75, width = 0.6) +
      facet_wrap(~ Species, scales = "free_y", ncol = 3) +
      theme_minimal() +
      labs(y = "Elevation (m a.s.l.)", x = "", fill = "") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        strip.text = element_text(face = "italic", size = 10),
        legend.position = "bottom"
      ) +
      scale_fill_manual(values = c("Current" = "#2ecc71", "SSP1-2.6" = "#3498db", "SSP5-8.5" = "#e74c3c"))
    
    ggsave(p_ele, filename = file.path(fig_dir, "Pinus_elevation_shift.tiff"), width = 14, height = 12, units = 'cm', dpi = 600)
  }
}

# ==============================================================================
# SECTION 5: ABSOLUTE PROTECTED HABITAT AREA (KM2)
# ==============================================================================

if (file.exists(vn_pas_path)) {
  vn_pas <- terra::rast(vn_pas_path)
  
  process_species_pa_abs <- function(sp_clean, sp_real) {
    purrr::map_df(scenarios, function(scen) {
      suffix <- if (scen == "current") "_current.tif" else paste0("_", scen, "_ensemble.tif")
      path_main <- file.path(base_path, paste0(sp_clean, suffix))
      
      if (!file.exists(path_main)) {
        alt_clean <- paste0("Pinus_", substring(sp_clean, 7))
        path_main <- file.path(base_path, paste0(alt_clean, suffix))
      }
      if (!file.exists(path_main)) {
        return(data.frame(species = sp_real, scenario = scen, area = NA_real_, Min_GCM = NA_real_, Max_GCM = NA_real_))
      }
      
      r_main <- terra::rast(path_main)
      vn_pas_sync <- terra::resample(vn_pas, r_main, method = "near")
      sp_pas <- r_main * vn_pas_sync
      area_val <- as.numeric(terra::global(sp_pas, fun = "sum", na.rm = TRUE))
      
      min_gcm <- NA_real_
      max_gcm <- NA_real_
      
      if (scen != "current") {
        gcm_paths <- file.path(base_path, paste0(sp_clean, "_", gcms, "_", scen, ".tif"))
        existing_paths <- gcm_paths[file.exists(gcm_paths)]
        if (length(existing_paths) == 0) {
          alt_clean <- paste0("Pinus_", substring(sp_clean, 7))
          gcm_paths <- file.path(base_path, paste0(alt_clean, "_", gcms, "_", scen, ".tif"))
          existing_paths <- gcm_paths[file.exists(gcm_paths)]
        }
        
        if (length(existing_paths) > 0) {
          gcm_areas <- purrr::map_dbl(existing_paths, function(p) {
            r_gcm <- terra::rast(p)
            as.numeric(terra::global(r_gcm * vn_pas_sync, fun = "sum", na.rm = TRUE))
          })
          gcm_areas_clean <- gcm_areas[!is.na(gcm_areas)]
          if (length(gcm_areas_clean) > 0) {
            min_gcm <- min(gcm_areas_clean)
            max_gcm <- max(gcm_areas_clean)
          }
        }
      }
      
      data.frame(
        species  = sp_real,
        scenario = scen,
        area     = area_val,
        Min_GCM  = min_gcm,
        Max_GCM  = max_gcm
      )
    })
  }
  
  pine_pa_abs_data <- purrr::map2_dfr(clean_names, species_names, process_species_pa_abs) %>%
    dplyr::mutate(scenario = dplyr::recode(scenario,
                                           "current" = "Current",
                                           "ssp126"  = "SSP1-2.6",
                                           "ssp585"  = "SSP5-8.5"))
  
  p_pa_abs <- ggplot(pine_pa_abs_data, aes(x = factor(species, levels = species_names), y = area, fill = scenario)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", alpha = 0.85) +
    geom_errorbar(aes(ymin = Min_GCM, ymax = Max_GCM),
                  position = position_dodge(width = 0.8),
                  width = 0.25,
                  na.rm = TRUE) +
    labs(x = "", y = expression(paste("Protected Habitat Area (", km^2, ")")), fill = "") +
    scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
    theme_classic() +
    theme(
      axis.text.x = element_text(face = "italic", angle = 30, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      legend.position = c(0.8, 0.85),
      legend.background = element_rect(fill = "transparent")
    ) +
    scale_fill_manual(values = c("Current" = "#2ecc71", "SSP1-2.6" = "#3498db", "SSP5-8.5" = "#e74c3c"))
  
  ggsave(p_pa_abs, filename = file.path(fig_dir, "Pinus_pas_absolute.jpg"), width = 12, height = 10, units = 'cm', dpi = 600)
}

cat("\nAll post-modeling analyses completed successfully!\n")




