
suppressPackageStartupMessages({
  library(raster)
  library(terra)
  library(dismo)
  library(openxlsx)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

cat("==============================================================================\n")
cat("Starting Sensitivity & Robustness Analysis (Schoener's D) for All 6 Species\n")
cat("==============================================================================\n\n")

latlong <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Base Directories
dir_base    <- "E:/OneDrive/2024/DTCSCL_2024/Modelling"
dir_env     <- file.path(dir_base, "Enviromental/current")
dir_data    <- file.path(dir_base, "Data")
dir_output  <- file.path(dir_base, "output_maps")
dir_out_rev <- file.path(dir_base, "output_maps_26_8")
dir_figures <- file.path(dir_base, "Figures")
dir_maps    <- file.path(dir_base, "maps")

if (!dir.exists(dir_out_rev)) dir.create(dir_out_rev, recursive = TRUE)
if (!dir.exists(dir_figures)) dir.create(dir_figures, recursive = TRUE)

vn_mask <- raster::raster(file.path(dir_maps, "VN_raster.tif"))
fallback_bg <- file.path(dir_data, "backg.csv")

# Load regional unmasked predictors
env_files <- list.files(dir_env, full.names = TRUE, pattern = "\\.tif$")
predictors_full <- raster::stack(env_files)

# Species Configurations
species_configs <- list(
  cernua     = list(name = "Pinus cernua",     occ_file = "Pinus cernua_thin.xlsx",     bg_file = "pinus_cernua_backg.csv",     env_vars = c('bio6','bio7','bio13','bio14','bio15','evergreen_forest','slope','soil','SOC'),           opt_rm = 1.0, opt_fc = "LQ"),
  dalatensis = list(name = "Pinus dalatensis", occ_file = "Pinus dalatensis_thin.xlsx", bg_file = "Pinus_dalatensis_backg.csv", env_vars = c('bio1','bio2','bio4','bio7','bio14','bio15','bio16','bio18','evergreen_forest','slope','soil','SOC'), opt_rm = 1.0, opt_fc = "LQH"),
  fenzeliana = list(name = "Pinus fenzeliana", occ_file = "Pinus fenzeliana_thin.xlsx", bg_file = "Pinus_fenzeliana_backg.csv", env_vars = c('bio1','bio2','bio3','bio13','bio17','bio18','evergreen_forest','slope','soil','SOC'),           opt_rm = 1.0, opt_fc = "LQH"),
  henryi     = list(name = "Pinus henryi",     occ_file = "Pinus henryi_thin.xlsx",     bg_file = "pinus_henryi_backg.csv",     env_vars = c('bio3','bio5','bio6','bio7','bio12','bio14','bio15','slope','evergreen_forest','soil','SOC'),           opt_rm = 2.0, opt_fc = "LQ"),
  krempfii   = list(name = "Pinus krempfii",   occ_file = "Pinus krempfii_thin.xlsx",   bg_file = "pinus_krempfii_backg.csv",   env_vars = c('bio2','bio7','bio8','bio13','bio14','bio18','slope','evergreen_forest','soil','SOC'),           opt_rm = 1.0, opt_fc = "LQH"),
  latteri    = list(name = "Pinus latteri",    occ_file = "Pinus latteri_thin.xlsx",    bg_file = "pinus_latteri_backg.csv",    env_vars = c('bio2','bio4','bio5','bio13','bio18','bio19','slope','evergreen_forest','soil','SOC'),           opt_rm = 1.0, opt_fc = "LQ")
)

# Helper function to read background data safely
read_bg_data <- function(bg_path, fallback_path) {
  if (file.exists(bg_path)) {
    bg <- read.csv(bg_path)
    if (ncol(bg) >= 3) {
      bg <- bg[, c(2, 3)]
    } else {
      bg <- bg[, c(1, 2)]
    }
  } else {
    bg <- read.csv(fallback_path)[, c(2, 3)]
  }
  colnames(bg) <- c("Longitude", "Latitude")
  return(bg)
}

# Function to compute Schoener's D and Warren's I between two raster predictions
calc_niche_overlap <- function(r1, r2) {
  p1 <- na.omit(raster::values(r1))
  p2 <- na.omit(raster::values(r2))
  
  p1 <- p1 / sum(p1)
  p2 <- p2 / sum(p2)
  
  D <- 1 - 0.5 * sum(abs(p1 - p2))
  hellinger <- sqrt(sum((sqrt(p1) - sqrt(p2))^2))
  I_val <- 1 - 0.5 * (hellinger^2)
  r_cor <- cor(p1, p2)
  
  return(c(Schoener_D = D, Warren_I = I_val, Pearson_r = r_cor))
}

# Helper function to generate Maxent arguments based on RM and FC string
get_maxent_args <- function(rm_val, fc_str) {
  fc_upper <- toupper(fc_str)
  c(
    paste0("betamultiplier=", rm_val),
    paste0("linear=", ifelse(grepl("L", fc_upper), "true", "false")),
    paste0("quadratic=", ifelse(grepl("Q", fc_upper), "true", "false")),
    paste0("hinge=", ifelse(grepl("H", fc_upper), "true", "false")),
    paste0("product=", ifelse(grepl("P", fc_upper), "true", "false")),
    paste0("threshold=", ifelse(grepl("T", fc_upper), "true", "false")),
    "threads=2"
  )
}

sensitivity_results_list <- list()

for (sp_key in names(species_configs)) {
  cfg <- species_configs[[sp_key]]
  cat(sprintf("\n>>> Running Sensitivity Overlap for %s (Selected: RM=%.1f, FC=%s) <<<\n", cfg$name, cfg$opt_rm, cfg$opt_fc))
  
  # Load occurrence and background data
  occ_df <- openxlsx::read.xlsx(file.path(dir_data, cfg$occ_file), sheet = 1)
  occ_coords <- occ_df[, c("Longitude", "Latitude")]
  
  sp_env <- raster::subset(predictors_full, cfg$env_vars)
  
  # Filter complete cases
  occ_env_vals <- raster::extract(sp_env, occ_coords)
  valid_occ_idx <- which(complete.cases(occ_env_vals))
  occ_valid <- occ_coords[valid_occ_idx, ]
  
  bg_coords <- read_bg_data(file.path(dir_data, cfg$bg_file), fallback_bg)
  bg_env_vals <- raster::extract(sp_env, bg_coords)
  valid_bg_idx <- which(complete.cases(bg_env_vals))
  bg_valid <- bg_coords[valid_bg_idx, ]
  
  occ_sp <- occ_valid
  coordinates(occ_sp) <- ~ Longitude + Latitude
  crs(occ_sp) <- latlong
  
  # Vietnam prediction stack
  sp_env_vn <- sp_env * vn_mask
  names(sp_env_vn) <- cfg$env_vars
  
  # 1. Fit Selected (Baseline) Model with its specific opt_rm and opt_fc
  base_args <- get_maxent_args(cfg$opt_rm, cfg$opt_fc)
  xm_base <- dismo::maxent(
    x = sp_env, p = occ_sp, a = bg_valid, factors = "soil",
    args = base_args
  )
  pred_base <- dismo::predict(sp_env_vn, xm_base)
  
  # 2. Fit Competing Parameter Variations (Only FC: LQ, LQH and RM: -0.5, +0.5)
  rm_low <- max(0.5, cfg$opt_rm - 0.5)
  rm_high <- cfg$opt_rm + 0.5
  
  variations <- list(
    "FC_LQ"        = list(name = "FC_LQ",        args = get_maxent_args(cfg$opt_rm, "LQ")),
    "FC_LQH"       = list(name = "FC_LQH",       args = get_maxent_args(cfg$opt_rm, "LQH")),
    "RM_minus_0.5" = list(name = "RM_minus_0.5", args = get_maxent_args(rm_low, cfg$opt_fc)),
    "RM_plus_0.5"  = list(name = "RM_plus_0.5",  args = get_maxent_args(rm_high, cfg$opt_fc))
  )
  
  for (v_name in names(variations)) {
    v_cfg <- variations[[v_name]]
    xm_alt <- dismo::maxent(x = sp_env, p = occ_sp, a = bg_valid, factors = "soil", args = v_cfg$args)
    pred_alt <- dismo::predict(sp_env_vn, xm_alt)
    
    overlap_metrics <- calc_niche_overlap(pred_base, pred_alt)
    
    cat(sprintf("  Comparison [%s]: Schoener's D = %.4f | Warren's I = %.4f | Pearson r = %.4f\n",
                v_name, overlap_metrics["Schoener_D"], overlap_metrics["Warren_I"], overlap_metrics["Pearson_r"]))
    
    sensitivity_results_list[[paste0(sp_key, "_", v_name)]] <- data.frame(
      Species                  = cfg$name,
      Selected_Configuration   = paste0("RM = ", cfg$opt_rm, ", FC = ", cfg$opt_fc),
      Alternative_Configuration= v_name,
      Schoener_D               = round(overlap_metrics["Schoener_D"], 4),
      Warren_I                 = round(overlap_metrics["Warren_I"], 4),
      Pearson_Correlation_r    = round(overlap_metrics["Pearson_r"], 4),
      stringsAsFactors         = FALSE
    )
  }
}

sens_df <- do.call(rbind, sensitivity_results_list)

# ------------------------------------------------------------------------------
# SAVE SENSITIVITY TABLE & PLOT
# ------------------------------------------------------------------------------
cat("\n========================================================================================================\n")
cat("                       SENSITIVITY ANALYSIS: SPATIAL NICHE OVERLAP (SCHOENER'S D)                       \n")
cat("========================================================================================================\n")
print(sens_df, row.names = FALSE)
cat("========================================================================================================\n\n")

sens_xlsx_path <- file.path(dir_data, "Sensitivity_Analysis_Schoener_D_All_6_Species_26_8.xlsx")
openxlsx::write.xlsx(sens_df, file = sens_xlsx_path, overwrite = TRUE)
cat(sprintf(">> Sensitivity analysis table successfully saved to:\n   %s\n\n", sens_xlsx_path))

sens_df <- read.xlsx(file.path(dir_data, "Sensitivity_Analysis_Schoener_D_All_6_Species_26_8.xlsx"))
# Generate Publication-ready Barplot
p_sens <- ggplot(sens_df, aes(x = Alternative_Configuration, y = Schoener_D)) +
  geom_bar(stat = "identity", width = 0.6, fill = "steelblue", color = "black", alpha = 0.85) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "firebrick", linewidth = 0.8) +
  annotate("text", x = 1.5, y = 0.88, label = " ", color = "firebrick", fontface = "italic", size = 3.5) +
  facet_wrap(~ Species, ncol = 3) +
  ylim(0, 1.05) +
  labs(
    #title = "Sensitivity Analysis of Model Predictions across Alternative Configurations",
    #ubtitle = "Spatial overlap (Schoener's D) between Selected Model and Competing Parameter Variations",
    x = "Alternative Parameter Variation",
    y = "Schoener's D (Spatial Niche Overlap)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 8.5),
    legend.position = "none",
    strip.text = element_text(face = "italic", size = 10),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey30")
  )

sens_fig_path <- file.path(dir_figures, "Sensitivity_Analysis_Overlap_All_Species_26_8.jpg")
ggsave(filename = sens_fig_path, plot = p_sens, width = 11, height = 7, dpi = 600)
cat(sprintf(">> Sensitivity plot saved to:\n   %s\n\n", sens_fig_path))
cat("Sensitivity analysis completed successfully for all 6 Pinus species!\n")
