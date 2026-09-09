# ==============================================================================
# Script: Modelling_Pinus henryi.R
# Purpose: Complete SDM Workflow for Pinus henryi in Vietnam:
#          1. Background point generation (buffer sampling)
#          2. Environmental variable multicollinearity evaluation
#          3. Model calibration & cross-validation with ENMeval
#          4. Final MaxEnt model fitting & projection (Current, 2070 SSP1-2.6 & SSP5-8.5)
#          5. Ensemble consensus forecasting (>= 3 GCMs)
#          6. Range change dynamics mapping
#          7. Variable contribution assessment
#          8. Multivariate Environmental Similarity Surface (MESS) analysis
# ==============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(raster)
  library(terra)
  library(dismo)
  library(maxnet)
  library(ENMeval)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggcorrplot)
  library(readxl)
  library(openxlsx)
})

options(java.parameters = "-Xmx32g")

# Global CRS and directories
latlong <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
dir_base       <- "E:/OneDrive/2024/DTCSCL_2024/Modelling"
dir_data       <- file.path(dir_base, "Data")
dir_env_curr   <- file.path(dir_base, "Enviromental/current")
dir_env_2070   <- file.path(dir_base, "Enviromental/2070")
dir_maps       <- file.path(dir_base, "maps")
dir_out_maps   <- file.path(dir_base, "output_maps")
dir_maxent_out <- file.path(dir_base, "Maxent_responses/Pinus_henryi")
dir_figures    <- file.path(dir_base, "Figures")

if (!dir.exists(dir_out_maps))   dir.create(dir_out_maps, recursive = TRUE)
if (!dir.exists(dir_maxent_out)) dir.create(dir_maxent_out, recursive = TRUE)
if (!dir.exists(dir_figures))    dir.create(dir_figures, recursive = TRUE)

selected_environment <- c('bio3','bio5','bio6','bio7','bio12','bio14','bio15','slope','evergreen_forest','soil','SOC')

# Vietnam mask raster
vn_raster <- raster::raster(file.path(dir_maps, "VN_raster.tif"))

# ==============================================================================
# SECTION 1: BACKGROUND POINT GENERATION
# ==============================================================================
cat(">> Generating background points within geographic buffer...\n")
elevation <- raster::raster(file.path(dir_env_curr, "elevation.tif"))
occ_raw <- openxlsx::read.xlsx(file.path(dir_data, "Pinus henryi_thin.xlsx"), sheet = 1)

pinus_henryi_sp <- occ_raw
coordinates(pinus_henryi_sp) <- ~ Longitude + Latitude
crs(pinus_henryi_sp) <- latlong

pinus_henryi_sf <- sf::st_as_sf(pinus_henryi_sp)
pinus_henryi_buffer <- sf::st_buffer(pinus_henryi_sf, dist = 100000) # 100 km buffer
elevation_buffer <- raster::mask(elevation, pinus_henryi_buffer)

set.seed(42)
backg_pts <- dismo::randomPoints(elevation_buffer, 10000)
write.csv(backg_pts, file = file.path(dir_data, "pinus_henryi_backg.csv"), row.names = FALSE)

# ==============================================================================
# SECTION 2: VARIABLE MULTICOLLINEARITY (CORRELATION MATRIX)
# ==============================================================================
cat(">> Calculating environmental variable correlation matrix...\n")
env_files_all <- list.files(dir_env_curr, full.names = TRUE, pattern = "\\.tif$")
predictors_all <- raster::stack(env_files_all)
bioclim_vars <- c('bio1','bio2','bio3','bio4','bio5','bio6','bio7','bio8','bio9',
                  'bio10','bio11','bio12','bio13','bio14','bio15','bio16','bio17','bio18','bio19', 'elevation', 'slope')
predictors_check <- raster::subset(predictors_all, bioclim_vars)

bg_check <- read.csv(file.path(dir_data, "pinus_henryi_backg.csv"))
if (ncol(bg_check) >= 3) bg_check <- bg_check[, c(2, 3)] else bg_check <- bg_check[, c(1, 2)]
colnames(bg_check) <- c("Longitude", "Latitude")

pt_hy_env <- raster::extract(predictors_check, bg_check, df = TRUE)
pt_hy_env <- pt_hy_env[!is.na(as.numeric(as.character(pt_hy_env$bio1))), ]
pt_hy_env <- pt_hy_env[!is.na(as.numeric(as.character(pt_hy_env$elevation))), ]
pt_hy_env <- pt_hy_env[!is.na(as.numeric(as.character(pt_hy_env$slope))), ]
env_cor <- data.frame(cor(pt_hy_env[, 2:22], method = "spearman"))

p_cor <- ggcorrplot(
  env_cor, method = "square", show.legend = TRUE, show.diag = TRUE, lab = TRUE,
  ggtheme = ggplot2::theme_classic(), colors = c("#6D9EC1", "white", "#E46726"),
  hc.order = FALSE, type = 'lower', legend.title = 'Coefficient',
  sig.level = TRUE, pch = 2, pch.col = 'red', digits = 1, title = ""
) + theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7), text = element_text(size = 7))

ggsave(p_cor, filename = file.path(dir_figures, "pinus_henryi_correlation_matrix.jpg"), width = 20, height = 20, units = 'cm', dpi = 300)

# ==============================================================================
# SECTION 3: MODEL CALIBRATION & PERFORMANCE EVALUATION (ENMeval)
# ==============================================================================
cat(">> Running ENMeval calibration and evaluation...\n")
pinus_henryi_occ <- openxlsx::read.xlsx(file.path(dir_data, "Pinus henryi_thin.xlsx"), sheet = 1) %>%
  dplyr::select(Longitude, Latitude)

backg_eval <- read.csv(file.path(dir_data, "pinus_henryi_backg.csv"))
if (ncol(backg_eval) >= 3) backg_eval <- backg_eval[, c(2, 3)] else backg_eval <- backg_eval[, c(1, 2)]
colnames(backg_eval) <- c("Longitude", "Latitude")

predictors_sel <- raster::subset(predictors_all, selected_environment)
study_extent <- raster::extent(min(backg_eval$Longitude), max(backg_eval$Longitude),
                               min(backg_eval$Latitude), max(backg_eval$Latitude))
predictors_crop <- raster::crop(predictors_sel, study_extent)

res <- ENMeval::ENMevaluate(
  occs = pinus_henryi_occ,
  envs = terra::rast(predictors_crop),
  bg = backg_eval,
  tune.args = list(fc = c("LQ", "LQH"), rm = seq(0.5, 5, 0.5)),
  partitions = "block",
  algorithm = "maxnet",
  numCores = 1
)

# Extract environmental data per fold and compute TSS and Kappa metrics
occ_env_eval <- as.data.frame(raster::extract(predictors_crop, pinus_henryi_occ))
bg_env_eval  <- as.data.frame(raster::extract(predictors_crop, backg_eval))
fold_ids <- sort(unique(res@occs.grp))
fold_data_list <- list()

for (k in fold_ids) {
  val_p   <- occ_env_eval[res@occs.grp == k, , drop = FALSE]
  val_a   <- bg_env_eval[res@bg.grp == k, , drop = FALSE]
  train_p <- occ_env_eval[res@occs.grp != k, , drop = FALSE]
  train_a <- bg_env_eval[res@bg.grp != k, , drop = FALSE]
  
  train_p <- train_p[complete.cases(train_p), , drop = FALSE]
  train_a <- train_a[complete.cases(train_a), , drop = FALSE]
  val_p   <- val_p[complete.cases(val_p), , drop = FALSE]
  val_a   <- val_a[complete.cases(val_a), , drop = FALSE]
  
  fold_data_list[[as.character(k)]] <- list(
    train_p = train_p, train_a = train_a,
    val_p = val_p, val_a = val_a,
    p_train_vec = c(rep(1, nrow(train_p)), rep(0, nrow(train_a))),
    data_train = rbind(train_p, train_a)
  )
}

res_df <- res@results
all_tss_avg <- numeric(nrow(res_df))
all_tss_sd  <- numeric(nrow(res_df))
all_kappa_avg <- numeric(nrow(res_df))
all_kappa_sd  <- numeric(nrow(res_df))
all_fold_eval_list <- list()

for (i in 1:nrow(res_df)) {
  fc_val <- as.character(res_df$fc[i])
  rm_val <- as.numeric(res_df$rm[i])
  curr_tune <- as.character(res_df$tune.args[i])
  
  f_args <- list(
    l = grepl("L", fc_val), q = grepl("Q", fc_val), h = grepl("H", fc_val),
    p = grepl("P", fc_val), t = grepl("T", fc_val)
  )
  f_classes <- paste(names(f_args)[unlist(f_args)], collapse = "")
  if (f_classes == "") f_classes <- "lq"
  
  f_tss   <- numeric(length(fold_ids))
  f_kappa <- numeric(length(fold_ids))
  
  for (idx in seq_along(fold_ids)) {
    k  <- fold_ids[idx]
    fd <- fold_data_list[[as.character(k)]]
    
    mod_k <- tryCatch({
      maxnet::maxnet(
        p = fd$p_train_vec, data = fd$data_train,
        f = maxnet::maxnet.formula(p = fd$p_train_vec, data = fd$data_train, classes = f_classes),
        regmult = rm_val
      )
    }, error = function(e) NULL)
    
    if (!is.null(mod_k)) {
      v_p <- tryCatch(as.numeric(predict(mod_k, fd$val_p, type = "cloglog")), error = function(e) NA)
      v_a <- tryCatch(as.numeric(predict(mod_k, fd$val_a, type = "cloglog")), error = function(e) NA)
      v_p <- na.omit(v_p)
      v_a <- na.omit(v_a)
      
      if (length(v_p) > 0 && length(v_a) > 0) {
        ev <- dismo::evaluate(p = v_p, a = v_a)
        f_tss[idx]   <- max(ev@TPR + ev@TNR - 1, na.rm = TRUE)
        f_kappa[idx] <- max(ev@kappa, na.rm = TRUE)
      } else {
        f_tss[idx] <- NA; f_kappa[idx] <- NA
      }
    } else {
      f_tss[idx] <- NA; f_kappa[idx] <- NA
    }
  }
  
  all_tss_avg[i]   <- mean(f_tss, na.rm = TRUE)
  all_tss_sd[i]    <- sd(f_tss, na.rm = TRUE)
  all_kappa_avg[i] <- mean(f_kappa, na.rm = TRUE)
  all_kappa_sd[i]  <- sd(f_kappa, na.rm = TRUE)
  
  all_fold_eval_list[[curr_tune]] <- data.frame(
    fold = as.character(fold_ids),
    TSS_max = as.numeric(f_tss),
    Kappa_max = as.numeric(f_kappa),
    stringsAsFactors = FALSE
  )
}

res@results$tss.val.avg   <- round(all_tss_avg, 4)
res@results$tss.val.sd    <- round(all_tss_sd, 4)
res@results$kappa.val.avg <- round(all_kappa_avg, 4)
res@results$kappa.val.sd  <- round(all_kappa_sd, 4)

best_tune <- res@results %>%
  filter(auc.val.avg == max(auc.val.avg, na.rm = TRUE)) %>%
  slice(1) %>%
  pull(tune.args)

cat(sprintf(">> Optimal parameter configuration: %s (AUC: %.4f, TSS: %.4f)\n",
            best_tune,
            res@results$auc.val.avg[res@results$tune.args == best_tune],
            res@results$tss.val.avg[res@results$tune.args == best_tune]))

# Fold summary
n_total_occ <- nrow(pinus_henryi_occ)
n_total_bg  <- nrow(backg_eval)
fold_counts_list <- list()

for (k in fold_ids) {
  n_test_occ  <- sum(res@occs.grp == k)
  n_train_occ <- n_total_occ - n_test_occ
  n_test_bg   <- sum(res@bg.grp == k)
  n_train_bg  <- n_total_bg - n_test_bg
  
  fold_counts_list[[as.character(k)]] <- data.frame(
    fold = as.character(k),
    train_occ = n_train_occ, test_occ = n_test_occ,
    train_bg = n_train_bg, test_bg = n_test_bg,
    stringsAsFactors = FALSE
  )
}
fold_counts <- do.call(rbind, fold_counts_list)
fold_res <- res@results.partitions %>% filter(tune.args == best_tune)
fold_res$fold <- as.character(fold_res$fold)

fold_summary <- fold_counts %>%
  left_join(fold_res, by = "fold") %>%
  left_join(all_fold_eval_list[[best_tune]], by = "fold") %>%
  dplyr::select(fold, train_occ, test_occ, train_bg, test_bg, tune.args, auc.val, TSS_max, Kappa_max, cbi.val, or.10p, or.mtp)

cat("\n=================== FOLD-LEVEL VALIDATION RESULTS ===================\n")
print(fold_summary)

# Save evaluation results
openxlsx::write.xlsx(res@results, file = file.path(dir_data, "pinus_henryi_model_performance_26_8.xlsx"), overwrite = TRUE)
openxlsx::write.xlsx(fold_summary, file = file.path(dir_data, "pinus_henryi_fold_performance_26_8.xlsx"), overwrite = TRUE)

# ==============================================================================
# SECTION 4: FINAL MODEL FITTING & PREDICTION (CURRENT & FUTURE)
# ==============================================================================
cat("\n>> Fitting final MaxEnt model...\n")
pinus_henryi_sp <- pinus_henryi_occ
coordinates(pinus_henryi_sp) <- ~ Longitude + Latitude
crs(pinus_henryi_sp) <- latlong

predictors_final <- raster::subset(predictors_all, selected_environment)

xm <- dismo::maxent(
  x = predictors_final,
  p = pinus_henryi_sp,
  a = backg_eval,
  factors = 'soil',
  path = dir_maxent_out,
  args = c(
    'betamultiplier=2',
    'linear=true',
    'quadratic=true',
    'hinge=true',
    'product=false',
    'threshold=false',
    'threads=2',
    'responsecurves=true',
    'jackknife=true',
    'askoverwrite=false',
    'autofeature=false'
  )
)

# Current prediction
cat(">> Projecting to current climate...\n")
pred_curr_crop <- (predictors_final * vn_raster)
names(pred_curr_crop) <- selected_environment

pred_pinus_henryi <- dismo::predict(pred_curr_crop, xm)
maxent_result <- read.csv(file.path(dir_maxent_out, "maxentResults.csv"), header = TRUE)
threshold_val <- maxent_result$Equal.training.sensitivity.and.specificity.Cloglog.threshold[1]
pinus_henryi_bin <- pred_pinus_henryi > threshold_val

raster::writeRaster(pred_pinus_henryi, filename = file.path(dir_out_maps, "pinus_henryi_current_raw.tif"), format = 'GTiff', overwrite = TRUE)
raster::writeRaster(pinus_henryi_bin, filename = file.path(dir_out_maps, "pinus_henryi_current.tif"), format = 'GTiff', overwrite = TRUE)

# Future projections: 2070 SSP1-2.6 & SSP5-8.5 across 5 GCMs
gcms <- c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")
scenarios <- c("ssp126", "ssp585")

for (scen in scenarios) {
  cat(sprintf(">> Projecting for 2070 %s across GCMs...\n", toupper(scen)))
  for (gcm in gcms) {
    gcm_files <- list.files(file.path(dir_env_2070, gcm, scen), full.names = TRUE, pattern = "\\.tif$")
    pred_2070 <- raster::stack(gcm_files)
    pred_2070 <- raster::subset(pred_2070, selected_environment)
    pred_2070_crop <- pred_2070 * vn_raster
    names(pred_2070_crop) <- selected_environment
    
    pred_gcm_raw <- dismo::predict(pred_2070_crop, xm)
    pred_gcm_bin <- pred_gcm_raw > threshold_val
    
    out_bin_path <- file.path(dir_out_maps, paste0("pinus_henryi_", gcm, "_", scen, ".tif"))
    out_raw_path <- file.path(dir_out_maps, paste0("pinus_henryi_", gcm, "_", scen, "_raw.tif"))
    
    raster::writeRaster(pred_gcm_bin, filename = out_bin_path, overwrite = TRUE)
    raster::writeRaster(pred_gcm_raw, filename = out_raw_path, overwrite = TRUE)
    cat(sprintf("   Finished %s - %s\n", gcm, scen))
  }
}

# ==============================================================================
# SECTION 5: CONSENSUS ENSEMBLE FORECASTING (>= 3 GCMS)
# ==============================================================================
cat(">> Generating multi-model ensemble consensus maps...\n")
for (scen in scenarios) {
  gcm_rasters <- lapply(gcms, function(g) {
    raster::raster(file.path(dir_out_maps, paste0("pinus_henryi_", g, "_", scen, ".tif")))
  })
  ensemble_r <- raster::overlay(raster::stack(gcm_rasters), fun = sum) >= 3
  out_ens_path <- file.path(dir_out_maps, paste0("pinus_henryi_", scen, "_ensemble.tif"))
  raster::writeRaster(ensemble_r, filename = out_ens_path, overwrite = TRUE)
  cat(sprintf("   Saved ensemble for %s: %s\n", scen, out_ens_path))
}

# ==============================================================================
# SECTION 6: RANGE DYNAMICS CLASSIFICATION
# ==============================================================================
cat(">> Calculating range dynamics rasters...\n")
p_curr <- terra::rast(file.path(dir_out_maps, "pinus_henryi_current.tif"))
p_126  <- terra::rast(file.path(dir_out_maps, "pinus_henryi_ssp126_ensemble.tif"))
p_585  <- terra::rast(file.path(dir_out_maps, "pinus_henryi_ssp585_ensemble.tif"))

p_dyn_126 <- p_curr + (p_126 * 2)
p_dyn_585 <- p_curr + (p_585 * 2)

terra::writeRaster(p_dyn_126, filename = file.path(dir_out_maps, "Pinus_henryi_change_ssp126.tif"), overwrite = TRUE)
terra::writeRaster(p_dyn_585, filename = file.path(dir_out_maps, "Pinus_henryi_change_ssp585.tif"), overwrite = TRUE)

# ==============================================================================
# SECTION 7: VARIABLE CONTRIBUTION PLOT
# ==============================================================================
cat(">> Plotting variable contributions...\n")
contrib_cols <- grep("\\.contribution$", names(maxent_result), value = TRUE)
plot_data <- maxent_result %>%
  dplyr::select(dplyr::all_of(contrib_cols)) %>%
  tidyr::pivot_longer(cols = dplyr::everything(), names_to = "Variable", values_to = "Contribution") %>%
  dplyr::mutate(Variable = gsub("\\.contribution", "", Variable)) %>%
  dplyr::mutate(Variable = dplyr::case_when(
    Variable == "evergreen_forest" ~ "EF",
    TRUE ~ paste0(toupper(substring(Variable, 1, 1)), substring(Variable, 2))
  )) %>%
  dplyr::group_by(Variable) %>%
  dplyr::summarise(Mean = mean(Contribution, na.rm = TRUE), SD = sd(Contribution, na.rm = TRUE)) %>%
  dplyr::arrange(dplyr::desc(Mean))

p_contrib <- ggplot(plot_data, aes(x = reorder(Variable, Mean), y = Mean)) +
  geom_bar(stat = "identity", color = "black", width = 0.7, fill = "#4682B4") +
  geom_errorbar(aes(ymin = pmax(0, Mean - SD), ymax = Mean + SD), width = 0.2, alpha = 0.5) +
  coord_flip() +
  labs(x = "", y = "Contribution (%)") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 10, face = "bold"),
    title = element_text(size = 12, face = "bold")
  )

ggsave(p_contrib, filename = file.path(dir_figures, "Pinus_henryi_variables.jpg"), width = 8, height = 6, units = 'cm', dpi = 600)

# ==============================================================================
# SECTION 8: MULTIVARIATE ENVIRONMENTAL SIMILARITY SURFACE (MESS)
# ==============================================================================
cat(">> Calculating Multivariate Environmental Similarity Surface (MESS)...\n")
env_backg_all <- rbind(backg_eval, pinus_henryi_occ)
env_backg_mat <- as.matrix(na.omit(raster::extract(predictors_final, env_backg_all)))

rasters_masked <- predictors_final * vn_raster
mess_curr <- dismo::mess(rasters_masked, env_backg_mat)
raster::writeRaster(mess_curr, filename = file.path(dir_out_maps, "Pinus henryi MESS current.tif"), format = "GTiff", overwrite = TRUE)

for (scen in scenarios) {
  for (gcm in gcms) {
    gcm_files <- list.files(file.path(dir_env_2070, gcm, scen), pattern = "\\.tif$", full.names = TRUE)
    gcm_stack <- raster::subset(raster::stack(gcm_files), selected_environment) * vn_raster
    mess_gcm <- dismo::mess(gcm_stack, env_backg_mat)
    out_mess_name <- sprintf("Pinus henryi MESS 2070 %s %s.tif", scen, gcm)
    raster::writeRaster(mess_gcm, filename = file.path(dir_out_maps, out_mess_name), format = "GTiff", overwrite = TRUE)
  }
}

# Export environmental values at presence points
pinus_henryi_env_vals <- as.matrix(na.omit(raster::extract(predictors_final, pinus_henryi_occ)))
openxlsx::write.xlsx(as.data.frame(pinus_henryi_env_vals), file = file.path(dir_data, "Pinus_henryi_env.xlsx"), overwrite = TRUE)

cat("\nPinus henryi modeling workflow completed successfully!\n")
