# ============================================================================
# 03_predict_soc_forest.R
#
# Purpose:
#   Use the trained Random Forest model for forest SOC to generate a
#   wall-to-wall SOC prediction map from individual predictor rasters
#   (one GeoTIFF per variable).
#
# Steps:
#   1) List all .tif files in the prediction folder
#   2) Build a mapping Predictor -> filename (auto or from CSV)
#   3) Load each raster, harmonize to a common geometry, and stack them
#   4) Export all cells to a CSV (cell index + predictors)
#   5) Load the RF model bundle (.RDS) from 01_master_soc_pipeline.R
#   6) Predict SOC in chunks (memory-safe)
#   7) Rebuild a GeoTIFF with the predicted SOC values
#
# Inputs (edit paths as needed):
#   - Folder with 14 predictor GeoTIFFs for forest (same CRS, similar extent)
#   - rf_cv10_results_best.rds: model bundle from 01_master_soc_pipeline.R (forest)
#
# Outputs:
#   - forest_predictors_<res>.csv
#   - SOC_pred_forest_<res>.tif
#
# Author: Miguel Valdez Vasquez
# ============================================================================

rm(list = ls()); gc()

suppressPackageStartupMessages({
  library(terra)
  library(data.table)
  library(ranger)
})

# -------------------------- USER CONFIG -------------------------------------

# Base directory where the forest predictor TIFFs live
base_dir <- "C:/forest/prediction"

# Resolution label for file names (e.g. "30m")
res_label <- "30m"

# Model bundle (forest) created by 01_master_soc_pipeline.R
model_rds <- file.path(base_dir, "rf_cv10_results_best.rds")

# Output CSV and SOC GeoTIFF
csv_out <- file.path(base_dir, paste0("forest_predictors_", res_label, ".csv"))
out_tif <- file.path(base_dir, paste0("SOC_pred_forest_", res_label, ".tif"))

# Mapping CSV (optional, used to pin each predictor to a specific .tif)
map_csv <- file.path(base_dir, "predictor_file_mapping.csv")

# Predictor names (must match forest model training headers)
predictors <- c(
  "Precipitation", "Temperature", "Elevation", "Soiltype",
  "Blue", "Green", "BSI", "GNDVI", "GSAVI", "NBR",
  "NDVI", "NMDI", "OSAVI", "SAVI"
)

# --------------------- 0) LIST TIFFS IN FOLDER -----------------------------

all_tifs <- list.files(base_dir, pattern = "\\.(tif|tiff)$", full.names = TRUE)

if (!length(all_tifs)) {
  stop("No .tif/.tiff files found in: ", base_dir)
}

cat("Found", length(all_tifs), "tif files in folder:\n")
print(basename(all_tifs))

# --------------------- 1) BUILD PREDICTOR -> FILE MAPPING ------------------

if (file.exists(map_csv)) {
  # Use existing mapping CSV
  cat("\nUsing existing mapping file:\n  ", map_csv, "\n", sep = "")
  map <- fread(map_csv)
  if (!all(c("Predictor", "Filename") %in% names(map))) {
    stop("Mapping file exists but is missing columns Predictor, Filename: ", map_csv)
  }

  miss_files <- map[!file.exists(file.path(base_dir, Filename))]
  if (nrow(miss_files)) {
    stop(
      "These mapped files do not exist:\n",
      paste(file.path(base_dir, miss_files$Filename), collapse = "\n")
    )
  }

  files_by_var <- setNames(
    file.path(base_dir, map$Filename),
    map$Predictor
  )

} else {
  # Try automatic matching by case-insensitive "contains"
  cat("\nNo mapping CSV found. Attempting automatic matching...\n")

  auto_pick <- function(var) {
    bn <- basename(all_tifs)
    hits <- which(grepl(var, bn, ignore.case = TRUE))
    if (length(hits)) return(bn[hits[1]])
    character(0)
  }

  picked <- sapply(predictors, auto_pick, USE.NAMES = FALSE)
  files_by_var <- setNames(file.path(base_dir, picked), predictors)

  unmatched <- predictors[!file.exists(files_by_var)]
  if (length(unmatched)) {
    cat("Could not auto-match these predictors: ",
        paste(unmatched, collapse = ", "), "\n")

    # Write template CSV for the user to fill
    templ <- data.table(
      Predictor = predictors,
      Suggested = basename(picked),
      Filename  = ifelse(predictors %in% unmatched, "", basename(picked))
    )
    fwrite(templ, map_csv)

    stop(
      "Wrote mapping template:\n  ", map_csv,
      "\nPlease edit the 'Filename' column so that each Predictor\n",
      "points to the correct .tif in this folder, save, and rerun this script."
    )
  }
}

cat("\n=== Final Predictor ¡÷ file mapping ===\n")
print(data.table(Predictor = names(files_by_var),
                 Filename  = basename(files_by_var)))

# --------------------- 2) LOAD & HARMONIZE RASTERS -------------------------

cat("\nLoading and harmonizing predictor rasters...\n")

tmpl <- rast(files_by_var[[1]])
if (any(is.na(res(tmpl)))) stop("Template raster has invalid resolution.")
crs_tmpl <- crs(tmpl)

r_list <- vector("list", length(predictors))
names(r_list) <- predictors

for (nm in predictors) {
  f <- files_by_var[[nm]]
  cat("  -", nm, ":", basename(f), "\n")
  r <- rast(f)

  # Reproject if CRS differs
  if (!identical(crs(r), crs_tmpl)) {
    r <- project(r, tmpl,
                 method = if (nm == "Soiltype") "near" else "bilinear")
  }

  # Resample if geometry differs
  if (!compareGeom(r, tmpl, stopOnError = FALSE, crs = FALSE)) {
    r <- resample(r, tmpl,
                  method = if (nm == "Soiltype") "near" else "b
