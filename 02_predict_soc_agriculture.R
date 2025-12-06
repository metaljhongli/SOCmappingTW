# ============================================================================
# 02_predict_soc_agriculture.R
#
# Purpose:
#   Use the trained Random Forest model for agricultural SOC to generate a
#   wall-to-wall SOC prediction map from gridded predictor rasters.
#
# Steps:
#   1) Read a 30 m SpatRaster stack of agriculture predictors
#   2) Export all cells to a CSV (cell index + predictors)
#   3) Load the best RF model bundle (.RDS) from cross-validation
#   4) Predict SOC in chunks (memory-safe)
#   5) Rebuild a GeoTIFF with the predicted SOC values
#
# Inputs (edit paths as needed):
#   - soc_agriculture.tif      : 14-band predictor stack (30 m) for agriculture
#   - rf_cv10_results_best*.rds: model bundle from 01_master_soc_pipeline.R
#
# Outputs:
#   - soc_predictors_30m.csv
#   - SOC_pred_30m_fromcsv.tif
#
# Author: Miguel Valdez Vasquez
# ============================================================================

rm(list = ls()); gc()

# -------------------------- USER CONFIG -------------------------------------

# Base directory for agriculture prediction
base_dir  <- "D:/Backup_Post_Doc_research/TARI/RFmodelSubmission3/agriculture/prediction/60mtrue"

# Predictor stack (SpatRaster with 14 layers for agriculture)
stack_tif <- file.path(base_dir, "soc_agriculture.tif")

# Where to write the intermediate CSV with predictors
csv_out   <- file.path(base_dir, "soc_predictors_60m.csv")

# Model bundle created by 01_master_soc_pipeline.R
# NOTE: adjust the file name if needed, for example:
#   "rf_cv10_results_best.rds"  or  "rf_cv10_results_best_fold.rds"
model_rds <- file.path(base_dir, "rf_cv10_results_best_fold.rds")

# Output SOC prediction raster
out_tif   <- file.path(base_dir, "SOC_pred_60m_fromcsv.tif")

# ------------------------ 1) RASTER -> CSV ----------------------------------

library(terra)
library(data.table)

cat("Reading agriculture predictor stack:\n  ", stack_tif, "\n", sep = "")
rs <- rast(stack_tif)
stopifnot(inherits(rs, "SpatRaster"))

if (nlyr(rs) != 14) {
  stop("Expected 14 predictor layers, but found: ", nlyr(rs))
}

# Give explicit names to the 14 predictor layers
# (must be consistent with how the RF model was trained)
names(rs) <- c(
  "elevation",    # digital elevation model
  "precipitation",
  "soiltype",
  "slope",
  "temperature",
  "green",
  "red",
  "gndvi",
  "gsavi",
  "ndvi",
  "savi",
  "evi",
  "bsi",
  "swir1"
)

# Optional: if 0 means "no soil class", set it to NA for soiltype
rs[["soiltype"]] <- classify(rs[["soiltype"]],
                             rcl = matrix(c(0, NA), ncol = 2, byrow = TRUE))

cat("Stack has", nlyr(rs), "layers with names:\n")
print(names(rs))

# Chunked export of the raster to CSV for memory safety
if (file.exists(csv_out)) file.remove(csv_out)

n_rows <- nrow(rs)
n_cols <- ncol(rs)
chunk  <- 2000L         # number of rows per chunk (tune if needed)

cat("Exporting predictors to CSV in chunks...\n")
readStart(rs); on.exit(readStop(rs), add = TRUE)
header_written <- FALSE

for (r0 in seq(1L, n_rows, by = chunk)) {
  nr <- min(chunk, n_rows - r0 + 1L)

  # vals: matrix (#cells_in_chunk x 14)
  vals  <- readValues(rs, row = r0, nrows = nr, mat = TRUE)

  # Cell indices for this block
  cells <- seq(
    cellFromRowCol(rs, r0, 1L),
    length.out = nr * n_cols
  )

  DT <- as.data.table(vals)
  setnames(DT, names(rs))
  DT[, cell := cells]
  setcolorder(DT, c("cell", names(rs)))

  fwrite(DT, csv_out, append = header_written)
  header_written <- TRUE

  if ((r0 %% (chunk * 25L)) == 1L) {
    cat(sprintf("  Rows %d..%d of %d written\n", r0, r0 + nr - 1L, n_rows))
  }
}

cat("Predictor CSV written:\n  ", csv_out, "\n\n", sep = "")

# ------------------------ 2) CSV -> PREDICT -> GeoTIFF ----------------------

rm(rs); gc()

library(ranger)

cat("Loading template raster for geometry:\n  ", stack_tif, "\n", sep = "")
templ <- rast(stack_tif)

cat("Loading RF model bundle:\n  ", model_rds, "\n", sep = "")
mdl_raw <- readRDS(model_rds)

# Your model bundle is expected to be a list with at least $best_model
stopifnot(is.list(mdl_raw), !is.null(mdl_raw$best_model))

mdl <- mdl_raw$best_model
if (!inherits(mdl, "ranger") && is.list(mdl) && !is.null(mdl$forest)) {
  class(mdl) <- "ranger"
}

# Predictor names expected by the model
needed <- NULL
if (!is.null(mdl$forest$independent.variable.names)) {
  needed <- mdl$forest$independent.variable.names
} else if (!is.null(mdl_raw$params$predictors_present)) {
  needed <- mdl_raw$params$predictors_present
}
if (is.null(needed)) {
  stop("Cannot determine predictor names from model (no independent.variable.names or params$predictors_present).")
}

cat("Model expects predictors:\n")
print(needed)

# Read CSV (predictors only)
cat("Reading predictor CSV for prediction:\n  ", csv_out, "\n", sep = "")
DT <- fread(csv_out, showProgress = TRUE)
stopifnot("cell" %in% names(DT))

# Soiltype ¡÷ factor if used
if ("soiltype" %in% needed) {
  cat("Converting 'soiltype' to factor...\n")
  DT[, soiltype := factor(as.integer(round(soiltype)))]
}

# Ensure we have all needed columns
missing_cols <- setdiff(needed, names(DT))
if (length(missing_cols) > 0) {
  stop("Missing predictors in CSV: ", paste(missing_cols, collapse = ", "))
}

X <- DT[, ..needed]
n <- nrow(X)
cat("Total rows to predict:", n, "\n")

# Predict in chunks
chunk_size <- 2e6L   # 2 million rows per chunk; adjust if memory is low
pred <- numeric(n)
i <- 1L

cat("Predicting SOC in chunks...\n")
while (i <= n) {
  j <- min(i + chunk_size - 1L, n)
  pr <- predict(mdl, data = X[i:j, , drop = FALSE])$predictions
  pred[i:j] <- as.numeric(pr)
  cat(sprintf("  Predicted rows %d..%d\n", i, j))
  i <- j + 1L
}

# Rebuild raster using cell index
cat("Rebuilding prediction raster and writing GeoTIFF...\n")
out <- rast(templ, nlyr = 1)
names(out) <- "SOC"

vals <- rep(NA_real_, ncell(templ))

if (max(DT$cell, na.rm = TRUE) > ncell(templ)) {
  stop("'cell' index exceeds template size in template raster.")
}

vals[DT$cell] <- pred
values(out)   <- vals

writeRaster(
  out,
  out_tif,
  overwrite = TRUE,
  wopt = list(gdal = "COMPRESS=LZW", datatype = "FLT4S")
)

cat("SOC prediction raster written:\n  ", out_tif, "\n", sep = "")
cat("Value range (SOC):\n")
print(global(out, "range", na.rm = TRUE))

na_cells <- freq(is.na(out), value = TRUE)
na_pct   <- if (is.null(na_cells)) 0 else 100 * na_cells[1, "count"] / ncell(out)
cat(sprintf("NA pixels: %.3f%%\n", na_pct))

cat("\nDone.\n")
