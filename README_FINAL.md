# SOC Modeling and Prediction for Agricultural and Forest Lands in Taiwan

This repository provides the complete, reproducible workflow used to train and apply Random Forest (RF) models for predicting soil organic carbon (SOC) across agricultural and forest landscapes in Taiwan. The workflow includes data preprocessing, 10-fold cross-validation, model evaluation, and spatial prediction at 30–100 m resolution.

All code is written in **R (≥ 4.0)** using open-source packages.

---

## ⭐ Overview

This repository contains three core scripts:

### **01_master_soc_pipeline.R**
- Trains the agriculture and forest SOC models using 10-fold cross-validation  
- Computes RMSE, MAE, R², and Lin’s CCC  
- Exports diagnostic plots, variable importance, and partial dependence plots  
- Saves the best model objects for each land-use type  

### **02_predict_soc_agriculture.R**
- Applies the trained agriculture RF model to a 14-band raster stack  
- Predictors are derived from **Landsat-5 TM** and **Landsat-8 OLI** imagery (30 m)  
- Generates a **30 m SOC prediction map**

### **03_predict_soc_forest.R**
- Applies the trained forest RF model to individual predictor GeoTIFFs  
- Aligns predictors to a common grid  
- Generates **30–100 m forest SOC maps**, depending on the spatial resolution of available predictors  
  (Landsat spectral indices, DEM, terrain variables, climate data)

The repository includes example datasets so reviewers can run the full workflow without needing access to government-restricted SOC or predictor data.

---

## 📁 Folder Structure

```
SOCmappingTW/
├─ README.md
├─ LICENSE
├─ 01_master_soc_pipeline.R
├─ 02_predict_soc_agriculture.R
├─ 03_predict_soc_forest.R
├─ SOCagric_20pct_dummy_for_master.csv
├─ SOCforest_dummy_for_master.csv
└─ (user-added folders, not included in repo)
   ├─ prediction_inputs/
   │   ├─ agriculture/   # 14-band 30 m raster stack
   │   └─ forest/        # individual predictor rasters (30–100 m)
   └─ outputs/
       ├─ agriculture/
       └─ forest/
```

**Note:** Prediction rasters (Landsat spectral bands, DEM, NDVI/EVI, climate layers) are not included due to file size and licensing restrictions. Users may substitute equivalent datasets.

---

## 🔧 Requirements

Install required R packages:

```r
install.packages(c(
  "data.table",
  "ranger",
  "terra",
  "ggplot2",
  "grid",
  "png"
))
```

---

## 🚀 How to Run the Workflow

### **1. Train the SOC Models**

Run the unified modeling pipeline:

```r
source("01_master_soc_pipeline.R")
```

The script:

- Loads the included test datasets:  
  `SOCagric_20pct_dummy_for_master.csv`  
  `SOCforest_dummy_for_master.csv`
- Runs 10-fold cross-validation  
- Computes RMSE, MAE, R², Lin’s CCC  
- Exports diagnostic plots  
- Saves best models:
  - `rf_cv10_results_best_fold.rds` (agriculture)
  - `rf_cv10_results_best.rds` (forest)

---

### **2. Predict SOC for Agriculture (30 m)**

Place in your working directory:

- `soc_agriculture.tif` — 14-band 30 m raster stack  
- `rf_cv10_results_best_fold.rds` — agriculture model

Then run:

```r
source("02_predict_soc_agriculture.R")
```

Output:

```
SOC_pred_agriculture_30m.tif
```

---

### **3. Predict SOC for Forest (30–100 m)**

Place the predictor GeoTIFFs in:

```
prediction_inputs/forest/
```

If names do not match automatically, the script will generate:

```
predictor_file_mapping.csv
```

Edit this file so each predictor matches the correct raster.

Then run:

```r
source("03_predict_soc_forest.R")
```

Output (depending on resolution):

```
SOC_pred_forest_30m.tif
SOC_pred_forest_100m.tif
```

---

## 🧪 Example Data (for Reviewers)

The repository includes fully functional example datasets:

- `SOCagric_20pct_dummy_for_master.csv`  
- `SOCforest_dummy_for_master.csv`

These contain:

- A **20% authorized subset of SOC field observations** (released with TARI approval)  
- Synthetic (dummy) predictor variables with the correct structure  
- All required columns for running the full modeling workflow

This ensures **complete reproducibility of the code** without exposing restricted government data.

---

## 🔒 Data Restrictions

The full SOC field dataset is owned by the  
**Taiwan Agricultural Research Institute (TARI)** and cannot be publicly released due to government data-use regulations.

A **20% authorized SOC subset** is provided in this repository with institutional permission.  
This subset preserves the structure and statistical behavior of the full dataset and is sufficient to run the complete workflow.

Satellite-based predictors (Landsat-5, Landsat-8, DEM, climate layers) are publicly available from NASA/USGS archives but are not included due to file size.  
Processing instructions are provided in the manuscript Methods section.

---

## 📄 License

This code is released under the **MIT License**.

---

## 📚 Citation

Syu, C.-H., Valdez-Vásquez, M., et al. (2025).  
*Soil organic carbon modeling across agricultural and forest landscapes in Taiwan using long-term field data, remote sensing, and Random Forest models.* Journal Name.

---

## 👤 Contact

For questions or collaborations:

**Chien-Hui Syu** – Taiwan Agricultural Research Institute  
**Miguel Valdez Vásquez** – National Central University  
Center for Space and Remote Sensing Research (CSRSR)
