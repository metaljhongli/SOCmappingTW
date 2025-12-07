# # SOC Modeling and Prediction for Agricultural and Forest Lands in Taiwan

This repository provides the complete, reproducible workflow used to train and apply Random Forest (RF) models for predicting soil organic carbon (SOC) in agricultural and forest landscapes across Taiwan. The workflow includes data cleaning, 10-fold cross-validation, model diagnostics, and spatial prediction at 30–100 m resolution.

All code is written in **R (≥ 4.0)** and uses only open-source packages.

---

## ⭐ Overview

This repository contains three scripts:

1. **01_master_soc_pipeline.R**  
   Trains two RF models (agriculture + forest) using 10-fold cross-validation.  
   Exports model diagnostics, partial dependence plots, variable importance,  
   and saves the best model bundle for each land-use type.

2. **02_predict_soc_agriculture.R**  
   Applies the trained agriculture model to a 14-band raster stack  
   and generates a 30 or 60 m SOC prediction map.

3. **03_predict_soc_forest.R**  
   Applies the trained forest model to individual predictor GeoTIFFs  
   (one per variable), aligns them to a common grid, and generates a  
   forest SOC map (30 or 100 m).

This repository also contains **example dummy data** to allow reviewers to  
run the scripts without needing access to government-restricted SOC datasets.

---

## 📁 Folder Structure

```
taiwan-soc-modeling/
├─ README.md
├─ LICENSE
├─ R/
│  ├─ 01_master_soc_pipeline.R
│  ├─ 02_predict_soc_agriculture.R
│  └─ 03_predict_soc_forest.R
├─ example_data/
SOC_data_forest
SOC_data_agriculture
└─ (user-added folders, not included in repo)
   ├─ prediction_inputs/
   │  ├─ agriculture/   # 14-band raster stack (60 m)
   │  └─ forest/        # 14 GeoTIFF predictors
   └─ outputs/
      ├─ agriculture/
      └─ forest/
```

**Note:**  
Prediction rasters (HLS, DEM, NDVI, etc.) are not included because of file size  
and licensing restrictions. Users may substitute any equivalent predictor dataset.

---

## 🔧 Requirements

Install required packages:

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
Run the unified training pipeline:

```r
source("R/01_master_soc_pipeline.R")
```

The script:

- Reads `SOCagric_20pct_dummy_for_master.csv` and `SOCforest_20pct_dummy_for_master.csv`  
- Performs 10-fold cross-validation  
- Computes RMSE, MAE, R², and Lin’s CCC  
- Exports diagnostic plots  
- Saves the best model objects as:  
  - `rf_cv10_results_best_fold.rds` (agriculture)  
  - `rf_cv10_results_best.rds` (forest)

Multicollinearity filtering is optional (`do_collinearity = TRUE/FALSE`).  
It is disabled by default to match the published results.

---

### **2. Predict SOC for Agriculture (30-60 m)**

Ensure the directory contains:

```
soc_agriculture.tif   # 14-band stack
rf_cv10_results_best_fold.rds
```

Then run:

```r
source("R/02_predict_soc_agriculture.R")
```

Output:

```
SOC_pred_60m_fromcsv.tif
```

---

### **3. Predict SOC for Forest (30–100 m)**

Place the 14 predictor GeoTIFFs in the folder specified by `base_dir`.

If automatic name matching fails, the script will generate:

```
predictor_file_mapping.csv
```

Edit this file so each predictor is linked to the correct GeoTIFF.

Run:

```r
source("R/03_predict_soc_forest.R")
```

Output:

```
SOC_pred_forest_100m.tif
```

(or change resolution label in the script)

---

## 🧪 Example Data (for reviewers)

Because SOC field samples and predictor rasters are government-restricted and cannot be publicly released, the repository includes:

SOC_samples_agriculture.csv
SOC_samples_forest.csv

These files contain small random datasets with the correct structure  
(14 predictors + cell index), allowing the scripts to be executed  
without real data.

Reviewers and users can inspect the code logic, pipelines, and outputs  
without requiring access to restricted datasets.

## 🔒 Data Restrictions

The SOC field observations used for model training are owned by the  
**Taiwan Agricultural Research Institute (TARI)** and are not publicly  
accessible due to government data-use policies.  

All remote sensing predictors (HLS, DEM, NDVI, climate) are publicly available  
from NASA/USGS repositories but are not included here due to large file sizes.  
Users may download equivalent data to reproduce the spatial predictions.

---

## 📄 License

This code is released under the **MIT License** and is free for academic  
and research use. See `LICENSE` for details.

---

## 📚 Citation

If you use this workflow, please cite the associated publication:

**Syu, C.H; Valdez-Vásquez, M., et al. (2025). Soil organic carbon modeling across agricultural and forest landscapes in Taiwan using long-term field data, remote sensing, and Random Forest models. *Journal Name*.**

---

## 👤 Contact

For questions or collaborations:

**Chien-Hui Syu, Miguel Valdez Vásquez**  
Taiwan Agriculture Research Institute, Taiwan
National Central University, Taiwan  
Center for Space and Remote Sensing Research (CSRSR)
