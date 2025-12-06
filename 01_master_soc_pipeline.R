# ==== SOC end-to-end: Train (CV), Plots, Save RDS, PNGa?・TIFF (with OPTIONAL collinearity) ====
# - Trains two RF models (Agriculture & Forest) with 10-fold CV using ranger
# - Saves metrics, VI, PDP panel (top-4 excluding Soiltype), scatter plots (4200x3600 @ 600dpi)
# - Saves best model bundle (.rds) per run
# - Batch converts any PNGs in your manuscript figures folder to TIFF @ 600 dpi (no magick required)
# - OPTIONAL: multicollinearity diagnostics/pruning (OFF by default)

rm(list = ls()); gc()

suppressPackageStartupMessages({
  library(data.table)
  library(ranger)
  library(ggplot2)
  library(grid)
})

# ------------------------- CONFIG -------------------------------------------
cfg <- list(
  # Agriculture block
  agri = list(
    name       = "agriculture",
    base_dir   = "D:/Backup_Post_Doc_research/TARI/RFmodelSubmission3/agriculture/final_run",
    csv        = "SOCagriculture_clean3.csv",
    predictors = c("Precipitation","Temperature","Elevation","Slope","Soiltype","Green","Red",
                   "SWIR1","BSI","EVI","GNDVI","GSAVI","NDVI","NMDI","SAVI"),
    exclude_pdp = "Soiltype",      # PDP panel will skip Soiltype
    elev_max    = 2500,             # Elevation PDP x-max (m)
    gndvi_rng   = c(-0.5, 1.0)      # GNDVI PDP range
  ),

  # Forest block
  forest = list(
    name       = "forest",
    base_dir   = "D:/Backup_Post_Doc_research/TARI/RFmodelSubmission3/forest",
    csv        = "SOCforest_clean.csv",
    predictors = c("SAVI","GSAVI","OSAVI","Blue","NDVI","Green","BSI","Soiltype",
                   "Precipitation","Temperature","NMDI","GNDVI","NBR","NIR"),
    exclude_pdp = "Soiltype",
    elev_max    = NA_real_,         # not forcing Elevation limits by default
    gndvi_rng   = c(-0.5, 1.0)
  ),

  # Global RF/CV settings
  rf = list(num_trees = 1000, min_node_size = 5, n_folds = 10, seed = 123),

  # PNGa?・TIFF convert (final step)
  convert_png_dir = "D:/Backup_Post_Doc_research/TARI/manuscript/submission3/figures_results",

  # OPTIONAL collinearity (kept FALSE to match your "no filtering" default)
  do_collinearity = FALSE,    # set TRUE to enable pruning
  corr_thresh     = 0.90      # drop one of any pair with |r| >= threshold
)

# ------------------------- METRIC HELPERS -----------------------------------
rmse <- function(o,p) sqrt(mean((o-p)^2))
mae  <- function(o,p) mean(abs(o-p))
r2   <- function(o,p){ ssr <- sum((o-p)^2); sst <- sum((o-mean(o))^2); 1 - ssr/sst }
ccc  <- function(o,p){ m1<-mean(o); m2<-mean(p); v1<-var(o); v2<-var(p); c12<-cov(o,p); 2*c12/(v1+v2+(m1-m2)^2) }

# ------------------------- PDP HELPERS --------------------------------------
pdp_numeric <- function(model, data, var, grid_n=101L, elev_max=NA_real_, gndvi_rng=NULL, predictors=NULL){
  X <- data[, predictors, with=FALSE]
  if (identical(var, "Elevation") && is.finite(elev_max)) {
    lo <- max(0, as.numeric(quantile(X[[var]], 0.01, na.rm = TRUE)))
    xs <- seq(lo, elev_max, length.out = grid_n)
  } else if (identical(var, "GNDVI") && !is.null(gndvi_rng)) {
    xs <- seq(gndvi_rng[1], gndvi_rng[2], length.out = grid_n)
  } else {
    xs <- unique(quantile(X[[var]], probs = seq(0.01, 0.99, length.out = grid_n), na.rm = TRUE))
  }
  out <- lapply(xs, function(v){
    Xtmp <- X; Xtmp[[var]] <- v
    pd <- mean(predict(model, data=Xtmp)$predictions, na.rm=TRUE)
    data.table(var=var, x=as.numeric(v), pd=pd)
  })
  rbindlist(out)
}

make_pdp_panel <- function(best_model, DT, vi, outdir, exclude_var=NULL, elev_max=NA_real_, gndvi_rng=NULL, predictors=NULL){
  vi_order <- vi$Variable
  vars <- head(setdiff(vi_order, exclude_var), 4L)
  DT_pd <- copy(DT); if ("Soiltype" %in% names(DT_pd)) DT_pd$Soiltype <- droplevels(DT_pd$Soiltype)

  pd_list <- list(); chosen <- character(0)
  for (v in vars) {
    if (is.factor(DT_pd[[v]])) {
      levs <- levels(DT_pd[[v]]); keep <- names(sort(table(DT_pd[[v]]), TRUE))[seq_len(min(12L, length(levs)))]
      pdv <- rbindlist(lapply(keep, function(lv){
        Xtmp <- DT_pd[, predictors, with=FALSE]; Xtmp[[v]] <- factor(lv, levels=levs)
        data.table(var=v, x=lv, pd=mean(predict(best_model, data=Xtmp)$predictions, na.rm=TRUE))
      }))
    } else {
      pdv <- pdp_numeric(best_model, DT_pd, v, 101L, elev_max, gndvi_rng, predictors)
    }
    if (nrow(pdv) > 1) { pd_list[[v]] <- pdv; chosen <- c(chosen, v) }
  }
  if (length(chosen) < 4L) {
    pool <- setdiff(vi_order, c(exclude_var, chosen))
    for (v in pool) {
      if (is.factor(DT_pd[[v]])) {
        levs <- levels(DT_pd[[v]]); keep <- names(sort(table(DT_pd[[v]]), TRUE))[seq_len(min(12L, length(levs)))]
        pdv <- rbindlist(lapply(keep, function(lv){
          Xtmp <- DT_pd[, predictors, with=FALSE]; Xtmp[[v]] <- factor(lv, levels=levs)
          data.table(var=v, x=lv, pd=mean(predict(best_model, data=Xtmp)$predictions, na.rm=TRUE))
        }))
      } else {
        pdv <- pdp_numeric(best_model, DT_pd, v, 101L, elev_max, gndvi_rng, predictors)
      }
      if (nrow(pdv) > 1) { pd_list[[v]] <- pdv; chosen <- c(chosen, v); if (length(chosen) >= 4L) break }
    }
  }
  vars <- head(chosen, 4L)
  labs_abcd <- paste0(letters[seq_along(vars)], ") ", vars)

  theme_compact <- theme_minimal(base_size = 8) +
    theme(
      legend.position = "none",
      plot.title  = element_text(hjust = 0.5, face = "bold", size = 9, margin = margin(b = 4), colour = "black"),
      axis.title.x= element_text(size = 8, margin = margin(t = 4), colour = "black"),
      axis.title.y= element_text(size = 8, margin = margin(r = 4), colour = "black"),
      axis.text   = element_text(size = 7, colour = "black"),
      panel.border= element_rect(colour = "black", fill = NA, linewidth = 0.6),
      plot.margin = margin(4, 6, 4, 6)
    )

  make_plot <- function(v, ttl){
    pdv <- pd_list[[v]]
    if (is.factor(DT_pd[[v]])) {
      ggplot(pdv, aes(x=x, y=pd, group=1)) +
        geom_point(size=0.8) + geom_line(linewidth=0.5) +
        labs(title=ttl, x=v, y="Partial dependence") +
        theme_compact + theme(axis.text.x = element_text(size=6, angle=90, vjust=0.5, hjust=1))
    } else {
      p <- ggplot(pdv, aes(x=x, y=pd)) + geom_line(linewidth=0.5) +
        labs(title=ttl, x=v, y="Partial dependence") + theme_compact
      if (identical(v, "GNDVI") && !is.null(gndvi_rng)) {
        p <- p + scale_x_continuous(limits = gndvi_rng, breaks = seq(gndvi_rng[1], gndvi_rng[2], by=0.2),
                                    labels = scales::label_number(accuracy=0.1), expand = expansion(mult=c(0.01,0.03)))
      } else if (v %in% c("Precipitation","Temperature","Elevation")) {
        if (identical(v, "Elevation") && is.finite(elev_max)) {
          p <- p + scale_x_continuous(limits=c(0, elev_max), breaks=seq(0, elev_max, by=500),
                                      labels = scales::label_number(accuracy=1), expand = expansion(mult=c(0.01,0.03)))
        } else {
          p <- p + scale_x_continuous(labels = scales::label_number(accuracy=1), expand = expansion(mult=c(0.01,0.03)))
        }
      } else {
        p <- p + scale_x_continuous(labels = scales::label_number(accuracy=0.01, big.mark=","), expand = expansion(mult=c(0.01,0.03)))
      }
      p + scale_y_continuous(labels = scales::label_number(accuracy=0.01, big.mark=","), expand = expansion(mult=c(0.02,0.05)))
    }
  }

  plots <- Map(make_plot, vars, labs_abcd)
  png(file.path(outdir, "pdp_top4_panel.png"), width=7, height=5, units="in", res=600)
  grid.newpage(); pushViewport(viewport(layout=grid.layout(2,2)))
  vp <- function(r,c) viewport(layout.pos.row=r, layout.pos.col=c)
  if (length(plots)>=1) print(plots[[1]], vp=vp(1,1))
  if (length(plots)>=2) print(plots[[2]], vp=vp(1,2))
  if (length(plots)>=3) print(plots[[3]], vp=vp(2,1))
  if (length(plots)>=4) print(plots[[4]], vp=vp(2,2))
  grid.rect(gp=gpar(col="black", fill=NA, lwd=1.5))
  dev.off()
}

# ------------------------- MAIN RUNNER --------------------------------------
run_pipeline <- function(block){
  # Output folder
  stamp  <- format(Sys.time(), "%Y%m%d_%H%M%S")
  outdir <- file.path(block$base_dir, paste0("rf_cv10_", block$name, "_outputs_", stamp))
  dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
  csv_path <- file.path(block$base_dir, block$csv)

  # Load + typing
  DT <- fread(csv_path)
  stopifnot("SOC" %in% names(DT))
  missing <- setdiff(block$predictors, names(DT))
  if (length(missing)) stop("Missing predictors in CSV: ", paste(missing, collapse=", "))
  DT <- DT[, c("SOC", block$predictors), with=FALSE]
  DT <- DT[complete.cases(DT)]
  if (!is.factor(DT$Soiltype)) { if (is.numeric(DT$Soiltype)) DT[, Soiltype := as.integer(round(Soiltype))]; DT[, Soiltype := factor(Soiltype)] }
  for (nm in setdiff(block$predictors, "Soiltype")) if (!is.numeric(DT[[nm]])) DT[[nm]] <- as.numeric(DT[[nm]])
  if (!is.numeric(DT$SOC)) DT[, SOC := as.numeric(SOC)]

  # ---------- OPTIONAL: multicollinearity diagnostics/pruning ---------------
  if (isTRUE(cfg$do_collinearity)) {
    num_vars <- setdiff(block$predictors, "Soiltype")
    M <- try(cor(DT[, ..num_vars], use = "pairwise.complete.obs"), silent = TRUE)
    if (!inherits(M, "try-error")) {
      keep <- character(0); drop <- character(0)
      for (v in colnames(M)) {
        if (v %in% drop) next
        keep <- c(keep, v)
        high <- setdiff(names(which(abs(M[v, ]) >= cfg$corr_thresh)), v)
        drop <- unique(c(drop, high))
      }
      pruned <- unique(keep)
      message("[", block$name, "] Collinearity pruning ON (|r|>=", cfg$corr_thresh,
              ") a?・ kept: ", paste(c(pruned, intersect("Soiltype", block$predictors)), collapse=", "),
              "; dropped: ", paste(setdiff(num_vars, pruned), collapse=", "))
      block$predictors <- c(pruned, intersect("Soiltype", block$predictors))
      # Retain only the selected columns in DT
      DT <- DT[, c("SOC", block$predictors), with=FALSE]
    } else {
      warning("[", block$name, "] Correlation matrix failed; skipping collinearity step.")
    }
  } else {
    message("[", block$name, "] Collinearity pruning OFF (as requested)")
  }

  # ---------- CV + model fit ------------------------------------------------
  n <- nrow(DT); set.seed(cfg$rf$seed)
  fold_id <- sample(rep(1:cfg$rf$n_folds, length.out=n))
  mtry <- floor(sqrt(length(block$predictors)))

  metrics <- data.table(Fold=integer(), n_train=integer(), n_test=integer(), R2=numeric(), RMSE=numeric(), MAE=numeric(), CCC=numeric())
  pred_rows <- vector("list", cfg$rf$n_folds)
  best_fold <- NA_integer_; best_rmse <- Inf; best_model <- NULL

  for (k in 1:cfg$rf$n_folds) {
    trn <- DT[fold_id != k]; tst <- DT[fold_id == k]
    rf <- ranger(
      SOC ~ ., data = trn[, c("SOC", block$predictors), with=FALSE],
      num.trees = cfg$rf$num_trees, mtry = mtry, min.node.size = cfg$rf$min_node_size,
      importance = "permutation", splitrule = "variance", write.forest = TRUE,
      respect.unordered.factors = "order", num.threads = max(1, parallel::detectCores()-1),
      seed = cfg$rf$seed + k
    )
    pr  <- predict(rf, data = tst[, block$predictors, with=FALSE])$predictions
    obs <- tst$SOC
    m <- data.table(Fold=k, n_train=nrow(trn), n_test=nrow(tst), R2=r2(obs,pr), RMSE=rmse(obs,pr), MAE=mae(obs,pr), CCC=ccc(obs,pr))
    metrics <- rbind(metrics, m)
    pred_rows[[k]] <- data.table(row_id=which(fold_id==k), Fold=k, Observed=obs, Predicted=pr)
    if (m$RMSE < best_rmse) { best_rmse <- m$RMSE; best_fold <- k; best_model <- rf }
    cat(sprintf("[%s] Fold %02d/%02d | R2=%.3f RMSE=%.3f MAE=%.3f CCC=%.3f\n", block$name, k, cfg$rf$n_folds, m$R2, m$RMSE, m$MAE, m$CCC))
  }

  preds_all <- rbindlist(pred_rows)
  overall <- data.table(Fold=0L, n_train=NA_integer_, n_test=nrow(preds_all),
                        R2=r2(preds_all$Observed, preds_all$Predicted),
                        RMSE=rmse(preds_all$Observed, preds_all$Predicted),
                        MAE=mae(preds_all$Observed, preds_all$Predicted),
                        CCC=ccc(preds_all$Observed, preds_all$Predicted))
  metrics_all <- rbind(metrics, overall)

  # ---------- Variable Importance ------------------------------------------
  vi <- data.table(Variable = names(ranger::importance(best_model)), Importance = as.numeric(ranger::importance(best_model)))
  data.table::setorder(vi, -Importance)
  fwrite(vi, file.path(outdir, "variable_importance.csv"))

  p_vi <- ggplot(vi, aes(x=reorder(Variable, Importance), y=Importance)) +
    geom_col(fill="grey70", color="black", linewidth=0.2) + coord_flip() +
    labs(x=NULL, y="Permutation importance",
         title=paste0(ifelse(block$name=="agriculture","a) ","b) "), "Variable importance (SOC ", tools::toTitleCase(block$name), ")")) +
    scale_y_continuous(expand = expansion(mult=c(0,0.02)), limits = c(0, NA)) +
    theme_minimal(base_size=10) +
    theme(plot.title=element_text(hjust=0.5, colour="black"),
          axis.title=element_text(colour="black"), axis.text=element_text(colour="black"),
          panel.border=element_rect(colour="black", fill=NA, linewidth=0.5), panel.grid.minor=element_blank())
  ggsave(file.path(outdir, paste0("variable_importance_", block$name, ".png")), p_vi, width=6, height=5, dpi=600)

  # ---------- PDP Panel -----------------------------------------------------
  make_pdp_panel(best_model, DT, vi, outdir,
                 exclude_var = block$exclude_pdp,
                 elev_max    = block$elev_max,
                 gndvi_rng   = block$gndvi_rng,
                 predictors  = block$predictors)

  # ---------- Scatter (overall + per-fold) ---------------------------------
  lim_all <- range(c(preds_all$Observed, preds_all$Predicted), na.rm = TRUE)
  theme_fullcanvas <- theme_minimal(base_size=16) +
    theme(plot.margin=margin(0,0,0,0), panel.border=element_rect(colour="black", fill=NA, linewidth=0.5),
          panel.grid.minor=element_blank(), axis.title=element_text(colour="black"), axis.text=element_text(colour="black"),
          plot.title=element_text(hjust=0.5, face="plain", size=16))
  # Overall metrics box (includes MAPE)
  lab_overall <- sprintf("RMSE = %.3f\nMAE = %.3f\nMAPE%% = %.2f\nR\u00B2 = %.3f\nLin's CCC = %.3f",
                         overall$RMSE, overall$MAE,
                         mean(abs((preds_all$Observed - preds_all$Predicted)/preds_all$Observed))*100,
                         overall$R2, overall$CCC)
  dx <- 0.04*diff(lim_all); dy <- 0.04*diff(lim_all)

  p_all <- ggplot(preds_all, aes(Observed, Predicted)) +
    geom_point(alpha=0.6, shape=4, stroke=0.5, size=1.2) +
    geom_abline(slope=1, intercept=0, linetype=2) +
    labs(title="Observed vs Predicted", x="Observed SOC", y="Predicted SOC") +
    scale_x_continuous(limits=lim_all, expand=c(0,0)) + scale_y_continuous(limits=lim_all, expand=c(0,0)) +
    annotate("label", x=lim_all[2]-dx, y=lim_all[1]+dy, hjust=1, vjust=0,
             label=lab_overall, label.size=0.3, size=4.5, fill="white", color="black") +
    theme_fullcanvas
  ggsave(file.path(outdir, "scatter_all_folds.png"), p_all, width=4200, height=3600, units="px", dpi=600, limitsize=FALSE)

  # Per-fold
  fold_ids <- sort(unique(preds_all$Fold))
  for (k in fold_ids) {
    dd <- preds_all[Fold==k]; m  <- metrics[Fold==k]
    lab_fold <- sprintf("RMSE = %.3f\nMAE = %.3f\nMAPE%% = %.2f\nR\u00B2 = %.3f\nLin's CCC = %.3f",
                        m$RMSE, m$MAE, mean(abs((dd$Observed - dd$Predicted)/dd$Observed))*100, m$R2, m$CCC)
    p <- ggplot(dd, aes(Observed, Predicted)) +
      geom_point(alpha=0.7, shape=4, stroke=0.5, size=1.2) +
      geom_abline(slope=1, intercept=0, linetype=2) +
      labs(title="Observed vs Predicted", x="Observed SOC", y="Predicted SOC") +
      scale_x_continuous(limits=lim_all, expand=c(0,0)) + scale_y_continuous(limits=lim_all, expand=c(0,0)) +
      annotate("label", x=lim_all[2]-dx, y=lim_all[1]+dy, hjust=1, vjust=0,
               label=lab_fold, label.size=0.3, size=4.5, fill="white", color="black") +
      theme_fullcanvas
    ggsave(file.path(outdir, sprintf("scatter_fold_%02d.png", k)), p, width=4200, height=3600, units="px", dpi=600, limitsize=FALSE)
  }

  # ---------- Save artifacts bundle ----------------------------------------
  fwrite(metrics_all, file.path(outdir, "metrics_cv10.csv"))
  fwrite(preds_all,   file.path(outdir, "predictions_cv10.csv"))
  results <- list(
    metrics    = metrics_all,
    best_fold  = best_fold,
    best_model = best_model,
    preds_all  = preds_all,
    vi         = vi,
    params     = list(K=cfg$rf$n_folds, R=1, num_trees=cfg$rf$num_trees, min_node_size=cfg$rf$min_node_size,
                      mtry=floor(sqrt(length(block$predictors))), target="SOC", predictors_present=block$predictors, seed=cfg$rf$seed)
  )
  saveRDS(results, file.path(outdir, "rf_cv10_results_best.rds"))

  writeLines(c(
    paste0("Random Forest (ranger) 10-fold CV outputs (SOC ", tools::toTitleCase(block$name), ")"),
    paste("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste("Data:", csv_path),
    paste("Predictors:", paste(block$predictors, collapse=", ")),
    sprintf("RF: trees=%d, mtry=%d, min.node.size=%d", cfg$rf$num_trees, floor(sqrt(length(block$predictors))), cfg$rf$min_node_size),
    "Metrics: RMSE, MAE, R2, Lin's CCC",
    "Artifacts:",
    " - metrics_cv10.csv (fold 1..10 + overall row Fold=0)",
    " - predictions_cv10.csv (Observed/Predicted with fold id)",
    " - variable_importance.csv / variable_importance_*.png",
    " - pdp_top4_panel.png (600 dpi)",
    " - scatter_all_folds.png and scatter_fold_XX.png (4200x3600 @ 600 dpi)",
    " - rf_cv10_results_best.rds (best model + metadata)"
  ), file.path(outdir, "README.txt"))

  cat("\n[", block$name, "] Done. Outputs in:\n  ", outdir, "\n", sep="")
}

# ------------------------- EXECUTE BOTH -------------------------------------
run_pipeline(cfg$agri)
run_pipeline(cfg$forest)

# ------------------------- PNG a?・ TIFF (600 dpi) ------------------------------
# Uses base R 'png' + TIFF device; writes to ./tiff subfolder next to PNGs
if (!requireNamespace("png", quietly = TRUE)) install.packages("png")
library(png); library(grid)

in_dir  <- cfg$convert_png_dir
out_dir <- file.path(in_dir, "tiff")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pngs <- list.files(in_dir, pattern = "\\.png$", full.names = TRUE, ignore.case = TRUE)
if (length(pngs)) {
  for (f in pngs) {
    img <- try(readPNG(f, native = TRUE), silent = TRUE)
    if (inherits(img, "try-error")) { message("Skip (read error): ", f); next }
    wh <- attr(img, "native")
    if (is.null(wh)) { arr <- readPNG(f); h <- dim(arr)[1]; w <- dim(arr)[2]; img <- as.raster(arr) } else { w <- wh[1]; h <- wh[2] }
    out_file <- file.path(out_dir, sub("\\.png$", ".tiff", basename(f), ignore.case = TRUE))
    grDevices::tiff(filename = out_file, width = w, height = h, units = "px", res = 600, compression = "lzw")
    grid.newpage(); pushViewport(viewport(width = unit(1, "npc"), height = unit(1, "npc")))
    grid.raster(img, interpolate = FALSE)
    dev.off()
    message("TIFF -> ", out_file, " [", w, "x", h, " px @ 600 dpi]")
  }
  message("PNGa?・TIFF conversion complete: ", out_dir)
} else {
  message("No PNGs found for conversion in: ", in_dir)
}
