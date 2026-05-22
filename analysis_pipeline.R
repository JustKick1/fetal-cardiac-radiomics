# ============================================================
# Feature Selection + Cross-Validation Pipeline
# ============================================================
# Repeated stratified 5-fold cross-validation for:
#   1) Three-class one-vs-all classification (Control / GroupA / GroupB)
#   2) Binary classification (GroupA vs GroupB) with radiomics, clinical,
#      and combined feature sets
#   3) Sensitivity analyses (LHO subgroup; postnatal-echo-confirmed subgroup)
# ============================================================

# ---- 0. Packages ----
required_pkgs <- c("readxl","dplyr","pROC","glmnet","ranger","xgboost","tidyr")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
}
library(readxl); library(dplyr); library(pROC)
library(glmnet); library(ranger); library(xgboost); library(tidyr)

RNGkind("L'Ecuyer-CMRG")
set.seed(42)

# ============================================================
# ---- 1. Paths and parameters ----
# ============================================================
data_path     <- "./data/新Total.xlsx"
icc_path      <- "./data/ICC_select.csv"
out_dir       <- "./output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

clinical_vars <- c("LV_GLS","MV.annulus","RV.LV","GA")
N_REPEATS     <- 10
N_FOLDS       <- 5
LHO_types     <- c("aortic coarctation","aortic arch hypoplasia","aortic stenosis")

# Group coding: 0 = Control, 1 = GroupA, 2 = GroupB
# Instrument coding: 0 = E10, 1 = E8

# ============================================================
# ---- 2. Load data ----
# ============================================================
cat("Loading data...\n")
data <- read_excel(data_path)
cat(sprintf("Data dimensions: %d rows x %d columns\n", nrow(data), ncol(data)))

ad_col_idx <- which(toupper(names(data)) == "AD")
if (length(ad_col_idx)==0) {
  ad_col_idx <- 30
  cat(sprintf("Column 'AD' not found; using column %d as starting index\n", ad_col_idx))
} else {
  cat(sprintf("Radiomic features start at column %d (%s)\n", ad_col_idx, names(data)[ad_col_idx]))
}
all_radio_in_data <- names(data)[ad_col_idx:ncol(data)]

if (file.exists(icc_path) && file.info(icc_path)$size > 100) {
  icc_res    <- read.csv(icc_path, stringsAsFactors=FALSE)
  radio_vars <- intersect(icc_res$Feature, all_radio_in_data)
  cat(sprintf("Features passing ICC threshold: %d\n", length(radio_vars)))
  if (length(radio_vars)==0) stop("No overlap between ICC feature names and data column names")
} else {
  cat("Warning: ICC file not found; using all radiomic features\n")
  radio_vars <- all_radio_in_data
}

missing_c <- setdiff(clinical_vars, names(data))
if (length(missing_c)>0) stop(paste("Missing clinical columns:", paste(missing_c, collapse=", ")))

data$Group <- as.integer(data$Group)
cat(sprintf("Group distribution - Control: %d, GroupA: %d, GroupB: %d\n",
            sum(data$Group==0), sum(data$Group==1), sum(data$Group==2)))

if ("Disease" %in% names(data)) {
  cat("\nDisease distribution in GroupB:\n")
  print(sort(table(data$Disease[data$Group==2], useNA="ifany"), decreasing=TRUE))
}

# ============================================================
# ---- 3. Helper functions ----
# ============================================================

# ---- 3a. Stratified fold assignment ----
make_folds <- function(df, k, seed_val) {
  set.seed(seed_val)
  ga_q <- unique(quantile(df$GA, probs=c(0,.25,.5,.75,1), na.rm=TRUE))
  ga_group <- if (length(ga_q)<3) {
    ifelse(df$GA >= median(df$GA, na.rm=TRUE), "GA_H", "GA_L")
  } else {
    as.character(cut(df$GA, breaks=ga_q, include.lowest=TRUE,
                     labels=paste0("Q", seq_len(length(ga_q)-1))))
  }
  instrument_group <- if ("Instrument" %in% names(df)) {
    ins <- as.character(df$Instrument); ins[is.na(ins)] <- "Unknown"; ins
  } else "NoInstrument"
  strata   <- paste(df$Group, ga_group, instrument_group, sep="_")
  fold_ids <- integer(nrow(df))
  for (s in sort(unique(strata))) {
    idx <- which(strata==s)
    fold_ids[idx] <- if (length(idx)<k) sample(1:k, length(idx), replace=TRUE) else
      sample(rep_len(1:k, length(idx)))
  }
  fold_ids
}

# ---- 3b. Within-fold feature selection (Wilcoxon + Spearman pruning) ----
select_features <- function(X_train, y_train, feat_names,
                            p_thresh=0.05, cor_thresh=0.8) {
  if (length(feat_names)==0) return(character(0))
  p_vals <- setNames(sapply(feat_names, function(f) {
    tryCatch(wilcox.test(X_train[[f]] ~ y_train)$p.value, error=function(e) 1.0)
  }), feat_names)
  selected <- names(p_vals)[!is.na(p_vals) & p_vals < p_thresh]
  if (length(selected)<=1) return(selected)
  cor_mat <- cor(X_train[,selected,drop=FALSE], method="spearman",
                 use="pairwise.complete.obs")
  to_drop <- c()
  for (i in seq_len(nrow(cor_mat)-1)) {
    for (j in seq(i+1, ncol(cor_mat))) {
      fi <- selected[i]; fj <- selected[j]
      if (fi %in% to_drop || fj %in% to_drop) next
      if (!is.na(cor_mat[i,j]) && abs(cor_mat[i,j]) >= cor_thresh) {
        if (p_vals[fi] <= p_vals[fj]) to_drop <- c(to_drop, fj)
        else                          to_drop <- c(to_drop, fi)
      }
    }
  }
  setdiff(selected, to_drop)
}

# ---- 3c. Model training (LASSO / RF / XGBoost) with inverse-frequency weighting ----
train_predict <- function(X_tr, y_tr, X_te, model_type, seed_val=1) {
  y_num <- as.integer(y_tr=="1")
  n_pos <- max(sum(y_num==1),1); n_neg <- max(sum(y_num==0),1)
  w     <- ifelse(y_num==1, n_neg/n_pos, 1.0)
  tryCatch({
    if (model_type=="lasso") {
      set.seed(seed_val)
      fit <- cv.glmnet(as.matrix(X_tr), y_tr, family="binomial",
                       alpha=1, weights=w, nfolds=5, type.measure="auc")
      as.numeric(predict(fit, as.matrix(X_te), s="lambda.min", type="response"))
    } else if (model_type=="rf") {
      df_tr <- cbind(as.data.frame(X_tr), .y=factor(y_tr, levels=c("0","1")))
      fit   <- ranger(.y~., data=df_tr, probability=TRUE,
                      num.trees=500, case.weights=w, seed=seed_val)
      predict(fit, as.data.frame(X_te))$predictions[,"1"]
    } else if (model_type=="xgb") {
      set.seed(seed_val)
      dtr <- xgb.DMatrix(as.matrix(X_tr), label=y_num)
      dte <- xgb.DMatrix(as.matrix(X_te))
      fit <- xgb.train(
        params=list(objective="binary:logistic", eval_metric="auc",
                    max_depth=3, eta=0.1, scale_pos_weight=n_neg/n_pos,
                    seed=seed_val),
        data=dtr, nrounds=100, verbose=0)
      predict(fit, dte)
    }
  }, error=function(e) rep(NA_real_, nrow(X_te)))
}

# ---- 3d. Performance metrics (AUC, sensitivity, specificity, accuracy) ----
compute_metrics <- function(pred, label) {
  if (all(is.na(pred)) || length(unique(label))<2)
    return(list(auc=NA_real_, sens=NA_real_, spec=NA_real_, acc=NA_real_))
  roc_o <- pROC::roc(response=label, predictor=pred, quiet=TRUE)
  co    <- pROC::coords(roc_o, "best",
                        ret=c("sensitivity","specificity","accuracy"),
                        best.method="youden")
  list(auc  = as.numeric(pROC::auc(roc_o)),
       sens = as.numeric(co$sensitivity[1]),
       spec = as.numeric(co$specificity[1]),
       acc  = as.numeric(co$accuracy[1]))
}

# ---- 3e. Main CV function ----
run_cv <- function(df_subset, label_vec, feat_pool,
                   n_repeats, n_folds, base_seed=42) {
  models      <- c("lasso","rf","xgb")
  result_rows <- list()
  feat_count  <- if (length(feat_pool)>0)
    setNames(integer(length(feat_pool)), feat_pool) else integer(0)
  n_samp    <- nrow(df_subset)
  oof_store <- setNames(
    lapply(models, function(m) matrix(NA_real_, n_samp, n_repeats)), models)
  
  for (rep_i in seq_len(n_repeats)) {
    fold_ids <- make_folds(df_subset, n_folds, seed_val=base_seed + rep_i*100)
    for (fold_j in seq_len(n_folds)) {
      tr_idx <- which(fold_ids!=fold_j); te_idx <- which(fold_ids==fold_j)
      y_tr   <- label_vec[tr_idx]; y_te <- label_vec[te_idx]
      if (length(unique(y_te))<2) next
      X_tr <- df_subset[tr_idx, feat_pool, drop=FALSE]
      X_te <- df_subset[te_idx, feat_pool, drop=FALSE]
      sel_feats <- select_features(X_tr, y_tr, feat_pool)
      if (length(sel_feats)>0 && length(feat_count)>0) {
        m2 <- intersect(sel_feats, names(feat_count))
        feat_count[m2] <- feat_count[m2] + 1L
      }
      if (length(sel_feats)==0) next
      X_tr_s <- X_tr[,sel_feats,drop=FALSE]; X_te_s <- X_te[,sel_feats,drop=FALSE]
      preds <- setNames(lapply(seq_along(models), function(mi) {
        model_seed <- base_seed + rep_i*10000 + fold_j*100 + mi
        train_predict(X_tr_s, y_tr, X_te_s, models[mi], seed_val=model_seed)
      }), models)
      for (m in models) {
        if (!all(is.na(preds[[m]])))
          oof_store[[m]][te_idx, rep_i] <- preds[[m]]
        met <- compute_metrics(preds[[m]], y_te)
        result_rows[[length(result_rows)+1]] <- data.frame(
          repeat_id=rep_i, fold_id=fold_j, model=m,
          auc=met$auc, sens=met$sens, spec=met$spec, acc=met$acc,
          n_feats=length(sel_feats), stringsAsFactors=FALSE)
      }
    }
  }
  oof_mean <- setNames(
    lapply(models, function(m) rowMeans(oof_store[[m]], na.rm=TRUE)), models)
  empty_m <- data.frame(repeat_id=integer(), fold_id=integer(), model=character(),
                        auc=numeric(), sens=numeric(), spec=numeric(), acc=numeric(),
                        n_feats=integer())
  list(
    metrics     = if (length(result_rows)>0) bind_rows(result_rows) else empty_m,
    feat_count  = feat_count,
    total_folds = n_repeats * n_folds,
    oof_preds   = oof_mean,
    oof_store   = oof_store,
    oof_labels  = as.integer(label_vec == "1")
  )
}

# ---- 3f. Summary helpers ----
feat_count_to_df <- function(feat_count, total_folds, extra_cols=list()) {
  df <- if (length(feat_count)==0) {
    data.frame(Feature=character(), Count=integer(), Frequency=numeric(),
               stringsAsFactors=FALSE)
  } else {
    data.frame(Feature=names(feat_count), Count=as.integer(feat_count),
               Frequency=round(feat_count/total_folds, 3), stringsAsFactors=FALSE)
  }
  for (nm in names(extra_cols)) df[[nm]] <- extra_cols[[nm]]
  df
}

summarise_cv <- function(metrics_df) {
  metrics_df %>%
    group_by(model) %>%
    summarise(AUC_mean  = round(mean(auc,  na.rm=TRUE), 3),
              AUC_sd    = round(sd(auc,    na.rm=TRUE), 3),
              Sens_mean = round(mean(sens, na.rm=TRUE), 3),
              Sens_sd   = round(sd(sens,   na.rm=TRUE), 3),
              Spec_mean = round(mean(spec, na.rm=TRUE), 3),
              Spec_sd   = round(sd(spec,   na.rm=TRUE), 3),
              Acc_mean  = round(mean(acc,  na.rm=TRUE), 3),
              Acc_sd    = round(sd(acc,    na.rm=TRUE), 3),
              n_folds   = n(), .groups="drop")
}

select_best_model <- function(metrics_df) {
  metrics_df %>%
    group_by(model) %>%
    summarise(AUC_mean=mean(auc, na.rm=TRUE), .groups="drop") %>%
    arrange(desc(AUC_mean)) %>% slice(1) %>% pull(model)
}

# ---- 3g. Repeat-level DeLong test ----
compute_delong_repeated <- function(oof_store1, oof_store2, labels_num,
                                    name1="Model1", name2="Model2") {
  n_rep    <- ncol(oof_store1)
  per_rep  <- lapply(seq_len(n_rep), function(r) {
    p1 <- oof_store1[,r]; p2 <- oof_store2[,r]
    valid <- !is.na(p1) & !is.na(p2) & !is.na(labels_num)
    if (sum(valid)<10 || length(unique(labels_num[valid]))<2)
      return(list(p=NA_real_, auc1=NA_real_, auc2=NA_real_))
    tryCatch({
      roc1 <- pROC::roc(labels_num[valid], p1[valid], quiet=TRUE)
      roc2 <- pROC::roc(labels_num[valid], p2[valid], quiet=TRUE)
      test <- pROC::roc.test(roc1, roc2, method="delong", paired=TRUE)
      list(p=test$p.value,
           auc1=as.numeric(pROC::auc(roc1)),
           auc2=as.numeric(pROC::auc(roc2)))
    }, error=function(e) list(p=NA_real_, auc1=NA_real_, auc2=NA_real_))
  })
  p_vals <- sapply(per_rep, `[[`, "p")
  auc1s  <- sapply(per_rep, `[[`, "auc1")
  auc2s  <- sapply(per_rep, `[[`, "auc2")
  data.frame(
    Comparison   = paste(name1, "vs", name2),
    AUC_1_mean   = round(mean(auc1s, na.rm=TRUE), 3),
    AUC_2_mean   = round(mean(auc2s, na.rm=TRUE), 3),
    p_median     = round(median(p_vals, na.rm=TRUE), 4),
    n_sig_0.05   = sum(p_vals < 0.05, na.rm=TRUE),
    n_valid_reps = sum(!is.na(p_vals)),
    stringsAsFactors=FALSE)
}

# ---- 3h. Calibration (Brier, slope, intercept, binned plot data) ----
compute_calibration <- function(pred, labels_num, n_bins=10, label="Model") {
  valid <- !is.na(pred) & !is.na(labels_num)
  p <- pred[valid]; y <- labels_num[valid]
  
  brier <- round(mean((p - y)^2), 4)
  
  logit_p <- pmax(-10, pmin(10, log(p / (1 - p + 1e-8))))
  cal_fit  <- tryCatch(
    glm(y ~ logit_p, family=binomial, control=list(maxit=100)),
    error=function(e) NULL)
  intercept <- if (!is.null(cal_fit)) round(coef(cal_fit)[1], 3) else NA_real_
  slope     <- if (!is.null(cal_fit)) round(coef(cal_fit)[2], 3) else NA_real_
  
  bin_idx <- cut(p, breaks=seq(0,1,length.out=n_bins+1),
                 include.lowest=TRUE, labels=FALSE)
  cal_df <- data.frame(pred=p, obs=y, bin=bin_idx) %>%
    group_by(bin) %>%
    summarise(mean_pred=round(mean(pred),4), mean_obs=round(mean(obs),4),
              n=n(), .groups="drop") %>%
    mutate(model=label, Brier=brier, Cal_intercept=intercept, Cal_slope=slope)
  
  list(brier=brier, intercept=intercept, slope=slope, cal_df=cal_df)
}

# ---- 3i. Batch calibration for all models within a task ----
calibrate_all_models <- function(res_obj, labels_num, task_label) {
  models <- c("lasso","rf","xgb")
  bind_rows(lapply(models, function(m) {
    pred <- res_obj$oof_preds[[m]]
    if (all(is.na(pred))) return(NULL)
    cal <- compute_calibration(pred, labels_num,
                               label=paste0(task_label,"_",m))
    data.frame(task=task_label, model=m,
               Brier=cal$brier, Cal_slope=cal$slope,
               Cal_intercept=cal$intercept,
               stringsAsFactors=FALSE)
  }))
}

# ============================================================
# ---- 4. Analysis 1: three-class one-vs-all (primary) ----
# ============================================================
cat("\n========== Analysis 1: three-class OVA (primary) ==========\n")
ova_results <- list(); ova_feat_freq <- list(); ova_oof <- list()

for (pos_class in 0:2) {
  label_name <- c("Control","GroupA","GroupB")[pos_class+1]
  cat(sprintf("  OVA: %s vs Rest...\n", label_name))
  label_ova <- ifelse(data$Group==pos_class, "1", "0")
  res <- run_cv(data, label_ova, radio_vars, N_REPEATS, N_FOLDS)
  ova_results[[label_name]]   <- res$metrics %>% mutate(OVA_class=label_name)
  ova_feat_freq[[label_name]] <- feat_count_to_df(res$feat_count, res$total_folds,
                                                  extra_cols=list(OVA_class=label_name))
  ova_oof[[label_name]]       <- list(preds=res$oof_preds, store=res$oof_store,
                                      labels=res$oof_labels)
}

ova_summary <- bind_rows(lapply(names(ova_results), function(cls)
  summarise_cv(ova_results[[cls]]) %>% mutate(OVA_class=cls)))
ova_feat_df <- bind_rows(ova_feat_freq) %>%
  filter(Count>0) %>% arrange(OVA_class, desc(Frequency))
ova_best_models <- setNames(
  lapply(names(ova_results), function(cls) select_best_model(ova_results[[cls]])),
  names(ova_results))

cat("\n--- OVA best models ---\n")
print(data.frame(OVA_class=names(ova_best_models), Best_model=unlist(ova_best_models)))
cat("\n--- OVA performance summary ---\n"); print(ova_summary, row.names=FALSE)

# OVA calibration: all models
ova_cal_summary <- bind_rows(lapply(names(ova_oof), function(cls) {
  res_tmp <- list(oof_preds=ova_oof[[cls]]$preds)
  calibrate_all_models(res_tmp, ova_oof[[cls]]$labels,
                       task_label=paste0("OVA_",cls))
})) %>% mutate(best=(model == unlist(ova_best_models)[match(
  sub("OVA_","",task), names(ova_best_models))]))

cat("\n--- OVA calibration (all models) ---\n")
print(ova_cal_summary, row.names=FALSE)

# OVA calibration bins (best model only, for plotting)
ova_cal_bins <- bind_rows(lapply(names(ova_oof), function(cls) {
  bm  <- ova_best_models[[cls]]
  cal <- compute_calibration(ova_oof[[cls]]$preds[[bm]],
                             ova_oof[[cls]]$labels,
                             label=paste0(cls,"_",bm))
  cal$cal_df %>% mutate(OVA_class=cls, best_model=bm)
}))

# Save OVA OOF predictions
ova_oof_df <- do.call(cbind, lapply(names(ova_oof), function(cls) {
  setNames(data.frame(ova_oof[[cls]]$preds[[ova_best_models[[cls]]]]),
           paste0("prob_",cls))
}))
ova_oof_df <- cbind(
  data.frame(No=data$No, Group=data$Group,
             true_group=c("Control","GroupA","GroupB")[data$Group+1]),
  ova_oof_df)
write.csv(ova_oof_df, file.path(out_dir,"3class_OVA_OOF_predictions.csv"), row.names=FALSE)

# ============================================================
# ---- 5. Analysis 2: binary GroupA vs GroupB (primary) ----
# ============================================================
cat("\n========== Analysis 2: binary GroupA vs GroupB (primary) ==========\n")
df_bin    <- data %>% filter(Group %in% c(1,2))
label_bin <- ifelse(df_bin$Group==2, "1", "0")
lbl_bin_num <- as.integer(label_bin=="1")
cat(sprintf("  GroupA: %d, GroupB: %d\n",
            sum(df_bin$Group==1), sum(df_bin$Group==2)))

cat("  Model 1: Radiomics...\n")
res_radio <- run_cv(df_bin, label_bin, radio_vars, N_REPEATS, N_FOLDS)
cat("  Model 2: Clinical...\n")
res_clin  <- run_cv(df_bin, label_bin, clinical_vars, N_REPEATS, N_FOLDS)
cat("  Model 3: Combined...\n")
res_comb  <- run_cv(df_bin, label_bin, c(radio_vars,clinical_vars), N_REPEATS, N_FOLDS)

bin_metrics <- bind_rows(
  res_radio$metrics %>% mutate(feature_set="Radiomics"),
  res_clin$metrics  %>% mutate(feature_set="Clinical"),
  res_comb$metrics  %>% mutate(feature_set="Combined"))
bin_summary <- bin_metrics %>%
  group_by(feature_set, model) %>%
  summarise(AUC_mean  = round(mean(auc,  na.rm=TRUE), 3),
            AUC_sd    = round(sd(auc,    na.rm=TRUE), 3),
            Sens_mean = round(mean(sens, na.rm=TRUE), 3),
            Sens_sd   = round(sd(sens,   na.rm=TRUE), 3),
            Spec_mean = round(mean(spec, na.rm=TRUE), 3),
            Spec_sd   = round(sd(spec,   na.rm=TRUE), 3),
            Acc_mean  = round(mean(acc,  na.rm=TRUE), 3),
            Acc_sd    = round(sd(acc,    na.rm=TRUE), 3),
            n_folds   = n(), .groups="drop")

best_radio <- select_best_model(res_radio$metrics)
best_clin  <- select_best_model(res_clin$metrics)
best_comb  <- select_best_model(res_comb$metrics)
cat("\n--- Binary performance summary ---\n"); print(bin_summary, row.names=FALSE)
cat(sprintf("Best models -> Radiomics:%s | Clinical:%s | Combined:%s\n",
            best_radio, best_clin, best_comb))

# DeLong (repeat-level)
cat("\n--- Binary DeLong tests (repeat-level, median p) ---\n")
delong_results <- bind_rows(list(
  compute_delong_repeated(
    res_radio$oof_store[[best_radio]], res_clin$oof_store[[best_clin]],
    lbl_bin_num, paste0("Radiomics-",best_radio), paste0("Clinical-",best_clin)),
  compute_delong_repeated(
    res_comb$oof_store[[best_comb]], res_radio$oof_store[[best_radio]],
    lbl_bin_num, paste0("Combined-",best_comb), paste0("Radiomics-",best_radio)),
  compute_delong_repeated(
    res_comb$oof_store[[best_comb]], res_clin$oof_store[[best_clin]],
    lbl_bin_num, paste0("Combined-",best_comb), paste0("Clinical-",best_clin))
))
print(delong_results, row.names=FALSE)

# Binary calibration: all 9 models
bin_cal_summary <- bind_rows(
  calibrate_all_models(res_radio, lbl_bin_num, "Radiomics"),
  calibrate_all_models(res_clin,  lbl_bin_num, "Clinical"),
  calibrate_all_models(res_comb,  lbl_bin_num, "Combined")
)
cat("\n--- Binary calibration (all models) ---\n")
print(bin_cal_summary, row.names=FALSE)

# Binary calibration bins (best model only, for plotting)
bin_cal_bins <- bind_rows(
  compute_calibration(res_radio$oof_preds[[best_radio]], lbl_bin_num,
                      label=paste0("Radiomics_",best_radio))$cal_df,
  compute_calibration(res_clin$oof_preds[[best_clin]],   lbl_bin_num,
                      label=paste0("Clinical_",best_clin))$cal_df,
  compute_calibration(res_comb$oof_preds[[best_comb]],   lbl_bin_num,
                      label=paste0("Combined_",best_comb))$cal_df
)

bin_oof_df <- data.frame(
  No             = df_bin$No,
  Group          = df_bin$Group,
  true_label     = lbl_bin_num,
  pred_Radiomics = res_radio$oof_preds[[best_radio]],
  pred_Clinical  = res_clin$oof_preds[[best_clin]],
  pred_Combined  = res_comb$oof_preds[[best_comb]])

bin_feat_freq <- feat_count_to_df(res_radio$feat_count, res_radio$total_folds) %>%
  filter(Count>0) %>% arrange(desc(Frequency))

bin_comb_feat_freq <- feat_count_to_df(res_comb$feat_count, res_comb$total_folds) %>%
  filter(Feature %in% radio_vars, Count > 0) %>%
  arrange(desc(Frequency))

# ============================================================
# ---- 6. Analysis 3: sensitivity analysis (GroupB restricted to LHO) ----
# ============================================================
cat("\n========== Analysis 3: sensitivity analysis (LHO subgroup) ==========\n")
data_sa <- data %>%
  filter(Group %in% c(0,1) |
           (Group==2 &
              tolower(trimws(as.character(Disease))) %in% tolower(LHO_types)))
cat(sprintf("  SA sample: Control=%d, GroupA=%d, GroupB(LHO)=%d (total %d)\n",
            sum(data_sa$Group==0),
            sum(data_sa$Group==1),
            sum(data_sa$Group==2), nrow(data_sa)))

# SA OVA
cat("\n--- SA three-class OVA ---\n")
sa_ova_results <- list(); sa_ova_oof <- list()
for (pos_class in 0:2) {
  label_name   <- c("Control","GroupA","GroupB")[pos_class+1]
  label_ova_sa <- ifelse(data_sa$Group==pos_class, "1", "0")
  cat(sprintf("  SA OVA: %s vs Rest...\n", label_name))
  res <- run_cv(data_sa, label_ova_sa, radio_vars, N_REPEATS, N_FOLDS)
  sa_ova_results[[label_name]] <- res$metrics %>% mutate(OVA_class=label_name)
  sa_ova_oof[[label_name]]     <- list(preds=res$oof_preds, store=res$oof_store,
                                       labels=res$oof_labels)
}
sa_ova_summary <- bind_rows(lapply(names(sa_ova_results), function(cls)
  summarise_cv(sa_ova_results[[cls]]) %>% mutate(OVA_class=cls)))
sa_ova_best <- setNames(
  lapply(names(sa_ova_results), function(cls) select_best_model(sa_ova_results[[cls]])),
  names(sa_ova_results))
cat("\n--- SA OVA performance summary ---\n"); print(sa_ova_summary, row.names=FALSE)

# SA OVA calibration: all models
sa_ova_cal_summary <- bind_rows(lapply(names(sa_ova_oof), function(cls) {
  res_tmp <- list(oof_preds=sa_ova_oof[[cls]]$preds)
  calibrate_all_models(res_tmp, sa_ova_oof[[cls]]$labels,
                       task_label=paste0("SA_OVA_",cls))
}))

# SA OVA calibration bins (best model only)
sa_ova_cal_bins <- bind_rows(lapply(names(sa_ova_oof), function(cls) {
  bm  <- sa_ova_best[[cls]]
  cal <- compute_calibration(sa_ova_oof[[cls]]$preds[[bm]],
                             sa_ova_oof[[cls]]$labels,
                             label=paste0("SA_",cls,"_",bm))
  cal$cal_df %>% mutate(OVA_class=cls, best_model=bm)
}))

# Save SA OVA OOF predictions
sa_ova_oof_df <- do.call(cbind, lapply(names(sa_ova_oof), function(cls) {
  setNames(data.frame(sa_ova_oof[[cls]]$preds[[sa_ova_best[[cls]]]]),
           paste0("prob_",cls))
}))
sa_ova_oof_df <- cbind(
  data.frame(No=data_sa$No, Group=data_sa$Group,
             true_group=c("Control","GroupA","GroupB")[data_sa$Group+1]),
  sa_ova_oof_df)
write.csv(sa_ova_oof_df, file.path(out_dir,"SA_3class_OVA_OOF_predictions.csv"),
          row.names=FALSE)

# SA binary
cat("\n--- SA binary GroupA vs GroupB ---\n")
df_sa_bin   <- data_sa %>% filter(Group %in% c(1,2))
label_sa    <- ifelse(df_sa_bin$Group==2, "1", "0")
lbl_sa_num  <- as.integer(label_sa=="1")
cat(sprintf("  GroupA: %d, GroupB(LHO): %d\n",
            sum(df_sa_bin$Group==1), sum(df_sa_bin$Group==2)))

cat("  SA Model 1: Radiomics...\n")
res_sa_r    <- run_cv(df_sa_bin, label_sa, radio_vars, N_REPEATS, N_FOLDS)
cat("  SA Model 2: Clinical...\n")
res_sa_c    <- run_cv(df_sa_bin, label_sa, clinical_vars, N_REPEATS, N_FOLDS)
cat("  SA Model 3: Combined...\n")
res_sa_comb <- run_cv(df_sa_bin, label_sa, c(radio_vars,clinical_vars), N_REPEATS, N_FOLDS)

sa_bin_metrics <- bind_rows(
  res_sa_r$metrics    %>% mutate(feature_set="Radiomics"),
  res_sa_c$metrics    %>% mutate(feature_set="Clinical"),
  res_sa_comb$metrics %>% mutate(feature_set="Combined"))
sa_bin_summary <- sa_bin_metrics %>%
  group_by(feature_set, model) %>%
  summarise(AUC_mean  = round(mean(auc,  na.rm=TRUE), 3),
            AUC_sd    = round(sd(auc,    na.rm=TRUE), 3),
            Sens_mean = round(mean(sens, na.rm=TRUE), 3),
            Sens_sd   = round(sd(sens,   na.rm=TRUE), 3),
            Spec_mean = round(mean(spec, na.rm=TRUE), 3),
            Spec_sd   = round(sd(spec,   na.rm=TRUE), 3),
            Acc_mean  = round(mean(acc,  na.rm=TRUE), 3),
            Acc_sd    = round(sd(acc,    na.rm=TRUE), 3),
            n_folds   = n(), .groups="drop")

sa_best_radio <- select_best_model(res_sa_r$metrics)
sa_best_clin  <- select_best_model(res_sa_c$metrics)
sa_best_comb  <- select_best_model(res_sa_comb$metrics)
cat("\n--- SA binary performance summary ---\n"); print(sa_bin_summary, row.names=FALSE)
cat(sprintf("SA best models -> Radiomics:%s | Clinical:%s | Combined:%s\n",
            sa_best_radio, sa_best_clin, sa_best_comb))

# SA DeLong
cat("\n--- SA DeLong tests (repeat-level) ---\n")
sa_delong <- bind_rows(list(
  compute_delong_repeated(
    res_sa_r$oof_store[[sa_best_radio]], res_sa_c$oof_store[[sa_best_clin]],
    lbl_sa_num, paste0("SA_Radiomics-",sa_best_radio), paste0("SA_Clinical-",sa_best_clin)),
  compute_delong_repeated(
    res_sa_comb$oof_store[[sa_best_comb]], res_sa_r$oof_store[[sa_best_radio]],
    lbl_sa_num, paste0("SA_Combined-",sa_best_comb), paste0("SA_Radiomics-",sa_best_radio)),
  compute_delong_repeated(
    res_sa_comb$oof_store[[sa_best_comb]], res_sa_c$oof_store[[sa_best_clin]],
    lbl_sa_num, paste0("SA_Combined-",sa_best_comb), paste0("SA_Clinical-",sa_best_clin))
))
print(sa_delong, row.names=FALSE)

# SA binary calibration: all 9 models
sa_bin_cal_summary <- bind_rows(
  calibrate_all_models(res_sa_r,    lbl_sa_num, "SA_Radiomics"),
  calibrate_all_models(res_sa_c,    lbl_sa_num, "SA_Clinical"),
  calibrate_all_models(res_sa_comb, lbl_sa_num, "SA_Combined")
)

# SA binary calibration bins (best model only)
sa_bin_cal_bins <- bind_rows(
  compute_calibration(res_sa_r$oof_preds[[sa_best_radio]], lbl_sa_num,
                      label=paste0("SA_Radiomics_",sa_best_radio))$cal_df,
  compute_calibration(res_sa_c$oof_preds[[sa_best_clin]],  lbl_sa_num,
                      label=paste0("SA_Clinical_",sa_best_clin))$cal_df,
  compute_calibration(res_sa_comb$oof_preds[[sa_best_comb]], lbl_sa_num,
                      label=paste0("SA_Combined_",sa_best_comb))$cal_df
)

sa_oof_df <- data.frame(
  No             = df_sa_bin$No,
  Group          = df_sa_bin$Group,
  true_label     = lbl_sa_num,
  pred_Radiomics = res_sa_r$oof_preds[[sa_best_radio]],
  pred_Clinical  = res_sa_c$oof_preds[[sa_best_clin]],
  pred_Combined  = res_sa_comb$oof_preds[[sa_best_comb]])

# ============================================================
# ---- 7. Primary vs SA comparison ----
# ============================================================
cat("\n========== Primary vs SA comparison (best single models) ==========\n")
compare_ova <- left_join(
  ova_summary %>% filter(model %in% unlist(ova_best_models)) %>%
    group_by(OVA_class) %>% slice_max(AUC_mean, n=1) %>%
    select(OVA_class, Best_model=model, Main_AUC=AUC_mean, Main_SD=AUC_sd),
  sa_ova_summary %>% filter(model %in% unlist(sa_ova_best)) %>%
    group_by(OVA_class) %>% slice_max(AUC_mean, n=1) %>%
    select(OVA_class, SA_Best_model=model, SA_AUC=AUC_mean, SA_SD=AUC_sd),
  by="OVA_class") %>% mutate(Delta=round(SA_AUC-Main_AUC, 3))
cat("\n--- Three-class OVA comparison ---\n"); print(compare_ova, row.names=FALSE)

make_bin_best_row <- function(summary_df, best_model, prefix="") {
  summary_df %>% filter(model==best_model) %>%
    select(feature_set,
           !!paste0(prefix,"Best_model"):=model,
           !!paste0(prefix,"AUC"):=AUC_mean,
           !!paste0(prefix,"SD"):=AUC_sd)
}
compare_bin <- left_join(
  bind_rows(make_bin_best_row(bin_summary%>%filter(feature_set=="Radiomics"),best_radio,"Main_"),
            make_bin_best_row(bin_summary%>%filter(feature_set=="Clinical"), best_clin, "Main_"),
            make_bin_best_row(bin_summary%>%filter(feature_set=="Combined"), best_comb, "Main_")),
  bind_rows(make_bin_best_row(sa_bin_summary%>%filter(feature_set=="Radiomics"),sa_best_radio,"SA_"),
            make_bin_best_row(sa_bin_summary%>%filter(feature_set=="Clinical"), sa_best_clin, "SA_"),
            make_bin_best_row(sa_bin_summary%>%filter(feature_set=="Combined"), sa_best_comb, "SA_")),
  by="feature_set") %>% mutate(Delta=round(SA_AUC-Main_AUC, 3))
cat("\n--- Binary GroupA vs GroupB comparison ---\n"); print(compare_bin, row.names=FALSE)

# ============================================================
# ---- 8. GA group difference test ----
# ============================================================
cat("\n========== GA group difference ==========\n")
print(data %>% group_by(Group) %>%
        summarise(n=n(), GA_mean=round(mean(GA,na.rm=T),1),
                  GA_sd=round(sd(GA,na.rm=T),1), GA_min=min(GA), GA_max=max(GA)))
cat("\nKruskal-Wallis test:\n"); print(kruskal.test(GA~Group, data=data))
cat("\nPairwise comparisons (Bonferroni):\n")
print(pairwise.wilcox.test(data$GA, data$Group, p.adjust.method="bonferroni"))

# ============================================================
# ---- 9. Write output files ----
# ============================================================
cat("\n========== Writing output files ==========\n")

# Three-class primary
write.csv(ova_summary,
          file.path(out_dir,"3class_OVA_summary.csv"),             row.names=FALSE)
write.csv(bind_rows(ova_results),
          file.path(out_dir,"3class_OVA_raw_folds.csv"),           row.names=FALSE)
write.csv(ova_feat_df,
          file.path(out_dir,"3class_OVA_feature_frequency.csv"),   row.names=FALSE)
write.csv(ova_cal_summary,
          file.path(out_dir,"3class_OVA_calibration_summary.csv"), row.names=FALSE)
write.csv(ova_cal_bins,
          file.path(out_dir,"3class_OVA_calibration_bins.csv"),    row.names=FALSE)

# Binary primary
write.csv(bin_summary,
          file.path(out_dir,"binary_AvsB_summary.csv"),             row.names=FALSE)
write.csv(bin_metrics,
          file.path(out_dir,"binary_AvsB_raw_folds.csv"),           row.names=FALSE)
write.csv(bin_feat_freq,
          file.path(out_dir,"binary_AvsB_feature_frequency.csv"),   row.names=FALSE)
write.csv(bin_comb_feat_freq,
          file.path(out_dir,"binary_AvsB_combined_radiomic_feature_frequency.csv"),
          row.names=FALSE)
write.csv(delong_results,
          file.path(out_dir,"binary_AvsB_DeLong.csv"),              row.names=FALSE)
write.csv(bin_cal_summary,
          file.path(out_dir,"binary_AvsB_calibration_summary.csv"), row.names=FALSE)
write.csv(bin_cal_bins,
          file.path(out_dir,"binary_AvsB_calibration_bins.csv"),    row.names=FALSE)
write.csv(bin_oof_df,
          file.path(out_dir,"binary_AvsB_OOF_predictions.csv"),     row.names=FALSE)

# Sensitivity analysis (LHO)
write.csv(sa_ova_summary,
          file.path(out_dir,"SA_3class_OVA_summary.csv"),             row.names=FALSE)
write.csv(sa_ova_cal_summary,
          file.path(out_dir,"SA_3class_OVA_calibration_summary.csv"), row.names=FALSE)
write.csv(sa_ova_cal_bins,
          file.path(out_dir,"SA_3class_OVA_calibration_bins.csv"),    row.names=FALSE)
write.csv(sa_bin_summary,
          file.path(out_dir,"SA_binary_AvsB_summary.csv"),            row.names=FALSE)
write.csv(sa_delong,
          file.path(out_dir,"SA_binary_AvsB_DeLong.csv"),             row.names=FALSE)
write.csv(sa_bin_cal_summary,
          file.path(out_dir,"SA_binary_AvsB_calibration_summary.csv"),row.names=FALSE)
write.csv(sa_bin_cal_bins,
          file.path(out_dir,"SA_binary_AvsB_calibration_bins.csv"),   row.names=FALSE)
write.csv(sa_oof_df,
          file.path(out_dir,"SA_binary_AvsB_OOF_predictions.csv"),    row.names=FALSE)

# Comparison tables
write.csv(compare_ova,
          file.path(out_dir,"comparison_OVA_main_vs_SA.csv"),       row.names=FALSE)
write.csv(compare_bin,
          file.path(out_dir,"comparison_binary_main_vs_SA.csv"),    row.names=FALSE)

cat("Primary outputs and sensitivity analysis files written to:", out_dir, "\n")

# ============================================================
# ---- 10. Analysis 4: postnatal-echo sensitivity analysis ----
# ============================================================
cat("\n========== Analysis 4: postnatal-echo sensitivity analysis ==========\n")

df_echo <- data %>%
  filter(Group %in% c(1,2),
         PostnatalEcho == 1)

cat(sprintf("  Echo SA sample: GroupA=%d, GroupB=%d (total %d)\n",
            sum(df_echo$Group==1),
            sum(df_echo$Group==2), nrow(df_echo)))

label_echo   <- ifelse(df_echo$Group==2, "1", "0")
lbl_echo_num <- as.integer(label_echo=="1")

cat("  Echo SA Model 1: Radiomics...\n")
res_echo_r    <- run_cv(df_echo, label_echo, radio_vars,   N_REPEATS, N_FOLDS)
cat("  Echo SA Model 2: Clinical...\n")
res_echo_c    <- run_cv(df_echo, label_echo, clinical_vars, N_REPEATS, N_FOLDS)
cat("  Echo SA Model 3: Combined...\n")
res_echo_comb <- run_cv(df_echo, label_echo,
                        c(radio_vars,clinical_vars), N_REPEATS, N_FOLDS)

echo_bin_metrics <- bind_rows(
  res_echo_r$metrics    %>% mutate(feature_set="Radiomics"),
  res_echo_c$metrics    %>% mutate(feature_set="Clinical"),
  res_echo_comb$metrics %>% mutate(feature_set="Combined"))

echo_bin_summary <- echo_bin_metrics %>%
  group_by(feature_set, model) %>%
  summarise(AUC_mean  = round(mean(auc,  na.rm=TRUE), 3),
            AUC_sd    = round(sd(auc,    na.rm=TRUE), 3),
            Sens_mean = round(mean(sens, na.rm=TRUE), 3),
            Sens_sd   = round(sd(sens,   na.rm=TRUE), 3),
            Spec_mean = round(mean(spec, na.rm=TRUE), 3),
            Spec_sd   = round(sd(spec,   na.rm=TRUE), 3),
            Acc_mean  = round(mean(acc,  na.rm=TRUE), 3),
            Acc_sd    = round(sd(acc,    na.rm=TRUE), 3),
            n_folds   = n(), .groups="drop")

echo_best_radio <- select_best_model(res_echo_r$metrics)
echo_best_clin  <- select_best_model(res_echo_c$metrics)
echo_best_comb  <- select_best_model(res_echo_comb$metrics)

cat("\n--- Echo SA performance summary ---\n")
print(echo_bin_summary, row.names=FALSE)

cat("\n--- Primary vs Echo SA comparison (best single models) ---\n")

compare_echo <- left_join(
  bind_rows(
    make_bin_best_row(bin_summary %>% filter(feature_set=="Radiomics"),
                      best_radio, "Main_"),
    make_bin_best_row(bin_summary %>% filter(feature_set=="Clinical"),
                      best_clin,  "Main_"),
    make_bin_best_row(bin_summary %>% filter(feature_set=="Combined"),
                      best_comb,  "Main_")),
  bind_rows(
    make_bin_best_row(echo_bin_summary %>% filter(feature_set=="Radiomics"),
                      echo_best_radio, "Echo_"),
    make_bin_best_row(echo_bin_summary %>% filter(feature_set=="Clinical"),
                      echo_best_clin,  "Echo_"),
    make_bin_best_row(echo_bin_summary %>% filter(feature_set=="Combined"),
                      echo_best_comb,  "Echo_")),
  by="feature_set") %>%
  mutate(Delta=round(Echo_AUC - Main_AUC, 3))

print(compare_echo, row.names=FALSE)

write.csv(echo_bin_summary,
          file.path(out_dir,"echo_SA_binary_AvsB_summary.csv"),  row.names=FALSE)
write.csv(compare_echo,
          file.path(out_dir,"echo_SA_vs_main_comparison.csv"),   row.names=FALSE)

# ============================================================
# ---- 11. Per-repeat OOF predictions (for ROC CI band plotting) ----
# ============================================================

# OVA
for (cls in names(ova_oof)) {
  mat <- ova_oof[[cls]]$store[[ova_best_models[[cls]]]]
  df  <- as.data.frame(mat)
  colnames(df) <- paste0("rep", 1:ncol(df))
  df$true_label <- ova_oof[[cls]]$labels
  write.csv(df, file.path(out_dir, sprintf("OVA_%s_repeat_preds.csv", cls)),
            row.names=FALSE)
}

# Binary
for (fs in c("radio","clin","comb")) {
  res_obj <- list(radio=res_radio, clin=res_clin, comb=res_comb)[[fs]]
  bm      <- list(radio=best_radio, clin=best_clin, comb=best_comb)[[fs]]
  fs_name <- c(radio="Radiomics", clin="Clinical", comb="Combined")[[fs]]
  mat     <- res_obj$oof_store[[bm]]
  df      <- as.data.frame(mat)
  colnames(df) <- paste0("rep", 1:ncol(df))
  df$true_label <- lbl_bin_num
  write.csv(df, file.path(out_dir, sprintf("Binary_%s_repeat_preds.csv", fs_name)),
            row.names=FALSE)
}

cat("\n=== Pipeline finished ===\n")