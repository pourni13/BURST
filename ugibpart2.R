# UGIB THESIS -- PART 2: MODEL TRAINING

.v246_required_pkgs <- c(
  "dplyr", "tidyr", "ggplot2",
  "randomForest", "glmnet", "pROC", "MASS", "mice",
  "xgboost", "fastshap", "shapviz", "Hmisc",
  "writexl"
)
.v246_missing <- .v246_required_pkgs[
  !vapply(.v246_required_pkgs,
          function(p) requireNamespace(p, quietly = TRUE),
          logical(1))]
if (length(.v246_missing) > 0L) {
  message(sprintf("V2.4.6: installing %d missing package(s): %s",
                  length(.v246_missing),
                  paste(.v246_missing, collapse = ", ")))
  install.packages(.v246_missing,
                   repos = "https://cloud.r-project.org",
                   dependencies = TRUE)
}
rm(.v246_required_pkgs, .v246_missing)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(randomForest); library(glmnet); library(pROC)
})
suppressWarnings({
  for (loc in c("C.UTF-8","en_US.UTF-8",
                "English_United States.utf8","English_United States.1252")) {
    if (nzchar(tryCatch(Sys.setlocale("LC_ALL", loc), error=function(e) ""))) break
  }
})
while (!is.null(dev.list())) dev.off()

PRIMARY    <- "high_risk"
SECONDARY  <- "major_tx"
TERTIARY   <- c("outcome_death","rebleeding")
ALL_OUTCOMES <- c(PRIMARY, SECONDARY, TERTIARY)
cat(sprintf("Primary outcome:    %s\n", PRIMARY))
cat(sprintf("Secondary outcome:  %s\n", SECONDARY))
cat(sprintf("Tertiary outcomes:  %s\n\n", paste(TERTIARY, collapse=", ")))

theme_ugib <- function(base=12)
  theme_minimal(base_size=base) +
  theme(plot.title    = element_text(face="bold", size=base+1, hjust=0),
        plot.subtitle = element_text(size=base-1, colour="grey40"),
        plot.caption  = element_text(size=base-3, colour="grey55"),
        axis.text     = element_text(size=base-1),
        panel.grid.minor = element_blank(), legend.position="bottom")

theme_pub <- function(base=12)
  theme_minimal(base_size=base) +
  theme(plot.title    = element_text(face="bold", size=base+1, hjust=0.5),
        plot.subtitle = element_text(size=base-1, colour="grey30", hjust=0.5),
        axis.title    = element_text(size=base),
        axis.text     = element_text(size=base-1, colour="grey20"),
        panel.grid.minor = element_blank(),
        panel.border  = element_rect(colour="grey80", fill=NA, linewidth=0.5),
        legend.position="right", legend.text=element_text(size=base-2))

COL_BLUE <- "#2E75B6"; COL_RED <- "#C0392B"; COL_GREEN <- "#27AE60"
COL_AMB  <- "#E67E22"; COL_GREY <- "#7F7F7F"; COL_PURP <- "#8E44AD"
pct_fmt <- function(x) paste0(round(x*100), "%")

plots <- list()

cat("\n============================================================\n")
cat("LOADING PART 1 OUTPUTS\n")
cat("============================================================\n\n")

df_full <- readRDS("ugib_df_full.rds")
df_ml   <- readRDS("ugib_df_ml.rds")
feats   <- readRDS("final_feature_set.rds")
roc_sc  <- readRDS("roc_scores.rds")

for (d in c("df_full","df_ml")) {
  df_tmp <- get(d)
  if (all(c("hb_ed","hb_entry") %in% names(df_tmp)) & !"hb_mean" %in% names(df_tmp))
    df_tmp$hb_mean  <- (df_tmp$hb_ed  + df_tmp$hb_entry)  / 2
  if (all(c("hct_ed","hct_entry") %in% names(df_tmp)) & !"hct_mean" %in% names(df_tmp))
    df_tmp$hct_mean <- (df_tmp$hct_ed + df_tmp$hct_entry) / 2
  assign(d, df_tmp)
}

for (d in c("df_full","df_ml")) {
  df_tmp <- get(d)
  if ("major_transfusion" %in% names(df_tmp) & !"major_tx" %in% names(df_tmp))
    df_tmp$major_tx <- df_tmp$major_transfusion
  if ("high_risk_outcome" %in% names(df_tmp) & !"high_risk" %in% names(df_tmp))
    df_tmp$high_risk <- df_tmp$high_risk_outcome
  assign(d, df_tmp)
}

n <- nrow(df_full)
cat(sprintf("Loaded: %d patients | %d confirmed features\n", n, length(feats)))
cat("Features: ", paste(feats, collapse=", "), "\n\n")

# SECTION 7: MODEL TRAINING

display_names <- c(
  age = "Age (years)",
  symptom_duration = "Symptom duration (d)",
  sbp = "Systolic BP",
  dbp = "Diastolic BP",
  hr = "Heart rate",
  altered_mental = "Altered mental status",
  syncope = "Syncope",
  hematemesis_red = "Hematemesis (bright red)",
  hematemesis_coffee = "Hematemesis (coffee-grounds)",
  melena = "Melena",
  hematochezia = "Hematochezia",
  rectal_exam = "Rectal exam (positive)",
  gi_bleed_recent = "Recent GI bleed (<=12m)",
  gi_bleed_remote = "Remote GI bleed (>12m)",
  gi_bleed_lt3m = "GI bleed <3m",
  gi_bleed_3to6m = "GI bleed 3-6m",
  gi_bleed_gt6m = "GI bleed 6-12m",
  gi_bleed_1to5y = "GI bleed 1-5y",
  gi_bleed_gt5y = "GI bleed >5y",
  hb_mean = "Hemoglobin",
  lactate_ed = "Lactate",
  wbc_entry = "White-cell count",
  plt_entry = "Platelet count",
  ri_pct = "Reticulocyte index",
  creatinine = "Creatinine",
  urea = "Urea",
  bun = "BUN",
  albumin = "Albumin",
  inr = "INR",
  aptt = "aPTT",
  pt_pct = "Prothrombin time (%)",
  sgot = "AST",
  sgpt = "ALT",
  tbil = "Total bilirubin",
  ferritin = "Ferritin",
  ldh = "LDH",
  drug_act = "Anticoagulant",
  drug_antiplt = "Antiplatelet",
  drug_ppi = "PPI",
  drug_nsaid = "NSAID",
  drug_steroid = "Steroid",
  liver_cirrhosis = "Liver cirrhosis",
  renal_severity_3 = "Renal disease (3-level)",
  heart_failure = "Heart failure",
  cad = "Coronary artery disease",
  afib = "Atrial fibrillation",
  hypertension = "Hypertension",
  copd = "COPD",
  asthma = "Asthma",
  dm_any = "Diabetes mellitus",
  active_cancer = "Active cancer",
  comorbidity_hematol = "Haematologic comorbidity",
  ibd = "Inflammatory bowel disease",
  ed_levosim = "ED vasopressor",
  hb_ed = "Hemoglobin (ABG)",
  hb_entry = "Hemoglobin (analyser, entry)",
  hb_d1 = "Hemoglobin (day 1)",
  hct_ed = "Hematocrit (ABG)",
  hct_entry = "Hematocrit (analyser, entry)",
  hct_d1 = "Hematocrit (day 1)",
  hct_mean = "Hematocrit (mean)",
  total_wb = "Whole blood units (total)",
  rbcs = "pRBC units",
  total_plt = "Platelet units (total)",
  total_plasma = "Plasma units (total)",
  wb_ed = "Whole blood (ED)",
  wb_wkup = "Whole blood (work-up)",
  time_to_endo_h = "Time to endoscopy (h)",
  ed_fluids_l = "Fluid resuscitation (L)",
  hosp_days = "Hospital stay (days)"
)
relabel_var <- function(v) {
  out <- display_names[as.character(v)]
  out[is.na(out)] <- as.character(v)[is.na(out)]
  unname(out)
}
relabel_gender <- function(g) {
  g <- as.character(g)
  ifelse(g %in% c("M","male","Male","1"), "Sex (male)",
         ifelse(g %in% c("F","female","Female","0"), "Sex (female)", "Sex"))
}
LOW_SIGNAL_DROP <- c("asthma", "comorbidity_hematol", "ibd", "ed_levosim")
filter_low_signal <- function(df_in, var_col = "variable") {
  df_in[!df_in[[var_col]] %in% LOW_SIGNAL_DROP, , drop = FALSE]
}

cat("============================================================\n")
cat("SECTION 7: MODEL TRAINING\n")
cat("============================================================\n\n")

ed_feats <- feats[feats %in% names(df_ml)]

nzv <- sapply(ed_feats, function(v){x<-df_ml[[v]]
if(is.numeric(x)) return(var(x,na.rm=T)<1e-4)
tbl<-table(x,useNA="no"); length(tbl)<2 | max(tbl)/sum(tbl)>.97})
if(any(nzv)){
  cat("NZV removed:", paste(ed_feats[nzv],collapse=", "), "\n")
  ed_feats <- ed_feats[!nzv]
}

df_ed <- df_ml |>
  dplyr::select(all_of(c("patient_id",ed_feats)),
                transfusion_class, major_tx, high_risk) |>
  filter(!is.na(transfusion_class))

cat(sprintf("Training: n=%d | features=%d\n",nrow(df_ed),length(ed_feats)))
cat("Class distribution (transfusion_class):\n"); print(table(df_ed$transfusion_class))
cat(sprintf("PRIMARY (%s) events: %d / %d\n",
            PRIMARY, sum(df_ed[[PRIMARY]]=="YES"), nrow(df_ed)))
cat(sprintf("SECONDARY (%s) events: %d / %d\n\n",
            SECONDARY, sum(df_ed[[SECONDARY]]=="YES"), nrow(df_ed)))

cl_levels  <- levels(df_ed$transfusion_class)
n_min      <- min(table(df_ed$transfusion_class))
samp_sizes <- setNames(rep(n_min, 3), cl_levels)
class_wts  <- setNames(c(1, 3, 3), cl_levels)

cat(sprintf("Stratified sampsize per class: %d\n", n_min))
cat("Class weights:", paste(names(class_wts),class_wts,sep="=",collapse=" / "), "\n\n")

fit_lr <- function(df_in, ed_feats, outcome) {
  X <- as.data.frame(lapply(df_in[, ed_feats, drop=FALSE], function(v) {
    if (is.numeric(v)) v else as.numeric(v == "YES")
  }))
  y <- as.numeric(df_in[[outcome]] == "YES")
  meds <- vapply(X, function(c) median(c, na.rm=TRUE), numeric(1))
  for (j in seq_along(X)) X[[j]][is.na(X[[j]])] <- meds[[j]]
  mu <- vapply(X, mean, numeric(1)); sg <- vapply(X, sd, numeric(1))
  sg[is.na(sg) | sg == 0] <- 1
  Xs <- as.matrix(sweep(sweep(X, 2, mu, "-"), 2, sg, "/"))
  n <- length(y); n0 <- sum(y == 0); n1 <- sum(y == 1)
  w <- ifelse(y == 1, n / (2*max(n1,1)), n / (2*max(n0,1)))
  fit <- tryCatch(
    suppressWarnings(cv.glmnet(Xs, y, family="binomial", alpha=0,
                               weights=w, nfolds=5, type.measure="auc")),
    error = function(e) NULL)
  if (is.null(fit)) {
    fit <- glmnet(Xs, y, family="binomial", alpha=0, weights=w, lambda=0.01)
  }
  list(fit=fit, meds=meds, mu=mu, sg=sg, X=Xs, y=y,
       is_cv=inherits(fit, "cv.glmnet"))
}
predict_lr <- function(model, df_in, ed_feats) {
  Xn <- as.data.frame(lapply(df_in[, ed_feats, drop=FALSE], function(v) {
    if (is.numeric(v)) v else as.numeric(v == "YES")
  }))
  for (j in seq_along(Xn)) Xn[[j]][is.na(Xn[[j]])] <- model$meds[[j]]
  Xs <- as.matrix(sweep(sweep(Xn, 2, model$mu, "-"), 2, model$sg, "/"))
  if (model$is_cv) as.numeric(predict(model$fit, newx=Xs, s="lambda.min", type="response"))
  else as.numeric(predict(model$fit, newx=Xs, type="response"))
}

df_hr <- df_ml |>
  dplyr::select(all_of(ed_feats), high_risk) |>
  filter(!is.na(high_risk))

cat("==== PRIMARY OUTCOME: high_risk ====\n")

set.seed(2024)
hr_min <- min(table(df_hr$high_risk))
rf_hr <- randomForest(
  reformulate(ed_feats, "high_risk"),
  data=df_hr, ntree=500, importance=TRUE,
  strata=df_hr$high_risk,
  sampsize=c("NO"=hr_min, "YES"=hr_min))
oob_prob_hr <- rf_hr$votes[,"YES"] / rowSums(rf_hr$votes)
roc_ml_hr   <- roc(df_hr$high_risk, oob_prob_hr,
                   levels=c("NO","YES"), direction="<", quiet=TRUE)
cat(sprintf("7A. Binary RF OOB AUC (high_risk):    %.3f\n", auc(roc_ml_hr)))

set.seed(2024)
lr_hr <- fit_lr(df_hr, ed_feats, "high_risk")
prob_lr_hr <- predict_lr(lr_hr, df_hr, ed_feats)
roc_lr_hr  <- roc(df_hr$high_risk, prob_lr_hr,
                  levels=c("NO","YES"), direction="<", quiet=TRUE)
cat(sprintf("7B. Binary LR APPARENT AUC (high_risk): %.3f  (honest CV in Section 8)\n\n",
            auc(roc_lr_hr)))

cat("==== SECONDARY OUTCOME: major_tx ====\n")

set.seed(2024)
rf_mt <- randomForest(
  reformulate(ed_feats, "major_tx"),
  data=df_ed, ntree=500, importance=TRUE,
  strata=df_ed$major_tx,
  sampsize=c("NO"=sum(df_ed$major_tx=="YES"),
             "YES"=sum(df_ed$major_tx=="YES")))
oob_prob_mt <- rf_mt$votes[,"YES"] / rowSums(rf_mt$votes)
roc_ml_mt   <- roc(df_ed$major_tx, oob_prob_mt,
                   levels=c("NO","YES"), direction="<", quiet=TRUE)
cat(sprintf("7C. Binary RF OOB AUC (major_tx):    %.3f\n", auc(roc_ml_mt)))

set.seed(2024)
lr_mt <- fit_lr(df_ed, ed_feats, "major_tx")
prob_lr_mt <- predict_lr(lr_mt, df_ed, ed_feats)
roc_lr_mt  <- roc(df_ed$major_tx, prob_lr_mt,
                  levels=c("NO","YES"), direction="<", quiet=TRUE)
cat(sprintf("7D. Binary LR APPARENT AUC (major_tx): %.3f  (honest CV in Section 8)\n\n",
            auc(roc_lr_mt)))

cat("==== ORDINAL OUTCOME: transfusion_class (3 levels) ====\n")

set.seed(2024)
rf_tc <- randomForest(
  reformulate(ed_feats, "transfusion_class"),
  data=df_ed, ntree=500, importance=TRUE,
  mtry=floor(sqrt(length(ed_feats))),
  strata=df_ed$transfusion_class,
  sampsize=samp_sizes, classwt=class_wts)
oob_pred_tc <- rf_tc$predicted[!is.na(rf_tc$predicted)]
true_labs   <- df_ed$transfusion_class[!is.na(rf_tc$predicted)]
oob_acc_tc  <- round(mean(as.character(oob_pred_tc)==as.character(true_labs))*100,1)
cat(sprintf("7E. 3-class RF OOB accuracy: %.1f%%\n", oob_acc_tc))
cat("OOB confusion:\n"); print(rf_tc$confusion)
cat("\n")

# SECTION 8: 5x5 REPEATED STRATIFIED CROSS-VALIDATION
cat("============================================================\n")
cat("SECTION 8: 5x5 REPEATED STRATIFIED CV\n")
cat("============================================================\n\n")

N_REPEATS <- 5L
N_FOLDS   <- 5L

make_strat_folds <- function(y, k, seed) {
  set.seed(seed); folds <- integer(length(y))
  for (cls in unique(y)) {
    idx <- which(y == cls); idx <- sample(idx)
    folds[idx] <- cut(seq_along(idx), breaks=k, labels=FALSE)
  }
  folds
}

cv_binary <- function(df_in, ed_feats, outcome,
                      model_type = "RF",
                      n_repeats = N_REPEATS, n_folds = N_FOLDS,
                      seed = 2024) {
  y_full <- as.character(df_in[[outcome]])
  aucs   <- numeric(n_repeats)
  oof_all <- matrix(NA_real_, nrow=n_repeats, ncol=length(y_full))
  for (rep in seq_len(n_repeats)) {
    fld <- make_strat_folds(y_full, n_folds, seed + rep)
    oof <- numeric(length(y_full))
    for (k in seq_len(n_folds)) {
      tr <- which(fld != k); te <- which(fld == k)
      tr_df <- df_in[tr, , drop=FALSE]; te_df <- df_in[te, , drop=FALSE]
      if (length(unique(y_full[tr])) < 2) {
        oof[te] <- mean(y_full[tr] == "YES"); next
      }
      pred <- tryCatch({
        if (model_type == "RF") {
          n_min <- min(table(tr_df[[outcome]]))
          fit <- randomForest(reformulate(ed_feats, outcome),
                              data=tr_df, ntree=500,
                              strata=tr_df[[outcome]],
                              sampsize=c("NO"=n_min,"YES"=n_min))
          pp <- predict(fit, te_df, type="prob")
          pp[, "YES"]
        } else {
          m <- fit_lr(tr_df, ed_feats, outcome)
          predict_lr(m, te_df, ed_feats)
        }
      }, error = function(e) rep(mean(y_full[tr] == "YES"), length(te)))
      oof[te] <- pred
    }
    aucs[rep] <- tryCatch(
      as.numeric(auc(roc(y_full, oof, levels=c("NO","YES"),
                         direction="<", quiet=TRUE))),
      error = function(e) NA_real_)
    oof_all[rep, ] <- oof
  }
  list(aucs=aucs, oof=colMeans(oof_all), y=y_full)
}

cv_multiclass_acc <- function(df_in, ed_feats, outcome,
                              n_repeats = N_REPEATS, n_folds = N_FOLDS,
                              seed = 2024) {
  y_full <- as.character(df_in[[outcome]])
  acc_per_rep <- numeric(n_repeats)
  for (rep in seq_len(n_repeats)) {
    fld <- make_strat_folds(y_full, n_folds, seed + rep)
    correct <- 0L; total <- 0L
    for (k in seq_len(n_folds)) {
      tr <- which(fld != k); te <- which(fld == k)
      tr_df <- df_in[tr, , drop=FALSE]; te_df <- df_in[te, , drop=FALSE]
      if (length(unique(tr_df[[outcome]])) < 2) next
      n_min <- min(table(tr_df[[outcome]]))
      cl <- levels(tr_df[[outcome]])
      pred <- tryCatch({
        fit <- randomForest(reformulate(ed_feats, outcome),
                            data=tr_df, ntree=300,
                            mtry=floor(sqrt(length(ed_feats))),
                            strata=tr_df[[outcome]],
                            sampsize=setNames(rep(n_min, length(cl)), cl),
                            classwt=class_wts)
        as.character(predict(fit, te_df))
      }, error=function(e) rep(NA_character_, length(te)))
      ok <- !is.na(pred) & !is.na(y_full[te])
      correct <- correct + sum(pred[ok] == y_full[te][ok])
      total   <- total + sum(ok)
    }
    acc_per_rep[rep] <- if (total > 0) correct / total else NA_real_
  }
  acc_per_rep
}

cv_results <- list()

cat("--- 8A. transfusion_class (3-class RF) ---\n")
acc_tc <- cv_multiclass_acc(df_ed, ed_feats, "transfusion_class")
cv_results$transfusion_class <- list(model="RF", per_rep=acc_tc,
                                     mean=mean(acc_tc, na.rm=TRUE),
                                     sd=sd(acc_tc, na.rm=TRUE))
cat(sprintf("  5x5 CV accuracy: mean=%.1f%% sd=%.2f%%  per-repeat: %s\n\n",
            mean(acc_tc, na.rm=TRUE)*100, sd(acc_tc, na.rm=TRUE)*100,
            paste(sprintf("%.1f", acc_tc*100), collapse=", ")))

binary_outcomes <- list(major_tx = df_ed, high_risk = df_hr)
for (outcome in names(binary_outcomes)) {
  for (mtype in c("RF","LR")) {
    cat(sprintf("--- 8X. %s (%s) ---\n", outcome, mtype))
    cv <- cv_binary(binary_outcomes[[outcome]], ed_feats, outcome, model_type=mtype)
    key <- sprintf("%s_%s", outcome, mtype)
    cv_results[[key]] <- list(model=mtype, per_rep=cv$aucs,
                              mean=mean(cv$aucs, na.rm=TRUE),
                              sd=sd(cv$aucs, na.rm=TRUE),
                              oof=cv$oof, y=cv$y)
    cat(sprintf("  5x5 CV AUC: mean=%.3f sd=%.3f  per-repeat: %s\n\n",
                mean(cv$aucs, na.rm=TRUE), sd(cv$aucs, na.rm=TRUE),
                paste(sprintf("%.3f", cv$aucs), collapse=", ")))
  }
}

# SECTION 8F: Ridge-LR vs RF CV-comparable AUCs and paired DeLong
cat("\n--- Ridge-LR vs RF, 5x5 CV on PRIMARY ---\n")

.rf_aucs <- cv_results[["high_risk_RF"]]$per_rep
.lr_aucs <- cv_results[["high_risk_LR"]]$per_rep
.rf_oof  <- cv_results[["high_risk_RF"]]$oof
.lr_oof  <- cv_results[["high_risk_LR"]]$oof
.y_hr    <- cv_results[["high_risk_RF"]]$y

roc_rf_cv <- roc(.y_hr, .rf_oof, levels = c("NO","YES"), direction = "<", quiet = TRUE)
roc_lr_cv <- roc(.y_hr, .lr_oof, levels = c("NO","YES"), direction = "<", quiet = TRUE)
delong_rf_vs_lr <- tryCatch(
  roc.test(roc_rf_cv, roc_lr_cv, method = "delong", paired = TRUE),
  error = function(e) NULL)

cat(sprintf("  Ridge LR 5x5 CV AUC: mean=%.3f sd=%.3f range=%.3f-%.3f\n",
            mean(.lr_aucs, na.rm = TRUE), sd(.lr_aucs, na.rm = TRUE),
            min(.lr_aucs, na.rm = TRUE),  max(.lr_aucs, na.rm = TRUE)))
cat(sprintf("  RF       5x5 CV AUC: mean=%.3f sd=%.3f range=%.3f-%.3f\n",
            mean(.rf_aucs, na.rm = TRUE), sd(.rf_aucs, na.rm = TRUE),
            min(.rf_aucs, na.rm = TRUE),  max(.rf_aucs, na.rm = TRUE)))
cat(sprintf("  Held-out bagged AUC (RF): %.3f | (LR): %.3f\n",
            as.numeric(auc(roc_rf_cv)), as.numeric(auc(roc_lr_cv))))
if (!is.null(delong_rf_vs_lr)) {
  cat(sprintf("  Paired DeLong (RF vs LR, CV held-out): Z=%.3f p=%.4f\n",
              as.numeric(delong_rf_vs_lr$statistic),
              as.numeric(delong_rf_vs_lr$p.value)))
} else {
  cat("  Paired DeLong (RF vs LR, CV held-out): FAILED (see roc.test() error)\n")
}

saveRDS(list(
  rf_aucs = .rf_aucs, lr_aucs = .lr_aucs,
  rf_oof = .rf_oof, lr_oof = .lr_oof, y = .y_hr,
  roc_rf_cv = roc_rf_cv, roc_lr_cv = roc_lr_cv,
  delong_rf_vs_lr = delong_rf_vs_lr),
  "block_b_lr_vs_rf_cv.rds")
cat("Saved: block_b_lr_vs_rf_cv.rds\n")

cat(sprintf("\nRidge LR 5x5 CV AUC: mean=%.3f sd=%.3f\n",
            mean(.lr_aucs, na.rm = TRUE), sd(.lr_aucs, na.rm = TRUE)))
cat(sprintf("RF       5x5 CV AUC: mean=%.3f sd=%.3f\n",
            mean(.rf_aucs, na.rm = TRUE), sd(.rf_aucs, na.rm = TRUE)))
if (!is.null(delong_rf_vs_lr)) {
  cat(sprintf("Paired DeLong RF vs LR (CV held-out): p=%.4f\n",
              as.numeric(delong_rf_vs_lr$p.value)))
} else {
  cat("Paired DeLong RF vs LR (CV held-out): p=NA\n")
}

cv_long <- do.call(rbind, lapply(names(cv_results), function(k) {
  r <- cv_results[[k]]
  is_acc <- k == "transfusion_class"
  .cv_key_map <- c(
    "high_risk_RF"  = "High risk RF",
    "high_risk_LR"  = "High risk LR",
    "major_tx_RF"   = "Major transfusion RF",
    "major_tx_LR"   = "Major transfusion LR",
    "transfusion_class" = "Transfusion class")
  k_disp <- ifelse(k %in% names(.cv_key_map), .cv_key_map[k], k)
  data.frame(model_outcome = k_disp, repeat_id = seq_along(r$per_rep),
             metric = if (is_acc) r$per_rep else r$per_rep,
             metric_type = if (is_acc) "accuracy" else "AUC",
             stringsAsFactors = FALSE)
}))
plots[["p_cv_stability"]] <- ggplot(
  cv_long |> filter(metric_type == "AUC"),
  aes(x = model_outcome, y = metric, colour = model_outcome)) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.75) +
  stat_summary(fun = mean, geom = "point", shape = 95, size = 12, colour = "black") +
  scale_colour_manual(values = c("Major transfusion RF"=COL_BLUE, "Major transfusion LR"=COL_AMB,
                                 "High risk RF"=COL_GREEN, "High risk LR"=COL_RED)) +
  labs(title = "BURST CV stability: per-repeat AUCs",
       subtitle = sprintf("%d repeats x %d stratified folds; black bar = mean",
                          N_REPEATS, N_FOLDS),
       x = NULL, y = "AUC per repeat") +
  ylim(0, 1) + theme_ugib() + theme(legend.position = "none")
print(plots[["p_cv_stability"]])

saveRDS(cv_results, "cv_results.rds"); cat("\nSaved: cv_results.rds\n")

# SECTION 8G: NESTED (SELECTION-INCLUSIVE) CROSS-VALIDATION
RUN_NESTED_CV <- TRUE

if (RUN_NESTED_CV) {
  suppressPackageStartupMessages(library(rpart))
  cat("\n============================================================\n")
  cat("SECTION 8G: NESTED (SELECTION-INCLUSIVE) 5x5 CV\n")
  cat("============================================================\n")
  
  meta       <- readRDS("nested_cv_meta.rds")
  nest_cands <- meta$feat_cands
  nest_num   <- meta$vars_num
  nest_cat   <- meta$vars_cat
  PRIM       <- meta$PRIMARY
  
  B_INNER        <- 60L
  N_REPEATS_NEST <- 5L
  N_FOLDS_NEST   <- 5L
  
  nest_cands <- nest_cands[nest_cands %in% names(df_full)]
  nest_num   <- intersect(nest_num, nest_cands)
  nest_cat   <- intersect(nest_cat, nest_cands)
  df_nest    <- df_full[, c(nest_cands, PRIM), drop = FALSE]
  df_nest    <- df_nest[!is.na(df_nest[[PRIM]]), , drop = FALSE]
  df_nest[[PRIM]] <- factor(as.character(df_nest[[PRIM]]), levels = c("NO", "YES"))
  for (v in nest_cat) df_nest[[v]] <- factor(as.character(df_nest[[v]]))
  
  # fold-internal imputation: fit on TRAIN, apply to TRAIN and TEST
  .fit_imp <- function(dtr) {
    med <- setNames(lapply(nest_num, function(v) median(dtr[[v]], na.rm = TRUE)), nest_num)
    mod <- setNames(lapply(nest_cat, function(v) {
      tb <- sort(table(dtr[[v]]), decreasing = TRUE); if (length(tb)) names(tb)[1] else NA
    }), nest_cat)
    list(med = med, mod = mod)
  }
  .apply_imp <- function(d, im) {
    if ("lactate_ed" %in% names(d) && "lactate_missing" %in% names(d))
      d$lactate_missing <- as.integer(is.na(d$lactate_ed))
    for (v in nest_num) { x <- d[[v]]; x[is.na(x)] <- im$med[[v]]; d[[v]] <- x }
    for (v in nest_cat) {
      x <- as.character(d[[v]]); x[is.na(x)] <- im$mod[[v]]
      d[[v]] <- factor(x, levels = levels(df_nest[[v]]))
    }
    d
  }
  
  # balanced-RF OOB ROC on a feature set (used by the forward step)
  .rf_oob_roc <- function(dd, fs) {
    mn <- min(table(dd[[PRIM]]))
    set.seed(2024)
    rf <- randomForest(reformulate(fs, PRIM), data = dd, ntree = 500,
                       strata = dd[[PRIM]], sampsize = c("NO" = mn, "YES" = mn))
    p <- rf$votes[, "YES"] / rowSums(rf$votes)
    roc(dd[[PRIM]], p, levels = c("NO", "YES"), direction = "<", quiet = TRUE)
  }
  
  # Method 1: univariate nominal p<.05 (Wilcoxon / Fisher-chisq)
  .univ_nom <- function(dd) {
    out <- character(0)
    for (v in nest_num) {
      x <- dd[[v]]; mt <- dd[[PRIM]]; ok <- !is.na(x) & !is.na(mt)
      if (sum(ok) < 10 || length(unique(mt[ok])) < 2) next
      p <- tryCatch(suppressWarnings(
        wilcox.test(x[ok & mt == "YES"], x[ok & mt == "NO"])$p.value),
        error = function(e) NA)
      if (!is.na(p) && p < .05) out <- c(out, v)
    }
    for (v in nest_cat) {
      x <- dd[[v]]; mt <- dd[[PRIM]]; ok <- !is.na(x) & !is.na(mt)
      if (sum(ok) < 10) next
      tb <- table(x[ok], mt[ok]); if (nrow(tb) < 2 || ncol(tb) < 2) next
      p <- tryCatch(suppressWarnings(
        if (any(chisq.test(tb)$expected < 5))
          fisher.test(tb, simulate.p.value = TRUE, B = 2000)$p.value
        else chisq.test(tb)$p.value), error = function(e) NA)
      if (!is.na(p) && p < .05) out <- c(out, v)
    }
    out
  }
  
  # Method 3: RF top-15 Gini
  .rf_top15 <- function(dd, ntree = 500) {
    rf <- tryCatch(randomForest(reformulate(nest_cands, PRIM), data = dd,
                                ntree = ntree, mtry = floor(sqrt(length(nest_cands))),
                                importance = FALSE), error = function(e) NULL)
    if (is.null(rf)) return(character(0))
    g <- importance(rf)[, "MeanDecreaseGini"]
    names(sort(g, decreasing = TRUE))[1:min(15, length(g))]
  }
  
  # Method 4: LASSO lambda.1se variables
  .lasso_1se <- function(dd) {
    X <- tryCatch(model.matrix(reformulate(nest_cands), data = dd)[, -1, drop = FALSE],
                  error = function(e) NULL)
    y <- as.numeric(dd[[PRIM]] == "YES")
    if (is.null(X) || length(unique(y)) < 2 || nrow(X) < 20) return(character(0))
    cv <- tryCatch(cv.glmnet(X, y, family = "binomial", alpha = 1,
                             nfolds = 5, type.measure = "deviance"),
                   error = function(e) NULL)
    if (is.null(cv)) return(character(0))
    cf <- coef(cv, s = "lambda.1se"); tm <- rownames(cf)[as.numeric(cf) != 0]
    tm <- setdiff(tm, "(Intercept)")
    unique(unlist(lapply(nest_cands, function(v) if (any(grepl(paste0("^", v), tm))) v)))
  }
  
  # FULL headline-selection pipeline on an imputed training frame
  .select_headline <- function(dtr) {
    u <- .univ_nom(dtr)
    set.seed(2024)
    tr <- tryCatch(rpart(reformulate(nest_cands, PRIM), data = dtr, method = "class",
                         control = rpart.control(cp = 0, minsplit = 10,
                                                 maxdepth = 5, xval = 10)),
                   error = function(e) NULL)
    tree_top <- character(0)
    if (!is.null(tr)) {
      cpt <- tr$cptable; im <- which.min(cpt[, "xerror"])
      th  <- cpt[im, "xerror"] + cpt[im, "xstd"]
      cp1 <- cpt[which(cpt[, "xerror"] <= th)[1], "CP"]
      tf  <- prune(tr, cp = cp1)
      sp  <- unique(tf$frame$var[tf$frame$var != "<leaf>"])
      if (length(sp) >= 3 && !is.null(tf$variable.importance) &&
          length(tf$variable.importance) > 0)
        tree_top <- names(sort(tf$variable.importance, decreasing = TRUE))[
          1:min(15, length(tf$variable.importance))]
    }
    r <- .rf_top15(dtr, ntree = 500)
    l <- .lasso_1se(dtr)
    allc  <- unique(c(u, tree_top, r, l))
    votes <- vapply(allc, function(v)
      (v %in% u) + (v %in% tree_top) + (v %in% r) + (v %in% l), integer(1))
    confirmed <- sort(allc[votes >= 2L])
    if (!"hb_mean" %in% confirmed) confirmed <- sort(c(confirmed, "hb_mean"))
    
    nB <- nrow(dtr)
    sv <- matrix(0, nrow = length(nest_cands), ncol = B_INNER,
                 dimnames = list(nest_cands, NULL))
    for (b in seq_len(B_INNER)) {
      ib <- sample.int(nB, nB, replace = TRUE); db <- dtr[ib, , drop = FALSE]
      if (length(unique(db[[PRIM]])) < 2) next
      ub <- .univ_nom(db); rb <- .rf_top15(db, ntree = 300); lb <- .lasso_1se(db)
      su <- unique(c(ub, rb, lb))
      for (v in su)
        if (((v %in% ub) + (v %in% rb) + (v %in% lb)) >= 2 && v %in% rownames(sv))
          sv[v, b] <- 1
    }
    sp_stab <- rowSums(sv) / B_INNER
    stable  <- names(sp_stab)[sp_stab >= 0.60]
    if (length(stable) == 0) stable <- "hb_mean"
    
    tier1 <- "hb_mean"; tier2 <- setdiff(stable, tier1)
    tier3 <- setdiff(confirmed, c(tier1, tier2))
    base  <- unique(c(tier1, tier2))
    base_roc <- .rf_oob_roc(dtr, base)
    base_ci  <- ci.auc(base_roc); base_auc <- as.numeric(auc(base_roc))
    retained <- character(0)
    for (cand in tier3) {
      tr_roc <- .rf_oob_roc(dtr, c(base, cand)); tr_ci <- ci.auc(tr_roc)
      dp <- tryCatch(roc.test(base_roc, tr_roc, method = "delong", paired = TRUE)$p.value,
                     error = function(e) NA_real_)
      if (as.numeric(auc(tr_roc)) > base_auc &&
          ((!is.na(dp) && dp < 0.05) || tr_ci[1] > base_ci[3]))
        retained <- c(retained, cand)
    }
    sort(unique(c(base, retained)))
  }
  
  # OUTER 5x5 nested loop
  y_nest       <- as.character(df_nest[[PRIM]])
  nest_oof_all <- matrix(NA_real_, nrow = N_REPEATS_NEST, ncol = length(y_nest))
  nest_rep_auc <- numeric(N_REPEATS_NEST)
  feat_freq    <- setNames(integer(length(nest_cands)), nest_cands)
  n_fold_total <- 0L
  t0 <- Sys.time()
  for (rp in seq_len(N_REPEATS_NEST)) {
    fld <- make_strat_folds(y_nest, N_FOLDS_NEST, seed = 2024 + rp)
    oof <- rep(NA_real_, length(y_nest))
    for (k in seq_len(N_FOLDS_NEST)) {
      tr_i <- which(fld != k); te_i <- which(fld == k)
      dtr  <- df_nest[tr_i, , drop = FALSE]; dte <- df_nest[te_i, , drop = FALSE]
      im   <- .fit_imp(dtr)
      dtr  <- .apply_imp(dtr, im); dte <- .apply_imp(dte, im)
      fs   <- tryCatch(.select_headline(dtr), error = function(e) "hb_mean")
      hit  <- intersect(fs, names(feat_freq)); feat_freq[hit] <- feat_freq[hit] + 1L
      n_fold_total <- n_fold_total + 1L
      mn   <- min(table(dtr[[PRIM]]))
      pr   <- tryCatch({
        set.seed(2024)
        rf <- randomForest(reformulate(fs, PRIM), data = dtr, ntree = 500,
                           strata = dtr[[PRIM]], sampsize = c("NO" = mn, "YES" = mn))
        predict(rf, dte, type = "prob")[, "YES"]
      }, error = function(e) rep(mean(dtr[[PRIM]] == "YES"), nrow(dte)))
      oof[te_i] <- pr
    }
    nest_oof_all[rp, ] <- oof
    nest_rep_auc[rp]   <- tryCatch(as.numeric(auc(roc(y_nest, oof, levels = c("NO", "YES"),
                                                      direction = "<", quiet = TRUE))),
                                   error = function(e) NA_real_)
    cat(sprintf("  repeat %d/%d done | within-repeat pooled AUC = %.4f | elapsed %.1f min\n",
                rp, N_REPEATS_NEST, nest_rep_auc[rp],
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  
  nest_oof_mean   <- colMeans(nest_oof_all, na.rm = TRUE)
  nest_pooled_roc <- roc(y_nest, nest_oof_mean, levels = c("NO", "YES"),
                         direction = "<", quiet = TRUE)
  nest_pooled_auc <- as.numeric(auc(nest_pooled_roc))
  nest_pooled_ci  <- as.numeric(ci.auc(nest_pooled_roc))
  
  naive_oob <- as.numeric(auc(roc_ml_hr))
  naive_cv  <- tryCatch({
    cc <- cv_results[["high_risk_RF"]]
    v  <- if (!is.null(cc$per_rep)) cc$per_rep else if (!is.null(cc$aucs)) cc$aucs else NA_real_
    mean(v, na.rm = TRUE)
  }, error = function(e) NA_real_)
  nest_mean <- mean(nest_rep_auc, na.rm = TRUE)
  
  cat("\n---------------- NESTED-CV SUMMARY (primary: high_risk) ----------------\n")
  cat(sprintf("  Naive OOB AUC (fixed features, Sec 7A)      : %.4f\n", naive_oob))
  if (!is.na(naive_cv))
    cat(sprintf("  Naive 5x5 CV AUC (fixed features, Sec 8)    : %.4f\n", naive_cv))
  cat(sprintf("  NESTED 5x5 CV AUC (selection re-run/fold)   : %.4f  (sd %.4f, %d repeats)\n",
              nest_mean, sd(nest_rep_auc, na.rm = TRUE), N_REPEATS_NEST))
  cat(sprintf("  NESTED pooled AUC (all OOF)                 : %.4f  (95%% CI %.4f-%.4f)\n",
              nest_pooled_auc, nest_pooled_ci[1], nest_pooled_ci[3]))
  cat(sprintf("  SELECTION OPTIMISM (naive OOB - nested CV)  : %+.4f\n",
              naive_oob - nest_mean))
  cat(sprintf("\n  Feature-selection frequency across %d outer folds:\n", n_fold_total))
  ff <- sort(feat_freq[feat_freq > 0] / n_fold_total, decreasing = TRUE)
  for (nm in names(ff)) cat(sprintf("    %-16s %3.0f%%\n", nm, 100 * ff[nm]))
  
  saveRDS(list(nested_rep_auc = nest_rep_auc, nested_mean = nest_mean,
               nested_sd = sd(nest_rep_auc, na.rm = TRUE),
               nested_pooled_auc = nest_pooled_auc, nested_pooled_ci = nest_pooled_ci,
               naive_oob = naive_oob, naive_cv = naive_cv,
               selection_optimism = naive_oob - nest_mean,
               feat_freq = feat_freq, n_folds = n_fold_total,
               B_inner = B_INNER, n_repeats = N_REPEATS_NEST, n_folds_cv = N_FOLDS_NEST),
          "nested_cv_result.rds")
  cat("\nSaved: nested_cv_result.rds\n")
  
  cat("\n>>> PASTE-READY SENTENCE FOR THESIS (numbers auto-filled):\n")
  cat("    \"When the entire feature-selection pipeline (four-method vote, bootstrap\n")
  cat("     stability, and forward selection) was re-run inside every fold of the 5x5\n")
  cat("     repeated stratified cross-validation, with imputation fitted on the training\n")
  cat(sprintf("     fold only, the selection-inclusive AUC was %.3f (vs the fixed-feature OOB\n",
              nest_mean))
  cat(sprintf("     AUC of %.3f), an optimism of %.3f.\"\n",
              naive_oob, naive_oob - nest_mean))
} else {
  cat("\nSection 8G (nested CV) SKIPPED -- set RUN_NESTED_CV <- TRUE to run.\n")
}

# SECTION 8B: HB ABLATION
cat("\n============================================================\n")
cat("SECTION 8B: HB ABLATION\n")
cat("============================================================\n\n")

if ("hb_entry" %in% names(df_ml) && "hb_mean" %in% ed_feats) {
  ed_feats_hbentry <- c(setdiff(ed_feats, "hb_mean"), "hb_entry")
  ed_feats_hbentry <- ed_feats_hbentry[ed_feats_hbentry %in% names(df_ml)]
  
  df_hr_alt <- df_ml |>
    dplyr::select(all_of(ed_feats_hbentry), all_of(PRIMARY)) |>
    filter(!is.na(.data[[PRIMARY]]))
  df_hr_alt[[PRIMARY]] <- factor(df_hr_alt[[PRIMARY]], levels=c("NO","YES"))
  
  hr_min_alt <- min(table(df_hr_alt[[PRIMARY]]))
  set.seed(2024)
  rf_hr_alt <- randomForest(
    reformulate(ed_feats_hbentry, PRIMARY), data=df_hr_alt,
    ntree=500, strata=df_hr_alt[[PRIMARY]],
    sampsize=c("NO"=hr_min_alt, "YES"=hr_min_alt))
  
  oob_p_alt <- rf_hr_alt$votes[, "YES"] / rowSums(rf_hr_alt$votes)
  roc_alt   <- roc(df_hr_alt[[PRIMARY]], oob_p_alt,
                   levels=c("NO","YES"), direction="<", quiet=TRUE)
  ci_alt    <- ci.auc(roc_alt)
  
  delong_hb <- tryCatch(
    roc.test(roc_ml_hr, roc_alt, method="delong", paired=TRUE)$p.value,
    error=function(e) NA_real_)
  
  cat(sprintf("hb_mean  (headline):  AUC=%.4f (%.4f-%.4f)\n",
              as.numeric(auc(roc_ml_hr)),
              ci.auc(roc_ml_hr)[1], ci.auc(roc_ml_hr)[3]))
  cat(sprintf("hb_entry (ablation):  AUC=%.4f (%.4f-%.4f)\n",
              as.numeric(auc(roc_alt)), ci_alt[1], ci_alt[3]))
  cat(sprintf("DeLong test (paired): p=%.4f\n", delong_hb))
  cat(sprintf("Delta AUC: %+.4f -- %s\n",
              as.numeric(auc(roc_alt)) - as.numeric(auc(roc_ml_hr)),
              ifelse(!is.na(delong_hb) && delong_hb < 0.05,
                     "DIFFERENT", "EQUIVALENT (use hb_mean as headline per Bland-Altman)")))
  
  saveRDS(list(auc_hbmean = as.numeric(auc(roc_ml_hr)),
               ci_hbmean  = ci.auc(roc_ml_hr),
               auc_hbentry = as.numeric(auc(roc_alt)),
               ci_hbentry  = ci_alt,
               delong_p   = delong_hb),
          "hb_ablation.rds")
  cat("Saved: hb_ablation.rds\n")
} else {
  cat("hb_entry not in df_ml or hb_mean not in ed_feats -- ablation skipped\n")
}

# SECTION 8C: MICE-RF SENSITIVITY
cat("\n============================================================\n")
cat("SECTION 8C: MICE-RF SENSITIVITY\n")
cat("============================================================\n\n")

mice_path <- "ugib_mice_object.rds"
if (file.exists(mice_path)) {
  mice_obj <- readRDS(mice_path)
  m_imp    <- mice_obj$m
  cat(sprintf("Loaded MICE object with m=%d imputations.\n", m_imp))
  
  imp_aucs <- numeric(m_imp)
  for (k in seq_len(m_imp)) {
    df_k <- tryCatch(mice::complete(mice_obj, action=k),
                     error=function(e) NULL)
    if (is.null(df_k)) { imp_aucs[k] <- NA_real_; next }
    
    if (!PRIMARY %in% names(df_k) && nrow(df_k) == nrow(df_full)) {
      df_k[[PRIMARY]] <- df_full[[PRIMARY]]
    }
    if ("hb_mean" %in% ed_feats && !"hb_mean" %in% names(df_k) &&
        all(c("hb_ed","hb_entry") %in% names(df_k))) {
      df_k$hb_mean <- (df_k$hb_ed + df_k$hb_entry) / 2
    }
    
    feats_k <- intersect(ed_feats, names(df_k))
    if (!PRIMARY %in% names(df_k) || length(feats_k) < 2) {
      imp_aucs[k] <- NA_real_; next
    }
    df_k[[PRIMARY]] <- factor(df_k[[PRIMARY]], levels=c("NO","YES"))
    if (length(unique(df_k[[PRIMARY]])) < 2) {
      imp_aucs[k] <- NA_real_; next
    }
    cls_min_k <- min(table(df_k[[PRIMARY]]))
    set.seed(2024 + k)
    rf_k <- tryCatch(randomForest(
      reformulate(feats_k, PRIMARY), data=df_k,
      ntree=500, strata=df_k[[PRIMARY]],
      sampsize=c("NO"=cls_min_k, "YES"=cls_min_k)),
      error=function(e) NULL)
    if (is.null(rf_k)) { imp_aucs[k] <- NA_real_; next }
    oob_p_k <- rf_k$votes[, "YES"] / rowSums(rf_k$votes)
    r_k <- tryCatch(roc(df_k[[PRIMARY]], oob_p_k,
                        levels=c("NO","YES"), direction="<", quiet=TRUE),
                    error=function(e) NULL)
    imp_aucs[k] <- if (is.null(r_k)) NA_real_ else as.numeric(auc(r_k))
  }
  
  imp_auc_mean <- mean(imp_aucs, na.rm=TRUE)
  imp_auc_sd   <- sd(imp_aucs,   na.rm=TRUE)
  imp_auc_lo   <- imp_auc_mean - 1.96 * imp_auc_sd
  imp_auc_hi   <- imp_auc_mean + 1.96 * imp_auc_sd
  
  cat(sprintf("Per-imputation AUCs: range %.4f - %.4f (n=%d/%d converged)\n",
              min(imp_aucs, na.rm=TRUE), max(imp_aucs, na.rm=TRUE),
              sum(!is.na(imp_aucs)), m_imp))
  cat(sprintf("Pooled (mean +/- 1.96 SD): %.4f (%.4f-%.4f)\n",
              imp_auc_mean, imp_auc_lo, imp_auc_hi))
  cat(sprintf("Median-imputation headline AUC: %.4f\n", as.numeric(auc(roc_ml_hr))))
  cat(sprintf("Delta: %+.4f -- %s\n",
              imp_auc_mean - as.numeric(auc(roc_ml_hr)),
              ifelse(abs(imp_auc_mean - as.numeric(auc(roc_ml_hr))) < 0.01,
                     "EQUIVALENT", "DIFFERS by >0.01")))
  
  saveRDS(list(imp_aucs=imp_aucs, m=m_imp,
               auc_mean=imp_auc_mean, auc_sd=imp_auc_sd,
               auc_lo=imp_auc_lo, auc_hi=imp_auc_hi,
               headline_median_auc=as.numeric(auc(roc_ml_hr))),
          "mice_rf_sensitivity.rds")
  cat("Saved: mice_rf_sensitivity.rds\n")
} else {
  cat(sprintf("MICE object %s not found -- sensitivity skipped\n", mice_path))
}

# SECTION 9: MODEL EVALUATION
cat("============================================================\n")
cat("SECTION 9: MODEL EVALUATION\n")
cat("============================================================\n\n")

cat("--- 9A. Confusion matrix ---\n\n")
conf_mat <- table(Predicted=oob_pred_tc, Actual=true_labs)
print(conf_mat)

class_metrics <- do.call(rbind, lapply(cl_levels, function(cls){
  tp <- sum(oob_pred_tc==cls & true_labs==cls)
  fp <- sum(oob_pred_tc==cls & true_labs!=cls)
  fn <- sum(oob_pred_tc!=cls & true_labs==cls)
  tn <- sum(oob_pred_tc!=cls & true_labs!=cls)
  sens <- ifelse(tp+fn>0, tp/(tp+fn), NA)
  spec <- ifelse(tn+fp>0, tn/(tn+fp), NA)
  ppv  <- ifelse(tp+fp>0, tp/(tp+fp), NA)
  npv  <- ifelse(tn+fn>0, tn/(tn+fn), NA)
  f1   <- ifelse(!is.na(ppv)&!is.na(sens)&(ppv+sens)>0, 2*ppv*sens/(ppv+sens), NA)
  data.frame(Class=cls, TP=tp, FP=fp, FN=fn, TN=tn,
             Sensitivity=round(sens,3), Specificity=round(spec,3),
             PPV=round(ppv,3), NPV=round(npv,3), F1=round(f1,3))
}))
cat("\nPer-class metrics:\n"); print(class_metrics)

po <- mean(as.character(oob_pred_tc)==as.character(true_labs))
pe <- sum(sapply(cl_levels, function(cls)
  mean(oob_pred_tc==cls,na.rm=T) * mean(true_labs==cls,na.rm=T)))
kappa <- round((po-pe)/(1-pe), 3)
oob_acc <- round(po*100, 1)
cat(sprintf("\nOOB accuracy: %.1f%% | Cohen kappa: %.3f\n\n", oob_acc, kappa))

conf_df <- as.data.frame(conf_mat)
.conf_label_map <- c("0_units"="0 units", "1to2"="1-2 units", "3plus"="3 plus units")
conf_df$Actual    <- factor(.conf_label_map[as.character(conf_df$Actual)],    levels=c("0 units","1-2 units","3 plus units"))
conf_df$Predicted <- factor(.conf_label_map[as.character(conf_df$Predicted)], levels=c("3 plus units","1-2 units","0 units"))
plots[["p_conf"]] <- ggplot(conf_df, aes(x=Actual, y=Predicted, fill=Freq)) +
  geom_tile(colour="white", linewidth=0.8) +
  geom_text(aes(label=ifelse(Freq>0,as.character(Freq),"")),
            size=6, fontface="bold") +
  scale_fill_gradient(low="#EBF5FB", high=COL_BLUE, name="Count") +
  labs(title="Confusion matrix: BURST 3-class model (OOB predictions)",
       subtitle=sprintf("OOB accuracy: %.1f%% | Cohen kappa: %.3f | Balanced sampling",
                        oob_acc, kappa),
       x="Actual class", y="Predicted class") + theme_ugib()
print(plots[["p_conf"]])

cat("--- 9B. AUC one-vs-rest (3-class) ---\n\n")
oob_probs_tc <- rf_tc$votes / rowSums(rf_tc$votes)
roc_ovr <- list()
cols_ovr <- c(COL_BLUE, COL_RED, COL_GREEN)
for (cls in cl_levels) {
  bin_out <- factor(ifelse(df_ed$transfusion_class==cls,"YES","NO"), levels=c("NO","YES"))
  ok <- !is.na(oob_probs_tc[,cls])
  roc_ovr[[cls]] <- tryCatch(
    roc(bin_out[ok], oob_probs_tc[ok,cls], levels=c("NO","YES"), direction="<", quiet=TRUE),
    error=function(e) NULL)
  if(!is.null(roc_ovr[[cls]])){
    ci <- ci.auc(roc_ovr[[cls]])
    cat(sprintf("  %-12s AUC=%.3f (95%% CI %.3f-%.3f)\n",
                cls, auc(roc_ovr[[cls]]), ci[1], ci[3]))
  }
}
cat(sprintf("\n  major_tx binary AUC=%.3f\n\n", auc(roc_ml_mt)))

par(mar=c(5,5,4,2))
plot(roc_ovr[[1]], col=cols_ovr[1], lwd=2.5, legacy.axes = TRUE,
     main="BURST AUC one-vs-rest: 3-class transfusion model",
     xlim = c(1, 0), ylim = c(0, 1), asp = NA,
     xaxs = "i", yaxs = "i")
for(i in 2:3) if(!is.null(roc_ovr[[i]]))
  plot(roc_ovr[[i]], col=cols_ovr[i], lwd=2.5, legacy.axes = TRUE, add=TRUE)
legend("bottomright",bty="n",lwd=2.5,col=cols_ovr,
       legend=sapply(seq_along(cl_levels), function(i)
         sprintf("%s AUC=%.3f", cl_levels[i],
                 ifelse(is.null(roc_ovr[[i]]),NA,auc(roc_ovr[[i]])))))

cat("--- 9C. Hosmer-Lemeshow calibration ---\n\n")
hl_df <- data.frame(prob=oob_prob_mt,
                    actual=as.numeric(df_ed$major_tx=="YES")) |>
  filter(!is.na(prob)) |> arrange(prob) |>
  mutate(grp=cut(seq_len(n()), breaks=quantile(seq_len(n()),
                                               probs=seq(0,1,.1)), include.lowest=TRUE, labels=FALSE))

hl_g <- hl_df |> group_by(grp) |>
  summarise(n=n(), obs=sum(actual), exp=sum(prob), .groups="drop") |>
  mutate(obs_no=n-obs, exp_no=n-exp,
         obs_rate=obs/n, exp_rate=exp/n,
         se=sqrt(obs_rate*(1-obs_rate)/n))

hl_stat <- sum((hl_g$obs-hl_g$exp)^2/hl_g$exp +
                 (hl_g$obs_no-hl_g$exp_no)^2/hl_g$exp_no)
hl_p <- 1 - pchisq(hl_stat, df=8)
cat(sprintf("Hosmer-Lemeshow: chi2=%.3f df=8 p=%.4f\n", hl_stat, hl_p))
cat(ifelse(hl_p>0.05, "Calibration: ADEQUATE\n\n", "Calibration: POOR\n\n"))

plots[["p_cal"]] <- ggplot(hl_g, aes(x=exp_rate, y=obs_rate)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey50", linewidth=0.9) +
  geom_point(size=4, colour=COL_BLUE) +
  geom_errorbar(aes(ymin=pmax(obs_rate-1.96*se,0),
                    ymax=pmin(obs_rate+1.96*se,1)), width=0.01, colour=COL_BLUE, linewidth=0.8) +
  annotate("text",x=0.08,y=0.85,label=sprintf("HL p = %.3f",hl_p),size=4,colour="grey30") +
  scale_x_continuous(limits=c(0,1), labels=pct_fmt) +
  scale_y_continuous(limits=c(0,1), labels=pct_fmt) +
  labs(title="Calibration: BURST on major transfusion",
       subtitle=sprintf("Hosmer-Lemeshow chi2=%.2f, df=8, p=%.4f", hl_stat, hl_p),
       x="Mean predicted probability", y="Observed event rate") +
  theme_ugib() + theme(legend.position="none")
print(plots[["p_cal"]])

cat("--- 9D. Brier score ---\n\n")
brier   <- round(mean((oob_prob_mt - as.numeric(df_ed$major_tx=="YES"))^2, na.rm=T), 4)
prev_mt <- mean(df_ed$major_tx=="YES", na.rm=T)
brier_ref <- prev_mt*(1-prev_mt)
bss     <- round(1-brier/brier_ref, 4)
cat(sprintf("Brier: %.4f | Naive: %.4f | Skill score: %.4f\n\n", brier, brier_ref, bss))

saveRDS(list(class_metrics=class_metrics, roc_ovr=roc_ovr,
             roc_mt=roc_ml_mt, roc_hr=roc_ml_hr,
             oob_acc=oob_acc, kappa=kappa,
             hl_stat=hl_stat, hl_p=hl_p,
             brier=brier, bss=bss), "model_evaluation.rds")

cat("\n--- 9F. Multi-outcome (death, rebleeding) -- V2.4.5 descriptive-only ---\n")

score_auc_safe_9f <- function(score_vec, outcome_vec) {
  ok <- !is.na(score_vec) & !is.na(outcome_vec)
  if (length(unique(outcome_vec[ok])) < 2L) return(NA_real_)
  as.numeric(auc(roc(outcome_vec[ok], score_vec[ok],
                     levels=levels(factor(outcome_vec[ok])),
                     direction="<", quiet=TRUE)))
}

multi_outcome_results <- list()
for (oc in c("outcome_death","rebleeding")) {
  if (!oc %in% names(df_full)) {
    cat(sprintf("  %s: column not in df_full -- skipping\n", oc)); next
  }
  oc_vec <- df_full[[oc]]
  n_evt  <- sum(oc_vec == "YES", na.rm=TRUE)
  n_tot  <- sum(!is.na(oc_vec))
  cat(sprintf("  %s: %d / %d events (%.1f%%)\n",
              oc, n_evt, n_tot, 100*n_evt/n_tot))
  cat(sprintf("    n_events=%d below threshold for stable 5x5 CV. Reporting descriptive only.\n",
              n_evt))
  cat(sprintf("    Comparator-score AUCs (no CV; for thesis lit-benchmarking):\n"))
  out_oc <- list(n_events = n_evt, n_total = n_tot,
                 event_rate_pct = round(100 * n_evt / n_tot, 1),
                 status = "descriptive_only_insufficient_for_CV")
  for (sc_n in c("gbs_score","aims65_score","pre_rockall",
                 "full_rockall","canuka")) {
    if (!sc_n %in% names(df_full)) next
    a <- score_auc_safe_9f(as.numeric(df_full[[sc_n]]), oc_vec)
    cat(sprintf("      %-14s AUC = %s\n", sc_n,
                ifelse(is.na(a), "NA", sprintf("%.3f", a))))
    out_oc[[sc_n]] <- a
  }
  multi_outcome_results[[oc]] <- out_oc
  cat("\n")
}
saveRDS(multi_outcome_results, "multi_outcome_results.rds")
cat("Saved: multi_outcome_results.rds (V2.4.5 descriptive-only)\n\n")

# SECTION 9G: PRIMARY (high_risk) CALIBRATION
cat("\n============================================================\n")
cat("SECTION 9G: PRIMARY CALIBRATION\n")
cat("============================================================\n\n")

hl_df_pri <- data.frame(prob = oob_prob_hr,
                        actual = as.numeric(df_hr$high_risk == "YES")) |>
  filter(!is.na(prob)) |> arrange(prob) |>
  mutate(grp = cut(seq_len(n()),
                   breaks = quantile(seq_len(n()), probs = seq(0, 1, .1)),
                   include.lowest = TRUE, labels = FALSE))

hl_g_pri <- hl_df_pri |> group_by(grp) |>
  summarise(n = n(), obs = sum(actual), exp = sum(prob), .groups = "drop") |>
  mutate(obs_no   = n - obs,
         exp_no   = n - exp,
         obs_rate = obs / n,
         exp_rate = exp / n,
         se       = sqrt(obs_rate * (1 - obs_rate) / n))

hl_stat_pri <- sum((hl_g_pri$obs - hl_g_pri$exp)^2 / hl_g_pri$exp +
                     (hl_g_pri$obs_no - hl_g_pri$exp_no)^2 / hl_g_pri$exp_no)
hl_p_pri <- 1 - pchisq(hl_stat_pri, df = 8)
cat(sprintf("Hosmer-Lemeshow (PRIMARY high_risk): chi2=%.3f df=8 p=%.4f\n",
            hl_stat_pri, hl_p_pri))
cat(ifelse(hl_p_pri > 0.05, "Calibration: ADEQUATE\n\n", "Calibration: POOR\n\n"))

brier_pri     <- round(mean((oob_prob_hr -
                               as.numeric(df_hr$high_risk == "YES"))^2, na.rm = TRUE), 4)
prev_hr       <- mean(df_hr$high_risk == "YES", na.rm = TRUE)
brier_ref_pri <- prev_hr * (1 - prev_hr)
bss_pri       <- round(1 - brier_pri / brier_ref_pri, 4)
cat(sprintf("Brier (PRIMARY): %.4f | Naive: %.4f | Skill score: %.4f\n\n",
            brier_pri, brier_ref_pri, bss_pri))

calib_lr <- tryCatch({
  eps <- 1e-6
  p_clip <- pmin(pmax(oob_prob_hr, eps), 1 - eps)
  logit_p <- log(p_clip / (1 - p_clip))
  y_num   <- as.numeric(df_hr$high_risk == "YES")
  ok      <- !is.na(logit_p) & !is.na(y_num)
  fit_int <- glm(y_num[ok] ~ offset(logit_p[ok]), family = binomial)
  fit_slp <- glm(y_num[ok] ~ logit_p[ok],         family = binomial)
  list(
    intercept_calib = round(coef(fit_int)[1], 4),
    intercept_se    = round(summary(fit_int)$coefficients[1, "Std. Error"], 4),
    slope           = round(coef(fit_slp)[2], 4),
    slope_se        = round(summary(fit_slp)$coefficients[2, "Std. Error"], 4),
    slope_ci_lo     = round(coef(fit_slp)[2] - 1.96 * summary(fit_slp)$coefficients[2, "Std. Error"], 4),
    slope_ci_hi     = round(coef(fit_slp)[2] + 1.96 * summary(fit_slp)$coefficients[2, "Std. Error"], 4))
}, error = function(e) {
  list(intercept_calib=NA_real_, intercept_se=NA_real_,
       slope=NA_real_, slope_se=NA_real_,
       slope_ci_lo=NA_real_, slope_ci_hi=NA_real_)
})

cat("Calibration-in-the-large + slope (Steyerberg 2009):\n")
cat(sprintf("  Intercept (slope fixed at 1): %.4f (SE %.4f)  [0 = perfect]\n",
            calib_lr$intercept_calib, calib_lr$intercept_se))
cat(sprintf("  Slope:                        %.4f (95%% CI %.4f-%.4f)  [1 = perfect]\n",
            calib_lr$slope, calib_lr$slope_ci_lo, calib_lr$slope_ci_hi))
if (!is.na(calib_lr$slope)) {
  if (calib_lr$slope < 0.85) {
    cat("  -> Slope <0.85: predictions likely TOO EXTREME (overfit). Consider shrinkage.\n")
  } else if (calib_lr$slope > 1.15) {
    cat("  -> Slope >1.15: predictions TOO CONSERVATIVE.\n")
  } else {
    cat("  -> Slope within 0.85-1.15: calibration acceptable.\n")
  }
}
cat("\n")

plots[["p_cal_pri"]] <- ggplot(hl_g_pri, aes(x = exp_rate, y = obs_rate)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey50", linewidth = 0.9) +
  geom_point(size = 4, colour = COL_BLUE) +
  geom_errorbar(aes(ymin = pmax(obs_rate - 1.96 * se, 0),
                    ymax = pmin(obs_rate + 1.96 * se, 1)),
                width = 0.01, colour = COL_BLUE, linewidth = 0.8) +
  annotate("text", x = 0.08, y = 0.92,
           label = sprintf("HL p = %.3f", hl_p_pri),
           size = 4, colour = "grey30") +
  annotate("text", x = 0.08, y = 0.85,
           label = sprintf("Brier = %.3f", brier_pri),
           size = 4, colour = "grey30") +
  scale_x_continuous(limits = c(0, 1), labels = pct_fmt) +
  scale_y_continuous(limits = c(0, 1), labels = pct_fmt) +
  labs(title    = "Calibration: BURST on primary outcome",
       subtitle = sprintf("Hosmer-Lemeshow chi2=%.2f, df=8, p=%.4f",
                          hl_stat_pri, hl_p_pri),
       x = "Mean predicted probability", y = "Observed event rate") +
  theme_ugib() + theme(legend.position = "none")
print(plots[["p_cal_pri"]])

saveRDS(list(hl_stat = hl_stat_pri, hl_p = hl_p_pri,
             brier = brier_pri, brier_ref = brier_ref_pri, bss = bss_pri,
             intercept_calib = calib_lr$intercept_calib,
             slope = calib_lr$slope,
             slope_ci_lo = calib_lr$slope_ci_lo,
             slope_ci_hi = calib_lr$slope_ci_hi,
             hl_groups = hl_g_pri),
        "calibration_primary.rds")
cat("Saved: calibration_primary.rds\n\n")

# SECTION 9H: PRIMARY THRESHOLD-SENSITIVITY TABLE
cat("\n============================================================\n")
cat("SECTION 9H: PRIMARY THRESHOLD SENSITIVITY TABLE\n")
cat("============================================================\n\n")

thr_grid <- c(0.30, 0.40, 0.50, 0.60, 0.70)
y_pri    <- factor(df_hr$high_risk, levels = c("NO", "YES"))
threshold_metrics <- do.call(rbind, lapply(thr_grid, function(t) {
  pred <- factor(ifelse(oob_prob_hr >= t, "YES", "NO"),
                 levels = c("NO", "YES"))
  tp <- sum(pred == "YES" & y_pri == "YES", na.rm = TRUE)
  fp <- sum(pred == "YES" & y_pri == "NO",  na.rm = TRUE)
  fn <- sum(pred == "NO"  & y_pri == "YES", na.rm = TRUE)
  tn <- sum(pred == "NO"  & y_pri == "NO",  na.rm = TRUE)
  sens <- if (tp + fn > 0) tp / (tp + fn) else NA_real_
  spec <- if (tn + fp > 0) tn / (tn + fp) else NA_real_
  ppv  <- if (tp + fp > 0) tp / (tp + fp) else NA_real_
  npv  <- if (tn + fn > 0) tn / (tn + fn) else NA_real_
  f1   <- if (!is.na(ppv) && !is.na(sens) && (ppv + sens) > 0)
    2 * ppv * sens / (ppv + sens) else NA_real_
  data.frame(
    threshold   = t,
    n_pos       = tp + fp,
    pct_pos     = round(100 * (tp + fp) / max(tp + fp + tn + fn, 1), 1),
    sensitivity = round(sens, 3),
    specificity = round(spec, 3),
    PPV         = round(ppv,  3),
    NPV         = round(npv,  3),
    F1          = round(f1,   3),
    Youden_J    = round(sens + spec - 1, 3),
    stringsAsFactors = FALSE
  )
}))
cat("Threshold sensitivity (PRIMARY: high_risk, OOB probabilities):\n")
print(threshold_metrics, row.names = FALSE)

youden_opt <- coords(roc_ml_hr, x = "best", best.method = "youden",
                     ret = c("threshold", "sensitivity", "specificity"),
                     transpose = FALSE)
cat(sprintf("\nOptimal Youden threshold: %.3f (sens=%.3f, spec=%.3f)\n",
            as.numeric(youden_opt["threshold"]),
            as.numeric(youden_opt["sensitivity"]),
            as.numeric(youden_opt["specificity"])))

event_probs <- oob_prob_hr[y_pri == "YES"]
thr_sens100 <- min(event_probs, na.rm = TRUE) - 1e-6
pred100 <- factor(ifelse(oob_prob_hr >= thr_sens100, "YES", "NO"),
                  levels = c("NO", "YES"))
tp100 <- sum(pred100 == "YES" & y_pri == "YES", na.rm = TRUE)
fp100 <- sum(pred100 == "YES" & y_pri == "NO",  na.rm = TRUE)
fn100 <- sum(pred100 == "NO"  & y_pri == "YES", na.rm = TRUE)
tn100 <- sum(pred100 == "NO"  & y_pri == "NO",  na.rm = TRUE)
sens100 <- if (tp100+fn100>0) tp100/(tp100+fn100) else NA_real_
spec100 <- if (tn100+fp100>0) tn100/(tn100+fp100) else NA_real_
ppv100  <- if (tp100+fp100>0) tp100/(tp100+fp100) else NA_real_
npv100  <- if (tn100+fn100>0) tn100/(tn100+fn100) else NA_real_
f1_100  <- if (!is.na(ppv100) && !is.na(sens100) && (ppv100+sens100)>0)
  2*ppv100*sens100/(ppv100+sens100) else NA_real_
sens100_row <- data.frame(
  threshold   = round(thr_sens100, 4),
  n_pos       = tp100 + fp100,
  pct_pos     = round(100 * (tp100 + fp100) / max(tp100+fp100+tn100+fn100, 1), 1),
  sensitivity = round(sens100, 3),
  specificity = round(spec100, 3),
  PPV         = round(ppv100, 3),
  NPV         = round(npv100, 3),
  F1          = round(f1_100, 3),
  Youden_J    = round(sens100 + spec100 - 1, 3),
  stringsAsFactors = FALSE)
threshold_metrics_v245 <- rbind(sens100_row, threshold_metrics)
threshold_metrics_v245 <- threshold_metrics_v245[order(threshold_metrics_v245$threshold), ]
rownames(threshold_metrics_v245) <- NULL
cat(sprintf("\n--- 9H V2.4.5 B3: 100%%-sensitivity threshold (per Shung 2020) ---\n"))
cat(sprintf("  Threshold at sens=1.00:  %.4f\n", thr_sens100))
cat(sprintf("  Specificity at sens=1.00: %.3f (%.1f%% of cohort safe-dischargeable)\n",
            spec100, 100 * (tn100 + fn100) / (tp100+fp100+tn100+fn100)))
cat(sprintf("  Reference: Shung 2020 reported spec=0.26 at sens=1.00 in external validation.\n"))
cat("\nUpdated threshold-sensitivity table (sens=1.00 row added):\n")
print(threshold_metrics_v245, row.names = FALSE)

saveRDS(list(grid_metrics      = threshold_metrics_v245,
             youden_opt       = youden_opt,
             sens100_thr      = thr_sens100,
             sens100_spec     = spec100,
             sens100_npv      = npv100,
             sens100_ppv      = ppv100,
             sens100_n_below  = tn100 + fn100,
             sens100_pct_below = 100 * (tn100 + fn100) /
               (tp100+fp100+tn100+fn100)),
        "threshold_metrics_primary.rds")
cat("\nSaved: threshold_metrics_primary.rds (V2.4.5 with sens=1.00 row)\n\n")

# SECTION 9H.1: Operating-point table
cat("\n--- Operating-point table at thresholds 0.30 / 0.50 / 0.70 ---\n")
.thr_block_a <- c(0.30, 0.50, 0.70)
operating_points <- do.call(rbind, lapply(.thr_block_a, function(t) {
  pred <- factor(ifelse(oob_prob_hr >= t, "YES", "NO"), levels = c("NO","YES"))
  obs  <- df_hr$high_risk
  tp <- sum(pred == "YES" & obs == "YES")
  fp <- sum(pred == "YES" & obs == "NO")
  fn <- sum(pred == "NO"  & obs == "YES")
  tn <- sum(pred == "NO"  & obs == "NO")
  data.frame(
    threshold = t,
    sens      = tp / (tp + fn),
    spec      = tn / (tn + fp),
    ppv       = tp / (tp + fp),
    npv       = tn / (tn + fn),
    n_above   = tp + fp,
    n_below   = tn + fn,
    stringsAsFactors = FALSE
  )
}))
print(round(operating_points, 3), row.names = FALSE)
saveRDS(operating_points, "operating_points.rds")
cat("Saved: operating_points.rds\n")

cat("\n--- Operating-point summary ---\n")
print(round(operating_points, 3), row.names = FALSE)

# SECTION 9I: LEARNING CURVE
cat("\n============================================================\n")
cat("SECTION 9I: LEARNING CURVE\n")
cat("============================================================\n\n")

lc_sizes <- c(40, 60, 80, 100, 130)
lc_sizes <- lc_sizes[lc_sizes <= nrow(df_hr)]
n_reps   <- 20

lc_results <- list()
for (n_sub in lc_sizes) {
  aucs <- numeric(n_reps)
  for (r in seq_len(n_reps)) {
    set.seed(2024 + r)
    idx_yes   <- which(df_hr$high_risk == "YES")
    idx_no    <- which(df_hr$high_risk == "NO")
    p_yes     <- length(idx_yes) / nrow(df_hr)
    n_yes_sub <- max(2, round(n_sub * p_yes))
    n_no_sub  <- n_sub - n_yes_sub
    n_yes_sub <- min(n_yes_sub, length(idx_yes))
    n_no_sub  <- min(n_no_sub,  length(idx_no))
    samp      <- c(sample(idx_yes, n_yes_sub),
                   sample(idx_no,  n_no_sub))
    d_sub     <- df_hr[samp, , drop = FALSE]
    if (length(unique(d_sub$high_risk)) < 2) { aucs[r] <- NA_real_; next }
    hr_min_sub <- min(table(d_sub$high_risk))
    rf_sub <- tryCatch(randomForest(
      reformulate(ed_feats, "high_risk"), data = d_sub, ntree = 500,
      strata = d_sub$high_risk,
      sampsize = c("NO" = hr_min_sub, "YES" = hr_min_sub)),
      error = function(e) NULL)
    if (is.null(rf_sub)) { aucs[r] <- NA_real_; next }
    p_oob <- rf_sub$votes[, "YES"] / rowSums(rf_sub$votes)
    r_sub <- tryCatch(roc(d_sub$high_risk, p_oob,
                          levels = c("NO", "YES"),
                          direction = "<", quiet = TRUE),
                      error = function(e) NULL)
    aucs[r] <- if (is.null(r_sub)) NA_real_ else as.numeric(auc(r_sub))
  }
  lc_results[[as.character(n_sub)]] <- aucs
  cat(sprintf("  n=%3d: mean AUC=%.3f  sd=%.3f  (%d/%d reps successful)\n",
              n_sub, mean(aucs, na.rm = TRUE),
              sd(aucs, na.rm = TRUE),
              sum(!is.na(aucs)), n_reps))
}

lc_df <- do.call(rbind, lapply(names(lc_results), function(k) {
  a <- lc_results[[k]]
  data.frame(n        = as.integer(k),
             mean_auc = mean(a, na.rm = TRUE),
             sd_auc   = sd(a,   na.rm = TRUE),
             lo       = mean(a, na.rm = TRUE) - sd(a, na.rm = TRUE),
             hi       = mean(a, na.rm = TRUE) + sd(a, na.rm = TRUE),
             stringsAsFactors = FALSE)
}))

plots[["p_learning_curve"]] <- ggplot(lc_df, aes(x = n, y = mean_auc)) +
  geom_ribbon(aes(ymin = pmax(lo, 0.5), ymax = pmin(hi, 1)),
              fill = COL_BLUE, alpha = 0.15) +
  geom_line(colour  = COL_BLUE, linewidth = 0.9) +
  geom_point(colour = COL_BLUE, size = 3) +
  geom_hline(yintercept = as.numeric(auc(roc_ml_hr)),
             linetype = "dashed", colour = COL_RED, linewidth = 0.7) +
  annotate("text", x = max(lc_df$n) * 0.55,
           y = as.numeric(auc(roc_ml_hr)) + 0.012,
           label = sprintf("Full-cohort OOB AUC = %.3f", auc(roc_ml_hr)),
           colour = COL_RED, size = 3.5) +
  scale_x_continuous(breaks = lc_df$n) +
  scale_y_continuous(limits = c(0.5, 1), breaks = seq(0.5, 1, 0.1)) +
  labs(title    = "BURST learning curve on primary outcome",
       subtitle = sprintf("Mean OOB AUC +/- 1 SD across %d stratified resamples per size",
                          n_reps),
       x = "Training cohort size (n)", y = "OOB AUC") +
  theme_ugib()
print(plots[["p_learning_curve"]])

saveRDS(list(sizes = lc_sizes, n_reps = n_reps,
             aucs = lc_results, summary = lc_df),
        "learning_curve.rds")
cat("Saved: learning_curve.rds\n\n")

# SECTION 9J: Age-subgroup AUCs
cat("\n============================================================\n")
cat("SECTION 9J: Age-subgroup AUCs\n")
cat("============================================================\n\n")

.df_hr_age <- df_ml |>
  dplyr::select(all_of(unique(c("age", ed_feats))), high_risk) |>
  filter(!is.na(high_risk))
.df_ed_age <- df_ml |>
  dplyr::select(all_of(unique(c("age", "patient_id", ed_feats))),
                transfusion_class, major_tx, high_risk) |>
  filter(!is.na(transfusion_class))

stopifnot(nrow(.df_hr_age) == nrow(df_hr))
stopifnot(nrow(.df_ed_age) == nrow(df_ed))

.age_hr <- .df_hr_age$age
.age_ed <- .df_ed_age$age

.grp_hr <- ifelse(.age_hr < 70, "under70", "ge70")
.grp_ed <- ifelse(.age_ed < 70, "under70", "ge70")

cat("--- Subgroup sizes ---\n")
cat(sprintf("  <70: n=%d (spec 50) | high_risk events=%d (spec 24) | major_tx events=%d (spec 10)\n",
            sum(.grp_hr == "under70"),
            sum(.grp_hr == "under70" & df_hr$high_risk == "YES"),
            sum(.grp_ed == "under70" & df_ed$major_tx == "YES")))
cat(sprintf("  >=70: n=%d (spec 83) | high_risk events=%d (spec 54) | major_tx events=%d (spec 18)\n",
            sum(.grp_hr == "ge70"),
            sum(.grp_hr == "ge70" & df_hr$high_risk == "YES"),
            sum(.grp_ed == "ge70" & df_ed$major_tx == "YES")))

cat("\n--- 9J.1 PRIMARY (high_risk) subgroup AUCs ---\n")
subgroup_primary <- do.call(rbind, lapply(c("under70", "ge70"), function(g) {
  idx <- which(.grp_hr == g)
  r <- tryCatch(
    roc(df_hr$high_risk[idx], oob_prob_hr[idx],
        levels = c("NO","YES"), direction = "<", quiet = TRUE),
    error = function(e) NULL)
  if (is.null(r)) return(data.frame(group = g, n = length(idx),
                                    events = NA, auc = NA, lci = NA, uci = NA))
  ci <- pROC::ci.auc(r, method = "delong")
  data.frame(group = g, n = length(idx),
             events = sum(df_hr$high_risk[idx] == "YES"),
             auc = as.numeric(auc(r)),
             lci = as.numeric(ci[1]),
             uci = as.numeric(ci[3]))
}))
for (i in seq_len(nrow(subgroup_primary))) {
  r <- subgroup_primary[i, ]
  cat(sprintf("  %s: n=%d events=%d AUC=%.3f (95%% CI %.3f-%.3f)\n",
              r$group, r$n, r$events, r$auc, r$lci, r$uci))
}

cat("\n--- 9J.2 SECONDARY (major_tx) subgroup AUCs ---\n")
subgroup_secondary <- do.call(rbind, lapply(c("under70", "ge70"), function(g) {
  idx <- which(.grp_ed == g)
  r <- tryCatch(
    roc(df_ed$major_tx[idx], oob_prob_mt[idx],
        levels = c("NO","YES"), direction = "<", quiet = TRUE),
    error = function(e) NULL)
  if (is.null(r)) return(data.frame(group = g, n = length(idx),
                                    events = NA, auc = NA, lci = NA, uci = NA))
  ci <- pROC::ci.auc(r, method = "delong")
  data.frame(group = g, n = length(idx),
             events = sum(df_ed$major_tx[idx] == "YES"),
             auc = as.numeric(auc(r)),
             lci = as.numeric(ci[1]),
             uci = as.numeric(ci[3]))
}))
for (i in seq_len(nrow(subgroup_secondary))) {
  r <- subgroup_secondary[i, ]
  cat(sprintf("  MT %s: n=%d events=%d AUC=%.3f (95%% CI %.3f-%.3f)\n",
              r$group, r$n, r$events, r$auc, r$lci, r$uci))
}

saveRDS(list(primary   = subgroup_primary,
             secondary = subgroup_secondary,
             age_cut   = 70L,
             method    = "DeLong CIs on RF OOB probabilities"),
        "age_subgroup_aucs.rds")
cat("\nSaved: age_subgroup_aucs.rds\n")

cat("\nSubgroup AUCs (primary):\n"); print(subgroup_primary,  row.names = FALSE)
cat("\nSubgroup AUCs (secondary):\n"); print(subgroup_secondary, row.names = FALSE)

# SECTION 10: HEAD-TO-HEAD VS TRADITIONAL SCORES
cat("============================================================\n")
cat("SECTION 10: HEAD-TO-HEAD COMPARISON\n")
cat("============================================================\n\n")

score_names <- c("GBS","AIMS65","Pre-Rockall","Full Rockall","CANUKA")
score_cols  <- c(COL_GREY, COL_RED, COL_GREEN, COL_AMB, "#8E44AD")

get_score_rocs <- function(outcome_var) {
  by_oc <- roc_sc$by_outcome
  if (is.null(by_oc) || !outcome_var %in% names(by_oc)) {
    return(list(roc_sc$gbs, roc_sc$aims65,
                roc_sc$pre_rockall, roc_sc$full_rockall, roc_sc$canuka))
  }
  list(by_oc[[outcome_var]][["GBS"]],
       by_oc[[outcome_var]][["AIMS65"]],
       by_oc[[outcome_var]][["Pre-Rockall"]],
       by_oc[[outcome_var]][["Full Rockall"]],
       by_oc[[outcome_var]][["CANUKA"]])
}

norm_score <- function(x) { r <- range(x, na.rm=TRUE); (x - r[1]) / (r[2] - r[1]) }

head_to_head <- function(outcome_var, oob_prob, roc_ml, model_label, df_in) {
  cat(sprintf("\n=========================================================\n"))
  cat(sprintf("HEAD-TO-HEAD ON %s [%s]\n", outcome_var,
              if (outcome_var == PRIMARY) "PRIMARY" else "SECONDARY"))
  cat(sprintf("=========================================================\n"))
  
  roc_list_sc <- get_score_rocs(outcome_var)
  
  cat(sprintf("\n--- 10A. DeLong tests: %s vs each score (outcome: %s) ---\n",
              model_label, outcome_var))
  cat(sprintf("ML model AUC = %.3f\n", auc(roc_ml)))
  delong_df <- data.frame(Score=character(), AUC_ML=numeric(),
                          AUC_Score=numeric(), Delta=numeric(), p_value=numeric(),
                          Sig=character(), stringsAsFactors=FALSE)
  for (i in seq_along(score_names)) {
    roc_sc_i <- roc_list_sc[[i]]
    if (is.null(roc_sc_i)) next
    method <- if (length(roc_ml$response) == length(roc_sc_i$response))
      "delong" else "bootstrap"
    test <- tryCatch(
      roc.test(roc_ml, roc_sc_i, method=method,
               boot.n=ifelse(method=="bootstrap", 1000, NULL)),
      error=function(e) NULL)
    if (is.null(test)) next
    p  <- round(test$p.value, 4)
    sg <- ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*","ns")))
    cat(sprintf("  ML(%.3f) vs %-14s(%.3f): delta=%+.3f p=%.4f %s\n",
                auc(roc_ml), score_names[i], auc(roc_sc_i),
                auc(roc_ml)-auc(roc_sc_i), p, sg))
    delong_df <- rbind(delong_df, data.frame(
      Score=score_names[i], AUC_ML=round(auc(roc_ml),3),
      AUC_Score=round(auc(roc_sc_i),3),
      Delta=round(auc(roc_ml)-auc(roc_sc_i),3),
      p_value=p, Sig=sg, stringsAsFactors=FALSE))
  }
  
  outcome_bin <- as.numeric(df_in[[outcome_var]] == "YES")
  score_probs <- list(
    GBS         = norm_score(df_in$gbs_score),
    AIMS65      = norm_score(df_in$aims65_score),
    Pre_Rockall = norm_score(df_in$pre_rockall),
    Rockall     = norm_score(df_in$full_rockall),
    CANUKA      = norm_score(df_in$canuka))
  
  cat(sprintf("\n--- 10B. Categorical NRI vs each score (outcome: %s) ---\n", outcome_var))
  cat("    Thresholds: 0.30 (conservative admit), 0.50 (balanced), 0.70 (escalation)\n")
  cat("    Per Pencina & Steyerberg Stat Med 2014;33:3415; Hmisc::improveProb\n")
  
  nri_thresholds <- c(0.30, 0.50, 0.70)
  nri_results <- list()
  cat(sprintf("\n  %-12s %s\n",
              "Comparator",
              "Threshold | NRI(events) NRI(non-events) NRI(total)"))
  cat(sprintf("  %s\n", paste(rep("-", 70), collapse = "")))
  for (sc in names(score_probs)) {
    p_new <- oob_prob; p_ref <- score_probs[[sc]]
    ok <- !is.na(p_new) & !is.na(p_ref) & !is.na(outcome_bin)
    y <- outcome_bin[ok]; pn <- p_new[ok]; pr <- p_ref[ok]
    nri_results[[sc]] <- list()
    for (t in nri_thresholds) {
      cls_new <- as.integer(pn >= t)
      cls_ref <- as.integer(pr >= t)
      ev  <- y == 1; nev <- y == 0
      up_ev   <- sum(cls_new == 1 & cls_ref == 0 & ev)
      down_ev <- sum(cls_new == 0 & cls_ref == 1 & ev)
      n_ev <- sum(ev)
      nri_ev <- if (n_ev > 0) (up_ev - down_ev) / n_ev else NA_real_
      down_nev <- sum(cls_new == 0 & cls_ref == 1 & nev)
      up_nev   <- sum(cls_new == 1 & cls_ref == 0 & nev)
      n_nev <- sum(nev)
      nri_nev <- if (n_nev > 0) (down_nev - up_nev) / n_nev else NA_real_
      nri_tot <- nri_ev + nri_nev
      nri_results[[sc]][[as.character(t)]] <- list(
        threshold = t,
        events    = round(nri_ev, 3),
        nonevents = round(nri_nev, 3),
        total     = round(nri_tot, 3),
        up_ev = up_ev, down_ev = down_ev,
        up_nev = up_nev, down_nev = down_nev,
        n_ev = n_ev, n_nev = n_nev)
      cat(sprintf("  %-12s    %.2f    | %+.3f       %+.3f          %+.3f\n",
                  sc, t, nri_ev, nri_nev, nri_tot))
    }
  }
  cat("\n  Note: continuous NRI omitted in V2.4.6 per Pencina-Steyerberg 2014 critique.\n")
  
  cat(sprintf("\n--- 10C. IDI vs each score (outcome: %s) ---\n", outcome_var))
  idi_results <- list()
  for (sc in names(score_probs)) {
    p_new <- oob_prob; p_ref <- score_probs[[sc]]
    ok <- !is.na(p_new) & !is.na(p_ref) & !is.na(outcome_bin)
    y <- outcome_bin[ok]; pn <- p_new[ok]; pr <- p_ref[ok]
    is_new <- mean(pn[y==1]) - mean(pn[y==0])
    is_ref <- mean(pr[y==1]) - mean(pr[y==0])
    idi <- round(is_new - is_ref, 4)
    idi_results[[sc]] <- idi
    cat(sprintf("  IDI vs %-8s: IS_new=%.4f IS_ref=%.4f IDI=%+.4f %s\n",
                sc, is_new, is_ref, idi,
                ifelse(idi > 0, "(ML better)", "(score better)")))
  }
  
  cat(sprintf("\n--- 10D. DCA (outcome: %s) ---\n", outcome_var))
  dca_fn <- function(prob, y, thr) {
    ok <- !is.na(prob) & !is.na(y); nn <- sum(ok); pp <- prob[ok]; yy <- y[ok]
    sapply(thr, function(pt) {
      tp <- sum(pp >= pt & yy == 1); fp <- sum(pp >= pt & yy == 0)
      tp/nn - (pt/(1 - pt)) * fp/nn })
  }
  thr  <- seq(0.05, 0.60, 0.01)
  prev <- mean(outcome_bin, na.rm=TRUE)
  dca_df <- bind_rows(
    data.frame(x=thr, nb=dca_fn(oob_prob, outcome_bin, thr), Model="BURST"),
    data.frame(x=thr, nb=dca_fn(score_probs$GBS,         outcome_bin, thr), Model="GBS"),
    data.frame(x=thr, nb=dca_fn(score_probs$AIMS65,      outcome_bin, thr), Model="AIMS65"),
    data.frame(x=thr, nb=dca_fn(score_probs$Pre_Rockall, outcome_bin, thr), Model="Pre-Rockall"),
    data.frame(x=thr, nb=dca_fn(score_probs$Rockall,     outcome_bin, thr), Model="Full Rockall"),
    data.frame(x=thr, nb=dca_fn(score_probs$CANUKA,      outcome_bin, thr), Model="CANUKA"),
    data.frame(x=thr, nb=sapply(thr, function(pt) prev - (pt/(1-pt)) * (1-prev)),
               Model="Treat all"))
  
  plot_key <- paste0("p_dca_", outcome_var)
  dca_df$Model <- factor(dca_df$Model,
                         levels = c("GBS","AIMS65","Pre-Rockall","Full Rockall",
                                    "CANUKA","Treat all","BURST"))
  .nb_max <- max(dca_df$nb[is.finite(dca_df$nb)], na.rm = TRUE)
  .nb_min <- min(dca_df$nb[is.finite(dca_df$nb)], na.rm = TRUE)
  .ylim_lo <- min(-0.05, .nb_min)
  .ylim_hi <- max(0.45,  .nb_max * 1.05)
  plots[[plot_key]] <<- ggplot(dca_df,
                               aes(x=x*100, y=nb, colour=Model, linetype=Model)) +
    geom_line(linewidth=1) +
    geom_hline(yintercept=0, linetype="dotted", colour="grey50") +
    scale_colour_manual(values=c("BURST"=COL_BLUE, "GBS"=COL_GREY,
                                 "AIMS65"=COL_RED, "Pre-Rockall"=COL_GREEN, "Full Rockall"=COL_AMB,
                                 "CANUKA"="#8E44AD", "Treat all"="black")) +
    scale_linetype_manual(values=c("BURST"="solid", "GBS"="dashed",
                                   "AIMS65"="dashed", "Pre-Rockall"="dashed", "Full Rockall"="dashed",
                                   "CANUKA"="dashed", "Treat all"="dotted")) +
    coord_cartesian(ylim = c(.ylim_lo, .ylim_hi)) +
    labs(title=ifelse(outcome_var == PRIMARY,
                      "Decision curve analysis on primary outcome",
                      "Decision curve analysis on major transfusion"),
         subtitle="Net benefit vs threshold | model above 'treat all' = clinical utility",
         x="Threshold probability (%)", y="Net benefit",
         colour=NULL, linetype=NULL) +
    theme_ugib()
  print(plots[[plot_key]])
  
  leg_lab <- sprintf("%s (AUC=%.3f)", model_label, auc(roc_ml))
  leg_col <- COL_BLUE; leg_lwd <- 3
  .sc_auc <- sapply(roc_list_sc, function(r) if (is.null(r)) NA_real_ else as.numeric(auc(r)))
  for (i in order(.sc_auc, decreasing = TRUE))
    if (!is.null(roc_list_sc[[i]])) {
      leg_lab <- c(leg_lab, sprintf("%-12s  AUC=%.3f",
                                    score_names[i], auc(roc_list_sc[[i]])))
      leg_col <- c(leg_col, score_cols[i]); leg_lwd <- c(leg_lwd, 2)
    }
  par(mar=c(5, 5, 4, 2))
  plot(roc_ml, col=COL_BLUE, lwd=3,
       main = ifelse(outcome_var == PRIMARY,
                     "ROC: BURST vs traditional scores on primary outcome",
                     "ROC: BURST vs traditional scores on major transfusion"),
       legacy.axes = TRUE,
       xlim = c(1, 0), ylim = c(0, 1), asp = NA,
       xaxs = "i", yaxs = "i")
  for (i in seq_along(roc_list_sc))
    if (!is.null(roc_list_sc[[i]]))
      plot(roc_list_sc[[i]], col=score_cols[i], lwd=2, legacy.axes = TRUE, add=TRUE)
  legend("bottomright", bty="n", lwd=leg_lwd, col=leg_col, legend=leg_lab)
  
  list(delong=delong_df, nri=nri_results, idi=idi_results,
       dca_df=dca_df, roc_ml=roc_ml, roc_list_sc=roc_list_sc)
}

h2h_primary   <- head_to_head(PRIMARY,   oob_prob_hr, roc_ml_hr,
                              "BURST", df_full)
h2h_secondary <- head_to_head(SECONDARY, oob_prob_mt, roc_ml_mt,
                              "BURST",  df_full)

delong_df    <- h2h_primary$delong
nri_results  <- h2h_primary$nri
idi_results  <- h2h_primary$idi
dca_df       <- h2h_primary$dca_df
roc_list_sc  <- h2h_primary$roc_list_sc

cat("\n--- 10E. Indirect external validation (literature benchmarking) ---\n")
cat("    Riley 2016 BMJ 353:i3140; Debray 2015 J Clin Epidemiol 68:279\n\n")

score_auc <- function(score_vec, outcome_vec) {
  ok <- !is.na(score_vec) & !is.na(outcome_vec)
  if (length(unique(outcome_vec[ok])) < 2) return(NA_real_)
  as.numeric(auc(roc(outcome_vec[ok], score_vec[ok],
                     levels=levels(factor(outcome_vec[ok])),
                     direction="<", quiet=TRUE)))
}

gbs_v <- as.numeric(df_full$gbs_score)
aim_v <- as.numeric(df_full$aims65_score)
pre_v <- as.numeric(df_full$pre_rockall)
ful_v <- as.numeric(df_full$full_rockall)
can_v <- as.numeric(df_full$canuka)

ext_val <- data.frame(
  Score   = c("GBS","AIMS65","Full Rockall","CANUKA"),
  Outcome = c("Intervention","Mortality","Mortality","Composite"),
  Cohort  = c("Stanley 2017 (BMJ, n=3,012)",
              "Saltzman 2011 (GIE, n=29,222)",
              "Vreeburg 1999 (Gut, n=951)",
              "Oakland 2019 (CGH, n=10,639)"),
  Lit_AUC = c(0.86, 0.77, 0.81, 0.77),
  Our_AUC = c(score_auc(gbs_v, df_full$high_risk),
              score_auc(aim_v, df_full$outcome_death),
              score_auc(ful_v, df_full$outcome_death),
              score_auc(can_v, df_full$high_risk)),
  stringsAsFactors = FALSE)
ext_val$Delta       <- round(ext_val$Our_AUC - ext_val$Lit_AUC, 3)
ext_val$Within_pm05 <- abs(ext_val$Delta) < 0.05

cat(sprintf("%-13s %-13s %-32s  %7s  %7s  %7s  %s\n",
            "Score","Outcome","Cohort","Lit AUC","Our AUC","Delta","|D|<.05?"))
cat(paste(rep("-", 100), collapse=""), "\n")
for (i in seq_len(nrow(ext_val))) {
  cat(sprintf("%-13s %-13s %-32s  %7.2f  %7.3f  %+7.3f  %s\n",
              ext_val$Score[i], ext_val$Outcome[i], ext_val$Cohort[i],
              ext_val$Lit_AUC[i], ext_val$Our_AUC[i],
              ext_val$Delta[i],
              ifelse(ext_val$Within_pm05[i], "yes", "no")))
}
cat(sprintf("\n%d / %d audit-verified benchmarks within +-0.05 of published AUCs.\n",
            sum(ext_val$Within_pm05, na.rm=TRUE), nrow(ext_val)))
if (any(!ext_val$Within_pm05 & ext_val$Delta > 0, na.rm=TRUE)) {
  out_idx <- which(!ext_val$Within_pm05 & ext_val$Delta > 0)
  cat(sprintf("Outlier(s) in better-than-published direction: %s\n",
              paste(sprintf("%s/%s", ext_val$Score[out_idx],
                            ext_val$Outcome[out_idx]),
                    collapse=", ")))
}

ext_val$row_label <- factor(
  sprintf("%s -- %s\n%s", ext_val$Score, ext_val$Outcome, ext_val$Cohort),
  levels = rev(sprintf("%s -- %s\n%s", ext_val$Score, ext_val$Outcome, ext_val$Cohort)))
ev_long <- rbind(
  data.frame(row_label=ext_val$row_label, AUC=ext_val$Lit_AUC,
             source="Published", stringsAsFactors=FALSE),
  data.frame(row_label=ext_val$row_label, AUC=ext_val$Our_AUC,
             source=sprintf("Our cohort (n=%d)", nrow(df_full)),
             stringsAsFactors=FALSE))
plots[["p_extval"]] <- ggplot(ev_long,
                              aes(x=AUC, y=row_label, colour=source)) +
  geom_vline(xintercept = c(0.7,0.8,0.9), colour="grey90", linewidth=.4) +
  geom_point(position=position_dodge(width=0.45), size=3.2) +
  scale_colour_manual(values=c("Published"=COL_GREY,
                               setNames(COL_BLUE,
                                        sprintf("Our cohort (n=%d)", nrow(df_full))))) +
  xlim(0.5, 1.0) +
  labs(title="Indirect external validation: our AUCs vs published cohorts",
       subtitle="Per Riley 2016 BMJ 353:i3140 framework for indirect external validation",
       x="AUC", y=NULL, colour=NULL) +
  theme_ugib() + theme(panel.grid.major.y=element_blank(),
                       axis.text.y=element_text(lineheight=1.0, size=8))
print(plots[["p_extval"]])

boot_p_paired <- function(y, s_ml, s_sc, B=1000, seed=2024) {
  ok <- !is.na(s_ml) & !is.na(s_sc) & !is.na(y)
  y <- y[ok]; s_ml <- s_ml[ok]; s_sc <- s_sc[ok]
  if (length(unique(y)) < 2) return(NA_real_)
  set.seed(seed); diffs <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample.int(length(y), length(y), replace=TRUE)
    if (length(unique(y[idx])) < 2) { diffs[b] <- NA_real_; next }
    a1 <- as.numeric(auc(roc(y[idx], s_ml[idx], quiet=TRUE, direction="<")))
    a2 <- as.numeric(auc(roc(y[idx], s_sc[idx], quiet=TRUE, direction="<")))
    diffs[b] <- a1 - a2
  }
  diffs <- diffs[!is.na(diffs)]
  if (length(diffs) < 2 || sd(diffs) == 0) return(NA_real_)
  obs <- as.numeric(auc(roc(y, s_ml, quiet=TRUE, direction="<"))) -
    as.numeric(auc(roc(y, s_sc, quiet=TRUE, direction="<")))
  2 * (1 - pnorm(abs(obs / sd(diffs))))
}
sc_aucs <- sapply(roc_list_sc, function(r) if (is.null(r)) NA else as.numeric(auc(r)))
best_sc_idx <- which.max(sc_aucs)
best_sc_name <- score_names[best_sc_idx]
best_sc_roc  <- roc_list_sc[[best_sc_idx]]
y_num <- as.numeric(df_full[[PRIMARY]] == "YES")
ml_probs_for_pair <- if (length(oob_prob_hr) == length(y_num)) oob_prob_hr else
  rep(NA_real_, length(y_num))
best_sc_col <- switch(best_sc_name,
                      "GBS"          = df_full$gbs_score,
                      "AIMS65"       = df_full$aims65_score,
                      "Pre-Rockall"  = df_full$pre_rockall,
                      "Full Rockall" = df_full$full_rockall,
                      rep(NA_real_, nrow(df_full)))
sc_probs_for_pair <- as.numeric(best_sc_col)
p_pair <- boot_p_paired(y_num, ml_probs_for_pair, sc_probs_for_pair)
cat(sprintf("\nPaired bootstrap p (ML vs best score = %s) on %s: p = %s\n",
            best_sc_name, PRIMARY,
            ifelse(is.na(p_pair), "NA", sprintf("%.3f", p_pair))))

saveRDS(list(delong=delong_df, nri=nri_results, idi=idi_results,
             dca_df=dca_df, roc_ml=roc_ml_hr,
             h2h_primary=h2h_primary, h2h_secondary=h2h_secondary,
             ext_val=ext_val, paired_p=p_pair), "headtohead_results.rds")

cat("\n============================================================\n")
cat("SECTION 10F: XGBoost SENSITIVITY COMPARATOR\n")
cat("(manual stratified 5-fold CV; xgboost-3.x-safe)\n")
cat("============================================================\n\n")

xgb_block_status <- tryCatch({
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    cat("xgboost not installed -- skipping §10F\n")
    "SKIPPED"
  } else {
    X_hr_df <- df_hr[, ed_feats, drop = FALSE]
    for (v in ed_feats) {
      if (is.numeric(X_hr_df[[v]]) && any(is.na(X_hr_df[[v]]))) {
        X_hr_df[[v]][is.na(X_hr_df[[v]])] <- median(X_hr_df[[v]], na.rm = TRUE)
      }
    }
    X_hr <- model.matrix(~ . - 1, data = X_hr_df)
    y_hr <- as.numeric(df_hr[[PRIMARY]] == "YES")
    
    cat(sprintf("XGBoost input: n=%d, predictors after dummy coding=%d, events=%d\n",
                nrow(X_hr), ncol(X_hr), sum(y_hr)))
    
    xgb_params <- list(
      objective        = "binary:logistic",
      eval_metric      = "auc",
      eta              = 0.05,
      max_depth        = 4,
      min_child_weight = 5,
      subsample        = 0.7,
      colsample_bytree = 0.7)
    XGB_NROUNDS <- 60L
    
    set.seed(2024)
    nfold <- 5L
    idx_pos <- which(y_hr == 1)
    idx_neg <- which(y_hr == 0)
    fold_pos <- sample(rep(seq_len(nfold), length.out = length(idx_pos)))
    fold_neg <- sample(rep(seq_len(nfold), length.out = length(idx_neg)))
    fold_id <- integer(length(y_hr))
    fold_id[idx_pos] <- fold_pos
    fold_id[idx_neg] <- fold_neg
    
    oof_pred  <- rep(NA_real_, length(y_hr))
    fold_aucs <- rep(NA_real_, nfold)
    fold_ok   <- logical(nfold)
    
    cat(sprintf("Running stratified %d-fold CV (manual loop, fixed %d rounds):\n",
                nfold, XGB_NROUNDS))
    for (k in seq_len(nfold)) {
      fold_result <- tryCatch({
        train_rows <- which(fold_id != k)
        test_rows  <- which(fold_id == k)
        
        if (length(test_rows) < 1L || length(train_rows) < 1L) {
          stop("empty fold")
        }
        
        dtrain_k <- xgboost::xgb.DMatrix(data  = X_hr[train_rows, , drop = FALSE],
                                         label = y_hr[train_rows])
        dtest_k  <- xgboost::xgb.DMatrix(data  = X_hr[test_rows, , drop = FALSE],
                                         label = y_hr[test_rows])
        
        fit_k <- xgboost::xgb.train(
          params  = xgb_params,
          data    = dtrain_k,
          nrounds = XGB_NROUNDS,
          evals   = list(train = dtrain_k, test = dtest_k),
          verbose = 0)
        
        pred_k <- as.numeric(predict(fit_k,
                                     newdata = X_hr[test_rows, , drop = FALSE]))
        
        if (length(pred_k) != length(test_rows) || all(is.na(pred_k))) {
          stop(sprintf("predict() returned %d values for %d test rows",
                       length(pred_k), length(test_rows)))
        }
        
        fold_auc_k <- tryCatch(
          as.numeric(pROC::auc(pROC::roc(y_hr[test_rows], pred_k, quiet = TRUE))),
          error = function(e) NA_real_)
        
        list(ok = TRUE, pred = pred_k, test_rows = test_rows,
             auc = fold_auc_k, n_train = length(train_rows),
             n_test = length(test_rows), n_events = sum(y_hr[test_rows]))
      }, error = function(e) {
        cat(sprintf("  Fold %d FAILED: %s\n", k, conditionMessage(e)))
        list(ok = FALSE)
      })
      
      if (isTRUE(fold_result$ok)) {
        oof_pred[fold_result$test_rows] <- fold_result$pred
        fold_aucs[k] <- fold_result$auc
        fold_ok[k]   <- TRUE
        cat(sprintf("  Fold %d: n_train=%d n_test=%d events_test=%d AUC=%.3f\n",
                    k, fold_result$n_train, fold_result$n_test,
                    fold_result$n_events, fold_result$auc))
      }
    }
    
    n_ok   <- sum(fold_ok)
    n_pred <- sum(!is.na(oof_pred))
    cat(sprintf("\n%d/%d folds succeeded; %d patients have OOF predictions\n",
                n_ok, nfold, n_pred))
    
    if (n_pred < 0.5 * length(y_hr)) {
      cat("Insufficient OOF predictions - skipping AUC / DeLong reporting\n")
      saveRDS(list(
        method        = "manual_stratified_5fold_cv_v248",
        status        = "INSUFFICIENT_OOF",
        fold_ok       = fold_ok,
        fold_aucs     = fold_aucs,
        n_pred        = n_pred,
        params        = xgb_params,
        nrounds       = XGB_NROUNDS),
        "xgboost_sensitivity.rds")
      cat("\nSaved: xgboost_sensitivity.rds (insufficient OOF)\n")
      "PARTIAL"
    } else {
      ok_idx  <- which(!is.na(oof_pred))
      roc_xgb <- pROC::roc(y_hr[ok_idx], oof_pred[ok_idx],
                           levels = c(0, 1), direction = "<", quiet = TRUE)
      auc_xgb <- as.numeric(pROC::auc(roc_xgb))
      ci_xgb  <- as.numeric(pROC::ci.auc(roc_xgb))
      auc_rf  <- as.numeric(pROC::auc(roc_ml_hr))
      ci_rf   <- as.numeric(pROC::ci.auc(roc_ml_hr))
      
      cat(sprintf("\nXGBoost (manual 5-fold CV, %d rounds):\n", XGB_NROUNDS))
      cat(sprintf("  RF (OOB):       AUC = %.4f (95%% CI %.4f-%.4f)\n",
                  auc_rf, ci_rf[1], ci_rf[3]))
      cat(sprintf("  XGBoost (5xCV): AUC = %.4f (95%% CI %.4f-%.4f)\n",
                  auc_xgb, ci_xgb[1], ci_xgb[3]))
      cat(sprintf("  Delta (RF - XGB): %+.4f\n", auc_rf - auc_xgb))
      
      use_paired <- (n_pred == length(y_hr))
      .roc_rf_for_xgb <- pROC::roc(y_hr[ok_idx], oob_prob_hr[ok_idx],
                                   levels = c(0, 1), direction = "<", quiet = TRUE)
      delong_xgb <- tryCatch(
        pROC::roc.test(.roc_rf_for_xgb, roc_xgb, method = "delong", paired = use_paired),
        error = function(e) tryCatch(
          pROC::roc.test(.roc_rf_for_xgb, roc_xgb, method = "delong", paired = FALSE),
          error = function(e2) tryCatch(
            pROC::roc.test(.roc_rf_for_xgb, roc_xgb, method = "bootstrap",
                           boot.n = 1000, paired = FALSE),
            error = function(e3) NULL)))
      
      if (!is.null(delong_xgb)) {
        cat(sprintf("  %s DeLong RF vs XGBoost: p = %.4f\n",
                    ifelse(use_paired, "Paired", "Unpaired"),
                    delong_xgb$p.value))
        cat(sprintf("  Conclusion: %s\n",
                    ifelse(delong_xgb$p.value > 0.05,
                           "EQUIVALENT (RF and XGBoost statistically indistinguishable at n=133)",
                           "DIFFERENT")))
      } else {
        cat("  DeLong test failed to compute\n")
      }
      
      saveRDS(list(
        method         = "manual_stratified_5fold_cv_v248",
        status         = "OK",
        auc_rf         = auc_rf,
        auc_xgb        = auc_xgb,
        ci_xgb         = ci_xgb,
        fold_aucs      = fold_aucs,
        fold_ok        = fold_ok,
        n_pred         = n_pred,
        delong_p       = if (!is.null(delong_xgb)) delong_xgb$p.value else NA_real_,
        delong_paired  = use_paired,
        params         = xgb_params,
        nrounds        = XGB_NROUNDS,
        oof_pred       = oof_pred,
        y_hr           = y_hr),
        "xgboost_sensitivity.rds")
      cat("\nSaved: xgboost_sensitivity.rds (cross-validated, no in-sample fallback)\n")
      "OK"
    }
  }
}, error = function(e) {
  cat(sprintf("\n§10F XGBoost: error %s\n",
              conditionMessage(e)))
  
  "FAILED"
})
cat(sprintf("\n§10F status: %s\n", xgb_block_status))

cat("\n============================================================\n")
cat("SECTION 10G: ri_pct EXCLUSION SENSITIVITY\n")
cat("============================================================\n\n")

ed_feats_5 <- setdiff(ed_feats, "ri_pct")
cat(sprintf("6-feature headline: %s\n",
            paste(ed_feats, collapse = ", ")))
cat(sprintf("5-feature subset:   %s\n\n",
            paste(ed_feats_5, collapse = ", ")))

df_hr_5 <- df_hr |>
  dplyr::select(all_of(ed_feats_5), all_of(PRIMARY))
df_hr_5[[PRIMARY]] <- factor(df_hr_5[[PRIMARY]], levels = c("NO", "YES"))

hr_min_5 <- min(table(df_hr_5[[PRIMARY]]))
set.seed(2024)
rf_hr_5 <- randomForest(
  reformulate(ed_feats_5, PRIMARY),
  data     = df_hr_5,
  ntree    = 500,
  strata   = df_hr_5[[PRIMARY]],
  sampsize = c("NO" = hr_min_5, "YES" = hr_min_5))

oob_p_5 <- rf_hr_5$votes[, "YES"] / rowSums(rf_hr_5$votes)
roc_5   <- roc(df_hr_5[[PRIMARY]], oob_p_5,
               levels = c("NO", "YES"), direction = "<", quiet = TRUE)
auc_5   <- as.numeric(auc(roc_5))
ci_5    <- ci.auc(roc_5)

auc_6 <- as.numeric(auc(roc_ml_hr))

delong_5v6 <- tryCatch(
  roc.test(roc_ml_hr, roc_5, method = "delong", paired = TRUE),
  error = function(e) NULL)

cat("Discrimination comparison (PRIMARY: high_risk):\n")
cat(sprintf("  6-feature headline (with ri_pct):    AUC = %.4f (95%% CI %.4f-%.4f)\n",
            auc_6, ci.auc(roc_ml_hr)[1], ci.auc(roc_ml_hr)[3]))
cat(sprintf("  5-feature subset   (without ri_pct): AUC = %.4f (95%% CI %.4f-%.4f)\n",
            auc_5, ci_5[1], ci_5[3]))
cat(sprintf("  Delta AUC (6-feat - 5-feat):         %+.4f\n", auc_6 - auc_5))

if (!is.null(delong_5v6)) {
  cat(sprintf("  Paired DeLong test:                  p = %.4f\n",
              delong_5v6$p.value))
  cat(sprintf("  Conclusion: %s\n",
              ifelse(delong_5v6$p.value > 0.05,
                     "EQUIVALENT discrimination (drop ri_pct safely; supports 5-feature alt)",
                     "DIFFERENT discrimination (retain ri_pct)")))
}

if (exists("cv_binary")) {
  cv_5 <- cv_binary(df_hr_5, ed_feats_5, PRIMARY, model_type = "RF")
  cat(sprintf("\n  5-feature RF 5x5 CV AUC: mean = %.4f (sd = %.4f)\n",
              mean(cv_5$aucs, na.rm = TRUE),
              sd(cv_5$aucs, na.rm = TRUE)))
  cat(sprintf("    per-repeat: %s\n",
              paste(sprintf("%.3f", cv_5$aucs), collapse = ", ")))
  cat(sprintf("    (V2.4.5 6-feature headline 5x5 CV AUC: 0.914 sd 0.008)\n"))
}

saveRDS(list(
  ed_feats_6 = ed_feats,
  ed_feats_5 = ed_feats_5,
  auc_6      = auc_6,
  auc_5      = auc_5,
  ci_5       = ci_5,
  delong_p   = if (!is.null(delong_5v6)) delong_5v6$p.value else NA_real_,
  cv_5_mean  = if (exists("cv_5")) mean(cv_5$aucs, na.rm = TRUE) else NA_real_,
  cv_5_sd    = if (exists("cv_5")) sd(cv_5$aucs, na.rm = TRUE) else NA_real_),
  "ri_pct_sensitivity.rds")
cat("\nSaved: ri_pct_sensitivity.rds\n")

cat("\n============================================================\n")
cat("SECTION 10H: Full 5-feature evaluation\n")
cat("============================================================\n\n")

cat("--- 10H.1. 5-feature calibration ---\n")

y_pri_5 <- as.numeric(df_hr_5[[PRIMARY]] == "YES")
hl_decile_breaks_5 <- unique(quantile(oob_p_5, probs = seq(0, 1, 0.1), na.rm = TRUE))
if (length(hl_decile_breaks_5) >= 3) {
  hl_groups_5 <- cut(oob_p_5, breaks = hl_decile_breaks_5,
                     include.lowest = TRUE, labels = FALSE)
  hl_obs_5 <- tapply(y_pri_5, hl_groups_5, sum, na.rm = TRUE)
  hl_exp_5 <- tapply(oob_p_5, hl_groups_5, sum, na.rm = TRUE)
  hl_n_5   <- tapply(y_pri_5, hl_groups_5, length)
  hl_chi2_5 <- sum((hl_obs_5 - hl_exp_5)^2 /
                     (hl_exp_5 * (1 - hl_exp_5 / hl_n_5)), na.rm = TRUE)
  hl_df_5 <- length(hl_obs_5) - 2
  hl_p_5  <- if (hl_df_5 > 0) pchisq(hl_chi2_5, df = hl_df_5, lower.tail = FALSE) else NA_real_
} else {
  hl_chi2_5 <- NA_real_; hl_df_5 <- NA_integer_; hl_p_5 <- NA_real_
}

oob_p_5_clip <- pmin(pmax(oob_p_5, 1e-4), 1 - 1e-4)
logit_p_5 <- log(oob_p_5_clip / (1 - oob_p_5_clip))
calib_int_fit_5 <- glm(y_pri_5 ~ offset(logit_p_5), family = binomial)
calib_int_5 <- coef(calib_int_fit_5)[1]
calib_int_se_5 <- summary(calib_int_fit_5)$coefficients[1, 2]
calib_slope_fit_5 <- glm(y_pri_5 ~ logit_p_5, family = binomial)
calib_slope_5 <- coef(calib_slope_fit_5)[2]
calib_slope_ci_5 <- confint.default(calib_slope_fit_5)["logit_p_5", ]

brier_5 <- mean((oob_p_5 - y_pri_5)^2)
prev_5 <- mean(y_pri_5)
brier_ref_5 <- mean((prev_5 - y_pri_5)^2)
bss_5 <- 1 - brier_5 / brier_ref_5

cat(sprintf("  HL chi2 = %.2f (df = %d), p = %.4f\n", hl_chi2_5, hl_df_5, hl_p_5))
cat(sprintf("  Calibration intercept = %.3f (SE %.3f)\n", calib_int_5, calib_int_se_5))
cat(sprintf("  Calibration slope     = %.3f (95%% CI %.3f-%.3f)\n",
            calib_slope_5, calib_slope_ci_5[1], calib_slope_ci_5[2]))
cat(sprintf("  Brier = %.4f, BSS = %.3f\n", brier_5, bss_5))

cat("\n--- 10H.2. 5-feature threshold-sensitivity table ---\n")
thr_grid_5 <- c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80)
threshold_5_df <- do.call(rbind, lapply(thr_grid_5, function(t) {
  cls  <- as.integer(oob_p_5 >= t)
  tp <- sum(cls == 1 & y_pri_5 == 1)
  fp <- sum(cls == 1 & y_pri_5 == 0)
  fn <- sum(cls == 0 & y_pri_5 == 1)
  tn <- sum(cls == 0 & y_pri_5 == 0)
  sens <- if (tp + fn > 0) tp / (tp + fn) else NA_real_
  spec <- if (tn + fp > 0) tn / (tn + fp) else NA_real_
  ppv  <- if (tp + fp > 0) tp / (tp + fp) else NA_real_
  npv  <- if (tn + fn > 0) tn / (tn + fn) else NA_real_
  f1   <- if (!is.na(sens) && !is.na(ppv) && (sens + ppv) > 0)
    2 * sens * ppv / (sens + ppv) else NA_real_
  yj   <- if (!is.na(sens) && !is.na(spec)) sens + spec - 1 else NA_real_
  data.frame(threshold = t, sens = sens, spec = spec, ppv = ppv,
             npv = npv, f1 = f1, youden = yj)
}))
print(round(threshold_5_df, 3), row.names = FALSE)

ok_thr_5 <- threshold_5_df$sens >= 0.999
sens100_thr_5 <- if (any(ok_thr_5, na.rm = TRUE))
  max(threshold_5_df$threshold[ok_thr_5], na.rm = TRUE) else NA_real_
all_unique_p <- sort(unique(c(oob_p_5, 0)))
sens_at_p <- sapply(all_unique_p, function(t) {
  cls <- as.integer(oob_p_5 >= t)
  tp <- sum(cls == 1 & y_pri_5 == 1)
  fn <- sum(cls == 0 & y_pri_5 == 1)
  if (tp + fn > 0) tp / (tp + fn) else NA_real_
})
sens100_thr_5 <- max(all_unique_p[sens_at_p >= 0.999], na.rm = TRUE)
spec_at_sens100 <- {
  cls <- as.integer(oob_p_5 >= sens100_thr_5)
  tn <- sum(cls == 0 & y_pri_5 == 0)
  fp <- sum(cls == 1 & y_pri_5 == 0)
  if (tn + fp > 0) tn / (tn + fp) else NA_real_
}
cat(sprintf("\n  100%%-sensitivity threshold (5-feature): %.4f (specificity %.3f)\n",
            sens100_thr_5, spec_at_sens100))

cat("\n--- 10H.3. 5-feature head-to-head vs comparator scores (PRIMARY) ---\n")
roc_5feat <- roc(y_pri_5, oob_p_5,
                 levels = c(0, 1), direction = "<", quiet = TRUE)

norm_score_h <- function(x) {
  x <- as.numeric(x)
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) > 0) (x - rng[1]) / diff(rng) else rep(0.5, length(x))
}
score_probs_5 <- list(
  GBS         = norm_score_h(df_full$gbs_score),
  AIMS65      = norm_score_h(df_full$aims65_score),
  Pre_Rockall = norm_score_h(df_full$pre_rockall),
  Rockall     = norm_score_h(df_full$full_rockall),
  CANUKA      = norm_score_h(df_full$canuka))

cat("\n  DeLong tests (5-feature ML vs each score):\n")
delong_5_df <- data.frame(Score=character(), AUC_5=numeric(), AUC_score=numeric(),
                          Delta=numeric(), p_value=numeric(), stringsAsFactors=FALSE)
for (sc in names(score_probs_5)) {
  pr <- score_probs_5[[sc]][seq_along(y_pri_5)]
  ok <- !is.na(pr) & !is.na(y_pri_5)
  if (sum(ok) < 10) next
  .roc_sc_5feat <- roc(y_pri_5[ok], pr[ok], levels = c(0, 1), direction = "<", quiet = TRUE)
  rt <- tryCatch(roc.test(roc_5feat, .roc_sc_5feat, method = "delong", paired = TRUE),
                 error = function(e) NULL)
  auc_sc <- as.numeric(auc(.roc_sc_5feat))
  if (!is.null(rt)) {
    cat(sprintf("    ML5(%.3f) vs %-12s (%.3f): delta=%+.3f p=%.4f %s\n",
                auc_5, sc, auc_sc, auc_5 - auc_sc, rt$p.value,
                ifelse(rt$p.value < 0.05, "*", "ns")))
    delong_5_df <- rbind(delong_5_df, data.frame(
      Score = sc, AUC_5 = auc_5, AUC_score = auc_sc,
      Delta = auc_5 - auc_sc, p_value = rt$p.value))
  }
}

cat("\n  Categorical NRI (5-feature ML vs each score, primary outcome):\n")
nri_thresholds_5 <- c(0.30, 0.50, 0.70)
nri_5_results <- list()
cat(sprintf("    %-12s %s\n", "Comparator",
            "Threshold | NRI(events) NRI(non-events) NRI(total)"))
for (sc in names(score_probs_5)) {
  pr <- score_probs_5[[sc]][seq_along(y_pri_5)]
  ok <- !is.na(pr) & !is.na(y_pri_5)
  pn <- oob_p_5[ok]; pre <- pr[ok]; yo <- y_pri_5[ok]
  nri_5_results[[sc]] <- list()
  for (tt in nri_thresholds_5) {
    cls_n <- as.integer(pn  >= tt)
    cls_r <- as.integer(pre >= tt)
    ev <- yo == 1; nev <- yo == 0
    up_ev   <- sum(cls_n == 1 & cls_r == 0 & ev)
    down_ev <- sum(cls_n == 0 & cls_r == 1 & ev)
    nri_ev  <- if (sum(ev)  > 0) (up_ev   - down_ev) / sum(ev)  else NA_real_
    down_nev <- sum(cls_n == 0 & cls_r == 1 & nev)
    up_nev   <- sum(cls_n == 1 & cls_r == 0 & nev)
    nri_nev  <- if (sum(nev) > 0) (down_nev - up_nev) / sum(nev) else NA_real_
    nri_tot  <- nri_ev + nri_nev
    nri_5_results[[sc]][[as.character(tt)]] <- list(
      threshold = tt, events = round(nri_ev, 3),
      nonevents = round(nri_nev, 3), total = round(nri_tot, 3))
    cat(sprintf("    %-12s    %.2f    | %+.3f       %+.3f          %+.3f\n",
                sc, tt, nri_ev, nri_nev, nri_tot))
  }
}

cat("\n  IDI (5-feature ML vs each score, primary outcome):\n")
idi_5_results <- list()
for (sc in names(score_probs_5)) {
  pr <- score_probs_5[[sc]][seq_along(y_pri_5)]
  ok <- !is.na(pr) & !is.na(y_pri_5)
  pn <- oob_p_5[ok]; pre <- pr[ok]; yo <- y_pri_5[ok]
  is_n <- mean(pn[yo == 1]) - mean(pn[yo == 0])
  is_r <- mean(pre[yo == 1]) - mean(pre[yo == 0])
  idi_v <- round(is_n - is_r, 4)
  idi_5_results[[sc]] <- idi_v
  cat(sprintf("    IDI vs %-12s: IS_5feat=%.4f IS_score=%.4f IDI=%+.4f\n",
              sc, is_n, is_r, idi_v))
}

cat("\n--- 10H.4. 5-feature SHAP global importance ---\n")
if (requireNamespace("fastshap", quietly = TRUE)) {
  pred_wrapper_rf_5 <- function(object, newdata) {
    predict(object, newdata, type = "prob")[, "YES"]
  }
  X_shap_5 <- df_hr_5[, ed_feats_5, drop = FALSE]
  cat(sprintf("  Computing SHAP on 5-feature model: n=%d, features=%d, nsim=50...\n",
              nrow(X_shap_5), length(ed_feats_5)))
  set.seed(2024)
  shap_5 <- tryCatch(
    fastshap::explain(object = rf_hr_5, X = X_shap_5,
                      pred_wrapper = pred_wrapper_rf_5,
                      nsim = 50, parallel = FALSE),
    error = function(e) {cat(sprintf("  SHAP failed: %s\n", e$message)); NULL})
  if (!is.null(shap_5)) {
    shap_mat_5 <- as.matrix(shap_5)
    shap_global_5 <- data.frame(
      variable      = colnames(shap_mat_5),
      mean_abs_shap = colMeans(abs(shap_mat_5), na.rm = TRUE),
      sd_shap       = apply(shap_mat_5, 2, sd, na.rm = TRUE),
      stringsAsFactors = FALSE)
    shap_global_5 <- shap_global_5[order(-shap_global_5$mean_abs_shap), ]
    shap_global_5$importance_pct <- round(
      100 * shap_global_5$mean_abs_shap / max(shap_global_5$mean_abs_shap), 1)
    rownames(shap_global_5) <- NULL
    cat("  Global SHAP importance (5-feature):\n")
    print(shap_global_5, row.names = FALSE)
  }
} else {
  shap_global_5 <- NULL
  cat("  fastshap not installed -- skipping\n")
}

saveRDS(list(
  rf_hr_5            = rf_hr_5,
  ed_feats_5         = ed_feats_5,
  oob_p_5            = oob_p_5,
  roc_5feat          = roc_5feat,
  auc_5              = auc_5,
  ci_5               = ci_5,
  cv_5_mean          = if (exists("cv_5")) mean(cv_5$aucs) else NA_real_,
  cv_5_sd            = if (exists("cv_5")) sd(cv_5$aucs) else NA_real_,
  hl_chi2_5          = hl_chi2_5,
  hl_p_5             = hl_p_5,
  calib_int_5        = calib_int_5,
  calib_int_se_5     = calib_int_se_5,
  calib_slope_5      = calib_slope_5,
  calib_slope_ci_5   = calib_slope_ci_5,
  brier_5            = brier_5,
  bss_5              = bss_5,
  threshold_5_df     = threshold_5_df,
  sens100_thr_5      = sens100_thr_5,
  spec_at_sens100    = spec_at_sens100,
  delong_5_df        = delong_5_df,
  nri_5_results      = nri_5_results,
  idi_5_results      = idi_5_results,
  shap_global_5      = if (exists("shap_global_5")) shap_global_5 else NULL),
  "five_feature_full_evaluation.rds")
cat("\nSaved: five_feature_full_evaluation.rds\n")

# SECTION 11: EXPLAINABILITY
cat("\n============================================================\n")
cat("SECTION 11: EXPLAINABILITY\n")
cat("============================================================\n\n")

imp_mat <- importance(rf_hr)
imp_df  <- data.frame(variable=rownames(imp_mat),
                      gini=imp_mat[,"MeanDecreaseGini"],
                      acc =imp_mat[,if ("MeanDecreaseAccuracy" %in% colnames(imp_mat))
                        "MeanDecreaseAccuracy" else "MeanDecreaseGini"],
                      stringsAsFactors=FALSE) |>
  mutate(gini_pct=round(gini/max(gini)*100,1),
         acc_pct=round(acc/max(acc)*100,1)) |>
  arrange(desc(gini_pct))

cat(sprintf("Importance from PRIMARY RF (rf_hr, outcome=%s)\n", PRIMARY))
cat("Top 15 by Gini:\n")
print(imp_df[1:min(15,nrow(imp_df)), c("variable","gini_pct","acc_pct")])

top15_gini <- imp_df$variable[1:min(15,nrow(imp_df))]
top15_acc  <- imp_df |> arrange(desc(acc_pct)) |> pull(variable) |> head(15)
divergent  <- setdiff(top15_acc, top15_gini)
cat("\nHigh accuracy but lower Gini (rare-but-decisive):\n")
cat(paste(divergent, collapse=", "), "\n\n")

plots[["p_gini"]] <- ggplot(imp_df[1:min(20,nrow(imp_df)),] |>
                              mutate(variable_disp = relabel_var(variable)),
                            aes(x=gini_pct, y=reorder(variable_disp, gini_pct),
                                fill=variable %in% top15_acc)) +
  geom_col(width=0.75) +
  scale_fill_manual(values=c("TRUE"=COL_BLUE,"FALSE"="#7EB5D6"),
                    labels=c("TRUE"="Top-15 both metrics","FALSE"="Gini only"), name=NULL) +
  labs(title = "Variable importance (high risk BURST) -- Gini top 20",
       x="Gini importance (% of max)", y=NULL) + theme_ugib()
print(plots[["p_gini"]])

cat("\n--- Permutation importance (PRIMARY) -- V2.4.5 patched (Strobl 2007) ---\n")
cat("    Per Breiman 2001 / Strobl 2007. B=50 shuffles per variable.\n")
cat("    Method: shuffle, RE-FIT forest, use OOB votes (NOT predict on training).\n")
B_perm_imp <- 50
baseline_auc <- as.numeric(auc(roc_ml_hr))
cat(sprintf("    Baseline OOB AUC: %.4f | Total fits: %d x %d = %d (~3-5 min)\n\n",
            baseline_auc, length(ed_feats), B_perm_imp,
            length(ed_feats) * B_perm_imp))
perm_imp_results <- data.frame(variable=character(), drop_mean=numeric(),
                               drop_sd=numeric(), drop_pct=numeric(),
                               stringsAsFactors=FALSE)
t_perm_start <- Sys.time()
for (v in ed_feats) {
  drops <- numeric(B_perm_imp)
  for (b in seq_len(B_perm_imp)) {
    set.seed(2024 + b * 1000L + which(ed_feats == v))
    df_shuf <- df_hr
    df_shuf[[v]] <- sample(df_shuf[[v]])
    hr_min_b <- min(table(df_shuf$high_risk))
    rf_b <- tryCatch(
      randomForest(reformulate(ed_feats, "high_risk"),
                   data=df_shuf, ntree=500,
                   strata=df_shuf$high_risk,
                   sampsize=c("NO"=hr_min_b, "YES"=hr_min_b)),
      error=function(e) NULL)
    if (is.null(rf_b)) { drops[b] <- NA_real_; next }
    p_b <- rf_b$votes[, "YES"] / rowSums(rf_b$votes)
    r_b <- tryCatch(roc(df_shuf$high_risk, p_b,
                        levels=c("NO","YES"), direction="<", quiet=TRUE),
                    error=function(e) NULL)
    drops[b] <- if (!is.null(r_b)) baseline_auc - as.numeric(auc(r_b)) else NA_real_
  }
  perm_imp_results <- rbind(perm_imp_results, data.frame(
    variable  = v,
    drop_mean = round(mean(drops, na.rm=TRUE), 4),
    drop_sd   = round(sd(drops,   na.rm=TRUE), 4),
    drop_pct  = NA_real_,
    stringsAsFactors=FALSE))
  cat(sprintf("    %-15s drop=%+.4f (sd=%.4f)\n",
              v, mean(drops, na.rm=TRUE), sd(drops, na.rm=TRUE)))
}
cat(sprintf("\n    Elapsed: %.1f minutes\n",
            as.numeric(difftime(Sys.time(), t_perm_start, units="mins"))))
max_drop <- max(perm_imp_results$drop_mean, na.rm=TRUE)
if (max_drop > 0) {
  perm_imp_results$drop_pct <- round(perm_imp_results$drop_mean / max_drop * 100, 1)
}
perm_imp_results <- perm_imp_results[order(-perm_imp_results$drop_mean), ]
rownames(perm_imp_results) <- NULL
cat("Permutation importance (drop in OOB AUC when shuffled):\n")
print(perm_imp_results)

plots[["p_perm_imp"]] <- ggplot(perm_imp_results |>
                                  mutate(variable_disp = relabel_var(variable)),
                                aes(x=drop_mean, y=reorder(variable_disp, drop_mean))) +
  geom_col(fill=COL_BLUE, width=0.75) +
  geom_errorbar(aes(xmin=pmax(drop_mean - drop_sd, 0),
                    xmax=drop_mean + drop_sd),
                orientation="y", width=0.25, colour="grey40") +
  labs(title = "BURST permutation importance on primary outcome",
       subtitle=sprintf("Mean +/- SD over B=%d shuffles | drop in OOB AUC vs baseline %.3f",
                        B_perm_imp, baseline_auc),
       x="Mean drop in AUC when feature shuffled", y=NULL) +
  theme_ugib()
print(plots[["p_perm_imp"]])

saveRDS(perm_imp_results, "permutation_importance.rds")
cat("Saved: permutation_importance.rds\n\n")

num_vars_avail <- ed_feats[sapply(ed_feats, function(v) is.numeric(df_hr[[v]]))]
top_cont <- imp_df |> filter(variable %in% num_vars_avail) |>
  arrange(desc(acc_pct)) |> pull(variable) |> head(6)

df_work_exp <- df_hr |> dplyr::select(all_of(ed_feats), high_risk)

cat("Computing partial dependence (PRIMARY) for:", paste(top_cont,collapse=", "), "\n")
pdp_list <- lapply(top_cont, function(v){
  xseq <- seq(quantile(df_work_exp[[v]],.05,na.rm=T),
              quantile(df_work_exp[[v]],.95,na.rm=T), length.out=25)
  avg <- sapply(xseq, function(xv){
    dt <- df_work_exp; dt[[v]] <- xv
    mean(predict(rf_hr, newdata=dt, type="prob")[,"YES"], na.rm=T)})
  data.frame(variable=v, x=xseq, prob=avg)
})
pdp_df <- bind_rows(pdp_list)

pdp_df$variable <- factor(relabel_var(pdp_df$variable), levels = unique(relabel_var(pdp_df$variable)))
plots[["p_pdp"]] <- ggplot(pdp_df, aes(x=x, y=prob*100)) +
  geom_line(colour=COL_BLUE, linewidth=1.2) +
  geom_ribbon(aes(ymin=pmax(prob*100-4,0), ymax=pmin(prob*100+4,100)),
              fill=COL_BLUE, alpha=0.12) +
  facet_wrap(~variable, scales="free_x", ncol=3) +
  labs(title = "BURST partial dependence plots on primary outcome",
       x="Variable value", y="Avg predicted probability (%)") + theme_ugib()
print(plots[["p_pdp"]])

ice_var <- top_cont[1]
xseq_ice <- seq(quantile(df_work_exp[[ice_var]],.05,na.rm=T),
                quantile(df_work_exp[[ice_var]],.95,na.rm=T), length.out=20)

set.seed(42)
sample_rows <- sample(seq_len(nrow(df_work_exp)), min(30,nrow(df_work_exp)))
ice_rows <- lapply(sample_rows, function(i){
  probs <- sapply(xseq_ice, function(xv){
    dt <- df_work_exp[i,,drop=FALSE]; dt[[ice_var]] <- xv
    predict(rf_hr, newdata=dt, type="prob")[,"YES"]})
  data.frame(patient=i, x=xseq_ice, prob=probs,
             actual=as.character(df_work_exp$high_risk[i]))})
ice_df <- bind_rows(ice_rows)

plots[["p_ice"]] <- ggplot(ice_df, aes(x=x,y=prob*100,group=patient,colour=actual)) +
  geom_line(alpha=0.4,linewidth=0.7) +
  stat_summary(aes(group=1),fun=mean,geom="line",colour="black",linewidth=1.5) +
  scale_colour_manual(values=c("NO"=COL_GREY,"YES"=COL_RED),
                      labels=c("NO"="No intervention","YES"="Required intervention"), name="Actual") +
  labs(title = sprintf("ICE plot on primary outcome: %s vs high risk predicted probability", relabel_var(ice_var)),
       subtitle="Each line=one patient. Black=average (PDP). Coloured by actual.",
       x = relabel_var(ice_var), y="Predicted probability (%)") + theme_ugib()
print(plots[["p_ice"]])

imp_mat_mt <- importance(rf_mt)
imp_df_mt  <- data.frame(variable=rownames(imp_mat_mt),
                         gini=imp_mat_mt[,"MeanDecreaseGini"],
                         acc =imp_mat_mt[,if ("MeanDecreaseAccuracy" %in% colnames(imp_mat_mt))
                           "MeanDecreaseAccuracy" else "MeanDecreaseGini"],
                         stringsAsFactors=FALSE) |>
  mutate(gini_pct=round(gini/max(gini)*100,1),
         acc_pct=round(acc/max(acc)*100,1)) |>
  arrange(desc(gini_pct))

cat(sprintf("\n--- Importance from SECONDARY RF (rf_mt, outcome=%s) ---\n", SECONDARY))
cat("Top 10 by Gini (for comparison vs PRIMARY):\n")
print(imp_df_mt[1:min(10,nrow(imp_df_mt)), c("variable","gini_pct","acc_pct")])

top10_p <- imp_df$variable[1:min(10,nrow(imp_df))]
top10_s <- imp_df_mt$variable[1:min(10,nrow(imp_df_mt))]
both_top10 <- intersect(top10_p, top10_s)
cat(sprintf("Top-10 features agreeing across both outcomes: %d -- %s\n\n",
            length(both_top10), paste(both_top10, collapse=", ")))

saveRDS(list(imp_df=imp_df, imp_df_mt=imp_df_mt,
             pdp_df=pdp_df, top_cont=top_cont),
        "explainability_results.rds")

cat("\n============================================================\n")
cat("SECTION 11B: SHAP PER-PATIENT EXPLAINABILITY\n")
cat("============================================================\n\n")

if (!requireNamespace("fastshap", quietly = TRUE)) {
  cat("fastshap not installed -- skipping\n")
} else {
  pred_wrapper_rf <- function(object, newdata) {
    predict(object, newdata, type = "prob")[, "YES"]
  }
  
  X_shap <- df_hr[, ed_feats, drop = FALSE]
  cat(sprintf("Computing SHAP values: n=%d patients, %d features, nsim=50...\n",
              nrow(X_shap), length(ed_feats)))
  cat("(Estimated runtime: 2-3 minutes on Apple Silicon)\n\n")
  
  t_shap_start <- Sys.time()
  set.seed(2024)
  shap_rf <- tryCatch(
    fastshap::explain(
      object       = rf_hr,
      X            = X_shap,
      pred_wrapper = pred_wrapper_rf,
      nsim         = 50,
      parallel     = FALSE),
    error = function(e) {cat(sprintf("fastshap::explain failed: %s\n", e$message)); NULL})
  
  if (is.null(shap_rf)) {
    cat("SHAP computation skipped (see error above)\n")
  } else {
    t_shap_elapsed <- as.numeric(difftime(Sys.time(), t_shap_start, units = "mins"))
    cat(sprintf("SHAP computation complete in %.1f minutes\n\n", t_shap_elapsed))
    
    shap_mat <- as.matrix(shap_rf)
    shap_global <- data.frame(
      variable      = colnames(shap_mat),
      mean_abs_shap = colMeans(abs(shap_mat), na.rm = TRUE),
      sd_shap       = apply(shap_mat, 2, sd, na.rm = TRUE),
      stringsAsFactors = FALSE)
    shap_global <- shap_global[order(-shap_global$mean_abs_shap), ]
    shap_global$importance_pct <- round(
      100 * shap_global$mean_abs_shap / max(shap_global$mean_abs_shap), 1)
    rownames(shap_global) <- NULL
    
    cat("Global SHAP importance (mean absolute SHAP value per feature):\n")
    print(shap_global, row.names = FALSE)
    
    oob_prob <- rf_hr$votes[, "YES"] / rowSums(rf_hr$votes)
    y_pri    <- df_hr[[PRIMARY]]
    
    no_idx <- which(y_pri == "NO")
    pid_low <- no_idx[which.min(oob_prob[no_idx])]
    yes_idx <- which(y_pri == "YES")
    pid_high <- yes_idx[which.max(oob_prob[yes_idx])]
    pid_mid <- which.min(abs(oob_prob - 0.5))
    
    cat("\nExemplar patient SHAP breakdowns:\n")
    for (label_pid in list(
      list(label = "Low-risk true-negative", pid = pid_low),
      list(label = "High-risk true-positive", pid = pid_high),
      list(label = "Borderline (P~0.5)", pid = pid_mid))) {
      pid <- label_pid$pid
      cat(sprintf("\n  %-28s patient row %d | actual = %s | predicted P(YES) = %.3f\n",
                  label_pid$label, pid,
                  as.character(y_pri[pid]), oob_prob[pid]))
      shap_pat <- shap_mat[pid, ]
      ord <- order(-abs(shap_pat))
      for (j in ord) {
        feat_name <- names(shap_pat)[j]
        raw_val   <- X_shap[[feat_name]][pid]
        if (is.factor(raw_val)) {
          val_str <- as.character(raw_val)
        } else if (is.numeric(raw_val)) {
          val_str <- format(round(raw_val, 3), width = 8)
        } else {
          val_str <- format(as.character(raw_val), width = 8)
        }
        cat(sprintf("    %-15s value=%-10s SHAP=%+.4f\n",
                    feat_name, val_str, shap_pat[j]))
      }
    }
    
    saveRDS(list(
      shap_values_matrix = shap_mat,
      shap_global        = shap_global,
      exemplar_patients  = list(
        low_risk     = list(pid = pid_low,  prob = oob_prob[pid_low]),
        high_risk    = list(pid = pid_high, prob = oob_prob[pid_high]),
        borderline   = list(pid = pid_mid,  prob = oob_prob[pid_mid])),
      nsim               = 50,
      runtime_min        = t_shap_elapsed),
      "shap_results.rds")
    cat("\nSaved: shap_results.rds\n")
    
    if (requireNamespace("shapviz", quietly = TRUE)) {
      sv <- shapviz::shapviz(shap_rf, X = X_shap)
      tryCatch({
        .sv_relab <- sv; colnames(.sv_relab$S) <- relabel_var(colnames(.sv_relab$S)); colnames(.sv_relab$X) <- relabel_var(colnames(.sv_relab$X))
        plots[["p_shap_importance"]] <- shapviz::sv_importance(.sv_relab) +
          ggtitle("BURST global SHAP importance (mean |SHAP| per feature)") +
          theme_ugib()
        print(plots[["p_shap_importance"]])
        
        plots[["p_shap_dep_hb"]] <- shapviz::sv_dependence(.sv_relab, "Hemoglobin") +
          labs(title = "SHAP dependence: Hemoglobin")
        ggtitle("SHAP dependence: hb_mean") +
          theme_ugib()
        print(plots[["p_shap_dep_hb"]])
        
        cat("Saved SHAP plots to the standard plot stream.\n")
      }, error = function(e) {
        cat(sprintf("shapviz plotting skipped: %s\n", e$message))
      })
    }
  }
}

# SECTION 12: ORDINAL REGRESSION
cat("\n============================================================\n")
cat("SECTION 12: ORDINAL REGRESSION\n")
cat("============================================================\n\n")

imputed <- tryCatch(readRDS("ugib_mice_object.rds"),
                    error=function(e){cat("ugib_mice_object.rds not found -- skipping\n"); NULL})

or_df <- NULL
if(!is.null(imputed)){
  if(!requireNamespace("MASS",quietly=T)) stop("Install MASS")
  if(!requireNamespace("mice",quietly=T)) stop("Install mice")
  
  polr_feats <- imp_df |>
    filter(variable %in% ed_feats) |>
    arrange(desc(acc_pct)) |>
    pull(variable) |> head(7)
  
  d_preview <- mice::complete(imputed, 1)
  if(all(c("hb_ed","hb_entry") %in% names(d_preview)))
    d_preview$hb_mean <- (d_preview$hb_ed + d_preview$hb_entry)/2
  feats_used    <- polr_feats[polr_feats %in% names(d_preview) &
                                sapply(polr_feats, function(v) is.numeric(d_preview[[v]]))]
  feats_dropped <- setdiff(polr_feats, feats_used)
  n_events <- sum(df_full$transfusion_class != "0_units", na.rm=TRUE)
  
  cat(sprintf("Ordinal predictors (%d numeric of %d candidates): %s\n",
              length(feats_used), length(polr_feats),
              paste(feats_used, collapse=", ")))
  if (length(feats_dropped) > 0) {
    cat(sprintf("  Dropped (non-numeric, factor predictors not fit): %s\n",
                paste(feats_dropped, collapse=", ")))
  }
  cat(sprintf("EPV: %d events / %d = %.1f (interpretive only)\n\n",
              n_events, length(feats_used), n_events / length(feats_used)))
  
  m_imp <- imputed$m
  fits  <- vector("list", m_imp)
  for(k in seq_len(m_imp)){
    d_k <- mice::complete(imputed, k)
    if(all(c("hb_ed","hb_entry") %in% names(d_k)))
      d_k$hb_mean <- (d_k$hb_ed + d_k$hb_entry)/2
    d_k$transfusion_class <- df_full$transfusion_class
    feats_ok <- polr_feats[polr_feats %in% names(d_k) &
                             sapply(polr_feats, function(v) is.numeric(d_k[[v]]))]
    if(length(feats_ok)<2) next
    frm_k <- as.formula(paste("transfusion_class ~", paste(feats_ok,collapse="+")))
    d_k <- na.omit(d_k[, c(feats_ok, "transfusion_class")])
    fits[[k]] <- tryCatch(
      suppressWarnings(MASS::polr(frm_k, data=d_k, method="logistic", Hess=TRUE)),
      error=function(e) NULL)
  }
  
  valid <- fits[!sapply(fits, is.null)]
  m_val <- length(valid)
  cat(sprintf("Converged: %d / %d\n\n", m_val, m_imp))
  
  if(m_val>=2){
    all_coefs <- lapply(valid, coef)
    all_se2   <- lapply(valid, function(f) {
      vc <- tryCatch(vcov(f), error=function(e) NULL)
      if (is.null(vc)) return(rep(NA_real_, length(coef(f))))
      cn <- names(coef(f))
      d  <- diag(vc)
      d[cn]
    })
    var_names <- names(all_coefs[[1]])
    Qbar <- rowMeans(do.call(cbind, all_coefs))
    Ubar <- rowMeans(do.call(cbind, all_se2))
    B    <- apply(do.call(cbind, all_coefs), 1, var)
    T_v  <- Ubar + (1+1/m_val)*B
    SE   <- sqrt(T_v)
    z    <- Qbar/SE
    pval <- 2*(1-pnorm(abs(z)))
    
    or_df <- data.frame(
      predictor=var_names, log_OR=round(Qbar,4),
      OR=round(exp(Qbar),3),
      CI_low=round(exp(Qbar-1.96*SE),3),
      CI_high=round(exp(Qbar+1.96*SE),3),
      p_value=round(pval,4),
      sig=case_when(pval<0.001~"***",pval<0.01~"**",pval<0.05~"*",
                    pval<0.1~".",TRUE~"ns"),
      stringsAsFactors=FALSE) |>
      filter(!grepl("^[0-9]", predictor))
    
    cat("Pooled ordinal regression (Rubin's rules):\n"); print(or_df)
    
    or_df$predictor <- relabel_var(or_df$predictor)
    plots[["p_or"]] <- ggplot(or_df |>
                                mutate(predictor=factor(predictor,levels=rev(predictor))),
                              aes(x=OR, y=predictor, colour=sig %in% c("*","**","***"))) +
      geom_point(size=3.5) +
      geom_errorbar(aes(xmin=CI_low,xmax=CI_high),
                    orientation="y", width=0.25, linewidth=0.9) +
      geom_vline(xintercept=1,linetype="dashed",colour="grey50",linewidth=0.8) +
      scale_colour_manual(values=c("TRUE"=COL_BLUE,"FALSE"=COL_GREY),
                          labels=c("TRUE"="p < 0.05","FALSE"="p >= 0.05"), name=NULL) +
      scale_x_log10() +
      labs(title="BURST ordinal logistic regression (pooled, Rubin's rules)",
           subtitle=sprintf("Outcome: transfusion class | m=%d MICE datasets",m_val),
           x="Odds ratio (log scale, 95% CI)", y=NULL) + theme_ugib()
    print(plots[["p_or"]])
    
    saveRDS(list(or_df=or_df, m_val=m_val, polr_feats=polr_feats),
            "ordinal_results.rds")
  } else cat("Insufficient convergence\n")
} else cat("Skipping ordinal regression\n")

# SECTION 12B: ORDINAL PREDICTION PERFORMANCE
cat("\n============================================================\n")
cat("SECTION 12B: ORDINAL PREDICTION\n")
cat("============================================================\n\n")

if (!requireNamespace("MASS", quietly=TRUE)) {
  cat("MASS not installed -- skipping Section 12B\n")
} else {
  ord_feats <- ed_feats[sapply(ed_feats, function(v) is.numeric(df_ed[[v]]))]
  if (exists("imp_df")) {
    ord_feats <- imp_df |> filter(variable %in% ord_feats) |>
      arrange(desc(acc_pct)) |> pull(variable) |> head(7)
  } else {
    ord_feats <- head(ord_feats, 7)
  }
  cat(sprintf("Ordinal predictors (%d): %s\n", length(ord_feats),
              paste(ord_feats, collapse=", ")))
  
  med_impute <- function(X) {
    meds <- vapply(X, function(c) median(c, na.rm=TRUE), numeric(1))
    list(X = as.data.frame(lapply(seq_along(X), function(j) {
      v <- X[[j]]; v[is.na(v)] <- meds[[j]]; v })) |>
        setNames(names(X)),
      meds = meds)
  }
  
  cv_ordinal <- function(df_in, feats, outcome,
                         n_repeats = N_REPEATS, n_folds = N_FOLDS, seed = 2024) {
    y_full  <- df_in[[outcome]]
    cl      <- levels(y_full)
    confmat <- matrix(0L, nrow=length(cl), ncol=length(cl), dimnames=list(cl, cl))
    acc_per <- numeric(n_repeats); mae_per <- numeric(n_repeats)
    auc_per <- vector("list", length(cl)); names(auc_per) <- cl
    for (i in seq_along(auc_per)) auc_per[[i]] <- numeric(n_repeats)
    y_int <- as.integer(y_full)
    for (rep in seq_len(n_repeats)) {
      fld <- make_strat_folds(as.character(y_full), n_folds, seed + rep)
      pred_int <- integer(length(y_full))
      pred_pp  <- matrix(NA_real_, nrow=length(y_full), ncol=length(cl),
                         dimnames=list(NULL, cl))
      for (k in seq_len(n_folds)) {
        tr <- which(fld != k); te <- which(fld == k)
        Xtr <- df_in[tr, feats, drop=FALSE]; Xte <- df_in[te, feats, drop=FALSE]
        imp <- med_impute(Xtr)
        Xtr_i <- imp$X
        Xte_i <- as.data.frame(lapply(seq_along(Xte), function(j) {
          v <- Xte[[j]]; v[is.na(v)] <- imp$meds[[j]]; v })) |> setNames(feats)
        ytr <- y_full[tr]
        if (length(unique(ytr)) < length(cl)) {
          maj <- names(which.max(table(ytr)))
          pred_int[te] <- as.integer(factor(maj, levels=cl))
          pred_pp[te, ] <- 0
          pred_pp[te, maj] <- 1
          next
        }
        fit <- tryCatch(suppressWarnings(
          MASS::polr(reformulate(feats, "ytr"),
                     data = cbind(Xtr_i, ytr=ytr), Hess=TRUE, method="logistic")
        ), error=function(e) NULL)
        if (is.null(fit)) {
          maj <- names(which.max(table(ytr)))
          pred_int[te] <- as.integer(factor(maj, levels=cl))
          pred_pp[te, maj] <- 1; next
        }
        pp_te <- predict(fit, newdata=Xte_i, type="probs")
        if (is.null(dim(pp_te))) pp_te <- matrix(pp_te, nrow=1)
        pred_pp[te, colnames(pp_te)] <- pp_te
        pred_int[te] <- max.col(pp_te, ties.method="first")
      }
      ok <- !is.na(pred_int) & pred_int > 0
      tab <- table(factor(y_int[ok], levels=seq_along(cl)),
                   factor(pred_int[ok], levels=seq_along(cl)))
      confmat <- confmat + as.matrix(tab)
      acc_per[rep] <- mean(pred_int[ok] == y_int[ok])
      mae_per[rep] <- mean(abs(pred_int[ok] - y_int[ok]))
      for (i in seq_along(cl)) {
        y_bin <- as.integer(y_int == i)
        s_bin <- pred_pp[, cl[i]]
        ok2 <- !is.na(s_bin) & !is.na(y_bin)
        if (length(unique(y_bin[ok2])) < 2) {
          auc_per[[cl[i]]][rep] <- NA_real_
        } else {
          auc_per[[cl[i]]][rep] <-
            as.numeric(auc(roc(y_bin[ok2], s_bin[ok2], quiet=TRUE, direction="<")))
        }
      }
    }
    list(acc=acc_per, mae=mae_per, auc_per_class=auc_per, confmat=confmat)
  }
  
  ord_pred <- cv_ordinal(df_ed, ord_feats, "transfusion_class")
  cat(sprintf("\n5x5 CV ordinal accuracy: %.1f%% (sd=%.2f%%)\n",
              mean(ord_pred$acc, na.rm=TRUE)*100,
              sd(ord_pred$acc, na.rm=TRUE)*100))
  cat(sprintf("Mean absolute class error: %.3f (sd=%.3f)\n",
              mean(ord_pred$mae, na.rm=TRUE), sd(ord_pred$mae, na.rm=TRUE)))
  cat("Per-class one-vs-rest AUC (mean across repeats):\n")
  for (cl_i in names(ord_pred$auc_per_class)) {
    a <- ord_pred$auc_per_class[[cl_i]]
    cat(sprintf("  %-10s  AUC=%.3f (sd=%.3f)\n",
                cl_i, mean(a, na.rm=TRUE), sd(a, na.rm=TRUE)))
  }
  cat("Confusion (rows=true, cols=pred, summed across all CV folds & repeats):\n")
  print(ord_pred$confmat)
  
  cls_v <- rownames(ord_pred$confmat)
  total <- sum(ord_pred$confmat)
  cat("\nPer-class sensitivity / specificity:\n")
  for (i in seq_along(cls_v)) {
    tp <- ord_pred$confmat[i, i]
    fn <- sum(ord_pred$confmat[i, ]) - tp
    fp <- sum(ord_pred$confmat[, i]) - tp
    tn <- total - tp - fn - fp
    sens <- if ((tp+fn)>0) tp/(tp+fn) else NA
    spec <- if ((tn+fp)>0) tn/(tn+fp) else NA
    cat(sprintf("  %-10s  sens=%.2f  spec=%.2f\n", cls_v[i], sens, spec))
  }
  
  cm_df <- as.data.frame(as.table(ord_pred$confmat))
  names(cm_df) <- c("True","Predicted","Count")
  .ord_cm_map <- c("0_units"="0 units", "1to2"="1-2 units", "3plus"="3 plus units")
  cm_df$Predicted <- factor(.ord_cm_map[as.character(cm_df$Predicted)], levels=c("0 units","1-2 units","3 plus units"))
  cm_df$True      <- factor(.ord_cm_map[as.character(cm_df$True)],      levels=c("3 plus units","1-2 units","0 units"))
  plots[["p_ord_cm"]] <- ggplot(cm_df, aes(x=Predicted, y=True, fill=Count)) +
    geom_tile(colour="white", linewidth=0.5) +
    geom_text(aes(label=Count), colour="white", fontface="bold", size=4.5) +
    scale_fill_gradient(low="#7EB5D6", high=COL_BLUE) +
    labs(title="BURST ordinal prediction of transfusion class -- pooled confusion matrix",
         subtitle=sprintf("MASS::polr, 5x5 CV; cells = sum across %d folds x %d repeats",
                          N_FOLDS, N_REPEATS),
         x="Predicted class", y="True class") +
    theme_ugib() + theme(panel.grid=element_blank())
  print(plots[["p_ord_cm"]])
  
  saveRDS(ord_pred, "ordinal_prediction.rds")
  cat("\nSaved: ordinal_prediction.rds\n")
}

# SECTION 14: SENSITIVITY ANALYSES
cat("\n============================================================\n")
cat("SECTION 14: SENSITIVITY ANALYSES\n")
cat("============================================================\n\n")

cat("--- 14A. V1 sensitivity: high_risk new (no rebleed) vs old (with rebleed) ---\n")
n_flip <- sum(df_full$high_risk != df_full$high_risk_old, na.rm=TRUE)
cat(sprintf("Patients flipping with redefinition: %d\n", n_flip))
cat(sprintf("Events: NEW=%d  OLD=%d\n",
            sum(df_full$high_risk == "YES", na.rm=TRUE),
            sum(df_full$high_risk_old == "YES", na.rm=TRUE)))

scores_for_sens <- list(GBS=as.numeric(df_full$gbs_score),
                        AIMS65=as.numeric(df_full$aims65_score),
                        Pre_Rockall=as.numeric(df_full$pre_rockall),
                        Full_Rockall=as.numeric(df_full$full_rockall),
                        CANUKA=as.numeric(df_full$canuka))
sens_v1_rows <- list()
cat(sprintf("\n%-14s %8s %8s %8s\n","Score","AUC NEW","AUC OLD","Delta"))
cat(paste(rep("-", 50), collapse=""), "\n")
for (sc_n in names(scores_for_sens)) {
  s <- scores_for_sens[[sc_n]]
  a_new <- score_auc(s, df_full$high_risk)
  a_old <- score_auc(s, df_full$high_risk_old)
  cat(sprintf("%-14s %8.3f %8.3f %+8.4f\n", sc_n, a_new, a_old, a_new - a_old))
  sens_v1_rows[[length(sens_v1_rows)+1]] <- data.frame(
    Score=sc_n, AUC_new=a_new, AUC_old=a_old, Delta=a_new-a_old,
    stringsAsFactors=FALSE)
}
sens_v1_df <- do.call(rbind, sens_v1_rows)
cat("\nVerdict: redefinition is statistically inert if all |delta| < 0.02\n")

sens_v1_long <- rbind(
  data.frame(Score=sens_v1_df$Score, AUC=sens_v1_df$AUC_new,
             Defn="Rebleed excluded", stringsAsFactors=FALSE),
  data.frame(Score=sens_v1_df$Score, AUC=sens_v1_df$AUC_old,
             Defn="Rebleed included", stringsAsFactors=FALSE))
plots[["p_sens_v1"]] <- ggplot(sens_v1_long,
                               aes(x=AUC, y=Score, colour=Defn)) +
  geom_point(position=position_dodge(width=0.4), size=3.5) +
  scale_colour_manual(values=c("Rebleed included"=COL_RED,
                               "Rebleed excluded"=COL_BLUE)) +
  xlim(0.4, 1.0) +
  labs(title="Sensitivity analysis: Comparator-score AUCs under high risk UGIB (with vs without rebleeding)",
       subtitle=sprintf("%d patients flip with rebleed dropped", n_flip),
       x="AUC", y=NULL, colour=NULL) +
  theme_ugib() + theme(panel.grid.major.y=element_blank())
print(plots[["p_sens_v1"]])

cat("\n--- 14B. V2 sensitivity: BLEED_CONFIRMED == YES subset ---\n")
if ("bleed_confirmed" %in% names(df_full)) {
  bc <- df_full$bleed_confirmed
} else if ("BLEED_CONFIRMED" %in% names(df_full)) {
  bc <- df_full$BLEED_CONFIRMED
} else {
  bc <- NULL
}
sens_v2_df <- NULL
if (!is.null(bc)) {
  mask <- bc == "YES" & !is.na(bc)
  df_sub <- df_full[mask, , drop=FALSE]
  cat(sprintf("Subset n=%d (vs %d in primary)\n", nrow(df_sub), nrow(df_full)))
  sens_v2_rows <- list()
  cat(sprintf("\n%-14s %-14s %8s\n","Score","Outcome","AUC"))
  cat(paste(rep("-", 45), collapse=""), "\n")
  for (sc_n in names(scores_for_sens)) {
    s <- scores_for_sens[[sc_n]][mask]
    for (oc in c("major_tx","high_risk","rebleeding")) {
      if (!oc %in% names(df_sub)) next
      a <- score_auc(s, df_sub[[oc]])
      cat(sprintf("%-14s %-14s %8.3f\n", sc_n, oc, a))
      sens_v2_rows[[length(sens_v2_rows)+1]] <- data.frame(
        Score=sc_n, Outcome=oc, AUC=a, stringsAsFactors=FALSE)
    }
  }
  sens_v2_df <- do.call(rbind, sens_v2_rows)
} else {
  cat("BLEED_CONFIRMED column not found -- skipping V2 sensitivity.\n")
}

saveRDS(list(sens_v1=sens_v1_df, sens_v2=sens_v2_df, n_flip=n_flip),
        "sensitivity_results.rds")
cat("\nSaved: sensitivity_results.rds\n")

# SECTION 13: TABLE 1 + TABLE 2 + PUBLICATION FIGURES
cat("\n============================================================\n")
cat("SECTION 13: TABLE 1 & PUBLICATION OUTPUTS\n")
cat("============================================================\n\n")

cont_row <- function(var, label, dig=1) {
  x<-df_full[[var]]; tc<-df_full$transfusion_class
  ok<-!is.na(x)&!is.na(tc); if(sum(ok)<5) return(NULL)
  fmt <- function(v) sprintf(paste0("%.",dig,"f (%.",dig,"f)"),
                             median(v,na.rm=T), IQR(v,na.rm=T))
  kw_p <- tryCatch(kruskal.test(x[ok]~tc[ok])$p.value, error=function(e) NA)
  data.frame(Variable=label, Overall=fmt(x[ok]),
             Class_0=fmt(x[ok&tc=="0_units"]), Class_1to2=fmt(x[ok&tc=="1to2"]),
             Class_3plus=fmt(x[ok&tc=="3plus"]),
             p=ifelse(is.na(kw_p),"NA",ifelse(kw_p<0.001,"<0.001",sprintf("%.3f",kw_p))),
             stringsAsFactors=FALSE)
}

cat_row <- function(var, label) {
  x  <- df_full[[var]]; tc <- df_full$transfusion_class
  x  <- ifelse(toupper(as.character(x)) %in% c("M","MALE"),   "YES",
               ifelse(toupper(as.character(x)) %in% c("F","FEMALE"), "NO", as.character(x)))
  ok <- !is.na(x) & !is.na(tc) & toupper(x) %in% c("YES","NO")
  if(sum(ok) < 5) return(NULL)
  x_ok  <- toupper(x[ok])
  tc_ok <- tc[ok]
  fmt <- function(idx) {
    n_y <- sum(x_ok[idx] == "YES", na.rm=TRUE)
    n_t <- sum(idx)
    if(n_t == 0) "0 (NA)" else sprintf("%d (%.1f%%)", n_y, n_y/n_t*100)
  }
  tbl <- table(x_ok, tc_ok)
  use_f <- any(suppressWarnings(tryCatch(chisq.test(tbl)$expected,
                                         error=function(e) matrix(0))) < 5)
  pv <- tryCatch(
    if(use_f) fisher.test(tbl, simulate.p.value=TRUE, B=1e4)$p.value
    else chisq.test(tbl)$p.value,
    error=function(e) NA)
  data.frame(Variable=label,
             Overall     = fmt(rep(TRUE, length(x_ok))),
             Class_0     = fmt(tc_ok == "0_units"),
             Class_1to2  = fmt(tc_ok == "1to2"),
             Class_3plus = fmt(tc_ok == "3plus"),
             p = ifelse(is.na(pv), "NA",
                        ifelse(pv < 0.001, "<0.001", sprintf("%.3f", pv))),
             stringsAsFactors=FALSE)
}

safe_row <- function(fn,...) tryCatch(fn(...), error=function(e) NULL)

table1 <- bind_rows(Filter(Negate(is.null), list(
  data.frame(Variable="n", Overall=sprintf("%d",n),
             Class_0     = sprintf("%d (%.1f%%)", sum(df_full$transfusion_class=="0_units"),
                                   100*mean(df_full$transfusion_class=="0_units")),
             Class_1to2  = sprintf("%d (%.1f%%)", sum(df_full$transfusion_class=="1to2"),
                                   100*mean(df_full$transfusion_class=="1to2")),
             Class_3plus = sprintf("%d (%.1f%%)", sum(df_full$transfusion_class=="3plus"),
                                   100*mean(df_full$transfusion_class=="3plus")),
             p="", stringsAsFactors=FALSE),
  safe_row(cont_row,"age","Age (years)",1),
  safe_row(cat_row,"gender","Male sex"),
  safe_row(cont_row,"hosp_days","Hospital stay (days)",1),
  data.frame(Variable="--- Presentation ---",Overall="",Class_0="",
             Class_1to2="",Class_3plus="",p="",stringsAsFactors=FALSE),
  safe_row(cat_row,"hemodynamic_instab","Hemodynamic instability"),
  safe_row(cat_row,"hematemesis_red","Hematemesis (bright red)"),
  safe_row(cat_row,"hematemesis_coffee","Hematemesis (coffee grounds)"),
  safe_row(cat_row,"melena","Melena"),
  safe_row(cat_row,"hematochezia","Hematochezia"),
  safe_row(cat_row,"syncope","Syncope"),
  safe_row(cat_row,"altered_mental","Altered mental status"),
  data.frame(Variable="--- Laboratory ---",Overall="",Class_0="",
             Class_1to2="",Class_3plus="",p="",stringsAsFactors=FALSE),
  safe_row(cont_row,"hb_mean","Haemoglobin mean (g/dL)",2),
  safe_row(cont_row,"albumin","Albumin (g/dL)",2),
  safe_row(cont_row,"urea","Urea (mg/dL)",1),
  safe_row(cont_row,"creatinine","Creatinine (mg/dL)",2),
  safe_row(cont_row,"inr","INR",2),
  safe_row(cont_row,"aptt","aPTT (sec)",1),
  safe_row(cont_row,"plt_entry","Platelets (x10^3/uL)",1),
  safe_row(cont_row,"wbc_entry","WBC (x10^3/uL)",1),
  safe_row(cont_row,"lactate_ed","Lactate ED (mmol/L)",2),
  data.frame(Variable="--- Scores ---",Overall="",Class_0="",
             Class_1to2="",Class_3plus="",p="",stringsAsFactors=FALSE),
  safe_row(cont_row,"gbs_score","GBS",1),
  safe_row(cont_row,"aims65_score","AIMS65",0),
  safe_row(cont_row,"pre_rockall","Pre-Rockall",1),
  safe_row(cont_row,"full_rockall","Full Rockall",1),
  data.frame(Variable="--- Comorbidities ---",Overall="",Class_0="",
             Class_1to2="",Class_3plus="",p="",stringsAsFactors=FALSE),
  safe_row(cat_row,"hypertension","Hypertension"),
  safe_row(cat_row,"heart_failure","Heart failure"),
  safe_row(cat_row,"afib","Atrial fibrillation"),
  safe_row(cat_row,"liver_cirrhosis","Liver cirrhosis"),
  safe_row(cat_row,"renal_ckd","CKD"),
  safe_row(cat_row,"renal_aki","AKI"),
  safe_row(cat_row,"cancer_gi","GI malignancy"),
  safe_row(cat_row,"cancer_extra","Extraintestinal malignancy"),
  safe_row(cat_row,"copd","COPD"),
  safe_row(cat_row,"dm_controlled","Diabetes (controlled)"),
  safe_row(cat_row,"dm_uncontrolled","Diabetes (uncontrolled)")
)))

colnames(table1) <- c("Variable","Overall","0 units","1-2 units",">=3 units","p-value")
cat(sprintf("Table 1: %d rows\n", nrow(table1)))
print(table1)

fmt_auc <- function(roc_obj) {
  if (is.null(roc_obj)) return(c(auc=NA_real_, ci_str=NA_character_))
  ci <- ci.auc(roc_obj)
  c(auc = round(as.numeric(auc(roc_obj)), 3),
    ci_str = sprintf("%.3f-%.3f", ci[1], ci[3]))
}

roc_sc_p <- get_score_rocs(PRIMARY)
roc_sc_s <- get_score_rocs(SECONDARY)

table2_rows <- list()
add_row <- function(model, type, roc_p, roc_s) {
  fp <- fmt_auc(roc_p); fs <- fmt_auc(roc_s)
  table2_rows[[length(table2_rows)+1]] <<- data.frame(
    Model = model, Type = type,
    AUC_primary = fp["auc"], CI_primary = fp["ci_str"],
    AUC_secondary = fs["auc"], CI_secondary = fs["ci_str"],
    stringsAsFactors=FALSE)
}
add_row("BURST",   "ML",
        roc_ml_hr, roc_ml_mt)
add_row("BURST (ridge LR)",   "ML",
        roc_lr_hr, roc_lr_mt)
add_row("GBS",                "Score", roc_sc_p[[1]], roc_sc_s[[1]])
add_row("AIMS65",             "Score", roc_sc_p[[2]], roc_sc_s[[2]])
add_row("Pre-Rockall",        "Score", roc_sc_p[[3]], roc_sc_s[[3]])
add_row("Full Rockall",       "Score", roc_sc_p[[4]], roc_sc_s[[4]])

table2 <- do.call(rbind, table2_rows)
names(table2) <- c("Model","Type",
                   sprintf("AUC vs %s (PRIMARY)",   PRIMARY),
                   sprintf("CI95 vs %s (PRIMARY)",  PRIMARY),
                   sprintf("AUC vs %s (SECONDARY)", SECONDARY),
                   sprintf("CI95 vs %s (SECONDARY)",SECONDARY))
cat(sprintf("\nTable 2: model performance vs PRIMARY (%s) and SECONDARY (%s)\n",
            PRIMARY, SECONDARY))
print(table2)

roc_list_sc_pub <- roc_sc_p
score_aucs_pub  <- sapply(roc_list_sc_pub,
                          function(r) if (is.null(r)) NA else as.numeric(auc(r)))
leg_lab_pub <- sprintf("BURST  AUC=%.3f", as.numeric(auc(roc_ml_hr)))
leg_col_pub <- COL_BLUE; leg_lwd_pub <- 3
for (i in order(score_aucs_pub, decreasing = TRUE))
  if (!is.null(roc_list_sc_pub[[i]])) {
    leg_lab_pub <- c(leg_lab_pub, sprintf("%-12s  AUC=%.3f",
                                          score_names[i], as.numeric(auc(roc_list_sc_pub[[i]]))))
    leg_col_pub <- c(leg_col_pub, score_cols[i]); leg_lwd_pub <- c(leg_lwd_pub, 2)
  }

par(mar=c(5,5,4,2), font.lab=2, cex.lab=1.1)
plot(roc_ml_hr, col=COL_BLUE, lwd=3,
     main = "BURST vs traditional score on high risk upper gastrointestinal bleeding",
     legacy.axes = TRUE,
     xlim = c(1, 0), ylim = c(0, 1), asp = NA,
     xaxs = "i", yaxs = "i",
     font.main=2, cex.main=1.1)
for(i in seq_along(roc_list_sc_pub))
  if(!is.null(roc_list_sc_pub[[i]]))
    plot(roc_list_sc_pub[[i]], col=score_cols[i], lwd=2, legacy.axes = TRUE, add=TRUE)
legend("bottomright", bty="n", lwd=leg_lwd_pub, col=leg_col_pub,
       legend=leg_lab_pub, cex=0.9)

plots[["p_pub_imp"]] <- ggplot(imp_df[1:min(15,nrow(imp_df)),] |>
                                 mutate(variable_disp = relabel_var(variable)),
                               aes(x=gini_pct, y=reorder(variable_disp, gini_pct),
                                   fill=variable %in% intersect(top15_gini, top15_acc))) +
  geom_col(width=0.75, colour="white") +
  scale_fill_manual(values=c("TRUE"=COL_BLUE,"FALSE"="#7EB5D6"),
                    labels=c("TRUE"="Top-15 both","FALSE"="Gini only"), name=NULL) +
  labs(title = "BURST variable importance on high-risk upper gastrointestinal bleeding (RF)",
       x="Mean Decrease Gini (% of max)", y=NULL) +
  theme_pub() + theme(legend.position="bottom")
print(plots[["p_pub_imp"]])

plots[["p_pub_cal"]] <- ggplot(hl_g_pri, aes(x=exp_rate,y=obs_rate)) +
  geom_abline(slope=1,intercept=0,linetype="dashed",colour="grey50",linewidth=1) +
  geom_point(size=4,colour=COL_BLUE) +
  geom_errorbar(aes(ymin=pmax(obs_rate-1.96*se,0),
                    ymax=pmin(obs_rate+1.96*se,1)), width=0.015,colour=COL_BLUE,linewidth=0.9) +
  annotate("text",x=0.08,y=0.85,label=sprintf("HL p = %.3f",hl_p_pri),size=4,colour="grey30") +
  scale_x_continuous(limits=c(0,1),labels=pct_fmt) +
  scale_y_continuous(limits=c(0,1),labels=pct_fmt) +
  labs(title="Calibration: BURST on primary outcome",
       subtitle=sprintf("Hosmer-Lemeshow | %s (PRIMARY)", PRIMARY),
       x="Mean predicted probability", y="Observed event rate") +
  theme_pub() + theme(legend.position="none")
print(plots[["p_pub_cal"]])

dca_key <- paste0("p_dca_", PRIMARY)
if (dca_key %in% names(plots)) {
  plots[["p_pub_dca"]] <- plots[[dca_key]] +
    labs(title = "BURST decision curve analysis on high-risk upper gastrointestinal bleeding") +
    theme_pub()
  print(plots[["p_pub_dca"]])
}

cat("\n============================================================\n")
cat("SAVING ALL OUTPUTS\n")
cat("============================================================\n\n")

safe_get <- function(obj_name)
  tryCatch(get(obj_name), error=function(e){
    cat(sprintf("  WARNING: %s not found -> NULL\n",obj_name)); NULL})

saveRDS(list(rf_tc=safe_get("rf_tc"), rf_mt=safe_get("rf_mt"),
             rf_hr=safe_get("rf_hr"), roc_mt=safe_get("roc_ml_mt"),
             roc_hr=safe_get("roc_ml_hr"), features=ed_feats),
        "ed_model_weighted.rds")

if(requireNamespace("writexl",quietly=T)){
  library(writexl)
  nri_long <- tryCatch({
    rows_l <- list()
    for (sc in names(nri_results)) {
      for (t_name in names(nri_results[[sc]])) {
        rec <- nri_results[[sc]][[t_name]]
        rows_l[[length(rows_l) + 1L]] <- data.frame(
          Score        = sc,
          Threshold    = as.numeric(t_name),
          NRI_events   = rec$events,
          NRI_nonevents = rec$nonevents,
          NRI_total    = rec$total,
          stringsAsFactors = FALSE)
      }
    }
    do.call(rbind, rows_l)
  }, error = function(e) data.frame(note = sprintf("NRI build error: %s", conditionMessage(e))))
  
  write_xlsx(list(
    "Table 1 Characteristics" = table1,
    "Table 2 Performance"     = table2,
    "OR table"  = if(!is.null(or_df)) or_df else data.frame(note="Not fitted"),
    "DeLong"    = delong_df,
    "NRI categorical" = nri_long,
    "Importance" = imp_df
  ), "thesis_tables.xlsx")
  cat("Saved: thesis_tables.xlsx\n")
} else cat("writexl not installed -- skipping Excel export\n")

cat("\n--- Exporting named PNGs at 300 dpi ---\n")

.fig_map <- list(
  list("p_cv_stability",      "4.25",   "CV stability per-repeat AUCs"),
  list("p_cal_pri",           "4.27",   "Calibration plot primary"),
  list("p_pub_cal",           "4.27b",  "Calibration with Brier annotation"),
  list("p_cal",               "4.28",   "Calibration plot secondary"),
  list("p_learning_curve",    "4.26",   "Learning curve primary"),
  list("p_dca_high_risk",     "4.31",   "DCA primary"),
  list("p_dca_major_tx",      "4.32",   "DCA secondary"),
  list("p_gini",              "4.33",   "Gini importance headline RF"),
  list("p_perm_imp",          "4.34",   "Permutation importance B=50"),
  list("p_pdp",               "4.37",   "Partial dependence top continuous"),
  list("p_ice",               "4.38",   "ICE plot hb_mean"),
  list("p_extval",            "4.49",   "Indirect external validation forest")
)

.png_out <- "final_figures"
if (!dir.exists(.png_out)) dir.create(.png_out, recursive = TRUE)

.manifest <- data.frame(
  filename = character(), fig_no = character(),
  caption  = character(), source = character(),
  stringsAsFactors = FALSE)

for (.m in .fig_map) {
  .key <- .m[[1]]; .fno <- .m[[2]]; .cap <- .m[[3]]
  if (is.null(plots[[.key]])) {
    cat(sprintf("  SKIP   plot[[%s]] missing -- Fig %s not exported\n", .key, .fno))
    next
  }
  .fn <- file.path(.png_out, sprintf("final_fig_%s_%s.png", .fno, .key))
  tryCatch({
    ggsave(.fn, plot = plots[[.key]], width = 10, height = 7, dpi = 300, bg = "white")
    cat(sprintf("  saved  %s  -> Fig %s | %s\n", basename(.fn), .fno, .cap))
    .manifest <- rbind(.manifest, data.frame(
      filename = basename(.fn), fig_no = .fno, caption = .cap,
      source = .key, stringsAsFactors = FALSE))
  }, error = function(e) {
    cat(sprintf("  FAIL   %s: %s\n", .fn, conditionMessage(e)))
  })
}

.fn_roc_pri <- file.path(.png_out, "final_fig_4.24_roc_ml_vs_scores_primary.png")
png(.fn_roc_pri, width = 10*300, height = 7*300, res = 300, bg = "white")
.roc_list_sc_pri <- get_score_rocs(PRIMARY)
par(mar = c(5, 5, 4, 2))
plot(roc_ml_hr, col = COL_BLUE, lwd = 3,
     main = sprintf("BURST vs traditional scores [PRIMARY: %s]", PRIMARY),
     legacy.axes = TRUE,
     xlim = c(1, 0), ylim = c(0, 1), asp = NA,
     xaxs = "i", yaxs = "i")
.leg_lab <- sprintf("BURST  AUC=%.3f", as.numeric(auc(roc_ml_hr)))
.leg_col <- COL_BLUE; .leg_lwd <- 3
.sc_auc_pri <- sapply(.roc_list_sc_pri, function(r) if (is.null(r)) NA_real_ else as.numeric(auc(r)))
for (i in order(.sc_auc_pri, decreasing = TRUE))
  if (!is.null(.roc_list_sc_pri[[i]])) {
    plot(.roc_list_sc_pri[[i]], col = score_cols[i], lwd = 2, legacy.axes = TRUE, add = TRUE)
    .leg_lab <- c(.leg_lab, sprintf("%-12s  AUC=%.3f",
                                    score_names[i],
                                    as.numeric(auc(.roc_list_sc_pri[[i]]))))
    .leg_col <- c(.leg_col, score_cols[i]); .leg_lwd <- c(.leg_lwd, 2)
  }
legend("bottomright", bty = "n", lwd = .leg_lwd, col = .leg_col,
       legend = .leg_lab, cex = 0.9)
dev.off()
cat(sprintf("  saved  %s  -> Fig 4.24 | ROC ML vs scores primary\n",
            basename(.fn_roc_pri)))
.manifest <- rbind(.manifest, data.frame(
  filename = basename(.fn_roc_pri), fig_no = "4.24",
  caption = "ROC overlay ML vs scores primary",
  source = "base-R roc_ml_hr + get_score_rocs(PRIMARY)",
  stringsAsFactors = FALSE))

.fn_roc_sec <- file.path(.png_out, "final_fig_4.30_roc_ml_vs_scores_secondary.png")
png(.fn_roc_sec, width = 10*300, height = 7*300, res = 300, bg = "white")
.roc_list_sc_sec <- get_score_rocs(SECONDARY)
par(mar = c(5, 5, 4, 2))
plot(roc_ml_mt, col = COL_BLUE, lwd = 3,
     main = sprintf("BURST vs traditional scores [SECONDARY: %s]", SECONDARY),
     legacy.axes = TRUE,
     xlim = c(1, 0), ylim = c(0, 1), asp = NA,
     xaxs = "i", yaxs = "i")
.leg_lab <- sprintf("BURST  AUC=%.3f", as.numeric(auc(roc_ml_mt)))
.leg_col <- COL_BLUE; .leg_lwd <- 3
.sc_auc_sec <- sapply(.roc_list_sc_sec, function(r) if (is.null(r)) NA_real_ else as.numeric(auc(r)))
for (i in order(.sc_auc_sec, decreasing = TRUE))
  if (!is.null(.roc_list_sc_sec[[i]])) {
    plot(.roc_list_sc_sec[[i]], col = score_cols[i], lwd = 2, legacy.axes = TRUE, add = TRUE)
    .leg_lab <- c(.leg_lab, sprintf("%-12s  AUC=%.3f",
                                    score_names[i],
                                    as.numeric(auc(.roc_list_sc_sec[[i]]))))
    .leg_col <- c(.leg_col, score_cols[i]); .leg_lwd <- c(.leg_lwd, 2)
  }
legend("bottomright", bty = "n", lwd = .leg_lwd, col = .leg_col,
       legend = .leg_lab, cex = 0.9)
dev.off()
cat(sprintf("  saved  %s  -> Fig 4.30 | ROC ML vs scores secondary\n",
            basename(.fn_roc_sec)))
.manifest <- rbind(.manifest, data.frame(
  filename = basename(.fn_roc_sec), fig_no = "4.30",
  caption = "ROC overlay ML vs scores secondary",
  source = "base-R roc_ml_mt + get_score_rocs(SECONDARY)",
  stringsAsFactors = FALSE))

write.csv(.manifest, file.path(.png_out, "final_figure_manifest.csv"),
          row.names = FALSE)
cat(sprintf("\nManifest: %s (%d figures exported)\n",
            file.path(.png_out, "final_figure_manifest.csv"),
            nrow(.manifest)))

pdf("ugib_part2.pdf", width=13, height=9)
for(nm in names(plots)) tryCatch(print(plots[[nm]]), error=function(e) NULL)
par(mar=c(5,5,4,2))
plot(roc_ovr[[1]],col=cols_ovr[1],lwd=2.5,
     main="BURST AUC one-vs-rest: 3-class model",
     legacy.axes = TRUE,
     xlim = c(1, 0), ylim = c(0, 1), asp = NA,
     xaxs = "i", yaxs = "i")
for(i in 2:3) if(!is.null(roc_ovr[[i]]))
  plot(roc_ovr[[i]],col=cols_ovr[i],lwd=2.5,add=T)
legend("bottomright",bty="n",lwd=2.5,col=cols_ovr,
       legend=sapply(seq_along(cl_levels),function(i)
         sprintf("%s AUC=%.3f", c("0_units"="0 units","1to2"="1-2 units","3plus"="3 plus units")[cl_levels[i]],
                 ifelse(is.null(roc_ovr[[i]]),NA,auc(roc_ovr[[i]])))))
par(mar=c(5,5,4,2))
roc_list_sc_pdf <- get_score_rocs(PRIMARY)
leg_lab_pdf <- sprintf("BURST  AUC=%.3f", as.numeric(auc(roc_ml_hr)))
leg_col_pdf <- COL_BLUE; leg_lwd_pdf <- 3
.sc_auc_pdf <- sapply(roc_list_sc_pdf, function(r) if (is.null(r)) NA_real_ else as.numeric(auc(r)))
for (i in order(.sc_auc_pdf, decreasing = TRUE))
  if (!is.null(roc_list_sc_pdf[[i]])) {
    leg_lab_pdf <- c(leg_lab_pdf, sprintf("%-12s  AUC=%.3f",
                                          score_names[i], as.numeric(auc(roc_list_sc_pdf[[i]]))))
    leg_col_pdf <- c(leg_col_pdf, score_cols[i])
    leg_lwd_pdf <- c(leg_lwd_pdf, 2)
  }
plot(roc_ml_hr, col=COL_BLUE, lwd=3,
     main=sprintf("BURST vs traditional scores [PRIMARY: %s]", PRIMARY),
     legacy.axes = TRUE,
     xlim = c(1, 0), ylim = c(0, 1), asp = NA,
     xaxs = "i", yaxs = "i")
for (i in seq_along(roc_list_sc_pdf))
  if (!is.null(roc_list_sc_pdf[[i]]))
    plot(roc_list_sc_pdf[[i]], col=score_cols[i], lwd=2, legacy.axes = TRUE, add=TRUE)
legend("bottomright", bty="n", lwd=leg_lwd_pdf, col=leg_col_pdf,
       legend=leg_lab_pdf, cex=0.9)
dev.off()


cat("\n=== PART 2 COMPLETE ===\n")
cat("Output files:\n")
cat("  ugib_part2.pdf           -- all plots\n")
cat("  thesis_tables.xlsx       -- Table 1, 2, OR, DeLong, NRI, importance\n")
cat("  ed_model_weighted.rds    -- trained RF models\n")
cat("  model_evaluation.rds     -- metrics\n")
cat("  headtohead_results.rds   -- DeLong, NRI, IDI, DCA\n")
cat("  explainability_results.rds -- importance, PDP\n")
cat("  ordinal_results.rds      -- OR table (if MICE)\n")

cat("\n============================================================\n")
cat("APPENDIX: REPRODUCIBILITY INFORMATION (V2.5.8)\n")
cat("============================================================\n\n")
cat(sprintf("R version: %s\n", R.version.string))
cat(sprintf("Platform:  %s\n", R.version$platform))
cat(sprintf("Run date:  %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))

cat("Random seeds used:\n")
cat("  set.seed(2024) -- main analyses (RF training, CV, bootstrap, perm-imp base)\n")
cat("  set.seed(42)   -- ICE plot patient sample, sensitivity analyses\n")
cat("  Per-feature/permutation seeds: 2024 + b*1000 + which(ed_feats == v)\n\n")

cat("Cross-validation: 5 repeats x 5 stratified folds\n")
cat("RF: ntree=500, sampsize-balanced; LR: ridge (alpha=0), 5-fold lambda CV\n\n")

.si <- sessionInfo()
print(.si)
writeLines(capture.output(.si), "sessionInfo_v258.txt")
cat("\nSaved: sessionInfo_v258.txt\n")

cat("\nSource code: eosfeatureselect_v258.R + ugib_part2_v258.R\n")
cat("============================================================\n")
