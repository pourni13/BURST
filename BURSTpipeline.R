# =============================================================================
# BURST -- Bedside UGIB Risk Stratification and Transfusion model
# Random Forest risk model for upper GI bleeding (UGIB)
# One pipeline, three tasks, one shared feature engine:
#   TASK 1  high_risk          (primary, binary)   -> BURST-HighRisk
#   TASK 2  major_tx           (>=3 pRBC, binary)  -> BURST-MajorTx
#   TASK 3  transfusion_class  (0 / 1-2 / 3+ )     -> BURST-Class
# Shared upstream: data prep -> candidate pool -> feature selection.
# Per-task heads: balanced Random Forest (BURST) + ridge-LR comparator.
# =============================================================================
# Pournaras G., Malousi A., Haidich A.B., Kontopoulou T.
# =============================================================================
# STEP 0: CONFIGURATION & PACKAGES

DATA_PATH   <- "Desktop/UGIB_thesis/RstudioWD/UGIBfinal.xlsx"
SEED        <- 2024
NTREE       <- 500L
N_REPEATS   <- 5L
N_FOLDS     <- 5L
B_STABILITY <- 100L
STAB_THRESH <- 0.60
VOTE_THRESH <- 2L
CLASS_WTS3  <- c(1, 3, 3)
RUN_NESTED_CV <- TRUE
RUN_SHAP      <- TRUE

PRIMARY   <- "high_risk"
SECONDARY <- "major_tx"
ORDINAL   <- "transfusion_class"
CLASS_LEVELS <- c("0_units", "1to2", "3plus")

.pkgs <- c("readxl","dplyr","tidyr","rpart","randomForest","glmnet","pROC",
           "fastshap","shapviz")
.miss <- .pkgs[!vapply(.pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(.miss)) install.packages(.miss, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr)
  library(rpart);  library(randomForest); library(glmnet); library(pROC)
})
suppressWarnings(for (loc in c("C.UTF-8","en_US.UTF-8","C")) {
  if (nzchar(tryCatch(Sys.setlocale("LC_ALL", loc), error = function(e) ""))) break
})

# STEP 1: LOAD & CLEAN DATA

clean_colnames <- function(nms) {
  nms <- gsub("\u03bc", "", nms, fixed = TRUE)
  nms <- gsub("\u03b3", "", nms, fixed = TRUE)
  nms <- gsub("\xce\xbc", "", nms, useBytes = TRUE)
  nms <- gsub("\xce\xb3", "", nms, useBytes = TRUE)
  nms <- iconv(nms, from = "UTF-8", to = "ASCII//TRANSLIT", sub = "")
  nms <- gsub("/", "_per_", nms); nms <- gsub("[()]", "", nms)
  nms <- gsub("%", "", nms);      nms <- gsub(" ", "_", nms)
  nms <- gsub("-", "_", nms);     nms <- gsub("_+", "_", nms)
  nms <- gsub("_$", "", nms);     nms <- gsub("^_", "", nms)
  nms
}

load_raw <- function(path = DATA_PATH) {
  df_raw <- read_excel(path, sheet = 1)
  names(df_raw) <- clean_colnames(names(df_raw))
  df_raw |>
    mutate(across(where(is.character),
                  ~ ifelse(. %in% c("N/A","NA",""), NA_character_, .))) |>
    dplyr::select(-any_of("Endoscopy_date_hours"))
}

rename_vars <- function(df) {
  df |> rename(
    patient_id = Number, age = AgeAge, gender = Gender,
    sbp = Systolic_BP_mmHg, dbp = Diastolic_BP_mmHg, hr = Heart_Rate_bpm,
    symptom_duration = Duration_of_symptoms, ed_levosim = ED_Levo,
    hemodynamic_instab = Hemodynamic_Instability,
    altered_mental = Alteration_Mental_Status, syncope = Syncope,
    hematemesis_red = Hematemesis_Bright_Red,
    hematemesis_coffee = Hematemesis_Coffee_Grounds,
    melena = Melena, hematochezia = Hematochezia, rectal_exam = Digital_Rectal_Exam,
    gi_bleed_lt3m = History_GI_Bleeding_under3months,
    gi_bleed_3to6m = History_GI_Bleeding_3to6months,
    gi_bleed_gt6m = History_GI_Bleeding_6to12months,
    gi_bleed_1to5y = History_GI_Bleeding_1to5years,
    gi_bleed_gt5y = History_GI_Bleeding_over5years,
    hb_ed = Hb_ED_mg_per_dl, hct_ed = HCT_ED, lactate_ed = ED_Lac_mmol_per_L,
    hb_entry = Hb_Entry_mg_per_dl, hct_entry = HCT_Entry, rbcs = RBCs, ri_pct = RI,
    wbc_entry = WBC_Entry_x103_per_L, plt_entry = PLT_Entry_x103_per_L,
    creatinine = Creatinine_Entry_mg_per_dl, urea = Urea_Entry_mg_per_dl, bun = BUN,
    sgot = SGOT_IU_per_L, sgpt = SGPT_IU_per_L, ldh = LDH_IU_per_L,
    tbil = TBIL, ferritin = Fer_ng_per_ml, albumin = Albumin_g_per_dl,
    inr = INR, aptt = aPTT, pt_pct = PT,
    drug_act = Drugs_ACTs, drug_antiplt = Drugs_antiPLT, drug_ppi = Drugs_PPI,
    drug_nsaid = Drugs_NSAID, drug_steroid = Drugs_Steroids,
    liver_cirrhosis = Comorbidities_Hepatic_Liver_Cirrhosis,
    heart_failure = Comorbidities_Cardiac_Cardiac_Insufficiency,
    cad = Heart_CAD, afib = Heart_Atrial_Fibrilation,
    hypertension = Comorbidities_Arterial_Hypertension,
    copd = Comorbidities_COPD, asthma = Comorbidities_Pulmonary_Asthma,
    ibd = Comorbidities_GI_IBD, comorbidity_hematol = Comorbidities_Hematological,
    cancer_gi = Comorbidities_Cancer_GI, cancer_extra = Cancer_Extraintestinal,
    cancer_under6m = Cancer_under6months, cancer_over6m = Cancer_over6months,
    cancer_gt5y = Cancer_gt5years,
    renal_ckd = Comorbidities_Renal_CKD, renal_aki = Renal_AKD,
    renal_aki_on_ckd = Comorbidities_Renal_AKD_CKD,
    dm_controlled = Comorbidities_Diabetes_Melitus_Regulated,
    dm_uncontrolled = Comorbidities_Metabolic_Diabetes_Melitus_Dysregulated,
    endo_clips = Workup_Endoscopy_Clips, endo_adrenaline = Workup_Endoscopy_Adrenaline,
    endo_clips_adr = Workup_Endoscopy_Clips_Adrenaline, endo_apc = Workup_Endoscopy_APC,
    forrest_ia = Peptic_Ulcers_Forrest_IA, forrest_ib = Peptic_Ulcers_Forest_IB,
    forrest_iia = Peptic_Ulcers_Forest_IIA, forrest_iib = Peptic_Ulcers_Forest_IIB,
    forrest_iic = Peptic_Ulcers_Forest_IIC, forrest_iii = Peptic_Ulcers_Forrest_III,
    varices = Portal_hypertension_Varices,
    portal_gastropathy = Portal_Hypertension_Portal_Hypertensive_Gastropathy,
    endo_malignancy = Endo_Malignancy, mallory_weiss = Laceration_Mallory_Weiss,
    esoph_laceration = Laceration_oesophagus,
    erosive_gastritis = Erosive_Gastritis, erosive_esophagitis = Erosive_Esophagitis,
    angioectasia = Vascular_Anomalies_Angioectasias,
    dieulafoy = Vascular_Anomalies_Dieulafoy_lesion,
    aortoenteric = Aortoenteric_Fistula, atrophic_gastritis = Atrophic_Gastritis,
    total_wb = Total_Transfusion_Whole_Blood,
    outcome_death = Outcome_Death, rebleeding = Rebleeding_InHosp,
    gbs_score = GBS, aims65_score_xl = AIMS65_Score,
    pre_rockall_xl = Pre_Rockall, full_rockall_xl = Full_Rockall_Computed,
    canuka_xl = CANUKA)
}

coerce_types <- function(df) {
  num_cols <- c("age","sbp","dbp","hr","symptom_duration",
                "hb_ed","hct_ed","lactate_ed","hb_entry","hct_entry","rbcs","ri_pct",
                "wbc_entry","plt_entry","creatinine","urea","bun","sgot","sgpt","ldh",
                "tbil","ferritin","albumin","inr","aptt","pt_pct","total_wb",
                "gbs_score","aims65_score_xl","pre_rockall_xl","full_rockall_xl","canuka_xl")
  yn_cols <- c("ed_levosim","hemodynamic_instab","altered_mental","syncope",
               "hematemesis_red","hematemesis_coffee","melena","hematochezia",
               "gi_bleed_lt3m","gi_bleed_3to6m","gi_bleed_gt6m","gi_bleed_1to5y","gi_bleed_gt5y",
               "drug_act","drug_antiplt","drug_ppi","drug_nsaid","drug_steroid",
               "liver_cirrhosis","heart_failure","cad","afib","hypertension","copd",
               "asthma","ibd","comorbidity_hematol","cancer_gi","cancer_extra",
               "cancer_under6m","cancer_over6m","cancer_gt5y","renal_ckd","renal_aki",
               "renal_aki_on_ckd","dm_controlled","dm_uncontrolled",
               "endo_clips","endo_adrenaline","endo_clips_adr","endo_apc",
               "forrest_ia","forrest_ib","forrest_iia","forrest_iib","forrest_iic","forrest_iii",
               "varices","portal_gastropathy","endo_malignancy","mallory_weiss",
               "esoph_laceration","erosive_gastritis","erosive_esophagitis","angioectasia",
               "dieulafoy","aortoenteric","atrophic_gastritis",
               "outcome_death","rebleeding")
  df <- df |> mutate(across(all_of(num_cols[num_cols %in% names(df)]),
                            ~ suppressWarnings(as.numeric(.))))
  df <- df |> mutate(across(all_of(yn_cols[yn_cols %in% names(df)]),
                            ~ factor(case_when(toupper(.) == "YES" ~ "YES",
                                               toupper(.) == "NO"  ~ "NO",
                                               TRUE ~ NA_character_),
                                     levels = c("NO","YES"))))
  df |> mutate(
    gender = factor(gender, levels = c("F","M")),
    rectal_exam = factor(ifelse(rectal_exam %in% c("NEG","POS"), rectal_exam, NA),
                         levels = c("NEG","POS")))
}

# STEP 2: DERIVE PREDICTORS & OUTCOMES

derive_features <- function(df) {
  df <- df |> mutate(
    hb_mean  = (hb_ed + hb_entry) / 2,
    hct_mean = (hct_ed + hct_entry) / 2,
    gi_bleed_recent = factor(case_when(
      is.na(gi_bleed_lt3m) & is.na(gi_bleed_3to6m) & is.na(gi_bleed_gt6m) ~ NA_character_,
      gi_bleed_lt3m == "YES" | gi_bleed_3to6m == "YES" | gi_bleed_gt6m == "YES" ~ "YES",
      TRUE ~ "NO"), levels = c("NO","YES")),
    gi_bleed_remote = factor(case_when(
      is.na(gi_bleed_1to5y) & is.na(gi_bleed_gt5y) ~ NA_character_,
      gi_bleed_1to5y == "YES" | gi_bleed_gt5y == "YES" ~ "YES",
      TRUE ~ "NO"), levels = c("NO","YES")),
    dm_any = factor(case_when(
      is.na(dm_controlled) & is.na(dm_uncontrolled) ~ NA_character_,
      dm_controlled == "YES" | dm_uncontrolled == "YES" ~ "YES",
      TRUE ~ "NO"), levels = c("NO","YES")),
    renal_severity_3 = factor(case_when(
      is.na(renal_ckd) & is.na(renal_aki) & is.na(renal_aki_on_ckd) ~ NA_character_,
      renal_aki == "YES" | renal_aki_on_ckd == "YES" ~ "Acute",
      renal_ckd == "YES" ~ "Chronic",
      TRUE ~ "None"), levels = c("None","Chronic","Acute")),
    active_cancer = factor(case_when(
      is.na(cancer_gi) & is.na(cancer_extra) ~ NA_character_,
      (is.na(cancer_gi)   | cancer_gi   != "YES") &
        (is.na(cancer_extra)| cancer_extra!= "YES") ~ "NO",
      !is.na(cancer_gt5y) & cancer_gt5y == "YES" &
        !is.na(cancer_under6m) & cancer_under6m == "NO" &
        !is.na(cancer_over6m)  & cancer_over6m  == "NO" ~ "NO",
      !is.na(cancer_under6m) & cancer_under6m == "NO" &
        !is.na(cancer_over6m)  & cancer_over6m  == "NO" &
        !is.na(cancer_gt5y)    & cancer_gt5y    == "NO" ~ "NO",
      TRUE ~ "YES"), levels = c("NO","YES")))
  df |> mutate(
    transfusion_class = factor(case_when(
      total_wb == 0 ~ "0_units", total_wb <= 2 ~ "1to2", TRUE ~ "3plus"),
      levels = CLASS_LEVELS, ordered = TRUE),
    major_tx = factor(ifelse(total_wb >= 3, "YES", "NO"), levels = c("NO","YES")),
    endo_haemostasis = (endo_clips == "YES" | endo_adrenaline == "YES" |
                          endo_clips_adr == "YES" | endo_apc == "YES") & !is.na(endo_clips),
    high_risk = factor(ifelse(outcome_death == "YES" | total_wb > 0 | endo_haemostasis,
                              "YES","NO"), levels = c("NO","YES")),
    lactate_missing = as.integer(is.na(lactate_ed)))
}

# STEP 3: DERIVE TRADITIONAL SCORES (for head-to-head benchmarking)

derive_scores <- function(df) {
  is_yes <- function(x) !is.na(x) & x == "YES"
  active_yes <- !is.na(df$active_cancer) & df$active_cancer == "YES"
  df <- df |> mutate(
    pre_rockall =
      case_when(age < 60 ~ 0L, age < 80 ~ 1L, age >= 80 ~ 2L, TRUE ~ NA_integer_) +
      case_when(sbp >= 100 & hr < 100 ~ 0L, sbp >= 100 & hr >= 100 ~ 1L,
                sbp < 100 ~ 2L, TRUE ~ NA_integer_) +
      case_when(renal_ckd == "YES" | renal_aki == "YES" | renal_aki_on_ckd == "YES" |
                  liver_cirrhosis == "YES" | active_yes ~ 3L,
                heart_failure == "YES" | cad == "YES" ~ 2L, TRUE ~ 0L),
    rockall_diagnosis = case_when(
      mallory_weiss == "YES" ~ 0L, endo_malignancy == "YES" ~ 2L,
      forrest_ia == "YES" | forrest_ib == "YES" | forrest_iia == "YES" |
        forrest_iib == "YES" | forrest_iic == "YES" | forrest_iii == "YES" |
        varices == "YES" | portal_gastropathy == "YES" | esoph_laceration == "YES" |
        erosive_gastritis == "YES" | erosive_esophagitis == "YES" | angioectasia == "YES" |
        dieulafoy == "YES" | aortoenteric == "YES" | atrophic_gastritis == "YES" ~ 1L,
      TRUE ~ 0L),
    rockall_stigmata = case_when(
      forrest_ia == "YES" | forrest_ib == "YES" |
        forrest_iia == "YES" | forrest_iib == "YES" ~ 2L,
      TRUE ~ 0L),
    full_rockall = pre_rockall + rockall_diagnosis + rockall_stigmata,
    aims65_score = case_when(
      is.na(albumin) | is.na(inr) | is.na(sbp) | is.na(age) ~ NA_integer_,
      TRUE ~ as.integer((albumin < 3) + (inr > 1.5) + is_yes(altered_mental) +
                          (sbp <= 90) + (age > 65))),
    hb_for_gbs = ifelse(!is.na(hb_ed) & !is.na(hb_entry), (hb_ed + hb_entry)/2,
                        ifelse(!is.na(hb_entry), hb_entry, hb_ed)),
    bun_mmol = bun * 0.357,
    gbs_bun = case_when(is.na(bun_mmol) ~ NA_integer_, bun_mmol < 6.5 ~ 0L,
                        bun_mmol < 8 ~ 2L, bun_mmol < 10 ~ 3L, bun_mmol < 25 ~ 4L, TRUE ~ 6L),
    gbs_hb = case_when(
      is.na(hb_for_gbs) | is.na(gender) ~ NA_integer_,
      gender == "M" & hb_for_gbs >= 13 ~ 0L, gender == "M" & hb_for_gbs >= 12 ~ 1L,
      gender == "M" & hb_for_gbs >= 10 ~ 3L, gender == "M" ~ 6L,
      gender == "F" & hb_for_gbs >= 12 ~ 0L, gender == "F" & hb_for_gbs >= 10 ~ 1L, TRUE ~ 6L),
    gbs_sbp = case_when(is.na(sbp) ~ NA_integer_, sbp >= 110 ~ 0L, sbp >= 100 ~ 1L,
                        sbp >= 90 ~ 2L, TRUE ~ 3L),
    gbs_hr = as.integer(!is.na(hr) & hr >= 100),
    gbs_melena = as.integer(is_yes(melena)), gbs_syncope = as.integer(is_yes(syncope)) * 2L,
    gbs_liver = as.integer(is_yes(liver_cirrhosis)) * 2L,
    gbs_hf = as.integer(is_yes(heart_failure)) * 2L,
    gbs_computed = gbs_bun + gbs_hb + gbs_sbp + gbs_hr + gbs_melena +
      gbs_syncope + gbs_liver + gbs_hf)
  if ("gbs_score" %in% names(df))
    df$gbs_score <- ifelse(is.na(df$gbs_score), df$gbs_computed, df$gbs_score)
  bun_mmol_int <- as.integer(round(df$bun * 0.357 * 10))
  df <- df |> mutate(
    canuka_age = case_when(is.na(age) ~ NA_integer_, age < 50 ~ 0L, age < 65 ~ 1L, TRUE ~ 2L),
    canuka_melena = ifelse(is.na(melena), NA_integer_, as.integer(is_yes(melena))),
    canuka_hematem = case_when(
      is.na(hematemesis_red) & is.na(hematemesis_coffee) ~ NA_integer_,
      is_yes(hematemesis_red) | is_yes(hematemesis_coffee) ~ 1L, TRUE ~ 0L),
    canuka_syncope = ifelse(is.na(syncope), NA_integer_, as.integer(is_yes(syncope))),
    canuka_hb = case_when(is.na(hb_mean) ~ NA_integer_, hb_mean > 12 ~ 0L,
                          hb_mean > 10 ~ 1L, hb_mean > 8 ~ 2L, TRUE ~ 3L),
    canuka_bun = case_when(is.na(bun) ~ NA_integer_, bun_mmol_int < 50 ~ 0L,
                           bun_mmol_int < 100 ~ 1L, bun_mmol_int < 150 ~ 2L, TRUE ~ 3L),
    canuka_liver = ifelse(is_yes(liver_cirrhosis), 2L, 0L),
    canuka_malig = ifelse(active_yes, 2L, 0L),
    canuka_hr = case_when(is.na(hr) ~ NA_integer_, hr < 100 ~ 0L, hr < 125 ~ 1L, TRUE ~ 2L),
    canuka_sbp = case_when(is.na(sbp) ~ NA_integer_, sbp >= 120 ~ 0L, sbp >= 100 ~ 1L,
                           sbp >= 80 ~ 2L, TRUE ~ 3L))
  canuka_cols <- c("canuka_age","canuka_melena","canuka_hematem","canuka_syncope",
                   "canuka_hb","canuka_bun","canuka_liver","canuka_malig","canuka_hr","canuka_sbp")
  df$canuka <- rowSums(df[, canuka_cols])
  df
}

prepare_data <- function(path = DATA_PATH) {
  load_raw(path) |> rename_vars() |> coerce_types() |> derive_features() |> derive_scores()
}

# STEP 4: CANDIDATE POOL & SIMPLE IMPUTATION

candidate_pool <- function(df) {
  pool <- c("age","gender","symptom_duration","sbp","dbp","hr","hemodynamic_instab",
            "ed_levosim","altered_mental","syncope","hematemesis_red","hematemesis_coffee",
            "melena","hematochezia","rectal_exam","gi_bleed_recent","gi_bleed_remote",
            "hb_mean","hct_mean","lactate_ed","wbc_entry","plt_entry","ri_pct","creatinine",
            "urea","bun","albumin","inr","aptt","pt_pct","sgot","sgpt","tbil","ferritin",
            "ldh","drug_act","drug_antiplt","drug_ppi","drug_nsaid","drug_steroid",
            "liver_cirrhosis","renal_severity_3","heart_failure","cad","afib","hypertension",
            "copd","asthma","dm_any","active_cancer","comorbidity_hematol","ibd")
  pool <- pool[pool %in% names(df)]
  pool <- pool[!pool %in% c("hct_mean","bun","hemodynamic_instab","ed_levosim")]
  c(pool, "lactate_missing")
}

impute_simple <- function(df, vars) {
  for (v in vars) {
    x <- df[[v]]
    if (is.numeric(x)) {
      x[is.na(x)] <- median(x, na.rm = TRUE)
    } else {
      md <- names(sort(table(x), decreasing = TRUE))[1]
      if (!is.null(md)) x[is.na(x)] <- md
    }
    df[[v]] <- x
  }
  df
}

# STEP 5: FEATURE SELECTION (4-METHOD VOTE -> STABILITY -> TIER FORWARD)

vars_by_type <- function(df, cands) list(
  num = cands[vapply(cands, function(v) is.numeric(df[[v]]), logical(1))],
  cat = cands[vapply(cands, function(v) is.factor(df[[v]]) || is.character(df[[v]]),
                     logical(1))])

univ_nominal <- function(df, cands, outcome, fisher_B = 10000L) {
  vt  <- vars_by_type(df, cands); out <- character(0); mt <- df[[outcome]]
  for (v in vt$num) {
    x <- df[[v]]; ok <- !is.na(x) & !is.na(mt)
    if (sum(ok) < 10 || length(unique(mt[ok])) < 2) next
    p <- tryCatch(suppressWarnings(
      wilcox.test(x[ok & mt == "YES"], x[ok & mt == "NO"])$p.value),
      error = function(e) NA)
    if (!is.na(p) && p < .05) out <- c(out, v)
  }
  for (v in vt$cat) {
    x <- df[[v]]; ok <- !is.na(x) & !is.na(mt)
    if (sum(ok) < 10) next
    tb <- table(x[ok], mt[ok]); if (nrow(tb) < 2 || ncol(tb) < 2) next
    p <- tryCatch(suppressWarnings(
      if (any(chisq.test(tb)$expected < 5))
        fisher.test(tb, simulate.p.value = TRUE, B = fisher_B)$p.value
      else chisq.test(tb)$p.value), error = function(e) NA)
    if (!is.na(p) && p < .05) out <- c(out, v)
  }
  out
}

tree_top <- function(df, cands, outcome) {
  set.seed(SEED)
  tr <- tryCatch(rpart(reformulate(cands, outcome), data = df, method = "class",
                       control = rpart.control(cp = 0, minsplit = 10, maxdepth = 5, xval = 10)),
                 error = function(e) NULL)
  if (is.null(tr)) return(character(0))
  cpt <- tr$cptable; im <- which.min(cpt[, "xerror"])
  th  <- cpt[im, "xerror"] + cpt[im, "xstd"]
  cp1 <- cpt[which(cpt[, "xerror"] <= th)[1], "CP"]
  tf  <- prune(tr, cp = cp1)
  sp  <- unique(tf$frame$var[tf$frame$var != "<leaf>"])
  if (length(sp) < 3 || is.null(tf$variable.importance)) return(character(0))
  names(sort(tf$variable.importance, decreasing = TRUE))[1:min(15, length(tf$variable.importance))]
}

rf_top15 <- function(df, cands, outcome, ntree = NTREE) {
  set.seed(SEED)
  rf <- tryCatch(randomForest(reformulate(cands, outcome), data = df, ntree = ntree,
                              mtry = floor(sqrt(length(cands))), importance = FALSE),
                 error = function(e) NULL)
  if (is.null(rf)) return(character(0))
  g <- importance(rf)[, "MeanDecreaseGini"]
  names(sort(g, decreasing = TRUE))[1:min(15, length(g))]
}

lasso_1se <- function(df, cands, outcome) {
  X <- tryCatch(model.matrix(reformulate(cands), data = df)[, -1, drop = FALSE],
                error = function(e) NULL)
  y <- as.numeric(df[[outcome]] == "YES")
  if (is.null(X) || length(unique(y)) < 2 || nrow(X) < 20) return(character(0))
  set.seed(SEED)
  cv <- tryCatch(cv.glmnet(X, y, family = "binomial", alpha = 1, nfolds = 5,
                           type.measure = "deviance"), error = function(e) NULL)
  if (is.null(cv)) return(character(0))
  cf <- coef(cv, s = "lambda.1se"); tm <- setdiff(rownames(cf)[as.numeric(cf) != 0], "(Intercept)")
  unique(unlist(lapply(cands, function(v) if (any(grepl(paste0("^", v), tm))) v)))
}

balanced_rf_roc <- function(df, fs, outcome) {
  df[[outcome]] <- factor(as.character(df[[outcome]]), levels = c("NO","YES"))
  mn <- min(table(df[[outcome]]))
  set.seed(SEED)
  rf <- randomForest(reformulate(fs, outcome), data = df, ntree = NTREE,
                     strata = df[[outcome]], sampsize = c("NO" = mn, "YES" = mn))
  p  <- rf$votes[, "YES"] / rowSums(rf$votes)
  r  <- roc(df[[outcome]], p, levels = c("NO","YES"), direction = "<", quiet = TRUE)
  list(rf = rf, prob = p, roc = r, auc = as.numeric(auc(r)), ci = as.numeric(ci.auc(r)))
}

vote_confirm <- function(df_imp, df_univ, cands, outcome) {
  u <- univ_nominal(df_univ, cands, outcome, fisher_B = 10000L)
  t <- tree_top(df_imp, cands, outcome)
  r <- rf_top15(df_imp, cands, outcome)
  l <- lasso_1se(df_imp, cands, outcome)
  allc  <- unique(c(u, t, r, l))
  votes <- vapply(allc, function(v) (v %in% u) + (v %in% t) + (v %in% r) + (v %in% l), integer(1))
  confirmed <- sort(allc[votes >= VOTE_THRESH])
  if (!"hb_mean" %in% confirmed) confirmed <- sort(c(confirmed, "hb_mean"))
  confirmed
}

stability_select <- function(df_imp, cands, outcome, B = B_STABILITY) {
  df_work <- df_imp[!is.na(df_imp[[outcome]]), , drop = FALSE]
  sv <- matrix(0, nrow = length(cands), ncol = B, dimnames = list(cands, NULL))
  n  <- nrow(df_work)
  set.seed(SEED)
  for (b in seq_len(B)) {
    db <- df_work[sample.int(n, n, replace = TRUE), , drop = FALSE]
    if (length(unique(db[[outcome]])) < 2) next
    ub <- univ_nominal(db, cands, outcome, fisher_B = 2000L)
    rb <- rf_top15(db, cands, outcome, ntree = 300)
    lb <- lasso_1se(db, cands, outcome)
    for (v in unique(c(ub, rb, lb)))
      if (((v %in% ub) + (v %in% rb) + (v %in% lb)) >= 2 && v %in% rownames(sv))
        sv[v, b] <- 1
  }
  sel_prob <- rowSums(sv) / B
  names(sel_prob)[sel_prob >= STAB_THRESH]
}

forward_select <- function(df_imp, confirmed, stable, outcome) {
  tier1 <- "hb_mean"
  tier2 <- setdiff(stable, tier1)
  tier3 <- setdiff(confirmed, c(tier1, tier2))
  base  <- unique(c(tier1, tier2))
  base_fit <- balanced_rf_roc(df_imp, base, outcome)
  retained <- character(0)
  for (cand in tier3) {
    test_fit <- balanced_rf_roc(df_imp, c(base, cand), outcome)
    dp <- tryCatch(roc.test(base_fit$roc, test_fit$roc, method = "delong",
                            paired = TRUE)$p.value, error = function(e) NA_real_)
    if (test_fit$auc > base_fit$auc &&
        ((!is.na(dp) && dp < 0.05) || test_fit$ci[1] > base_fit$ci[3]))
      retained <- c(retained, cand)
  }
  sort(unique(c(base, retained)))
}

select_features <- function(df_imp, cands, outcome, df_univ = df_imp) {
  confirmed <- vote_confirm(df_imp, df_univ, cands, outcome)
  stable    <- stability_select(df_imp, cands, outcome)
  if (length(stable) == 0) stable <- "hb_mean"
  list(confirmed = confirmed, stable = stable,
       final = forward_select(df_imp, confirmed, stable, outcome))
}

nzv_filter <- function(df, feats) {
  nzv <- sapply(feats, function(v) {
    x <- df[[v]]
    if (is.numeric(x)) return(var(x, na.rm = TRUE) < 1e-4)
    tbl <- table(x, useNA = "no"); length(tbl) < 2 || max(tbl) / sum(tbl) > .97
  })
  feats[!nzv]
}

# STEP 6: BURST TRAINING HEADS (BINARY RF + ORDINAL 3-CLASS RF + RIDGE-LR)

train_burst_binary <- function(df, feats, outcome) balanced_rf_roc(df, feats, outcome)

train_burst_class <- function(df, feats, outcome = ORDINAL) {
  df[[outcome]] <- factor(as.character(df[[outcome]]), levels = CLASS_LEVELS, ordered = TRUE)
  cl <- levels(df[[outcome]]); n_min <- min(table(df[[outcome]]))
  set.seed(SEED)
  randomForest(reformulate(feats, outcome), data = df, ntree = NTREE,
               mtry = floor(sqrt(length(feats))), strata = df[[outcome]],
               sampsize = setNames(rep(n_min, length(cl)), cl),
               classwt = setNames(CLASS_WTS3, cl))
}

fit_lr <- function(df, feats, outcome) {
  X <- as.data.frame(lapply(df[, feats, drop = FALSE],
                            function(v) if (is.numeric(v)) v else as.numeric(v == "YES")))
  y <- as.numeric(df[[outcome]] == "YES")
  meds <- vapply(X, function(c) median(c, na.rm = TRUE), numeric(1))
  for (j in seq_along(X)) X[[j]][is.na(X[[j]])] <- meds[[j]]
  mu <- vapply(X, mean, numeric(1)); sg <- vapply(X, sd, numeric(1)); sg[is.na(sg) | sg == 0] <- 1
  Xs <- as.matrix(sweep(sweep(X, 2, mu, "-"), 2, sg, "/"))
  n0 <- sum(y == 0); n1 <- sum(y == 1)
  w  <- ifelse(y == 1, length(y) / (2 * max(n1, 1)), length(y) / (2 * max(n0, 1)))
  fit <- tryCatch(suppressWarnings(cv.glmnet(Xs, y, family = "binomial", alpha = 0,
                                             weights = w, nfolds = 5, type.measure = "auc")),
                  error = function(e) glmnet(Xs, y, family = "binomial", alpha = 0,
                                             weights = w, lambda = 0.01))
  list(fit = fit, meds = meds, mu = mu, sg = sg, is_cv = inherits(fit, "cv.glmnet"))
}

predict_lr <- function(model, df, feats) {
  Xn <- as.data.frame(lapply(df[, feats, drop = FALSE],
                             function(v) if (is.numeric(v)) v else as.numeric(v == "YES")))
  for (j in seq_along(Xn)) Xn[[j]][is.na(Xn[[j]])] <- model$meds[[j]]
  Xs <- as.matrix(sweep(sweep(Xn, 2, model$mu, "-"), 2, model$sg, "/"))
  if (model$is_cv) as.numeric(predict(model$fit, newx = Xs, s = "lambda.min", type = "response"))
  else            as.numeric(predict(model$fit, newx = Xs, type = "response"))
}

# STEP 7: CROSS-VALIDATION (BINARY 5x5 / MULTICLASS / NESTED SELECTION-INCLUSIVE)

make_strat_folds <- function(y, k, seed) {
  set.seed(seed); folds <- integer(length(y))
  for (cls in unique(y)) {
    idx <- sample(which(y == cls))
    folds[idx] <- cut(seq_along(idx), breaks = k, labels = FALSE)
  }
  folds
}

cv_binary <- function(df, feats, outcome, model_type = "RF",
                      n_repeats = N_REPEATS, n_folds = N_FOLDS) {
  y_full <- as.character(df[[outcome]]); aucs <- numeric(n_repeats)
  oof_all <- matrix(NA_real_, nrow = n_repeats, ncol = length(y_full))
  for (rep in seq_len(n_repeats)) {
    fld <- make_strat_folds(y_full, n_folds, SEED + rep); oof <- numeric(length(y_full))
    for (k in seq_len(n_folds)) {
      tr <- which(fld != k); te <- which(fld == k)
      if (length(unique(y_full[tr])) < 2) { oof[te] <- mean(y_full[tr] == "YES"); next }
      oof[te] <- tryCatch({
        if (model_type == "RF") {
          n_min <- min(table(df[[outcome]][tr]))
          fit <- randomForest(reformulate(feats, outcome), data = df[tr, ], ntree = NTREE,
                              strata = df[[outcome]][tr], sampsize = c("NO" = n_min, "YES" = n_min))
          predict(fit, df[te, ], type = "prob")[, "YES"]
        } else {
          predict_lr(fit_lr(df[tr, ], feats, outcome), df[te, ], feats)
        }
      }, error = function(e) rep(mean(y_full[tr] == "YES"), length(te)))
    }
    aucs[rep] <- tryCatch(as.numeric(auc(roc(y_full, oof, levels = c("NO","YES"),
                                             direction = "<", quiet = TRUE))), error = function(e) NA_real_)
    oof_all[rep, ] <- oof
  }
  list(aucs = aucs, oof = colMeans(oof_all), y = y_full,
       mean = mean(aucs, na.rm = TRUE), sd = sd(aucs, na.rm = TRUE))
}

cv_multiclass_acc <- function(df, feats, outcome = ORDINAL,
                              n_repeats = N_REPEATS, n_folds = N_FOLDS) {
  y_full <- as.character(df[[outcome]]); acc <- numeric(n_repeats)
  for (rep in seq_len(n_repeats)) {
    fld <- make_strat_folds(y_full, n_folds, SEED + rep); correct <- 0L; total <- 0L
    for (k in seq_len(n_folds)) {
      tr <- which(fld != k); te <- which(fld == k)
      tr_df <- df[tr, , drop = FALSE]
      if (length(unique(tr_df[[outcome]])) < 2) next
      n_min <- min(table(tr_df[[outcome]])); cl <- levels(tr_df[[outcome]])
      pred <- tryCatch({
        fit <- randomForest(reformulate(feats, outcome), data = tr_df, ntree = 300,
                            mtry = floor(sqrt(length(feats))), strata = tr_df[[outcome]],
                            sampsize = setNames(rep(n_min, length(cl)), cl),
                            classwt = setNames(CLASS_WTS3, cl))
        as.character(predict(fit, df[te, ]))
      }, error = function(e) rep(NA_character_, length(te)))
      ok <- !is.na(pred) & !is.na(y_full[te])
      correct <- correct + sum(pred[ok] == y_full[te][ok]); total <- total + sum(ok)
    }
    acc[rep] <- if (total > 0) correct / total else NA_real_
  }
  list(per_rep = acc, mean = mean(acc, na.rm = TRUE), sd = sd(acc, na.rm = TRUE))
}

nested_cv <- function(df_full, cands, outcome,
                      n_repeats = N_REPEATS, n_folds = N_FOLDS) {
  df_n <- df_full[, c(cands, outcome), drop = FALSE]
  df_n <- df_n[!is.na(df_n[[outcome]]), , drop = FALSE]
  df_n[[outcome]] <- factor(as.character(df_n[[outcome]]), levels = c("NO","YES"))
  y <- as.character(df_n[[outcome]]); oof_all <- matrix(NA_real_, n_repeats, length(y))
  rep_auc <- numeric(n_repeats)
  for (rp in seq_len(n_repeats)) {
    fld <- make_strat_folds(y, n_folds, SEED + rp); oof <- rep(NA_real_, length(y))
    for (k in seq_len(n_folds)) {
      tr <- which(fld != k); te <- which(fld == k)
      dtr <- impute_simple(df_n[tr, , drop = FALSE], cands)
      dte <- impute_simple(df_n[te, , drop = FALSE], cands)
      sel <- tryCatch(select_features(dtr, cands, outcome, df_univ = dtr)$final,
                      error = function(e) "hb_mean")
      mn  <- min(table(dtr[[outcome]]))
      oof[te] <- tryCatch({
        set.seed(SEED)
        rf <- randomForest(reformulate(sel, outcome), data = dtr, ntree = NTREE,
                           strata = dtr[[outcome]], sampsize = c("NO" = mn, "YES" = mn))
        predict(rf, dte, type = "prob")[, "YES"]
      }, error = function(e) rep(mean(dtr[[outcome]] == "YES"), length(te)))
    }
    oof_all[rp, ] <- oof
    rep_auc[rp] <- tryCatch(as.numeric(auc(roc(y, oof, levels = c("NO","YES"),
                                               direction = "<", quiet = TRUE))), error = function(e) NA_real_)
  }
  pooled <- roc(y, colMeans(oof_all, na.rm = TRUE), levels = c("NO","YES"),
                direction = "<", quiet = TRUE)
  list(rep_auc = rep_auc, mean = mean(rep_auc, na.rm = TRUE),
       pooled_auc = as.numeric(auc(pooled)), pooled_ci = as.numeric(ci.auc(pooled)))
}

# STEP 8: EVALUATION (BINARY DISCRIMINATION/CALIBRATION; ORDINAL CONFUSION/AUC)

calibration_stats <- function(y, p) {
  yb <- as.numeric(y == "YES"); pc <- pmin(pmax(p, 1e-6), 1 - 1e-6); lp <- log(pc / (1 - pc))
  intercept <- coef(glm(yb ~ offset(lp), family = binomial))[1]
  slope     <- coef(glm(yb ~ lp, family = binomial))[2]
  brier <- mean((p - yb)^2); ref <- mean(yb) * (1 - mean(yb))
  list(intercept = as.numeric(intercept), slope = as.numeric(slope),
       brier = brier, brier_skill = 1 - brier / ref)
}

operating_points <- function(roc_obj, y, p) {
  yo <- coords(roc_obj, x = "best", best.method = "youden",
               ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
  list(youden = yo, sens100_threshold = min(p[y == "YES"]) - 1e-6)
}

evaluate_binary <- function(fit, df, outcome) {
  list(auc = fit$auc, ci = fit$ci,
       calibration = calibration_stats(df[[outcome]], fit$prob),
       operating = operating_points(fit$roc, as.character(df[[outcome]]), fit$prob))
}

evaluate_class <- function(rf, df, outcome = ORDINAL) {
  pred <- rf$predicted; ok <- !is.na(pred)
  truth <- df[[outcome]][ok]; pred <- pred[ok]
  conf <- table(Predicted = pred, Actual = truth)
  acc  <- mean(as.character(pred) == as.character(truth))
  po <- acc; pe <- sum(rowSums(conf) * colSums(conf)) / sum(conf)^2
  kappa <- (po - pe) / (1 - pe)
  probs <- rf$votes / rowSums(rf$votes)
  ovr <- lapply(CLASS_LEVELS, function(cls) {
    bin <- factor(ifelse(as.character(df[[outcome]]) == cls, "YES", "NO"), levels = c("NO","YES"))
    keep <- !is.na(bin) & !is.na(probs[, cls])
    if (length(unique(bin[keep])) < 2) return(NA_real_)
    as.numeric(auc(roc(bin[keep], probs[keep, cls], levels = c("NO","YES"),
                       direction = "<", quiet = TRUE)))
  })
  list(confusion = conf, accuracy = acc, kappa = kappa,
       ovr_auc = setNames(unlist(ovr), CLASS_LEVELS))
}

# STEP 9: HEAD-TO-HEAD VS TRADITIONAL SCORES (binary tasks)

head_to_head <- function(df, outcome, ml_prob) {
  scores <- c(GBS = "gbs_computed", AIMS65 = "aims65_score",
              `Pre-Rockall` = "pre_rockall", `Full Rockall` = "full_rockall",
              CANUKA = "canuka")
  y <- df[[outcome]]
  roc_ml <- roc(y, ml_prob, levels = c("NO","YES"), direction = "<", quiet = TRUE)
  out <- lapply(names(scores), function(nm) {
    col <- scores[[nm]]; if (!col %in% names(df)) return(NULL)
    sc <- df[[col]]; ok <- !is.na(sc) & !is.na(y)
    if (sum(ok) < 10 || length(unique(y[ok])) < 2) return(NULL)
    r  <- roc(y[ok], sc[ok], levels = c("NO","YES"), direction = "<", quiet = TRUE)
    dp <- if (length(roc_ml$response) == length(r$response))
      tryCatch(roc.test(roc_ml, r, method = "delong")$p.value, error = function(e) NA_real_)
    else
      tryCatch(roc.test(roc_ml, r, method = "bootstrap", boot.n = 1000)$p.value,
               error = function(e) NA_real_)
    data.frame(score = nm, auc_ml = as.numeric(auc(roc_ml)),
               auc_score = as.numeric(auc(r)),
               delta = as.numeric(auc(roc_ml)) - as.numeric(auc(r)), delong_p = dp)
  })
  do.call(rbind, Filter(Negate(is.null), out))
}

# STEP 10: EXPLAINABILITY (SHAP per binary task)

shap_explain <- function(rf, df, feats, nsim = 50) {
  if (!requireNamespace("fastshap", quietly = TRUE)) return(NULL)
  X <- df[, feats, drop = FALSE]
  pred_wrapper <- function(object, newdata) predict(object, newdata, type = "prob")[, "YES"]
  set.seed(SEED)
  sh <- fastshap::explain(object = rf, X = X, pred_wrapper = pred_wrapper,
                          nsim = nsim, parallel = FALSE)
  list(shap = sh, global = sort(colMeans(abs(as.matrix(sh))), decreasing = TRUE))
}

# STEP 11: FINALIZE & PREDICT (save BURST models, score new patients)

train_fill <- function(df, feats) {
  setNames(lapply(feats, function(v) {
    x <- df[[v]]
    if (is.numeric(x))
      list(type = "numeric", value = median(x, na.rm = TRUE))
    else
      list(type = "factor", levels = levels(as.factor(x)),
           value = names(sort(table(x), decreasing = TRUE))[1])
  }), feats)
}

apply_fill <- function(newdata, fill) {
  for (v in names(fill)) {
    if (!v %in% names(newdata)) next
    sp <- fill[[v]]
    if (sp$type == "numeric") {
      x <- suppressWarnings(as.numeric(newdata[[v]])); x[is.na(x)] <- sp$value
      newdata[[v]] <- x
    } else {
      x <- as.character(newdata[[v]]); x[is.na(x)] <- sp$value
      newdata[[v]] <- factor(x, levels = sp$levels)
    }
  }
  newdata
}

finalize_burst <- function(df, feats, model, outcome, kind, path) {
  obj <- list(model = model, features = feats, outcome = outcome, kind = kind,
              fill = train_fill(df, feats))
  saveRDS(obj, path)
  obj
}

predict_burst <- function(burst, newdata) {
  X <- apply_fill(newdata[, burst$features, drop = FALSE], burst$fill)
  if (burst$kind == "binary") predict(burst$model, X, type = "prob")[, "YES"]
  else list(class = predict(burst$model, X),
            probs = predict(burst$model, X, type = "prob"))
}

# STEP 12: RUN PIPELINE (all three tasks)

run_binary_task <- function(df_imp, cands, feats, outcome, label, nested = RUN_NESTED_CV) {
  cat(sprintf("\n=== %s (binary: %s) | features: %s ===\n",
              label, outcome, paste(feats, collapse = ", ")))
  df_task <- droplevels(df_imp[!is.na(df_imp[[outcome]]), , drop = FALSE])
  df_task[[outcome]] <- factor(as.character(df_task[[outcome]]), levels = c("NO","YES"))
  fit   <- train_burst_binary(df_task, feats, outcome)
  cv_rf <- cv_binary(df_task, feats, outcome, "RF")
  cv_lr <- cv_binary(df_task, feats, outcome, "LR")
  nest  <- if (nested) nested_cv(df_imp, cands, outcome) else NULL
  eval  <- evaluate_binary(fit, df_task, outcome)
  h2h   <- head_to_head(df_task, outcome, fit$prob)
  shap  <- if (RUN_SHAP) shap_explain(fit$rf, df_task, feats) else NULL
  burst <- finalize_burst(df_task, feats, fit$rf, outcome, "binary",
                          paste0("burst_", outcome, ".rds"))
  cat(sprintf("  OOB AUC %.3f (%.3f-%.3f) | 5x5 CV: RF %.3f (sd %.3f), LR %.3f (sd %.3f)\n",
              fit$auc, fit$ci[1], fit$ci[3], cv_rf$mean, cv_rf$sd, cv_lr$mean, cv_lr$sd))
  cat(sprintf("  calibration slope %.2f | Brier %.3f (skill %.3f)\n",
              eval$calibration$slope, eval$calibration$brier, eval$calibration$brier_skill))
  if (!is.null(nest))
    cat(sprintf("  nested CV AUC %.3f | selection optimism %+.3f\n",
                nest$mean, fit$auc - nest$mean))
  list(fit = fit, cv_rf = cv_rf, cv_lr = cv_lr, nested = nest,
       evaluation = eval, head_to_head = h2h, shap = shap, burst = burst)
}

run_ordinal_task <- function(df_imp, feats, label) {
  cat(sprintf("\n=== %s (3-class: %s) | features: %s ===\n",
              label, ORDINAL, paste(feats, collapse = ", ")))
  df_task <- droplevels(df_imp[!is.na(df_imp[[ORDINAL]]), , drop = FALSE])
  df_task[[ORDINAL]] <- factor(as.character(df_task[[ORDINAL]]),
                               levels = CLASS_LEVELS, ordered = TRUE)
  rf   <- train_burst_class(df_task, feats, ORDINAL)
  cv   <- cv_multiclass_acc(df_task, feats, ORDINAL)
  eval <- evaluate_class(rf, df_task, ORDINAL)
  burst <- finalize_burst(df_task, feats, rf, ORDINAL, "class",
                          paste0("burst_", ORDINAL, ".rds"))
  cat(sprintf("  OOB accuracy %.1f%% | kappa %.3f | 5x5 CV accuracy %.1f%% (sd %.1f%%)\n",
              eval$accuracy * 100, eval$kappa, cv$mean * 100, cv$sd * 100))
  cat(sprintf("  one-vs-rest AUC: %s\n",
              paste(sprintf("%s=%.3f", names(eval$ovr_auc), eval$ovr_auc), collapse = " ")))
  list(rf = rf, cv = cv, evaluation = eval, burst = burst)
}

run_burst_pipeline <- function(path = DATA_PATH) {
  df     <- prepare_data(path)
  cands  <- candidate_pool(df)
  df_imp <- impute_simple(df, cands)
  cat(sprintf("Cohort: n=%d | candidate features: %d\n", nrow(df), length(cands)))
  
  sel_hr  <- select_features(df_imp, cands, PRIMARY,   df_univ = df)
  sel_mt  <- select_features(df_imp, cands, SECONDARY, df_univ = df)
  feats_hr <- nzv_filter(df_imp, sel_hr$final)
  feats_mt <- nzv_filter(df_imp, sel_mt$final)
  
  task1 <- run_binary_task(df_imp, cands, feats_hr, PRIMARY,   "TASK 1 BURST-HighRisk", nested = RUN_NESTED_CV)
  task2 <- run_binary_task(df_imp, cands, feats_mt, SECONDARY, "TASK 2 BURST-MajorTx",  nested = FALSE)
  task3 <- run_ordinal_task(df_imp, feats_hr, "TASK 3 BURST-Class")
  
  list(data = df, candidates = cands,
       features = list(high_risk = feats_hr, major_tx = feats_mt),
       task1_high_risk = task1, task2_major_tx = task2, task3_transfusion_class = task3)
}

if (sys.nframe() == 0 && identical(environment(), globalenv())) {
  burst <- run_burst_pipeline()
}

# STEP 13: USING BURST FOR PREDICTION ON A NEW PATIENT
#
# Running this script trains and saves one model per task (finalize_burst):
#   burst_high_risk.rds          TASK 1  P(high-risk composite)
#   burst_major_tx.rds           TASK 2  P(>=3 units pRBC)
#   burst_transfusion_class.rds  TASK 3  ordinal class (0 / 1-2 / 3+ units)
# Each object holds: $model (Random Forest), $features (the 6-feature contract),
# $fill (train-set median/mode used to impute missing inputs), $outcome, $kind.
#
# To score a NEW first-hour ED patient, build a one-row data.frame with the six
# features (NA is allowed and imputed from $fill) and call predict_burst():
#   binary task -> numeric P(YES); ordinal task -> list(class, probs).
# hb_mean is (hb_ed + hb_entry)/2; rectal_exam is a factor with levels NEG/POS.

burst_score <- function(newdata, task = c("high_risk","major_tx","transfusion_class")) {
  task  <- match.arg(task)
  obj   <- readRDS(sprintf("burst_%s.rds", task))
  out   <- predict_burst(obj, newdata)
  if (obj$kind == "binary")
    list(prob = as.numeric(out),
         band = cut(as.numeric(out), c(-Inf, .30, .50, .70, Inf),
                    c("safe-discharge","standard-care","escalation","critical")))
  else out
}

# Example (high-risk triage):
#   new_patient <- data.frame(
#     hb_mean     = 8.4,
#     urea        = 95,
#     albumin     = 3.1,
#     lactate_ed  = 3.2,
#     ri_pct      = 1.4,
#     rectal_exam = factor("POS", levels = c("NEG","POS")))
#   res <- burst_score(new_patient, "high_risk")
#   cat(sprintf("BURST P(high-risk) = %.3f -> %s\n", res$prob, res$band))
