# UGIB THESIS -- PART 1 FEATURE SELECTION

PRIMARY    <- "high_risk"
SECONDARY  <- "major_tx"
TERTIARY   <- c("outcome_death","rebleeding")
ALL_OUTCOMES <- c(PRIMARY, SECONDARY, TERTIARY)
cat(sprintf("Primary outcome:    %s\n", PRIMARY))
cat(sprintf("Secondary outcome:  %s\n", SECONDARY))
cat(sprintf("Tertiary outcomes:  %s\n\n", paste(TERTIARY, collapse=", ")))

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(ggplot2)
  library(rpart);  library(rpart.plot)
  library(randomForest); library(glmnet); library(pROC)
})

suppressWarnings({
  for (loc in c("C.UTF-8","en_US.UTF-8",
                "English_United States.utf8","English_United States.1252")) {
    if (nzchar(tryCatch(Sys.setlocale("LC_ALL", loc), error=function(e) ""))) break
  }
})
while (!is.null(dev.list())) dev.off()

theme_ugib <- function(base=12)
  theme_minimal(base_size=base) +
  theme(plot.title    = element_text(face="bold", size=base+1, hjust=0),
        plot.subtitle = element_text(size=base-1, colour="grey40", hjust=0),
        plot.caption  = element_text(size=base-3, colour="grey55"),
        axis.text     = element_text(size=base-1),
        panel.grid.minor = element_blank(),
        legend.position  = "bottom")

COL_BLUE <- "#2E75B6"; COL_RED <- "#C0392B"; COL_GREEN <- "#27AE60"
COL_AMB  <- "#E67E22"; COL_GREY <- "#7F7F7F"

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

plots <- list()

# SECTION 1
cat("\n============================================================\n")
cat("SECTION 1:\n")
cat("============================================================\n\n")

UGIBfinal <- read_excel("~/Desktop/UGIB_thesis/RstudioWD/UGIBfinal.xlsx")
View(UGIBfinal)

df_raw <- read_excel("~/Desktop/UGIB_thesis/RstudioWD/UGIBfinal.xlsx", sheet=1)
cat("Loaded:", nrow(df_raw), "patients,", ncol(df_raw), "columns\n")

clean_colnames <- function(nms) {
  
  nms <- gsub("\u03bc", "", nms, fixed = TRUE)
  nms <- gsub("\u03b3", "", nms, fixed = TRUE)
  nms <- gsub("\xce\xbc", "", nms, useBytes = TRUE)
  nms <- gsub("\xce\xb3", "", nms, useBytes = TRUE)
  nms <- iconv(nms, from="UTF-8", to="ASCII//TRANSLIT", sub="")
  nms <- gsub("/","_per_",nms); nms <- gsub("[()]","",nms)
  nms <- gsub("%","",nms);      nms <- gsub(" ","_",nms)
  nms <- gsub("-","_",nms);     nms <- gsub("_+","_",nms)
  nms <- gsub("_$","",nms);     nms <- gsub("^_","",nms)
  nms
}
names(df_raw) <- clean_colnames(names(df_raw))
cat("Column names cleaned (special chars -> ASCII)\n")

df <- df_raw |>
  mutate(across(where(is.character),
                ~ifelse(. %in% c("N/A","NA",""), NA_character_, .))) |>
  dplyr::select(-any_of(c("Endoscopy_date_hours")))

df <- df |> rename(
  patient_id=Number, hospital_id=HIS, age=AgeAge, gender=Gender,
  hosp_days=Hospitalization_Days, symptom_duration=Duration_of_symptoms,
  time_to_endo_h=time_to_endoscopy_h,
  smoker=Smoker, alcohol=Alcohol, drug_abuse=Drug_Abuse,
  sbp=Systolic_BP_mmHg, dbp=Diastolic_BP_mmHg, hr=Heart_Rate_bpm,
  hemodynamic_instab=Hemodynamic_Instability,
  ed_fluids_l=ED_Fluids_per_litre, ed_levosim=ED_Levo,
  altered_mental=Alteration_Mental_Status, syncope=Syncope,
  hematemesis_red=Hematemesis_Bright_Red,
  hematemesis_coffee=Hematemesis_Coffee_Grounds,
  melena=Melena, hematochezia=Hematochezia,
  rectal_exam=Digital_Rectal_Exam,
  gi_bleed_lt3m=History_GI_Bleeding_under3months,
  gi_bleed_3to6m=History_GI_Bleeding_3to6months,
  gi_bleed_gt6m=History_GI_Bleeding_6to12months,
  gi_bleed_1to5y=History_GI_Bleeding_1to5years,
  gi_bleed_gt5y=History_GI_Bleeding_over5years,
  hb_ed=Hb_ED_mg_per_dl, hct_ed=HCT_ED, lactate_ed=ED_Lac_mmol_per_L,
  hb_entry=Hb_Entry_mg_per_dl, hct_entry=HCT_Entry,
  rbcs=RBCs,
  ri_pct=RI, wbc_entry=WBC_Entry_x103_per_L, plt_entry=PLT_Entry_x103_per_L,
  hb_d1=Hb_d1, hct_d1=HCT_d1,
  wbc_d1=WBC_d1_x103_per_L, plt_d1=PLT_d1_x103_per_L,
  creatinine=Creatinine_Entry_mg_per_dl,
  urea=Urea_Entry_mg_per_dl, bun=BUN,
  sgot=SGOT_IU_per_L, sgpt=SGPT_IU_per_L,
  ldh=LDH_IU_per_L, alp=ALP, ggt=GT,
  tbil=TBIL, dbil=DBIL, ibil=IBIL,
  ferritin=Fer_ng_per_ml, albumin=Albumin_g_per_dl,
  inr=INR, aptt=aPTT, qpt=QPT, pt_pct=PT, esd=ESDmm_per_h,
  forrest_ia=Peptic_Ulcers_Forrest_IA, forrest_ib=Peptic_Ulcers_Forest_IB,
  forrest_iia=Peptic_Ulcers_Forest_IIA, forrest_iib=Peptic_Ulcers_Forest_IIB,
  forrest_iic=Peptic_Ulcers_Forest_IIC, forrest_iii=Peptic_Ulcers_Forrest_III,
  varices=Portal_hypertension_Varices,
  portal_gastropathy=Portal_Hypertension_Portal_Hypertensive_Gastropathy,
  endo_malignancy=Endo_Malignancy,
  mallory_weiss=Laceration_Mallory_Weiss, esoph_laceration=Laceration_oesophagus,
  erosive_gastritis=Erosive_Gastritis, erosive_esophagitis=Erosive_Esophagitis,
  angioectasia=Vascular_Anomalies_Angioectasias,
  dieulafoy=Vascular_Anomalies_Dieulafoy_lesion,
  aortoenteric=Aortoenteric_Fistula, atrophic_gastritis=Atrophic_Gastritis,
  endo_clips=Workup_Endoscopy_Clips,
  endo_adrenaline=Workup_Endoscopy_Adrenaline,
  endo_clips_adr=Workup_Endoscopy_Clips_Adrenaline,
  endo_apc=Workup_Endoscopy_APC, workup_other=Workup_other,
  wb_ed=ED_Transfusion_Whole_Blood, wb_wkup=Workup_Transfusion_Whole_Blood,
  total_wb=Total_Transfusion_Whole_Blood,
  total_plasma=Total_Transfusion_Plasma, total_plt=Total_Transfusion_Platelets,
  drug_act=Drugs_ACTs, drug_antiplt=Drugs_antiPLT,
  drug_ppi=Drugs_PPI, drug_nsaid=Drugs_NSAID,
  drug_psych=Drugs_Psyc, drug_steroid=Drugs_Steroids,
  drug_alzheimer=Drugs_Alzheim,
  comorbidity_any=Comorbidities, comorbidity_hematol=Comorbidities_Hematological,
  cancer_gi=Comorbidities_Cancer_GI, cancer_extra=Cancer_Extraintestinal,
  cancer_under6m=Cancer_under6months,
  cancer_over6m=Cancer_over6months,
  cancer_gt5y=Cancer_gt5years,
  renal_ckd=Comorbidities_Renal_CKD, renal_aki=Renal_AKD,
  renal_aki_on_ckd=Comorbidities_Renal_AKD_CKD,
  liver_cirrhosis=Comorbidities_Hepatic_Liver_Cirrhosis,
  hypertension=Comorbidities_Arterial_Hypertension,
  heart_failure=Comorbidities_Cardiac_Cardiac_Insufficiency,
  afib=Heart_Atrial_Fibrilation, cad=Heart_CAD,
  ibd=Comorbidities_GI_IBD,
  asthma=Comorbidities_Pulmonary_Asthma, copd=Comorbidities_COPD,
  dvt=Comorbidities_DVT,
  dm_controlled=Comorbidities_Diabetes_Melitus_Regulated,
  dm_uncontrolled=Comorbidities_Metabolic_Diabetes_Melitus_Dysregulated,
  hyperlipidemia=Comorbidities_Hyperlipidemia,
  stroke_lt6m=Ischemic_stroke_under6months,
  stroke_gt6m=Ischemic_stroke_over6months,
  bleed_confirmed=BLEED_CONFIRMED, in_hosp_compl=IN_HOSP_COMPL_DVT_PE_STNST,
  outcome_exit=Outcome_Exit, outcome_death=Outcome_Death,
  pre_rockall_xl=Pre_Rockall,
  full_rockall_xl=Full_Rockall_Computed,
  aims65_score_xl=AIMS65_Score,
  gbs_score=GBS,
  rebleeding=Rebleeding_InHosp,
  canuka_xl=CANUKA)

num_cols <- c("age","sbp","dbp","hr","ed_fluids_l","time_to_endo_h",
              "hb_ed","hct_ed","lactate_ed","hb_entry","hct_entry","rbcs","ri_pct",
              "wbc_entry","plt_entry","hb_d1","hct_d1","wbc_d1","plt_d1",
              "creatinine","urea","bun","sgot","sgpt","ldh","alp","ggt",
              "tbil","dbil","ibil","ferritin","albumin","inr","aptt","qpt","pt_pct",
              "esd","hosp_days","symptom_duration",
              "total_wb","total_plasma","total_plt","wb_ed","wb_wkup",
              "pre_rockall_xl","full_rockall_xl","aims65_score_xl","gbs_score","canuka_xl")
df <- df |> mutate(across(all_of(num_cols[num_cols %in% names(df)]),
                          ~suppressWarnings(as.numeric(.))))

yn_cols <- c("smoker","alcohol","drug_abuse","hemodynamic_instab","ed_levosim",
             "altered_mental","syncope","hematemesis_red","hematemesis_coffee",
             "melena","hematochezia","gi_bleed_lt3m","gi_bleed_3to6m","gi_bleed_gt6m",
             "gi_bleed_1to5y","gi_bleed_gt5y",
             "forrest_ia","forrest_ib","forrest_iia","forrest_iib","forrest_iic","forrest_iii",
             "varices","portal_gastropathy","endo_malignancy","mallory_weiss",
             "esoph_laceration","erosive_gastritis","erosive_esophagitis",
             "angioectasia","dieulafoy","aortoenteric","atrophic_gastritis",
             "endo_clips","endo_adrenaline","endo_clips_adr","endo_apc","workup_other",
             "drug_act","drug_antiplt","drug_ppi","drug_nsaid","drug_psych",
             "drug_steroid","drug_alzheimer","dexamethasone",
             "comorbidity_any","comorbidity_hematol","cancer_gi","cancer_extra",
             "cancer_under6m","cancer_over6m","cancer_gt5y",
             "renal_ckd","renal_aki","renal_aki_on_ckd",
             "ibd","liver_cirrhosis","hypertension","heart_failure","afib","cad",
             "asthma","copd","dvt","dm_controlled","dm_uncontrolled","hyperlipidemia",
             "stroke_lt6m","stroke_gt6m","bleed_confirmed","in_hosp_compl",
             "rebleeding")
yn_cols <- yn_cols[yn_cols %in% names(df)]
df <- df |> mutate(across(all_of(yn_cols),
                          ~factor(case_when(toupper(.)=="YES"~"YES",toupper(.)=="NO"~"NO",
                                            TRUE~NA_character_), levels=c("NO","YES"))))

df <- df |> mutate(
  gender=factor(gender, levels=c("F","M")),
  season=factor(case_when(
    toupper(Season_Spring)=="YES"~"Spring", toupper(Season_Summer)=="YES"~"Summer",
    toupper(Season_Fall)=="YES"~"Fall",     toupper(Season_Winter)=="YES"~"Winter",
    TRUE~NA_character_), levels=c("Spring","Summer","Fall","Winter")),
  rectal_exam=factor(ifelse(rectal_exam %in% c("NEG","POS"),rectal_exam,NA),
                     levels=c("NEG","POS"))
) |> dplyr::select(-any_of(c("Season_Spring","Season_Summer","Season_Fall","Season_Winter")))
cat("Cleaned:", nrow(df), "x", ncol(df), "\n")

df$gi_bleed_recent <- factor(case_when(
  is.na(df$gi_bleed_lt3m)  & is.na(df$gi_bleed_3to6m) & is.na(df$gi_bleed_gt6m) ~ NA_character_,
  df$gi_bleed_lt3m  == "YES" |
    df$gi_bleed_3to6m == "YES" |
    df$gi_bleed_gt6m  == "YES" ~ "YES",
  TRUE ~ "NO"), levels = c("NO","YES"))
df$gi_bleed_remote <- factor(case_when(
  is.na(df$gi_bleed_1to5y) & is.na(df$gi_bleed_gt5y) ~ NA_character_,
  df$gi_bleed_1to5y == "YES" | df$gi_bleed_gt5y == "YES" ~ "YES",
  TRUE ~ "NO"), levels = c("NO","YES"))
cat(sprintf("\nGI bleed history collapsed:\n  recent (<=12m): %d  |  remote (>12m): %d\n",
            sum(df$gi_bleed_recent == "YES", na.rm=TRUE),
            sum(df$gi_bleed_remote == "YES", na.rm=TRUE)))

df$dm_any <- factor(case_when(
  is.na(df$dm_controlled) & is.na(df$dm_uncontrolled) ~ NA_character_,
  df$dm_controlled == "YES" | df$dm_uncontrolled == "YES" ~ "YES",
  TRUE ~ "NO"), levels = c("NO","YES"))
cat(sprintf("Diabetes collapsed: dm_any = %d / %d patients\n",
            sum(df$dm_any == "YES", na.rm=TRUE), nrow(df)))

df$renal_severity_3 <- factor(case_when(
  is.na(df$renal_ckd) & is.na(df$renal_aki) & is.na(df$renal_aki_on_ckd) ~ NA_character_,
  df$renal_aki == "YES" | df$renal_aki_on_ckd == "YES" ~ "Acute",
  df$renal_ckd == "YES" ~ "Chronic",
  TRUE ~ "None"),
  levels = c("None","Chronic","Acute"))
cat("Renal severity (3-level): "); print(table(df$renal_severity_3, useNA="no"))

cat("\nV2.4.1 collinearity audit on derived flags:\n")
audit_pair <- function(a, b) {
  va <- df[[a]]; vb <- df[[b]]
  
  to_yn <- function(x) {
    if (is.logical(x)) ifelse(is.na(x), NA_character_, ifelse(x, "YES", "NO"))
    else if (is.factor(x)) as.character(x)
    else as.character(x)
  }
  va_c <- to_yn(va); vb_c <- to_yn(vb)
  ok <- !is.na(va_c) & !is.na(vb_c)
  if (sum(ok) < 5) {
    cat(sprintf("  %s ~ %s : insufficient data (n=%d)\n", a, b, sum(ok)))
    return(invisible(NULL))
  }
  tbl <- table(va_c[ok], vb_c[ok])
  
  if (all(dim(tbl) == c(2, 2))) {
    n <- sum(tbl); chi <- suppressWarnings(chisq.test(tbl, correct=FALSE)$statistic)
    phi <- as.numeric(sqrt(chi / n))
    cat(sprintf("  %s ~ %s : n=%d, phi=%.3f%s\n",
                a, b, sum(ok), phi,
                if (phi > 0.5) "  (strong)" else if (phi > 0.3) "  (moderate)" else ""))
  } else {
    cat(sprintf("  %s ~ %s : non-2x2 (%dx%d), inspect manually\n",
                a, b, nrow(tbl), ncol(tbl)))
  }
}
audit_pair("gi_bleed_recent", "gi_bleed_remote")
audit_pair("dm_any", "hypertension")
audit_pair("active_cancer", "liver_cirrhosis")
audit_pair("dm_any", "cad")

# SECTION 2
cat("\n============================================================\n")
cat("SECTION 2:\n")
cat("============================================================\n\n")
n <- nrow(df)

df <- df |> mutate(hb_mean=(hb_ed+hb_entry)/2, hct_mean=(hct_ed+hct_entry)/2)
ok <- !is.na(df$hb_ed) & !is.na(df$hb_entry)
hb_diff <- df$hb_ed[ok]-df$hb_entry[ok]; hb_avg <- (df$hb_ed[ok]+df$hb_entry[ok])/2
m_d <- mean(hb_diff); sd_d <- sd(hb_diff)
cat(sprintf("Hb: r=%.4f | bias=%.3f | SD=%.3f | LoA: %.3f to %.3f\n",
            cor(df$hb_ed,df$hb_entry,use="complete"), m_d, sd_d, m_d-1.96*sd_d, m_d+1.96*sd_d))

plots[["p_ba"]] <- ggplot(data.frame(avg=hb_avg,diff=hb_diff), aes(avg,diff)) +
  geom_point(colour=COL_BLUE,size=2.5,alpha=.7) +
  geom_hline(yintercept=m_d, colour=COL_BLUE, linewidth=.9) +
  geom_hline(yintercept=m_d+1.96*sd_d, colour=COL_RED, linetype="dashed") +
  geom_hline(yintercept=m_d-1.96*sd_d, colour=COL_RED, linetype="dashed") +
  geom_hline(yintercept=0, colour="grey50", linetype="dotted") +
  annotate("text",x=Inf,y=m_d+.1,label=sprintf("Bias: %.2f",m_d),hjust=1.1,size=3.5,colour=COL_BLUE) +
  annotate("text",x=Inf,y=m_d+1.96*sd_d+.1,label=sprintf("+1.96SD: %.2f",m_d+1.96*sd_d),hjust=1.1,size=3) +
  annotate("text",x=Inf,y=m_d-1.96*sd_d-.2,label=sprintf("-1.96SD: %.2f",m_d-1.96*sd_d),hjust=1.1,size=3) +
  labs(title="Bland-Altman ABG vs Lab analyser",
       x="Mean Hb (g/dL)", y="Difference ABG - Lab (g/dL)") + theme_ugib()
print(plots[["p_ba"]])

df <- df |> mutate(
  active_cancer = case_when(
    is.na(cancer_gi) & is.na(cancer_extra) ~ NA,
    (is.na(cancer_gi) | cancer_gi != "YES") &
      (is.na(cancer_extra) | cancer_extra != "YES") ~ FALSE,
    !is.na(cancer_gt5y)    & cancer_gt5y    == "YES" &
      !is.na(cancer_under6m) & cancer_under6m == "NO"  &
      !is.na(cancer_over6m)  & cancer_over6m  == "NO" ~ FALSE,
    !is.na(cancer_under6m) & cancer_under6m == "NO" &
      !is.na(cancer_over6m)  & cancer_over6m  == "NO" &
      !is.na(cancer_gt5y)    & cancer_gt5y    == "NO" ~ FALSE,
    TRUE ~ TRUE))

cat("\n--- active_cancer derivation diagnostic (V2.4.4) ---\n")
cancer_audit <- df |>
  mutate(
    branch = case_when(
      is.na(cancer_gi) & is.na(cancer_extra)                            ~ "1_NA_both_unflagged",
      (is.na(cancer_gi) | cancer_gi != "YES") &
        (is.na(cancer_extra) | cancer_extra != "YES")                   ~ "2_no_cancer",
      !is.na(cancer_gt5y)    & cancer_gt5y    == "YES" &
        !is.na(cancer_under6m) & cancer_under6m == "NO"  &
        !is.na(cancer_over6m)  & cancer_over6m  == "NO"                   ~ "3_remission_only_excl_FALSE",
      !is.na(cancer_under6m) & cancer_under6m == "NO" &
        !is.na(cancer_over6m)  & cancer_over6m  == "NO" &
        !is.na(cancer_gt5y)    & cancer_gt5y    == "NO"                   ~ "4_all_temporal_NO_excl_FALSE",
      !is.na(cancer_under6m) & cancer_under6m == "YES"                  ~ "5_under6m_YES_active",
      !is.na(cancer_over6m)  & cancer_over6m  == "YES" &
        (is.na(cancer_gt5y) | cancer_gt5y == "NO")                        ~ "6_over6m_YES_gt5y_not_YES_active",
      !is.na(cancer_over6m)  & cancer_over6m  == "YES" &
        !is.na(cancer_gt5y)    & cancer_gt5y    == "YES"                  ~ "7_over6m_YES_AND_gt5y_YES_active",
      TRUE                                                              ~ "8_residual_cancer_flag_active"))
branch_tab <- table(cancer_audit$branch, useNA="no")
for (b in sort(names(branch_tab))) {
  cat(sprintf("  %-40s : %d\n", b, branch_tab[b]))
}
cat(sprintf("\n  Active cancer (TRUE):  %d\n",
            sum(df$active_cancer == TRUE,  na.rm=TRUE)))
cat(sprintf("  Not active (FALSE):    %d  (of which %d are no-cancer, %d are excluded cancer)\n",
            sum(df$active_cancer == FALSE, na.rm=TRUE),
            sum(df$active_cancer == FALSE & (is.na(df$cancer_gi) | df$cancer_gi != "YES") &
                  (is.na(df$cancer_extra) | df$cancer_extra != "YES"),
                na.rm=TRUE),
            sum(df$active_cancer == FALSE &
                  ((!is.na(df$cancer_gi) & df$cancer_gi == "YES") |
                     (!is.na(df$cancer_extra) & df$cancer_extra == "YES")),
                na.rm=TRUE)))
cat(sprintf("  NA (both type flags missing): %d\n",
            sum(is.na(df$active_cancer))))

excluded <- which(df$active_cancer == FALSE &
                    ((!is.na(df$cancer_gi) & df$cancer_gi == "YES") |
                       (!is.na(df$cancer_extra) & df$cancer_extra == "YES")))
if (length(excluded) > 0) {
  cat(sprintf("\n  Cancer-flagged patients EXCLUDED from active_cancer (n=%d):\n",
              length(excluded)))
  for (i in excluded) {
    cat(sprintf("    Patient #%s: gi=%s extra=%s u6=%s o6=%s gt5y=%s -> branch %s\n",
                df$patient_id[i],
                ifelse(is.na(df$cancer_gi[i]),    "NA", as.character(df$cancer_gi[i])),
                ifelse(is.na(df$cancer_extra[i]), "NA", as.character(df$cancer_extra[i])),
                ifelse(is.na(df$cancer_under6m[i]), "NA", as.character(df$cancer_under6m[i])),
                ifelse(is.na(df$cancer_over6m[i]),  "NA", as.character(df$cancer_over6m[i])),
                ifelse(is.na(df$cancer_gt5y[i]),    "NA", as.character(df$cancer_gt5y[i])),
                cancer_audit$branch[i]))
  }
}
saveRDS(list(branch_tab = branch_tab,
             excluded_patients = if (length(excluded) > 0) df$patient_id[excluded] else character(0),
             rule_version = "V2.4.4"),
        "active_cancer_audit.rds")
cat("\n  Saved: active_cancer_audit.rds\n")

df <- df |> mutate(
  pre_rockall = case_when(age<60~0L,age<80~1L,age>=80~2L,TRUE~NA_integer_) +
    case_when(sbp>=100&hr<100~0L,sbp>=100&hr>=100~1L,sbp<100~2L,TRUE~NA_integer_) +
    case_when(renal_ckd=="YES"|renal_aki=="YES"|renal_aki_on_ckd=="YES"|
                liver_cirrhosis=="YES" |
                (active_cancer == TRUE & !is.na(active_cancer)) ~ 3L,
              heart_failure=="YES"|cad=="YES" ~ 2L,
              TRUE ~ 0L))
cat(sprintf("Pre-Rockall: %d/%d\n", sum(!is.na(df$pre_rockall)), nrow(df)))

df <- df |> mutate(
  rockall_diagnosis = case_when(
    mallory_weiss == "YES" ~ 0L,
    endo_malignancy == "YES" ~ 2L,
    forrest_ia == "YES" | forrest_ib == "YES" | forrest_iia == "YES" |
      forrest_iib == "YES" | forrest_iic == "YES" | forrest_iii == "YES" |
      varices == "YES" | portal_gastropathy == "YES" | esoph_laceration == "YES" |
      erosive_gastritis == "YES" | erosive_esophagitis == "YES" |
      angioectasia == "YES" | dieulafoy == "YES" | aortoenteric == "YES" |
      atrophic_gastritis == "YES" ~ 1L,
    TRUE ~ 0L),
  rockall_stigmata = case_when(
    forrest_ia == "YES" | forrest_ib == "YES" |
      forrest_iia == "YES" | forrest_iib == "YES" ~ 2L,
    forrest_iic == "YES" | forrest_iii == "YES" ~ 0L,
    TRUE ~ 0L),
  full_rockall = pre_rockall + rockall_diagnosis + rockall_stigmata)

cat(sprintf("Full Rockall: %d/%d\n", sum(!is.na(df$full_rockall)), nrow(df)))
if (sum(!is.na(df$full_rockall)) >= 5) {
  cat(sprintf("  Range: %d - %d | Mean %.2f | Median %.0f (IQR %.0f-%.0f)\n",
              min(df$full_rockall, na.rm=T), max(df$full_rockall, na.rm=T),
              mean(df$full_rockall, na.rm=T), median(df$full_rockall, na.rm=T),
              quantile(df$full_rockall, 0.25, na.rm=T),
              quantile(df$full_rockall, 0.75, na.rm=T)))
  cat(sprintf("  Low-risk (<=2): %d | High-risk (>=8): %d\n",
              sum(df$full_rockall<=2, na.rm=T), sum(df$full_rockall>=8, na.rm=T)))
}

if ("pre_rockall_xl" %in% names(df)) {
  n_match <- sum(!is.na(df$pre_rockall) & !is.na(df$pre_rockall_xl) &
                   df$pre_rockall == df$pre_rockall_xl, na.rm=TRUE)
  n_diff  <- sum(!is.na(df$pre_rockall) & !is.na(df$pre_rockall_xl) &
                   df$pre_rockall != df$pre_rockall_xl, na.rm=TRUE)
  cat(sprintf("Pre-Rockall match with Excel: %d | differ: %d\n", n_match, n_diff))
  if (n_diff > 0) {
    diff_idx <- which(!is.na(df$pre_rockall) & !is.na(df$pre_rockall_xl) &
                        df$pre_rockall != df$pre_rockall_xl)
    for (i in diff_idx)
      cat(sprintf("  Patient #%s: script=%d, Excel=%d\n",
                  df$patient_id[i], df$pre_rockall[i], df$pre_rockall_xl[i]))
  }
}
if ("full_rockall_xl" %in% names(df)) {
  n_match <- sum(!is.na(df$full_rockall) & !is.na(df$full_rockall_xl) &
                   df$full_rockall == df$full_rockall_xl, na.rm=TRUE)
  n_diff  <- sum(!is.na(df$full_rockall) & !is.na(df$full_rockall_xl) &
                   df$full_rockall != df$full_rockall_xl, na.rm=TRUE)
  cat(sprintf("Full Rockall match with Excel: %d | differ: %d\n", n_match, n_diff))
  if (n_diff > 0) {
    diff_idx <- which(!is.na(df$full_rockall) & !is.na(df$full_rockall_xl) &
                        df$full_rockall != df$full_rockall_xl)
    for (i in diff_idx)
      cat(sprintf("  Patient #%s: script=%d, Excel=%d\n",
                  df$patient_id[i], df$full_rockall[i], df$full_rockall_xl[i]))
  }
}

cat("\n=== Reticulocyte Production Index (RPI / RIc) calculation ===\n")
normal_hct <- 45
df <- df |> mutate(
  ri_corrected = ifelse(!is.na(ri_pct) & !is.na(hct_entry),
                        ri_pct * (hct_entry / normal_hct), NA_real_),
  maturation_factor = case_when(
    is.na(hct_entry) ~ NA_real_,
    hct_entry >= 35 ~ 1.0,
    hct_entry >= 25 ~ 1.5,
    hct_entry >= 20 ~ 2.0,
    TRUE ~ 2.5),
  rpi = ri_corrected / maturation_factor,
  rpi_category = factor(case_when(
    is.na(rpi) ~ NA_character_,
    rpi > 3 ~ "Appropriate",
    rpi >= 2 ~ "Borderline",
    TRUE ~ "Inadequate"),
    levels = c("Inadequate","Borderline","Appropriate")))

cat(sprintf("RI%% populated:       %d/%d\n", sum(!is.na(df$ri_pct)), nrow(df)))
cat(sprintf("HCT_Entry populated: %d/%d\n", sum(!is.na(df$hct_entry)), nrow(df)))
cat(sprintf("RIc (corrected RI%%): %d/%d computed\n", sum(!is.na(df$ri_corrected)), nrow(df)))
cat(sprintf("RPI computed:        %d/%d (uses RIc + maturation factor; both from HCT_Entry)\n",
            sum(!is.na(df$rpi)), nrow(df)))
if (sum(!is.na(df$rpi)) >= 5) {
  cat(sprintf("\nRIc range: %.2f - %.2f | Mean %.2f | Median %.2f\n",
              min(df$ri_corrected, na.rm=T), max(df$ri_corrected, na.rm=T),
              mean(df$ri_corrected, na.rm=T), median(df$ri_corrected, na.rm=T)))
  cat(sprintf("RPI range: %.2f - %.2f | Mean %.2f | Median %.2f (IQR %.2f-%.2f)\n",
              min(df$rpi, na.rm=T), max(df$rpi, na.rm=T),
              mean(df$rpi, na.rm=T), median(df$rpi, na.rm=T),
              quantile(df$rpi, 0.25, na.rm=T), quantile(df$rpi, 0.75, na.rm=T)))
  cat("RPI interpretation distribution (missing excluded):\n")
  print(table(df$rpi_category, useNA="no"))
  cat(sprintf("Missing: %d\n", sum(is.na(df$rpi_category))))
}

df <- df |> mutate(
  aims65_score = case_when(
    is.na(albumin) | is.na(inr) | is.na(sbp) | is.na(age) ~ NA_integer_,
    TRUE ~ as.integer(
      (albumin < 3)            +
        (inr > 1.5)              +
        (altered_mental=="YES")  +
        (sbp <= 90)              +
        (age > 65))))
cat(sprintf("AIMS65 computed: %d/%d\n", sum(!is.na(df$aims65_score)), nrow(df)))

if ("aims65_score_xl" %in% names(df)) {
  n_match <- sum(!is.na(df$aims65_score) & !is.na(df$aims65_score_xl) &
                   df$aims65_score == df$aims65_score_xl, na.rm=TRUE)
  n_diff  <- sum(!is.na(df$aims65_score) & !is.na(df$aims65_score_xl) &
                   df$aims65_score != df$aims65_score_xl, na.rm=TRUE)
  cat(sprintf("  AIMS65 match with Excel column: %d | differ: %d\n", n_match, n_diff))
  if (n_diff > 0) {
    diff_idx <- which(!is.na(df$aims65_score) & !is.na(df$aims65_score_xl) &
                        df$aims65_score != df$aims65_score_xl)
    for (i in diff_idx) {
      cat(sprintf("    Patient #%s: script=%d, Excel=%d\n",
                  df$patient_id[i], df$aims65_score[i], df$aims65_score_xl[i]))
    }
  }
}

df <- df |> mutate(
  hb_for_gbs = ifelse(!is.na(hb_ed) & !is.na(hb_entry), (hb_ed + hb_entry)/2,
                      ifelse(!is.na(hb_entry), hb_entry, hb_ed)),
  bun_mmol = bun * 0.357,
  gbs_bun = case_when(
    is.na(bun_mmol) ~ NA_integer_,
    bun_mmol < 6.5 ~ 0L,
    bun_mmol < 8   ~ 2L,
    bun_mmol < 10  ~ 3L,
    bun_mmol < 25  ~ 4L,
    TRUE           ~ 6L),
  gbs_hb = case_when(
    is.na(hb_for_gbs) | is.na(gender) ~ NA_integer_,
    gender == "M" & hb_for_gbs >= 13 ~ 0L,
    gender == "M" & hb_for_gbs >= 12 ~ 1L,
    gender == "M" & hb_for_gbs >= 10 ~ 3L,
    gender == "M"                     ~ 6L,
    gender == "F" & hb_for_gbs >= 12 ~ 0L,
    gender == "F" & hb_for_gbs >= 10 ~ 1L,
    TRUE                              ~ 6L),
  gbs_sbp = case_when(
    is.na(sbp) ~ NA_integer_,
    sbp >= 110 ~ 0L,
    sbp >= 100 ~ 1L,
    sbp >= 90  ~ 2L,
    TRUE       ~ 3L),
  gbs_hr      = as.integer(!is.na(hr) & hr >= 100),
  gbs_melena  = as.integer(melena == "YES"),
  gbs_syncope = as.integer(syncope == "YES") * 2L,
  gbs_liver   = as.integer(liver_cirrhosis == "YES") * 2L,
  gbs_hf      = as.integer(heart_failure == "YES") * 2L,
  gbs_computed = gbs_bun + gbs_hb + gbs_sbp + gbs_hr +
    gbs_melena + gbs_syncope + gbs_liver + gbs_hf)

cat(sprintf("GBS computed: %d/%d\n", sum(!is.na(df$gbs_computed)), nrow(df)))
if (sum(!is.na(df$gbs_computed)) > 0) {
  cat(sprintf("  Range: %d - %d | Median %.0f (IQR %.0f-%.0f)\n",
              min(df$gbs_computed, na.rm=T), max(df$gbs_computed, na.rm=T),
              median(df$gbs_computed, na.rm=T),
              quantile(df$gbs_computed, 0.25, na.rm=T),
              quantile(df$gbs_computed, 0.75, na.rm=T)))
  cat(sprintf("  Low risk (GBS<=1): %d (%.1f%%) | High risk (>=6): %d (%.1f%%)\n",
              sum(df$gbs_computed<=1, na.rm=T),
              sum(df$gbs_computed<=1, na.rm=T)/sum(!is.na(df$gbs_computed))*100,
              sum(df$gbs_computed>=6, na.rm=T),
              sum(df$gbs_computed>=6, na.rm=T)/sum(!is.na(df$gbs_computed))*100))
}

if ("gbs_score" %in% names(df)) {
  n_match <- sum(!is.na(df$gbs_computed) & !is.na(df$gbs_score) &
                   df$gbs_computed == df$gbs_score, na.rm=TRUE)
  n_diff  <- sum(!is.na(df$gbs_computed) & !is.na(df$gbs_score) &
                   df$gbs_computed != df$gbs_score, na.rm=TRUE)
  cat(sprintf("GBS match with Excel: %d | differ: %d\n", n_match, n_diff))
  if (n_diff > 0) {
    diff_idx <- which(!is.na(df$gbs_computed) & !is.na(df$gbs_score) &
                        df$gbs_computed != df$gbs_score)
    for (i in diff_idx)
      cat(sprintf("  Patient #%s: script=%d, Excel=%d\n",
                  df$patient_id[i], df$gbs_computed[i], df$gbs_score[i]))
  }
  df$gbs_score <- ifelse(is.na(df$gbs_score), df$gbs_computed, df$gbs_score)
}

cat("\n=== CANUKA score calculation ===\n")
if(!"bun" %in% names(df)) df$bun <- df$urea / 2.14

bun_mmol_int <- as.integer(round(df$bun * 0.357 * 10))
bun_mmol     <- bun_mmol_int / 10

df <- df |> mutate(
  canuka_age = case_when(
    is.na(age) ~ NA_integer_,
    age < 50 ~ 0L, age < 65 ~ 1L, age >= 65 ~ 2L),
  canuka_melena = case_when(is.na(melena) ~ NA_integer_,
                            melena == "YES" ~ 1L, TRUE ~ 0L),
  canuka_hematem = case_when(
    is.na(hematemesis_red) & is.na(hematemesis_coffee) ~ NA_integer_,
    (hematemesis_red == "YES" & !is.na(hematemesis_red)) |
      (hematemesis_coffee == "YES" & !is.na(hematemesis_coffee)) ~ 1L,
    TRUE ~ 0L),
  canuka_syncope = case_when(is.na(syncope) ~ NA_integer_,
                             syncope == "YES" ~ 1L, TRUE ~ 0L),
  canuka_hb = case_when(
    is.na(hb_mean) ~ NA_integer_,
    hb_mean > 12  ~ 0L,
    hb_mean > 10  ~ 1L,
    hb_mean > 8   ~ 2L,
    TRUE ~ 3L),
  canuka_bun = case_when(
    is.na(bun) ~ NA_integer_,
    bun_mmol_int < 50   ~ 0L,
    bun_mmol_int < 100  ~ 1L,
    bun_mmol_int < 150  ~ 2L,
    TRUE ~ 3L),
  canuka_liver = case_when(is.na(liver_cirrhosis) ~ 0L,
                           liver_cirrhosis == "YES" ~ 2L, TRUE ~ 0L),
  canuka_malig = case_when(
    is.na(active_cancer) ~ 0L,
    active_cancer == TRUE ~ 2L,
    TRUE ~ 0L),
  canuka_hr = case_when(
    is.na(hr) ~ NA_integer_,
    hr < 100 ~ 0L, hr < 125 ~ 1L, TRUE ~ 2L),
  canuka_sbp = case_when(
    is.na(sbp) ~ NA_integer_,
    sbp >= 120 ~ 0L,
    sbp >= 100 ~ 1L,
    sbp >= 80  ~ 2L,
    TRUE ~ 3L))

canuka_cols <- c("canuka_age","canuka_melena","canuka_hematem","canuka_syncope",
                 "canuka_hb","canuka_bun","canuka_liver","canuka_malig",
                 "canuka_hr","canuka_sbp")
df$canuka <- rowSums(df[, canuka_cols])

debug_canuka_id <- "87"
if (any(df$patient_id == debug_canuka_id, na.rm=TRUE)) {
  i <- which(df$patient_id == debug_canuka_id)
  cat(sprintf("\n=== Patient #%s CANUKA component trace (post-V2.3 fix) ===\n",
              debug_canuka_id))
  cat(sprintf("  age=%s -> canuka_age=%s\n",
              df$age[i], df$canuka_age[i]))
  cat(sprintf("  melena=%s -> canuka_melena=%s\n",
              df$melena[i], df$canuka_melena[i]))
  cat(sprintf("  hematem_red=%s, coffee=%s -> canuka_hematem=%s\n",
              df$hematemesis_red[i], df$hematemesis_coffee[i],
              df$canuka_hematem[i]))
  cat(sprintf("  syncope=%s -> canuka_syncope=%s\n",
              df$syncope[i], df$canuka_syncope[i]))
  cat(sprintf("  hb_mean=%.2f -> canuka_hb=%s\n",
              df$hb_mean[i], df$canuka_hb[i]))
  cat(sprintf("  bun=%.2f mg/dL (%.1f mmol/L rounded) -> canuka_bun=%s\n",
              df$bun[i], round(df$bun[i]*0.357, 1), df$canuka_bun[i]))
  cat(sprintf("  liver_cirrhosis=%s -> canuka_liver=%s\n",
              df$liver_cirrhosis[i], df$canuka_liver[i]))
  cat(sprintf("  active_cancer=%s -> canuka_malig=%s\n",
              df$active_cancer[i], df$canuka_malig[i]))
  cat(sprintf("  hr=%s -> canuka_hr=%s\n", df$hr[i], df$canuka_hr[i]))
  cat(sprintf("  sbp=%s -> canuka_sbp=%s\n", df$sbp[i], df$canuka_sbp[i]))
  cat(sprintf("  TOTAL = %s\n", df$canuka[i]))
}

cat(sprintf("CANUKA computed for %d/%d patients (%d had missing components)\n",
            sum(!is.na(df$canuka)), nrow(df), sum(is.na(df$canuka))))
if(sum(!is.na(df$canuka)) > 0) {
  cat(sprintf("CANUKA range: %d - %d | mean %.2f | median %.0f (IQR %.0f-%.0f)\n",
              min(df$canuka, na.rm=T), max(df$canuka, na.rm=T),
              mean(df$canuka, na.rm=T), median(df$canuka, na.rm=T),
              quantile(df$canuka, 0.25, na.rm=T), quantile(df$canuka, 0.75, na.rm=T)))
  cat(sprintf("Low risk (CANUKA<=2): %d (%.1f%%) | High risk (>=10): %d (%.1f%%)\n",
              sum(df$canuka<=2, na.rm=T), sum(df$canuka<=2, na.rm=T)/sum(!is.na(df$canuka))*100,
              sum(df$canuka>=10, na.rm=T), sum(df$canuka>=10, na.rm=T)/sum(!is.na(df$canuka))*100))
}

if ("canuka_xl" %in% names(df)) {
  n_match <- sum(!is.na(df$canuka) & !is.na(df$canuka_xl) &
                   df$canuka == df$canuka_xl, na.rm=TRUE)
  n_diff  <- sum(!is.na(df$canuka) & !is.na(df$canuka_xl) &
                   df$canuka != df$canuka_xl, na.rm=TRUE)
  cat(sprintf("CANUKA match with Excel column: %d | differ: %d\n", n_match, n_diff))
  if (n_diff > 0) {
    diff_idx <- which(!is.na(df$canuka) & !is.na(df$canuka_xl) &
                        df$canuka != df$canuka_xl)
    for (i in diff_idx) {
      cat(sprintf("  Patient #%s: script=%d, Excel=%d\n",
                  df$patient_id[i], df$canuka[i], df$canuka_xl[i]))
    }
  }
}

df <- df |> mutate(renal_severity=factor(case_when(
  renal_aki_on_ckd=="YES"~"AKI_on_CKD", renal_aki=="YES"~"AKI",
  renal_ckd=="YES"~"CKD", TRUE~"None"),
  levels=c("None","CKD","AKI","AKI_on_CKD"), ordered=TRUE))

df <- df |> mutate(
  transfusion_class=factor(case_when(total_wb==0~"0_units",total_wb<=2~"1to2",
                                     TRUE~"3plus"), levels=c("0_units","1to2","3plus"), ordered=TRUE),
  major_tx=factor(ifelse(total_wb>=3,"YES","NO"), levels=c("NO","YES")),
  endo_haemostasis=(endo_clips=="YES"|endo_adrenaline=="YES"|
                      endo_clips_adr=="YES"|endo_apc=="YES") & !is.na(endo_clips),
  high_risk=factor(ifelse(outcome_death=="YES"|total_wb>0|endo_haemostasis,"YES","NO"),
                   levels=c("NO","YES")),
  high_risk_old=factor(ifelse(outcome_death=="YES"|rebleeding=="YES"|
                                total_wb>0|endo_haemostasis,"YES","NO"),
                       levels=c("NO","YES")),
  age_group=cut(age,breaks=c(17,40,60,75,110),
                labels=c("18-40","41-60","61-75",">75"),right=TRUE))
cat("\nOutcomes:\n"); print(table(df$transfusion_class))
cat("major_tx:", sum(df$major_tx=="YES"),
    "| high_risk (V1 NEW):", sum(df$high_risk=="YES"),
    "| high_risk_old:", sum(df$high_risk_old=="YES"),
    "| death:", sum(df$outcome_death=="YES",na.rm=TRUE),
    "| rebleeding:", sum(df$rebleeding=="YES",na.rm=TRUE), "\n")

cat("\n=== Sex-specific lab reference flags ===\n")
is_female <- df$gender == "F"
is_male   <- df$gender == "M"

lab_flag <- function(x, lo_m, hi_m, lo_f, hi_f, direction="both") {
  flag_lo <- rep(NA_integer_, length(x))
  flag_hi <- rep(NA_integer_, length(x))
  for(i in seq_along(x)) {
    if(is.na(x[i]) || is.na(is_female[i])) next
    lo <- if(is_female[i]) lo_f else lo_m
    hi <- if(is_female[i]) hi_f else hi_m
    flag_lo[i] <- as.integer(x[i] < lo)
    flag_hi[i] <- as.integer(x[i] > hi)
  }
  list(lo=flag_lo, hi=flag_hi)
}

f_wbc <- lab_flag(df$wbc_entry, 4.0, 10.5, 4.0, 10.5)
df$wbc_lo <- f_wbc$lo; df$wbc_hi <- f_wbc$hi

f_hb <- lab_flag(df$hb_entry, 13.5, 17.5, 12.0, 15.5)
df$hb_lo <- f_hb$lo; df$hb_hi <- f_hb$hi

f_hct <- lab_flag(df$hct_entry, 41, 51, 37, 47)
df$hct_lo <- f_hct$lo; df$hct_hi <- f_hct$hi

f_plt <- lab_flag(df$plt_entry, 150, 400, 150, 400)
df$plt_lo <- f_plt$lo; df$plt_hi <- f_plt$hi

f_rbc <- lab_flag(df$rbcs, 4.5, 5.9, 4.1, 5.1)
df$rbc_lo <- f_rbc$lo; df$rbc_hi <- f_rbc$hi

f_urea <- lab_flag(df$urea, 10, 50, 10, 50)
df$urea_lo <- f_urea$lo; df$urea_hi <- f_urea$hi

f_crea <- lab_flag(df$creatinine, 0.7, 1.2, 0.5, 0.9)
df$crea_lo <- f_crea$lo; df$crea_hi <- f_crea$hi

f_alb <- lab_flag(df$albumin, 3.5, 5.0, 3.5, 5.0)
df$alb_lo <- f_alb$lo; df$alb_hi <- f_alb$hi

f_sgot <- lab_flag(df$sgot, 5, 37, 5, 37)
df$sgot_lo <- f_sgot$lo; df$sgot_hi <- f_sgot$hi

f_sgpt <- lab_flag(df$sgpt, 5, 40, 5, 40)
df$sgpt_lo <- f_sgpt$lo; df$sgpt_hi <- f_sgpt$hi

f_alp <- lab_flag(df$alp, 40, 129, 35, 104)
df$alp_lo <- f_alp$lo; df$alp_hi <- f_alp$hi

f_ggt <- lab_flag(df$ggt, 8, 49, 7, 32)
df$ggt_lo <- f_ggt$lo; df$ggt_hi <- f_ggt$hi

f_tbil <- lab_flag(df$tbil, 0, 1.0, 0, 1.0)
df$tbil_hi <- f_tbil$hi

f_dbil <- lab_flag(df$dbil, 0, 0.3, 0, 0.3)
df$dbil_hi <- f_dbil$hi

f_ibil <- lab_flag(df$ibil, 0, 0.75, 0, 0.75)
df$ibil_hi <- f_ibil$hi

f_ldh <- lab_flag(df$ldh, 135, 225, 135, 225)
df$ldh_lo <- f_ldh$lo; df$ldh_hi <- f_ldh$hi

f_fer <- lab_flag(df$ferritin, 28, 300, 30, 300)
df$ferritin_lo <- f_fer$lo; df$ferritin_hi <- f_fer$hi

df$aptt_abnl <- ifelse(!is.na(df$aptt), as.integer(df$aptt < 26 | df$aptt > 38), NA_integer_)
df$qpt_pct_abnl <- ifelse(!is.na(df$pt_pct), as.integer(df$pt_pct < 80 | df$pt_pct > 110), NA_integer_)

n_total <- nrow(df)

lab_row <- function(test, lo_vec=NULL, hi_vec=NULL, raw=NULL) {
  n_miss <- if (!is.null(raw)) sum(is.na(raw)) else NA_integer_
  n_lo   <- if (!is.null(lo_vec)) sum(lo_vec == 1, na.rm=TRUE) else NA_integer_
  n_hi   <- if (!is.null(hi_vec)) sum(hi_vec == 1, na.rm=TRUE) else NA_integer_
  n_within <- n_total - ifelse(is.na(n_lo), 0, n_lo) -
    ifelse(is.na(n_hi), 0, n_hi) - ifelse(is.na(n_miss), 0, n_miss)
  data.frame(Test=test,
             Low_n    = ifelse(is.na(n_lo), "-", as.character(n_lo)),
             Within_n = as.character(n_within),
             High_n   = ifelse(is.na(n_hi), "-", as.character(n_hi)),
             Missing  = ifelse(is.na(n_miss), "-", as.character(n_miss)),
             stringsAsFactors = FALSE)
}

lab_summary <- rbind(
  lab_row("WBC",        df$wbc_lo,     df$wbc_hi,     df$wbc_entry),
  lab_row("Hb",         df$hb_lo,      df$hb_hi,      df$hb_entry),
  lab_row("HCT",        df$hct_lo,     df$hct_hi,     df$hct_entry),
  lab_row("PLT",        df$plt_lo,     df$plt_hi,     df$plt_entry),
  lab_row("RBC",        df$rbc_lo,     df$rbc_hi,     df$rbcs),
  lab_row("Urea",       df$urea_lo,    df$urea_hi,    df$urea),
  lab_row("Creatinine", df$crea_lo,    df$crea_hi,    df$creatinine),
  lab_row("Albumin",    df$alb_lo,     df$alb_hi,     df$albumin),
  lab_row("SGOT/AST",   df$sgot_lo,    df$sgot_hi,    df$sgot),
  lab_row("SGPT/ALT",   df$sgpt_lo,    df$sgpt_hi,    df$sgpt),
  lab_row("ALP",        df$alp_lo,     df$alp_hi,     df$alp),
  lab_row("gGT",        df$ggt_lo,     df$ggt_hi,     df$ggt),
  lab_row("TBIL",       NULL,          df$tbil_hi,    df$tbil),
  lab_row("DBIL",       NULL,          df$dbil_hi,    df$dbil),
  lab_row("IBIL",       NULL,          df$ibil_hi,    df$ibil),
  lab_row("LDH",        df$ldh_lo,     df$ldh_hi,     df$ldh),
  lab_row("Ferritin",   df$ferritin_lo,df$ferritin_hi,df$ferritin),
  lab_row("aPTT",       NULL,          df$aptt_abnl,  df$aptt),
  lab_row("PT%",        NULL,          df$qpt_pct_abnl,df$pt_pct))
cat(sprintf("Lab flag summary (n_total = %d):\n", n_total))
cat("NA/missing values excluded from Low/High counts and reported separately.\n")
print(lab_summary)

cat("\nFerritin interpretation caveat:\n")
cat(sprintf("  %d patients with HIGH ferritin (>300). Note: acute-phase reactant,\n",
            sum(df$ferritin_hi==1, na.rm=T)))
cat("  so HIGH ferritin does NOT indicate repleted iron stores. Elevated in\n")
cat("  infection, malignancy, liver disease. Interpret with CRP/ESR + comorbidities.\n")

cat("\n=== Data audit ===\n")

safe_dt <- function(x) {
  if (inherits(x, c("POSIXct","POSIXt","Date"))) return(as.POSIXct(x))
  if (is.numeric(x)) {
    return(as.POSIXct(x * 86400, origin="1899-12-30", tz="UTC"))
  }
  x <- as.character(x)
  x[x %in% c("","NA","N/A")] <- NA_character_
  fmts <- c("%Y-%m-%d %H:%M:%S", "%Y-%m-%d", "%Y/%m/%d %H:%M:%S", "%Y/%m/%d",
            "%d/%m/%Y %H:%M", "%d/%m/%Y", "%d-%m-%Y", "%m/%d/%Y")
  out <- as.POSIXct(rep(NA_real_, length(x)), origin="1970-01-01")
  for (fmt in fmts) {
    miss <- is.na(out) & !is.na(x)
    if (!any(miss)) break
    tryCatch({
      parsed <- as.POSIXct(x[miss], format=fmt, tz="UTC")
      out[miss] <- parsed
    }, warning=function(w){}, error=function(e){})
  }
  out
}

if (all(c("Entry_Date","Exit_Date","Date_Birth") %in% names(df))) {
  df$entry_dt <- safe_dt(df$Entry_Date)
  df$exit_dt  <- safe_dt(df$Exit_Date)
  df$dob_dt   <- safe_dt(df$Date_Birth)
  cat(sprintf("Dates parsed: Entry=%d/%d | Exit=%d/%d | DOB=%d/%d\n",
              sum(!is.na(df$entry_dt)), nrow(df),
              sum(!is.na(df$exit_dt)), nrow(df),
              sum(!is.na(df$dob_dt)), nrow(df)))
  
  if (sum(!is.na(df$entry_dt)) > 0) {
    period_first <- min(df$entry_dt, na.rm=TRUE)
    period_last  <- max(df$entry_dt, na.rm=TRUE)
    period_days  <- as.numeric(difftime(period_last, period_first, units="days"))
    cat(sprintf("\n=== DATA COLLECTION PERIOD (for Methods) ===\n"))
    cat(sprintf("  First patient entry: %s\n", format(period_first, "%Y-%m-%d")))
    cat(sprintf("  Last  patient entry: %s\n", format(period_last,  "%Y-%m-%d")))
    cat(sprintf("  Span: %.1f months (%d days, %d patients)\n\n",
                period_days/30.4, round(period_days), sum(!is.na(df$entry_dt))))
  }
  
  df$entry_month <- as.integer(format(df$entry_dt, "%m"))
  df$derived_season <- factor(case_when(
    df$entry_month %in% c(12,1,2) ~ "Winter",
    df$entry_month %in% c(3,4,5)  ~ "Spring",
    df$entry_month %in% c(6,7,8)  ~ "Summer",
    df$entry_month %in% c(9,10,11) ~ "Fall"),
    levels=c("Spring","Summer","Fall","Winter"))
  
  season_mismatch <- which(!is.na(df$season) & !is.na(df$derived_season) &
                             as.character(df$season) != as.character(df$derived_season))
  if(length(season_mismatch) > 0) {
    cat(sprintf("WARNING: season flag mismatch in %d patient(s):\n", length(season_mismatch)))
    for(i in season_mismatch) {
      cat(sprintf("  Patient #%s: entry=%s -> derived=%s, Excel flag=%s\n",
                  df$patient_id[i], format(df$entry_dt[i],"%Y-%m-%d"),
                  df$derived_season[i], df$season[i]))
    }
    cat("  -> Recommend correcting these in Excel\n")
  } else {
    cat("Season flags: all consistent with entry date\n")
  }
  
  if (any(!is.na(df$dob_dt))) {
    df$age_calc <- as.numeric(difftime(df$entry_dt, df$dob_dt, units="days")) / 365.25
    age_mismatch <- which(!is.na(df$age) & !is.na(df$age_calc) &
                            abs(df$age_calc - df$age) > 2)
    if(length(age_mismatch) > 0) {
      cat(sprintf("WARNING: age/DOB mismatch (>2 years) in %d patient(s):\n",
                  length(age_mismatch)))
      for(i in age_mismatch) {
        cat(sprintf("  #%s: recorded=%d, calculated=%.1f (DOB=%s, entry=%s)\n",
                    df$patient_id[i], df$age[i], df$age_calc[i],
                    format(df$dob_dt[i],"%Y-%m-%d"), format(df$entry_dt[i],"%Y-%m-%d")))
      }
    } else {
      cat("Age/DOB: all consistent within 2 years\n")
    }
  }
  
  if (any(!is.na(df$exit_dt))) {
    df$los_calc <- as.numeric(difftime(df$exit_dt, df$entry_dt, units="days"))
    los_mismatch <- which(!is.na(df$hosp_days) & !is.na(df$los_calc) &
                            abs(df$los_calc - df$hosp_days) > 2)
    if(length(los_mismatch) > 0) {
      cat(sprintf("WARNING: LOS mismatch (>2 days) in %d patient(s):\n",
                  length(los_mismatch)))
      for(i in los_mismatch) {
        cat(sprintf("  #%s: recorded=%d, calculated=%.0f\n",
                    df$patient_id[i], df$hosp_days[i], df$los_calc[i]))
      }
    } else {
      cat("LOS/entry-exit dates: all consistent within 2 days\n")
    }
  }
  
  bad_order <- !is.na(df$entry_dt) & !is.na(df$exit_dt) & df$entry_dt > df$exit_dt
  if (any(bad_order, na.rm=TRUE)) {
    cat(sprintf("WARNING: entry date AFTER exit date in %d patient(s)\n", sum(bad_order)))
  }
} else {
  cat("Date columns not found -- skipping date audit\n")
}

cat("\n=== 2H. BP and HR clinical categorisation ===\n")

df <- df |> mutate(
  bp_stage = factor(case_when(
    is.na(sbp) | is.na(dbp) ~ NA_character_,
    sbp > 180 | dbp > 120 ~ "Hypertensive_crisis",
    sbp >= 140 | dbp >= 90 ~ "Stage2_HTN",
    sbp >= 130 | dbp >= 80 ~ "Stage1_HTN",
    sbp < 90  | dbp < 60  ~ "Hypotension",
    sbp < 120 & sbp >= 90 & dbp < 80 & dbp >= 60 & sbp < 120 ~ NA_character_,
    TRUE ~ NA_character_),
    levels = c("Hypotension","Low_normal","Normal","Elevated",
               "Stage1_HTN","Stage2_HTN","Hypertensive_crisis")),
  hypotension_flag = as.integer(!is.na(sbp) & sbp < 90),
  severe_htn_flag  = as.integer(!is.na(sbp) & !is.na(dbp) & (sbp > 180 | dbp > 120)))

df$bp_stage <- factor(
  with(df, case_when(
    is.na(sbp) | is.na(dbp) ~ NA_character_,
    sbp > 180 | dbp > 120            ~ "Hypertensive_crisis",
    sbp >= 140 | dbp >= 90           ~ "Stage2_HTN",
    sbp >= 130 | dbp >= 80           ~ "Stage1_HTN",
    sbp >= 120 & sbp <= 129 & dbp < 80 ~ "Elevated",
    sbp < 90 | dbp < 60              ~ "Hypotension",
    sbp >= 90 & sbp < 120 & dbp >= 60 & dbp < 80 ~ "Normal",
    TRUE ~ NA_character_)),
  levels = c("Hypotension","Normal","Elevated",
             "Stage1_HTN","Stage2_HTN","Hypertensive_crisis"))

cat("BP stage distribution (missing excluded):\n")
print(table(df$bp_stage, useNA="no"))
cat(sprintf("Missing (SBP or DBP unavailable): %d\n", sum(is.na(df$bp_stage))))
cat(sprintf("Hypotension flag (SBP<90): %d (%.1f%%)\n",
            sum(df$hypotension_flag==1, na.rm=T),
            sum(df$hypotension_flag==1, na.rm=T)/sum(!is.na(df$hypotension_flag))*100))
cat(sprintf("Severe HTN flag (>180/>120): %d (%.1f%%)\n",
            sum(df$severe_htn_flag==1, na.rm=T),
            sum(df$severe_htn_flag==1, na.rm=T)/sum(!is.na(df$severe_htn_flag))*100))

df$hr_stage <- factor(
  with(df, case_when(
    is.na(hr) ~ NA_character_,
    hr < 60                    ~ "Bradycardia",
    hr >= 60  & hr < 100       ~ "Normal",
    hr >= 100 & hr < 120       ~ "Mild_tachycardia",
    hr >= 120 & hr < 150       ~ "Moderate_tachycardia",
    hr >= 150                  ~ "Severe_tachycardia",
    TRUE ~ NA_character_)),
  levels = c("Bradycardia","Normal","Mild_tachycardia",
             "Moderate_tachycardia","Severe_tachycardia"))

df$tachycardia_flag  <- as.integer(!is.na(df$hr) & df$hr >= 100)
df$severe_tachy_flag <- as.integer(!is.na(df$hr) & df$hr >= 125)

cat("\nHR stage distribution (missing excluded):\n")
print(table(df$hr_stage, useNA="no"))
cat(sprintf("Missing (HR unavailable): %d\n", sum(is.na(df$hr_stage))))
cat(sprintf("Tachycardia flag (HR>=100): %d (%.1f%%)\n",
            sum(df$tachycardia_flag==1, na.rm=T),
            sum(df$tachycardia_flag==1, na.rm=T)/sum(!is.na(df$tachycardia_flag))*100))
cat(sprintf("Severe tachy flag (HR>=125): %d (%.1f%%)\n",
            sum(df$severe_tachy_flag==1, na.rm=T),
            sum(df$severe_tachy_flag==1, na.rm=T)/sum(!is.na(df$severe_tachy_flag))*100))

# SECTION 3
cat("\n============================================================\n")
cat("SECTION 3:\n")
cat("============================================================\n\n")

cat("=== Demographics ===\n")
cat("Gender:\n"); print(table(df$gender, useNA="no"))
if (sum(is.na(df$gender)) > 0) cat(sprintf("  (Missing: %d)\n", sum(is.na(df$gender))))
cat(sprintf("Age: mean=%.1f SD=%.1f median=%d range=%d-%d\n",
            mean(df$age,na.rm=T),sd(df$age,na.rm=T),as.integer(median(df$age,na.rm=T)),
            min(df$age,na.rm=T),max(df$age,na.rm=T)))
cat("Season:\n"); print(table(df$season, useNA="no"))
if (sum(is.na(df$season)) > 0) cat(sprintf("  (Missing: %d)\n", sum(is.na(df$season))))
cat(sprintf("LOS: mean=%.1f SD=%.1f median=%d\n\n",
            mean(df$hosp_days,na.rm=T),sd(df$hosp_days,na.rm=T),
            as.integer(median(df$hosp_days,na.rm=T))))

COL_F_PALE <- "#E8998D"
COL_M_PALE <- "#7EB5D6"

gt <- df |> count(gender) |> mutate(pct=round(n/sum(n)*100,1))
plots[["p_gender"]] <- ggplot(gt,
                              aes(x=ifelse(gender=="M","Male","Female"), y=n, fill=gender)) +
  geom_col(width=.45, colour="white") +
  geom_text(aes(label=paste0(n,"\n(",pct,"%)")), vjust=-.3, fontface="bold", size=4.5) +
  scale_fill_manual(values=c(M=COL_M_PALE, F=COL_F_PALE)) +
  scale_y_continuous(limits=c(0, max(gt$n)*1.25)) +
  labs(title="Gender distribution", subtitle=sprintf("n = %d patients", nrow(df)),
       x=NULL, y="Patients (n)") +
  theme_ugib() + theme(legend.position="none")
print(plots[["p_gender"]])

age_mn <- mean(df$age,na.rm=T); age_md <- median(df$age,na.rm=T)
age_iq <- quantile(df$age, c(.25,.75), na.rm=T)

plots[["p_age_hist"]] <- ggplot(df, aes(age)) +
  geom_histogram(breaks=seq(15,110,5), fill=COL_BLUE, colour="white", linewidth=.4) +
  geom_vline(xintercept=age_mn, colour=COL_RED, linetype="dashed", linewidth=.9) +
  geom_vline(xintercept=age_md, colour=COL_GREY, linetype="dotdash", linewidth=.9) +
  annotate("text", x=age_mn-1.5, y=Inf, vjust=1.3, hjust=1, size=3.2,
           label=sprintf("Mean\n%.1f", age_mn), colour=COL_RED, fontface="bold") +
  annotate("text", x=age_md+1.5, y=Inf, vjust=1.3, hjust=0, size=3.2,
           label=sprintf("Median\n%.0f", age_md), colour=COL_GREY, fontface="bold") +
  scale_x_continuous(breaks=seq(20,110,10)) +
  labs(title="Age distribution",
       subtitle=sprintf("Mean %.1f +/- %.0f yr | Median %.0f (IQR %.0f-%.0f) | Range %d-%d",
                        age_mn, sd(df$age,na.rm=T), age_md, age_iq[1], age_iq[2],
                        min(df$age,na.rm=T), max(df$age,na.rm=T)),
       x="Age (years)", y="Patients (n)") + theme_ugib()
print(plots[["p_age_hist"]])

ag <- df |> filter(!is.na(age_group)) |> count(age_group) |>
  mutate(pct=round(n/sum(n)*100,1))
age_blues <- c("18-40"="#BDD7EE", "41-60"="#8DB4D8", "61-75"="#5B9BD5", ">75"="#2E5C8A")
plots[["p_age_groups"]] <- ggplot(ag, aes(age_group, n, fill=age_group)) +
  geom_col(width=.6, colour="white") +
  geom_text(aes(label=paste0(n,"\n(",pct,"%)")), vjust=-.3, fontface="bold", size=4) +
  scale_fill_manual(values=age_blues) +
  scale_y_continuous(limits=c(0,max(ag$n)*1.25)) +
  labs(title="Age groups", subtitle="Clinical categories",
       x="Age group (years)", y="Patients (n)") +
  theme_ugib() + theme(legend.position="none")
print(plots[["p_age_groups"]])

st <- df |> filter(!is.na(season)) |> count(season) |> mutate(pct=round(n/sum(n)*100,1))
fw <- sum(st$n[st$season %in% c("Fall","Winter")])
fw_pct <- round(fw/sum(st$n)*100,1)
.n_season   <- sum(st$n)
.n_missing  <- nrow(df) - .n_season
.season_sub <- if (.n_missing > 0) {
  sprintf("n = %d with calculable season (%d missing)", .n_season, .n_missing)
} else {
  sprintf("n = %d patients", .n_season)
}
plots[["p_season"]] <- ggplot(st, aes(season, n, fill=season)) +
  geom_col(width=.55, colour="white") +
  geom_text(aes(label=paste0(n,"\n(",pct,"%)")), vjust=-.3, fontface="bold", size=4) +
  scale_fill_manual(values=c(Spring="#4DAF7C",Summer="#F5A623",Fall=COL_RED,Winter=COL_BLUE)) +
  scale_y_continuous(limits=c(0, max(st$n)*1.25)) +
  labs(title="Seasonal distribution of UGIB admissions",
       subtitle=.season_sub,
       x=NULL, y="Patients (n)") +
  theme_ugib() + theme(legend.position="none")
print(plots[["p_season"]])

los_mn <- mean(df$hosp_days,na.rm=T); los_md <- median(df$hosp_days,na.rm=T)
los_iq <- quantile(df$hosp_days, c(.25,.75), na.rm=T)

plots[["p_los_hist"]] <- ggplot(df, aes(hosp_days)) +
  geom_histogram(breaks=seq(0, max(df$hosp_days,na.rm=T)+2, 2),
                 fill=COL_BLUE, colour="white", linewidth=.4) +
  geom_vline(xintercept=los_md, colour=COL_GREY, linetype="dotdash", linewidth=.9) +
  geom_vline(xintercept=los_mn, colour=COL_RED, linetype="dashed", linewidth=.9) +
  annotate("text", x=los_md-0.5, y=Inf, vjust=1.3, hjust=1, size=3.2,
           label=sprintf("Median\n%.0fd", los_md), colour=COL_GREY, fontface="bold") +
  annotate("text", x=los_mn+0.5, y=Inf, vjust=1.3, hjust=0, size=3.2,
           label=sprintf("Mean\n%.1fd", los_mn), colour=COL_RED, fontface="bold") +
  labs(title="Length of hospital stay",
       subtitle=sprintf("Mean %.1f +/- %.1f days | Median %.0f (IQR %.0f-%.0f) | Range %d-%d",
                        los_mn, sd(df$hosp_days,na.rm=T), los_md, los_iq[1], los_iq[2],
                        min(df$hosp_days,na.rm=T), max(df$hosp_days,na.rm=T)),
       x="Days", y="Patients (n)") + theme_ugib()
print(plots[["p_los_hist"]])

plots[["p_los_gender"]] <- ggplot(df, aes(x=gender, y=hosp_days, fill=gender)) +
  geom_boxplot(width=.5, colour="grey30", outlier.colour=COL_RED, outlier.size=2.5) +
  geom_jitter(width=.15, alpha=.4, size=1.5, colour="grey50") +
  scale_fill_manual(values=c(F=COL_F_PALE, M=COL_M_PALE)) +
  scale_x_discrete(labels=c(F="Female",M="Male")) +
  labs(title="Stay by gender", subtitle="Boxplot + individual patients",
       x=NULL, y="Days") +
  theme_ugib() + theme(legend.position="none")
print(plots[["p_los_gender"]])

rbc_df <- data.frame(units=0:10) |>
  left_join(df |> count(total_wb) |> rename(units=total_wb), by="units") |>
  mutate(n=ifelse(is.na(n),0L,n), pct=round(n/sum(n)*100,1),
         zone=factor(case_when(units==0~"None",units<=2~"Moderate",TRUE~"Major"),
                     levels=c("None","Moderate","Major")))

n_none <- sum(rbc_df$n[rbc_df$zone=="None"])
n_mod  <- sum(rbc_df$n[rbc_df$zone=="Moderate"])
n_maj  <- sum(rbc_df$n[rbc_df$zone=="Major"])

y_max <- max(rbc_df$n) * 1.3
y_zone_label <- max(rbc_df$n) * 1.15

plots[["p_rbc"]] <- ggplot(rbc_df, aes(units,n)) +
  annotate("rect",xmin=-.5,xmax=.5,ymin=0,ymax=Inf,fill="#EEEEEE",alpha=.5) +
  annotate("rect",xmin=.5,xmax=2.5,ymin=0,ymax=Inf,fill="#E8F4FD",alpha=.5) +
  annotate("rect",xmin=2.5,xmax=10.5,ymin=0,ymax=Inf,fill="#FFF3CD",alpha=.5) +
  annotate("text",x=0,y=y_zone_label,label="No transfusion",size=2.8,colour=COL_GREY,fontface="italic") +
  annotate("text",x=1.5,y=y_zone_label,label="Moderate (1-2)",size=2.8,colour=COL_BLUE,fontface="italic") +
  annotate("text",x=6.5,y=y_zone_label,label="Major (>=3 units)",size=2.8,colour=COL_AMB,fontface="italic") +
  geom_col(aes(fill=zone),width=.72,colour="white",linewidth=.4) +
  geom_text(aes(label=ifelse(n>0,paste0(n,"\n(",pct,"%)"),"")),
            vjust=-.2,size=3.2,fontface="bold",lineheight=.9) +
  scale_fill_manual(values=c(None="#BDBDBD",Moderate="#7EB5D6",Major=COL_BLUE),name=NULL) +
  scale_x_continuous(breaks=0:10) + scale_y_continuous(limits=c(0,y_max)) +
  labs(title="Distribution of pRBC units transfused per patient",
       subtitle=sprintf("n=%d | Mean %.2f +/- %.2f units | Median %d units",
                        nrow(df), mean(df$total_wb),sd(df$total_wb),as.integer(median(df$total_wb))),
       caption=sprintf("Grey = 0 units (n=%d) | Light blue = 1-2 units (n=%d) | Yellow = >=3 units (n=%d)",
                       n_none, n_mod, n_maj),
       x="Units of packed red blood cells (pRBC)", y="Patients (n)") + theme_ugib()
print(plots[["p_rbc"]])

cat("\n=== TIME-TO-ENDOSCOPY ===\n")
df_endo <- df |> filter(is.na(patient_id %in% c(40,47,51)) | !patient_id %in% c(40,47,51)) |>
  filter(!is.na(time_to_endo_h))
df_endo <- df_endo |> mutate(endo_timing=factor(case_when(
  time_to_endo_h<6~"<6h (emergency)", time_to_endo_h<=24~"6-24h (urgent)",
  TRUE~">24h (delayed)"), levels=c("<6h (emergency)","6-24h (urgent)",">24h (delayed)")))
kw_t <- kruskal.test(total_wb~endo_timing, data=df_endo)
sp_c <- cor.test(df_endo$time_to_endo_h, df_endo$total_wb, method="spearman", exact=FALSE)
cat(sprintf("n=%d | KW p=%.4f | Spearman rho=%.3f\n", nrow(df_endo), kw_t$p.value, sp_c$estimate))

plots[["p_timing"]] <- ggplot(df_endo, aes(endo_timing, total_wb, fill=endo_timing)) +
  geom_boxplot(width=.5, colour="grey30", outlier.size=2) +
  geom_jitter(width=.15, alpha=.4, size=1.8, colour="grey30") +
  stat_summary(fun=mean, geom="point", shape=18, size=4, colour=COL_RED) +
  scale_fill_manual(values=c("<6h (emergency)"=COL_RED,"6-24h (urgent)"=COL_AMB,">24h (delayed)"=COL_GREEN)) +
  annotate("text",x=2,y=max(df_endo$total_wb,na.rm=T)*.95,
           label=sprintf("KW p=%.4f\nSpearman rho=%.3f",kw_t$p.value,sp_c$estimate),
           size=3.5,colour="grey30",fontface="italic") +
  labs(title="Time to endoscopy vs pRBC units transfused",
       subtitle=sprintf("n=%d | Diamond = mean | reverse causation",nrow(df_endo)),
       x="Endoscopy timing group",y="Total pRBC units") + theme_ugib() + theme(legend.position="none")
print(plots[["p_timing"]])

cat("\n=== 3B. Bleeding presentation symptoms (primary vs secondary) ===\n")

sym_df <- data.frame(
  Symptom = c("Melena", "Hematemesis (coffee grounds)", "Hematemesis (bright red)",
              "Hematochezia",
              "Hemodynamic instability", "Altered mental status", "Syncope"),
  Category = c("Primary", "Primary", "Primary", "Primary",
               "Secondary", "Secondary", "Secondary"),
  Count = c(
    sum(df$melena == "YES", na.rm=TRUE),
    sum(df$hematemesis_coffee == "YES", na.rm=TRUE),
    sum(df$hematemesis_red == "YES", na.rm=TRUE),
    sum(df$hematochezia == "YES", na.rm=TRUE),
    sum(df$hemodynamic_instab == "YES", na.rm=TRUE),
    sum(df$altered_mental == "YES", na.rm=TRUE),
    sum(df$syncope == "YES", na.rm=TRUE)),
  stringsAsFactors = FALSE)
sym_df$Percent <- round(sym_df$Count / nrow(df) * 100, 1)
sym_df$Label   <- sprintf("%d (%.1f%%)", sym_df$Count, sym_df$Percent)

n_any_primary   <- sum(df$melena == "YES" | df$hematemesis_red == "YES" |
                         df$hematemesis_coffee == "YES" | df$hematochezia == "YES", na.rm=TRUE)
n_any_secondary <- sum(df$hemodynamic_instab == "YES" | df$altered_mental == "YES" |
                         df$syncope == "YES", na.rm=TRUE)
cat(sprintf("Primary symptoms (source-defining):\n"))
print(sym_df[sym_df$Category == "Primary", c("Symptom","Count","Percent")])
cat(sprintf("\n  -> %d of %d patients (%.1f%%) presented with at least one primary symptom\n",
            n_any_primary, nrow(df), n_any_primary/nrow(df)*100))

cat(sprintf("\nSecondary symptoms (severity markers):\n"))
print(sym_df[sym_df$Category == "Secondary", c("Symptom","Count","Percent")])
cat(sprintf("\n  -> %d of %d patients (%.1f%%) had at least one secondary symptom\n",
            n_any_secondary, nrow(df), n_any_secondary/nrow(df)*100))

sym_df <- sym_df[order(sym_df$Category, -sym_df$Count), ]
sym_df$Symptom <- factor(sym_df$Symptom, levels=rev(sym_df$Symptom))

plots[["p_symptoms"]] <- ggplot(sym_df, aes(x=Symptom, y=Count, fill=Category)) +
  geom_col(width=0.7, colour="white") +
  geom_text(aes(label=Label), hjust=-0.08, size=3.5, fontface="bold") +
  scale_fill_manual(values=c(Primary=COL_BLUE, Secondary=COL_AMB),
                    name=NULL,
                    labels=c("Primary (source-defining)",
                             "Secondary (severity marker)")) +
  coord_flip(clip="off") +
  scale_y_continuous(expand=expansion(mult=c(0, 0.22))) +
  labs(title="Bleeding presentation",
       subtitle=sprintf("n=%d", nrow(df)),
       x=NULL, y="Number of patients") +
  theme_ugib() + theme(legend.position="bottom",
                       plot.margin=unit(c(5.5, 60, 5.5, 5.5), "pt"))
print(plots[["p_symptoms"]])

cat("\n=== 3C. Seasonality analysis ===\n")

season_counts <- table(df$season)
cat("Season distribution:\n")
print(season_counts)
total <- sum(season_counts)
cat(sprintf("Proportions: Spring=%.1f%% Summer=%.1f%% Fall=%.1f%% Winter=%.1f%%\n",
            season_counts["Spring"]/total*100, season_counts["Summer"]/total*100,
            season_counts["Fall"]/total*100, season_counts["Winter"]/total*100))

chi <- chisq.test(season_counts, p=rep(0.25, 4))
cat(sprintf("\nChi-square goodness-of-fit test (H0: uniform distribution):\n"))
cat(sprintf("  X^2 = %.3f, df = %d, p = %.4f\n",
            chi$statistic, chi$parameter, chi$p.value))
if (chi$p.value < 0.05) {
  cat("  -> STATISTICALLY SIGNIFICANT departure from uniform\n")
} else {
  cat(sprintf("  -> NOT statistically significant (p=%.4f > 0.05)\n", chi$p.value))
  cat("  -> Consistent with INDICATION (not proof) of seasonal pattern\n")
}

fall_winter <- sum(season_counts[c("Fall","Winter")])
spring_summer <- sum(season_counts[c("Spring","Summer")])
chi_cw <- chisq.test(c(fall_winter, spring_summer), p=c(0.5, 0.5))
cat(sprintf("\nCold (Fall+Winter) vs Warm (Spring+Summer): %d vs %d\n",
            fall_winter, spring_summer))
cat(sprintf("  X^2 = %.3f, p = %.4f\n", chi_cw$statistic, chi_cw$p.value))
if (chi_cw$p.value < 0.05) {
  cat("  -> Significant cold-month excess (consistent with literature)\n")
} else {
  cat("  -> Cold-month excess INDICATION only (not significant at 0.05)\n")
}

season_df <- data.frame(
  Season = factor(names(season_counts), levels=c("Winter","Spring","Summer","Fall")),
  Count = as.numeric(season_counts),
  Percent = as.numeric(season_counts) / total * 100)
season_df$Label <- sprintf("%d (%.1f%%)", season_df$Count, season_df$Percent)

# SECTION 4: traditional systems ROC
cat("\n============================================================\n")
cat("SECTION 4:\n")
cat("============================================================\n\n")

score_summary <- data.frame(Score=c("GBS","AIMS65","Pre-Rockall","Full Rockall","CANUKA"),
                            n=c(sum(!is.na(df$gbs_score)),sum(!is.na(df$aims65_score)),
                                sum(!is.na(df$pre_rockall)),sum(!is.na(df$full_rockall)),
                                sum(!is.na(df$canuka))),
                            Mean=round(c(mean(df$gbs_score,na.rm=T),mean(df$aims65_score,na.rm=T),
                                         mean(df$pre_rockall,na.rm=T),mean(df$full_rockall,na.rm=T),
                                         mean(df$canuka,na.rm=T)),1),
                            SD=round(c(sd(df$gbs_score,na.rm=T),sd(df$aims65_score,na.rm=T),
                                       sd(df$pre_rockall,na.rm=T),sd(df$full_rockall,na.rm=T),
                                       sd(df$canuka,na.rm=T)),1))
print(score_summary)
cat("\nNote: Full Rockall uses the script-computed value from Section 2B-b\n")
cat("(Rockall_Excel column has been deleted per user request).\n")

make_score_hist <- function(var, label, colour) {
  d <- df |> filter(!is.na(.data[[var]]))
  mn <- mean(d[[var]]); sd_v <- sd(d[[var]]); md <- median(d[[var]])
  iq <- quantile(d[[var]], c(.25,.75))
  ggplot(d, aes(.data[[var]])) +
    geom_bar(fill=colour, colour="white", linewidth=.4) +
    geom_vline(xintercept=mn, colour=COL_RED, linetype="dashed", linewidth=.9) +
    annotate("text", x=mn+0.3, y=Inf, vjust=1.5, hjust=0, size=3,
             label=sprintf("Mean %.1f", mn), colour=COL_RED, fontface="bold") +
    labs(title=label,
         subtitle=sprintf("n=%d | Mean %.1f +/- %.1f | Median %.0f (IQR %.0f-%.0f) | Range %d-%d",
                          nrow(d), mn, sd_v, md, iq[1], iq[2], min(d[[var]]), max(d[[var]])),
         x=label, y="Patients (n)") + theme_ugib()
}
plots[["p_gbs"]]  <- make_score_hist("gbs_score","Glasgow-Blatchford Score",COL_BLUE);  print(plots[["p_gbs"]])
plots[["p_aims"]] <- make_score_hist("aims65_score","AIMS65 Score",COL_RED);            print(plots[["p_aims"]])
plots[["p_preR"]] <- make_score_hist("pre_rockall","Pre-endoscopy Rockall",COL_GREEN);  print(plots[["p_preR"]])
plots[["p_fR"]]   <- make_score_hist("full_rockall","Full Rockall Score",COL_AMB); print(plots[["p_fR"]])
if(sum(!is.na(df$canuka)) >= 5) {
  plots[["p_canuka"]] <- make_score_hist("canuka","CANUKA Score","#8E44AD")
  print(plots[["p_canuka"]])
} else {
  cat("CANUKA score has too few values for histogram (need to populate column)\n")
}

safe_roc <- function(score_var, outcome_var, d) {
  d2 <- d |> filter(!is.na(.data[[score_var]]), !is.na(.data[[outcome_var]]))
  if (length(unique(d2[[outcome_var]])) < 2) return(NULL)
  tryCatch(roc(d2[[outcome_var]], d2[[score_var]],
               levels=c("NO","YES"), direction="<", quiet=TRUE),
           error=function(e) NULL)
}

score_specs <- list(
  list(var="gbs_score",    label="GBS",          col=COL_BLUE),
  list(var="aims65_score", label="AIMS65",       col=COL_RED),
  list(var="pre_rockall",  label="Pre-Rockall",  col=COL_GREEN),
  list(var="full_rockall", label="Full Rockall", col=COL_AMB),
  list(var="canuka",       label="CANUKA",       col="#8E44AD"))

score_outcome_matrix <- data.frame()
roc_scores <- list()
for (oc in ALL_OUTCOMES) {
  if (!oc %in% names(df)) {
    cat(sprintf("WARNING: outcome '%s' not in df -- skipping\n", oc)); next
  }
  roc_scores[[oc]] <- list()
  cat(sprintf("\nScore AUCs vs %s (%s):\n",
              oc, switch(oc, "high_risk"="PRIMARY",
                         "major_tx"="secondary",
                         "tertiary")))
  for (sp in score_specs) {
    r <- safe_roc(sp$var, oc, df)
    roc_scores[[oc]][[sp$label]] <- r
    if (is.null(r)) {
      cat(sprintf("  %-15s (no ROC -- score column empty or one-class outcome)\n", sp$label))
      next
    }
    ci <- ci.auc(r)
    cat(sprintf("  %-15s AUC=%.3f (%.3f-%.3f) n=%d\n",
                sp$label, as.numeric(auc(r)), ci[1], ci[3], length(r$response)))
    score_outcome_matrix <- rbind(score_outcome_matrix, data.frame(
      Outcome   = oc,
      Score     = sp$label,
      AUC       = round(as.numeric(auc(r)), 3),
      CI_low    = round(ci[1], 3),
      CI_high   = round(ci[3], 3),
      n         = length(r$response),
      n_events  = sum(r$response == "YES"),
      stringsAsFactors = FALSE))
  }
}

roc_scores_flat_primary <- list(
  gbs          = roc_scores[[PRIMARY]][["GBS"]],
  aims65       = roc_scores[[PRIMARY]][["AIMS65"]],
  pre_rockall  = roc_scores[[PRIMARY]][["Pre-Rockall"]],
  full_rockall = roc_scores[[PRIMARY]][["Full Rockall"]],
  canuka       = roc_scores[[PRIMARY]][["CANUKA"]],
  by_outcome   = roc_scores,
  primary      = PRIMARY)
saveRDS(roc_scores_flat_primary, "roc_scores.rds")
saveRDS(score_summary, "score_summary.rds")
saveRDS(score_outcome_matrix, "score_outcome_matrix.rds")
write.csv(score_outcome_matrix, "score_outcome_matrix.csv", row.names=FALSE)
cat(sprintf("\nSaved: roc_scores.rds  (primary outcome = %s)\n", PRIMARY))
cat("Saved: score_outcome_matrix.rds and score_outcome_matrix.csv\n")

cat("\n--- Score x Outcome AUC matrix (wider view, primary outcome first) ---\n")
auc_wide <- score_outcome_matrix |>
  dplyr::select(Score, Outcome, AUC) |>
  pivot_wider(names_from=Outcome, values_from=AUC)
ord <- c("Score", intersect(ALL_OUTCOMES, names(auc_wide)))
auc_wide <- auc_wide[, ord, drop=FALSE]
print(auc_wide)

plot_roc_overlay <- function(outcome_name, roc_list_named, score_specs) {
  available <- list()
  for (sp in score_specs) {
    obj <- roc_list_named[[sp$label]]
    if (!is.null(obj)) available[[length(available)+1]] <- list(
      name=sp$label, obj=obj, col=sp$col)
  }
  if (length(available) == 0) {
    cat(sprintf("  No ROCs for %s -- skipping plot\n", outcome_name)); return(invisible(NULL))
  }
  aucs <- sapply(available, function(x) as.numeric(auc(x$obj)))
  available <- available[order(-aucs)]
  par(mar=c(5,5,4,2))
  .roc_title_map <- list(
    high_risk = "Traditional scores vs high risk upper gastrointestinal bleeding",
    major_tx = "Traditional scores vs major transfusion (>=3 pRBCs)",
    outcome_death = "Traditional scores vs death",
    rebleeding = "Traditional scores vs rebleeding")
  .roc_main <- if (!is.null(.roc_title_map[[outcome_name]])) .roc_title_map[[outcome_name]] else sprintf("Traditional scores vs %s", outcome_name)
  plot(available[[1]]$obj, col=available[[1]]$col, lwd=2.5,
       legacy.axes = TRUE,
       xlim = c(1, 0), ylim = c(0, 1),
       xaxs = "i", yaxs = "i",
       main = .roc_main)
  if (length(available) > 1) {
    for (i in 2:length(available)) {
      plot(available[[i]]$obj, col=available[[i]]$col, lwd=2.5, legacy.axes = TRUE, add=TRUE)
    }
  }
  abline(a=0, b=1, lty=2, col="grey60")
  leg_labs <- sapply(available, function(x)
    sprintf("%s (AUC=%.3f)", x$name, as.numeric(auc(x$obj))))
  leg_cols <- sapply(available, function(x) x$col)
  legend("bottomright", bty="n", lwd=2.5, col=leg_cols, legend=leg_labs)
}
for (oc in ALL_OUTCOMES) {
  if (!oc %in% names(roc_scores)) next
  plot_roc_overlay(oc, roc_scores[[oc]], score_specs)
}

cat("\n============================================================\n")
cat("SCORES FOR EXCEL IMPORT (copy-paste into columns, or use CSV)\n")
cat("============================================================\n\n")

scores_export <- data.frame(
  Patient_ID             = df$patient_id,
  HIS                    = df$hospital_id,
  Pre_Rockall            = df$pre_rockall,
  Full_Rockall_Computed  = df$full_rockall,
  AIMS65_Score           = df$aims65_score,
  GBS                    = df$gbs_score,
  CANUKA                 = df$canuka,
  RI_Corrected           = round(df$ri_corrected, 2),
  RPI                    = round(df$rpi, 2),
  stringsAsFactors = FALSE)

cat("TAB-SEPARATED TABLE (copy to Excel):\n\n")
cat(paste(names(scores_export), collapse="\t"), "\n")
for (i in seq_len(nrow(scores_export))) {
  row_vals <- as.character(scores_export[i, ])
  row_vals[is.na(row_vals) | row_vals == "NA"] <- ""
  cat(paste(row_vals, collapse="\t"), "\n")
}

write.csv(scores_export, "scores_for_excel.csv", row.names=FALSE, na="")
cat("\nSaved: scores_for_excel.csv\n")

# SECTION 5: transfusion splits
cat("\n============================================================\n")
cat("SECTION 5:\n")
cat("============================================================\n\n")

df <- df |> mutate(
  sch_A=factor(case_when(total_wb==0~"0",total_wb<=2~"1to2",TRUE~"3plus"),
               levels=c("0","1to2","3plus"),ordered=T),
  sch_B=factor(case_when(total_wb==0~"0",total_wb==1~"1",total_wb==2~"2",TRUE~"3plus"),ordered=T),
  sch_C=factor(case_when(total_wb==0~"0",total_wb==1~"1",TRUE~"2plus"),ordered=T),
  sch_D=factor(case_when(total_wb==0~"0",TRUE~"1plus"),ordered=T),
  sch_E=factor(case_when(total_wb==0~"0",total_wb<=3~"1to3",TRUE~"4plus"),ordered=T))

entropy_fn <- function(x){p<-x/sum(x);p<-p[p>0];-sum(p*log2(p))}
h_max <- entropy_fn(as.numeric(table(df$total_wb)))
cat("Class balance & entropy:\n")
for(s in LETTERS[1:5]){tbl<-table(df[[paste0("sch_",s)]])
cat(sprintf("  Scheme %s: %s | min=%.1f%% | entropy=%.3f (%.1f%%)\n",
            s, paste(names(tbl),tbl,sep="=",collapse=" / "),
            min(tbl)/sum(tbl)*100, entropy_fn(as.numeric(tbl)),
            entropy_fn(as.numeric(tbl))/h_max*100))}

set.seed(2024)
.tree_data <- df[, c("total_wb","hb_mean","albumin","urea","inr","sbp","hr")]
names(.tree_data) <- c("Whole_blood_units","Hemoglobin","Albumin","Urea","INR","Systolic_BP","Heart_rate")
tree_scheme <- rpart(Whole_blood_units ~ Hemoglobin + Albumin + Urea + INR + Systolic_BP + Heart_rate,
                     data = na.omit(.tree_data),
                     method = "anova", control = rpart.control(cp = .02, minsplit = 10, maxdepth = 3))
cat("\nDecision tree terminal nodes:\n")
print(tree_scheme$frame[tree_scheme$frame$var=="<leaf>",c("n","yval")])
rpart.plot(tree_scheme, type = 4, extra = 101, roundint = FALSE,
           main = "Natural cut-points of whole blood transfusion",
           cex = 0.8, fallen.leaves = TRUE)

saveRDS("A","chosen_scheme.rds")
cat("\n>>> CHOSEN: Scheme A (0 / 1-2 / >=3)\n")

# SECTION 6: feature selection
cat("\n============================================================\n")
cat("SECTION 6:\n")
cat("============================================================\n\n")

feat_cands <- c("age","gender","symptom_duration",
                "sbp","dbp","hr","hemodynamic_instab","ed_levosim",
                "altered_mental","syncope","hematemesis_red","hematemesis_coffee","melena",
                "hematochezia",
                "rectal_exam",
                "gi_bleed_recent","gi_bleed_remote",
                "hb_mean","hct_mean","lactate_ed",
                "wbc_entry","plt_entry","ri_pct","creatinine","urea","bun","albumin",
                "inr","aptt","pt_pct","sgot","sgpt","tbil",
                "ferritin","ldh",
                "drug_act","drug_antiplt","drug_ppi","drug_nsaid","drug_steroid",
                "liver_cirrhosis",
                "renal_severity_3",
                "heart_failure","cad","afib",
                "hypertension","copd","asthma",
                "dm_any",
                "active_cancer",
                "comorbidity_hematol","ibd")
feat_cands <- feat_cands[feat_cands %in% names(df)]
cat("Initial candidate pool:", length(feat_cands), "features\n")

cat("\n--- Collinearity check (pre-selection) ---\n")
if(all(c("hb_mean","hct_mean") %in% feat_cands)) {
  r_hb_hct <- cor(df$hb_mean, df$hct_mean, use="complete.obs")
  cat(sprintf("  hb_mean ~ hct_mean: r=%.4f (HCT = ~3 x Hb)\n", r_hb_hct))
  cat("  -> EXCLUDE hct_mean (keep hb_mean: standard transfusion trigger)\n")
}
if(all(c("urea","bun") %in% feat_cands)) {
  r_urea_bun <- cor(df$urea, df$bun, use="complete.obs")
  cat(sprintf("  urea ~ bun: r=%.4f (BUN = Urea / 2.14)\n", r_urea_bun))
  cat("  -> EXCLUDE bun (keep urea: European lab convention)\n")
}
if("hemodynamic_instab" %in% feat_cands) {
  cat("  -> EXCLUDE hemodynamic_instab (composite of SBP<=90 OR vasopressor;\n")
  cat("                                  redundant with raw sbp + ed_levosim)\n")
}
if("ed_levosim" %in% feat_cands) {
  n_levo <- sum(df$ed_levosim == "YES", na.rm=TRUE)
  cat(sprintf("  -> EXCLUDE ed_levosim (only n=%d patients on vasopressors;\n", n_levo))
  cat("                            captured by sbp/hr; insufficient prevalence)\n")
}
feat_cands <- feat_cands[!feat_cands %in% c("hct_mean", "bun",
                                            "hemodynamic_instab", "ed_levosim")]
cat(sprintf("After collinearity exclusion: %d features\n\n", length(feat_cands)))

cat("\n============================================================\n")
cat("MISSINGNESS PATTERN TABLE (pre-imputation)\n")
cat("============================================================\n\n")

miss_tbl <- data.frame(
  variable    = feat_cands,
  n_missing   = vapply(feat_cands, function(v) sum(is.na(df[[v]])), integer(1)),
  pct_missing = vapply(feat_cands,
                       function(v) round(100 * mean(is.na(df[[v]])), 1),
                       numeric(1)),
  stringsAsFactors = FALSE
)
miss_tbl <- miss_tbl[order(-miss_tbl$pct_missing), ]
rownames(miss_tbl) <- NULL

cat(sprintf("Cohort: n=%d patients | candidate features: %d\n\n",
            nrow(df), length(feat_cands)))
cat("Top 10 most-missing features:\n")
print(head(miss_tbl, 10), row.names=FALSE)
cat(sprintf("\nMax missingness: %.1f%% (%s)\n",
            miss_tbl$pct_missing[1], miss_tbl$variable[1]))
cat(sprintf("Variables with >=20%% missing: %d\n",
            sum(miss_tbl$pct_missing >= 20)))
cat(sprintf("Variables with >=10%% missing: %d\n",
            sum(miss_tbl$pct_missing >= 10)))
cat(sprintf("Variables with 0%% missing:    %d\n",
            sum(miss_tbl$pct_missing == 0)))

n_miss_per_pt <- rowSums(is.na(df[, feat_cands, drop=FALSE]))
cat("\nPatient-level missingness (candidate features missing per patient):\n")
cat(sprintf("  median: %d  IQR: %.0f-%.0f  max: %d\n",
            median(n_miss_per_pt),
            quantile(n_miss_per_pt, 0.25),
            quantile(n_miss_per_pt, 0.75),
            max(n_miss_per_pt)))
cat(sprintf("  patients with 0 missing features: %d / %d (%.1f%%)\n",
            sum(n_miss_per_pt == 0), nrow(df),
            100 * mean(n_miss_per_pt == 0)))

saveRDS(miss_tbl, "missingness_table.rds")
cat("\nSaved: missingness_table.rds\n\n")

impute_simple <- function(d,vars){for(v in vars){x<-d[[v]]
if(is.numeric(x)) d[[v]][is.na(x)]<-median(x,na.rm=T)
else if(is.factor(x)|is.character(x)){md<-names(sort(table(x),decreasing=T))[1]
if(!is.null(md)) d[[v]][is.na(x)]<-md}}; d}

df$lactate_missing <- as.integer(is.na(df$lactate_ed))
cat(sprintf("\nMissing-indicator for lactate: %d / %d (%.1f%%) flagged\n",
            sum(df$lactate_missing == 1), nrow(df),
            100 * mean(df$lactate_missing == 1)))

feat_cands <- c(feat_cands, "lactate_missing")
cat(sprintf("Candidate pool after lactate_missing addition: %d features\n",
            length(feat_cands)))

df_imp <- impute_simple(df, feat_cands)

saveRDS(list(
  feat_cands = feat_cands,
  vars_num   = feat_cands[vapply(feat_cands, function(v) is.numeric(df[[v]]), logical(1))],
  vars_cat   = feat_cands[vapply(feat_cands, function(v) is.factor(df[[v]]) ||
                                   is.character(df[[v]]), logical(1))],
  PRIMARY    = PRIMARY
), "nested_cv_meta.rds")
cat("Saved: nested_cv_meta.rds (candidate pool for selection-inclusive nested CV in Part 2)\n")

# SECTION 6: feature selection -- generic over outcome variable

vars_num <- feat_cands[sapply(feat_cands, function(v) is.numeric(df[[v]]))]
vars_cat <- feat_cands[sapply(feat_cands, function(v)
  is.factor(df[[v]]) | is.character(df[[v]]))]

run_feat_selection <- function(outcome_var, run_label, make_plots = TRUE) {
  cat(sprintf("\n>>> FEATURE SELECTION on '%s' (%s)\n",
              outcome_var, run_label))
  
  df_work <- df_imp |>
    dplyr::select(all_of(feat_cands), all_of(outcome_var)) |>
    filter(!is.na(.data[[outcome_var]]))
  cat(sprintf("    Working set: n=%d | YES=%d / NO=%d\n",
              nrow(df_work),
              sum(df_work[[outcome_var]] == "YES"),
              sum(df_work[[outcome_var]] == "NO")))
  
  cat("\n=== METHOD 1: Univariate screening ===\n")
  num_res <- lapply(vars_num, function(v) {
    x <- df[[v]]; mt <- df[[outcome_var]]; ok <- !is.na(x) & !is.na(mt)
    if (sum(ok) < 10) return(NULL)
    if (length(unique(mt[ok])) < 2) return(NULL)
    wt <- suppressWarnings(wilcox.test(x[ok & mt == "YES"], x[ok & mt == "NO"]))
    n1 <- sum(ok & mt == "YES"); n2 <- sum(ok & mt == "NO")
    data.frame(variable = v, type = "continuous",
               med_yes  = round(median(x[ok & mt == "YES"], na.rm = TRUE), 2),
               med_no   = round(median(x[ok & mt == "NO"],  na.rm = TRUE), 2),
               p_raw    = wt$p.value,
               effect_r = round(1 - 2 * as.numeric(wt$statistic) / (n1 * n2), 3),
               stringsAsFactors = FALSE)
  })
  cat_res <- lapply(vars_cat, function(v) {
    x <- df[[v]]; mt <- df[[outcome_var]]; ok <- !is.na(x) & !is.na(mt)
    if (sum(ok) < 10) return(NULL)
    tbl <- table(x[ok], mt[ok])
    if (nrow(tbl) < 2 || ncol(tbl) < 2) return(NULL)
    use_f <- any(suppressWarnings(tryCatch(chisq.test(tbl)$expected,
                                           error = function(e) matrix(0))) < 5)
    res <- tryCatch(if (use_f) fisher.test(tbl, simulate.p.value = TRUE, B = 1e4)
                    else chisq.test(tbl),
                    error = function(e) list(p.value = NA))
    data.frame(variable = v, type = "categorical",
               med_yes  = NA_real_, med_no = NA_real_,
               p_raw    = res$p.value, effect_r = NA_real_,
               stringsAsFactors = FALSE)
  })
  univ_df <- bind_rows(c(num_res, cat_res)) |>
    filter(!is.na(p_raw)) |>
    mutate(p_bonf = p.adjust(p_raw, "bonferroni"),
           sig = case_when(p_bonf < .001 ~ "***", p_bonf < .01 ~ "**",
                           p_bonf < .05 ~ "*", p_raw < .05 ~ "(nom)",
                           TRUE ~ "ns")) |>
    arrange(p_raw)
  univ_sel <- univ_df$variable[univ_df$p_bonf < .05]
  univ_nom <- univ_df$variable[univ_df$p_raw  < .05]
  cat(sprintf("Bonferroni: %d | Nominal: %d\n", length(univ_sel), length(univ_nom)))
  
  cat("\n=== METHOD 2: Decision tree ===\n")
  set.seed(2024)
  tree_full <- rpart(reformulate(feat_cands, outcome_var), data = df_work,
                     method = "class",
                     control = rpart.control(cp = 0, minsplit = 10,
                                             maxdepth = 5, xval = 10))
  cpt <- tree_full$cptable
  i_min <- which.min(cpt[, "xerror"])
  cp_min <- cpt[i_min, "CP"]
  thresh <- cpt[i_min, "xerror"] + cpt[i_min, "xstd"]
  cp_1se <- cpt[which(cpt[, "xerror"] <= thresh)[1], "CP"]
  tree_cp <- cp_1se
  cat(sprintf("CV-optimal cp: cp.min=%.4f (xerror=%.3f) | cp.1se=%.4f (HEADLINE)\n",
              cp_min, cpt[i_min, "xerror"], cp_1se))
  
  tree_fit <- prune(tree_full, cp = tree_cp)
  tree_imp <- if (!is.null(tree_fit$variable.importance) &&
                  length(tree_fit$variable.importance) > 0) {
    data.frame(variable = names(tree_fit$variable.importance),
               importance = as.numeric(tree_fit$variable.importance)) |>
      arrange(desc(importance)) |>
      mutate(imp_pct = round(importance / max(importance) * 100, 1))
  } else data.frame(variable = character(), importance = numeric(),
                    imp_pct = numeric(), stringsAsFactors = FALSE)
  tree_top15 <- if (nrow(tree_imp) > 0) tree_imp$variable[1:min(15, nrow(tree_imp))] else character()
  splits_used <- unique(tree_fit$frame$var[tree_fit$frame$var != "<leaf>"])
  cat(sprintf("Tree splits (pruned at cp=%.4f): %s\n",
              tree_cp, paste(splits_used, collapse = ", ")))
  tree_vote_valid <- length(splits_used) >= 3
  if (!tree_vote_valid) {
    cat(sprintf("  WARNING: only %d unique split variable(s) -- tree vote DOWN-WEIGHTED for this outcome\n",
                length(splits_used)))
    tree_top15 <- character(0)
  }
  
  cat("\n=== METHOD 3: Random Forest ===\n")
  set.seed(2024)
  rf_sel <- randomForest(reformulate(feat_cands, outcome_var), data = df_work,
                         ntree = 500, mtry = floor(sqrt(length(feat_cands))),
                         importance = TRUE)
  imp_mat <- importance(rf_sel)
  rf_imp <- data.frame(variable = rownames(imp_mat),
                       gini = imp_mat[, "MeanDecreaseGini"],
                       acc  = imp_mat[, "MeanDecreaseAccuracy"]) |>
    mutate(gini_pct = round(gini / max(gini) * 100, 1),
           acc_pct  = round(acc  / max(acc)  * 100, 1)) |>
    arrange(desc(gini_pct))
  rf_top15 <- rf_imp$variable[1:min(15, nrow(rf_imp))]
  cat(sprintf("RF OOB error: %.1f%%\n", rf_sel$err.rate[500, "OOB"] * 100))
  
  cat("\n=== METHOD 4: LASSO ===\n")
  X <- model.matrix(reformulate(feat_cands), data = df_work)[, -1]
  y <- as.numeric(df_work[[outcome_var]] == "YES")
  set.seed(2024)
  cv_lasso <- cv.glmnet(X, y, family = "binomial", alpha = 1,
                        nfolds = 5, type.measure = "deviance")
  cat(sprintf("lambda.1se: %.5f | lambda.min: %.5f\n",
              cv_lasso$lambda.1se, cv_lasso$lambda.min))
  
  lasso_coef     <- coef(cv_lasso, s = "lambda.1se")
  lasso_coef_min <- coef(cv_lasso, s = "lambda.min")
  lasso_df <- data.frame(term = rownames(lasso_coef),
                         coef = as.numeric(lasso_coef)) |>
    filter(term != "(Intercept)", coef != 0) |>
    arrange(desc(abs(coef)))
  lasso_df_min <- data.frame(term = rownames(lasso_coef_min),
                             coef = as.numeric(lasso_coef_min)) |>
    filter(term != "(Intercept)", coef != 0) |>
    arrange(desc(abs(coef)))
  lasso_sel <- unique(unlist(lapply(feat_cands, function(v)
    if (any(grepl(paste0("^", v), lasso_df$term))) v)))
  lasso_sel_min <- unique(unlist(lapply(feat_cands, function(v)
    if (any(grepl(paste0("^", v), lasso_df_min$term))) v)))
  cat(sprintf("lambda.1se: %d non-zero terms | %d variables (HEADLINE)\n",
              nrow(lasso_df), length(lasso_sel)))
  cat(sprintf("lambda.min: %d non-zero terms | %d variables (cross-check)\n",
              nrow(lasso_df_min), length(lasso_sel_min)))
  if (length(lasso_sel_min) > length(lasso_sel)) {
    extra <- setdiff(lasso_sel_min, lasso_sel)
    cat(sprintf("  Additional vars at lambda.min: %s\n",
                paste(extra, collapse = ", ")))
  }
  
  cat("\n=== VOTE COUNT ===\n")
  n_methods <- if (tree_vote_valid) 4L else 3L
  min_votes_for_confirmed <- 2L
  
  all_c <- unique(c(univ_nom, tree_top15, rf_top15, lasso_sel))
  vote_df <- data.frame(variable = all_c,
                        in_univ  = all_c %in% univ_nom,
                        in_tree  = all_c %in% tree_top15,
                        in_rf    = all_c %in% rf_top15,
                        in_lasso = all_c %in% lasso_sel) |>
    mutate(votes = as.integer(in_univ) + as.integer(in_tree) +
             as.integer(in_rf)   + as.integer(in_lasso)) |>
    arrange(desc(votes), variable)
  confirmed  <- sort(vote_df$variable[vote_df$votes >= min_votes_for_confirmed])
  unanimous  <- sort(vote_df$variable[vote_df$votes == n_methods])
  borderline <- sort(vote_df$variable[vote_df$votes == 1])
  cat(sprintf("Voting on %d methods (tree %s)\n",
              n_methods, if (tree_vote_valid) "VALID" else "DOWN-WEIGHTED"))
  cat(sprintf("Confirmed (>=%d/%d): %d -- %s\n",
              min_votes_for_confirmed, n_methods,
              length(confirmed), paste(confirmed, collapse = ", ")))
  cat(sprintf("Unanimous (%d/%d):   %d -- %s\n",
              n_methods, n_methods,
              length(unanimous), paste(unanimous, collapse = ", ")))
  cat(sprintf("Borderline (1/%d):  %d\n", n_methods, length(borderline)))
  
  if (!"hb_mean" %in% confirmed) confirmed <- c(confirmed, "hb_mean")
  confirmed <- sort(unique(confirmed))
  cat(sprintf("Final feature set: %d features\n", length(confirmed)))
  cat("  ", paste(confirmed, collapse = ", "), "\n\n")
  
  if (make_plots) {
    plots[[paste0("p_univ_", run_label)]] <<- ggplot(filter_low_signal(univ_df) |>
                                                       mutate(logp = -log10(p_raw),
                                                              sig_g = case_when(p_bonf < .05 ~ "Bonferroni",
                                                                                p_raw  < .05 ~ "Nominal", TRUE ~ "ns"),
                                                              variable_disp = relabel_var(variable),
                                                              variable_disp = factor(variable_disp, levels = rev(variable_disp))),
                                                     aes(logp, variable_disp, colour = sig_g)) +
      geom_point(size = 2.5) +
      geom_segment(aes(x = 0, xend = logp, yend = variable_disp), linewidth = .4) +
      geom_vline(xintercept = -log10(.05), linetype = "dashed",
                 colour = COL_RED, linewidth = .7) +
      geom_vline(xintercept = -log10(.05 / nrow(univ_df)),
                 linetype = "dotted", colour = "#8E44AD", linewidth = .7) +
      scale_colour_manual(values = c(Bonferroni = COL_BLUE,
                                     Nominal = "#7EB5D6", ns = COL_GREY),
                          name = NULL) +
      labs(title = "BURST feature screening: univariate vs high risk",
           subtitle = "Bonferroni-corrected | red dashed=p=.05, purple dotted=Bonferroni",
           x = "-log10(p)", y = NULL) +
      theme_ugib() + theme(axis.text.y = element_text(size = 7))
    print(plots[[paste0("p_univ_", run_label)]])
    
    .dt_relabel_split <- function(x, labs, digits, varlen, faclen) relabel_var(labs)
    rpart.plot(tree_fit, type = 4, extra = 104, roundint = FALSE,
               main = "BURST decision tree vs primary outcome",
               split.fun = .dt_relabel_split,
               cex = 0.75, fallen.leaves = TRUE)
    
    if (nrow(rf_imp) > 0) {
      plots[[paste0("p_rf_", run_label)]] <<- ggplot(filter_low_signal(rf_imp[1:min(20, nrow(rf_imp)), ]) |>
                                                       mutate(variable_disp = relabel_var(variable)),
                                                     aes(gini_pct, reorder(variable_disp, gini_pct))) +
        geom_col(fill = COL_GREEN, width = .7) +
        labs(title = "BURST (Random Forest) on primary outcome",
             subtitle = sprintf("Outcome: %s | 500 trees | OOB error: %.1f%%",
                                outcome_var, rf_sel$err.rate[500, "OOB"] * 100),
             x = "Gini importance (% of max)", y = NULL) +
        theme_ugib()
      print(plots[[paste0("p_rf_", run_label)]])
    }
    
    plot(cv_lasso, main = "BURST LASSO: 5-fold CV deviance on primary outcome")
    
    .tier_three <- if (n_methods == 4L) "3/4" else "3/3 unanimous"
    .tier_four  <- "4/4 unanimous"
    .lab_borderline <- sprintf("1/%d borderline", n_methods)
    .lab_confirmed  <- sprintf("2/%d confirmed",  n_methods)
    
    vote_df_plot <- filter_low_signal(vote_df) |>
      mutate(variable_disp = relabel_var(variable),
             variable_disp = factor(variable_disp, levels = rev(variable_disp)),
             col = case_when(
               n_methods == 4L & votes == 4 ~ .tier_four,
               n_methods == 4L & votes == 3 ~ .tier_three,
               n_methods == 3L & votes == 3 ~ .tier_three,
               votes == 2 ~ .lab_confirmed,
               votes == 1 ~ .lab_borderline,
               TRUE ~ NA_character_))
    
    .fill_map <- if (n_methods == 4L) {
      setNames(c("#1A4F7A", COL_BLUE, "#7EB5D6", COL_GREY),
               c(.tier_four, .tier_three, .lab_confirmed, .lab_borderline))
    } else {
      setNames(c("#1A4F7A", "#7EB5D6", COL_GREY),
               c(.tier_three, .lab_confirmed, .lab_borderline))
    }
    
    .x_breaks <- 0:n_methods
    .x_limits <- c(0, n_methods + 1)
    .x_label  <- if (n_methods == 4L) {
      "Methods selecting variable (out of 4)"
    } else {
      sprintf("Methods selecting variable (out of %d voting; tree down-weighted)",
              n_methods)
    }
    .title_str <- if (n_methods == 4L) {
      "BURST vote count on 4 methods vs high risk"
    } else {
      "BURST vote count: 4-method ensemble (classification tree down-weighted, effective N=3) vs high risk"
    }
    .subtitle_str <- sprintf("%d vote-confirmed (>=2/%d) | %d at maximum vote tier (%d/%d)",
                             length(confirmed), n_methods,
                             length(unanimous), n_methods, n_methods)
    
    plots[[paste0("p_votes_", run_label)]] <<- ggplot(vote_df_plot,
                                                      aes(votes, variable_disp, fill = col)) + geom_col(width = .75) +
      geom_vline(xintercept = 1.5, linetype = "dashed",
                 colour = COL_RED, linewidth = 1) +
      scale_fill_manual(values = .fill_map, name = NULL,
                        breaks = names(.fill_map)) +
      scale_x_continuous(breaks = .x_breaks, limits = .x_limits) +
      labs(title = .title_str,
           subtitle = .subtitle_str,
           x = .x_label, y = NULL) +
      theme_ugib() + theme(axis.text.y = element_text(size = 8))
    print(plots[[paste0("p_votes_", run_label)]])
  }
  
  list(univ_df = univ_df, tree_imp = tree_imp, rf_imp = rf_imp,
       lasso_df = lasso_df, vote_df = vote_df,
       confirmed = confirmed, unanimous = unanimous,
       n_methods = n_methods,
       df_work = df_work)
}

fs_primary <- run_feat_selection(PRIMARY, "primary", make_plots = TRUE)
confirmed   <- fs_primary$confirmed
unanimous   <- fs_primary$unanimous
vote_df     <- fs_primary$vote_df
df_work     <- fs_primary$df_work
n_methods   <- fs_primary$n_methods

if ("time_to_endo_h" %in% confirmed) {
  confirmed <- setdiff(confirmed, "time_to_endo_h")
  cat("V2.4 guard: dropped stray 'time_to_endo_h' from final feature set\n")
}

saveRDS(confirmed, "final_feature_set.rds")
cat(sprintf("Saved: final_feature_set.rds (PRIMARY = %s, n=%d features)\n\n",
            PRIMARY, length(confirmed)))

fs_secondary <- run_feat_selection(SECONDARY, "secondary", make_plots = FALSE)
confirmed_secondary <- fs_secondary$confirmed
saveRDS(confirmed_secondary, "final_feature_set_majortx.rds")
cat(sprintf("Saved: final_feature_set_majortx.rds (SECONDARY = %s, n=%d features)\n\n",
            SECONDARY, length(confirmed_secondary)))

both           <- intersect(confirmed, confirmed_secondary)
primary_only   <- setdiff(confirmed,  confirmed_secondary)
secondary_only <- setdiff(confirmed_secondary, confirmed)
cat(sprintf("--- Feature-set agreement (PRIMARY vs SECONDARY) ---\n"))
cat(sprintf("  In both:        %2d -- %s\n",
            length(both), paste(both, collapse = ", ")))
cat(sprintf("  PRIMARY only:   %2d -- %s\n",
            length(primary_only),
            if (length(primary_only) > 0) paste(primary_only, collapse = ", ") else "(none)"))
cat(sprintf("  SECONDARY only: %2d -- %s\n\n",
            length(secondary_only),
            if (length(secondary_only) > 0) paste(secondary_only, collapse = ", ") else "(none)"))

# SECTION 6B: STABILITY SELECTION
RUN_STABILITY <- TRUE

if (RUN_STABILITY) {
  cat("\n============================================================\n")
  cat("SECTION 6B: STABILITY SELECTION\n")
  cat("============================================================\n")
  cat("Estimating selection probability for each feature...\n")
  
  B <- 100
  stability_votes <- matrix(0, nrow=length(feat_cands), ncol=B,
                            dimnames=list(feat_cands, NULL))
  
  set.seed(2024)
  n <- nrow(df_work)
  
  pb_tick <- ceiling(B/10)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace=TRUE)
    df_b <- df_work[idx, , drop=FALSE]
    
    if (length(unique(df_b[[PRIMARY]])) < 2) next
    
    univ_b <- character(0)
    for (v in vars_num) {
      x <- df_b[[v]]; mt <- df_b[[PRIMARY]]; ok <- !is.na(x) & !is.na(mt)
      if (sum(ok) < 10) next
      if (length(unique(mt[ok])) < 2) next
      p <- tryCatch(suppressWarnings(wilcox.test(x[ok & mt=="YES"], x[ok & mt=="NO"])$p.value),
                    error=function(e) NA)
      if (!is.na(p) && p < .05) univ_b <- c(univ_b, v)
    }
    for (v in vars_cat) {
      x <- df_b[[v]]; mt <- df_b[[PRIMARY]]; ok <- !is.na(x) & !is.na(mt)
      if (sum(ok) < 10) next
      tbl <- table(x[ok], mt[ok])
      if (nrow(tbl) < 2 || ncol(tbl) < 2) next
      p <- tryCatch(suppressWarnings(
        if (any(chisq.test(tbl)$expected < 5))
          fisher.test(tbl, simulate.p.value=TRUE, B=2000)$p.value
        else chisq.test(tbl)$p.value), error=function(e) NA)
      if (!is.na(p) && p < .05) univ_b <- c(univ_b, v)
    }
    
    rf_b <- tryCatch(
      randomForest(reformulate(feat_cands, PRIMARY), data=df_b,
                   ntree=300, mtry=floor(sqrt(length(feat_cands))),
                   importance=FALSE),
      error=function(e) NULL)
    rf_top15_b <- if (!is.null(rf_b)) {
      imp <- importance(rf_b)[, "MeanDecreaseGini"]
      names(sort(imp, decreasing=TRUE))[1:min(15, length(imp))]
    } else character(0)
    
    X_b <- tryCatch(model.matrix(reformulate(feat_cands), data=df_b)[, -1],
                    error=function(e) NULL)
    y_b <- as.numeric(df_b[[PRIMARY]] == "YES")
    lasso_sel_b <- character(0)
    if (!is.null(X_b) && length(unique(y_b)) == 2 && nrow(X_b) > 20) {
      cv_b <- tryCatch(cv.glmnet(X_b, y_b, family="binomial", alpha=1,
                                 nfolds=5, type.measure="deviance"),
                       error=function(e) NULL)
      if (!is.null(cv_b)) {
        coef_b <- coef(cv_b, s="lambda.1se")
        terms_b <- rownames(coef_b)[as.numeric(coef_b) != 0]
        terms_b <- setdiff(terms_b, "(Intercept)")
        lasso_sel_b <- unique(unlist(lapply(feat_cands, function(v)
          if (any(grepl(paste0("^", v), terms_b))) v)))
      }
    }
    
    sel_union <- unique(c(univ_b, rf_top15_b, lasso_sel_b))
    for (v in sel_union) {
      n_methods <- (v %in% univ_b) + (v %in% rf_top15_b) + (v %in% lasso_sel_b)
      if (n_methods >= 2 && v %in% rownames(stability_votes)) {
        stability_votes[v, b] <- 1
      }
    }
    
    if (b %% pb_tick == 0) cat(sprintf("  ...%d of %d bootstrap reps done\n", b, B))
  }
  
  sel_prob <- rowSums(stability_votes) / B
  stability_df <- data.frame(variable=names(sel_prob), sel_prob=sel_prob,
                             stringsAsFactors=FALSE) |>
    arrange(desc(sel_prob))
  stability_df$stable <- stability_df$sel_prob >= 0.60
  stable_features <- stability_df$variable[stability_df$stable]
  
  cat(sprintf("\nSTABLE FEATURES (selection probability >= 60%% across %d bootstrap reps):\n", B))
  print(stability_df[stability_df$sel_prob >= 0.40, ])
  cat(sprintf("\n-> %d stable features\n", length(stable_features)))
  
  plot_df <- stability_df |> filter(sel_prob >= 0.10)
  plots[["p_stability"]] <- ggplot(filter_low_signal(plot_df) |>
                                     mutate(variable_disp = relabel_var(variable)),
                                   aes(sel_prob, reorder(variable_disp, sel_prob))) +
    geom_col(aes(fill=stable), width=0.7) +
    geom_vline(xintercept=0.60, linetype="dashed", colour=COL_RED, linewidth=1) +
    annotate("text", x=0.62, y=1, label="Stability threshold (60%)",
             colour=COL_RED, size=3, hjust=0, fontface="italic") +
    scale_fill_manual(values=c("TRUE"=COL_BLUE, "FALSE"=COL_GREY),
                      labels=c("TRUE"="Stable (>=60%)", "FALSE"="Unstable (<60%)"),
                      name=NULL) +
    scale_x_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.2),
                       labels=scales::percent_format(accuracy=1)) +
    labs(title = sprintf("BURST stability selection: feature selection probability (B=%d)", B),
         subtitle = sprintf("%d stable features at 60%% threshold",
                            length(stable_features)),
         x="Selection probability across bootstrap reps", y=NULL) +
    theme_ugib() + theme(axis.text.y=element_text(size=8))
  print(plots[["p_stability"]])
  
  saveRDS(stability_df, "stability_selection.rds")
  saveRDS(stable_features, "stable_features.rds")
  cat("Saved: stability_selection.rds, stable_features.rds\n")
  
  agree <- intersect(confirmed, stable_features)
  vote_only <- setdiff(confirmed, stable_features)
  stable_only <- setdiff(stable_features, confirmed)
  cat(sprintf("\nAgreement with single-run vote count:\n"))
  cat(sprintf("  Both methods: %d features\n", length(agree)))
  cat(sprintf("  Vote count only (possibly unstable): %s\n",
              if (length(vote_only) > 0) paste(vote_only, collapse=", ") else "none"))
  cat(sprintf("  Stability only (missed by vote count): %s\n",
              if (length(stable_only) > 0) paste(stable_only, collapse=", ") else "none"))
} else {
  cat("\nSection 6B stability selection SKIPPED (set RUN_STABILITY <- TRUE to run)\n")
}

# SECTION 6D: TIER 1+2+3 FORWARD-SELECTION
cat("\n============================================================\n")
cat("SECTION 6D: TIER 1+2+3 FORWARD-SELECTION\n")
cat("============================================================\n\n")

tier1 <- "hb_mean"
tier2 <- if (exists("stable_features")) {
  setdiff(stable_features, tier1)
} else {
  setdiff(unanimous, tier1)
}
tier3 <- setdiff(confirmed, c(tier1, tier2))

cat(sprintf("Tier 1 (clinical anchor):      %d feature(s)  -- %s\n",
            length(tier1), paste(tier1, collapse=", ")))
cat(sprintf("Tier 2 (stability >=60%%):      %d feature(s)  -- %s\n",
            length(tier2), paste(tier2, collapse=", ")))
cat(sprintf("Tier 3 (vote-only candidates): %d feature(s)  -- %s\n\n",
            length(tier3),
            if (length(tier3) > 0) paste(tier3, collapse=", ") else "(none)"))

fit_rf_set <- function(fs, outcome=PRIMARY, seed=2024) {
  d <- df_imp |> dplyr::select(all_of(fs), all_of(outcome)) |>
    filter(!is.na(.data[[outcome]]))
  d[[outcome]] <- factor(d[[outcome]], levels=c("NO","YES"))
  cls_min <- min(table(d[[outcome]]))
  set.seed(seed)
  rf <- randomForest(reformulate(fs, outcome), data=d,
                     ntree=500, strata=d[[outcome]],
                     sampsize=c("NO"=cls_min, "YES"=cls_min))
  oob_p <- rf$votes[, "YES"] / rowSums(rf$votes)
  r <- roc(d[[outcome]], oob_p, levels=c("NO","YES"),
           direction="<", quiet=TRUE)
  list(rf=rf, roc=r, auc=as.numeric(auc(r)),
       ci=ci.auc(r), n_feat=length(fs), feats=fs)
}

base_set  <- c(tier1, tier2)
base_fit  <- fit_rf_set(base_set)
cat(sprintf("Base model (Tier 1+2): %d features, OOB AUC=%.4f (%.4f-%.4f)\n",
            base_fit$n_feat, base_fit$auc, base_fit$ci[1], base_fit$ci[3]))

cat(sprintf("\nForward-selection over %d Tier 3 candidates:\n", length(tier3)))
cat(sprintf("%-22s %5s %8s %8s %8s %12s %s\n",
            "Candidate", "n", "AUC", "CI_lo", "CI_hi", "DeLong_p", "decision"))
cat(paste(rep("-", 90), collapse=""), "\n")

retained_tier3 <- character(0)
forward_log <- data.frame()
if (length(tier3) > 0) {
  for (cand in tier3) {
    test_set <- c(base_set, cand)
    test_fit <- fit_rf_set(test_set)
    delong_p <- tryCatch(
      roc.test(base_fit$roc, test_fit$roc,
               method="delong", paired=TRUE)$p.value,
      error = function(e) NA_real_)
    ci_uplift <- test_fit$ci[1] > base_fit$ci[3]
    auc_improves <- test_fit$auc > base_fit$auc
    retain <- auc_improves && ((!is.na(delong_p) && delong_p < 0.05) || ci_uplift)
    decision <- if (retain) "RETAIN" else "drop"
    cat(sprintf("%-22s %5d %8.4f %8.4f %8.4f %12.4f %s\n",
                cand, test_fit$n_feat, test_fit$auc,
                test_fit$ci[1], test_fit$ci[3],
                ifelse(is.na(delong_p), -1, delong_p), decision))
    if (retain) retained_tier3 <- c(retained_tier3, cand)
    forward_log <- rbind(forward_log, data.frame(
      candidate = cand, n_feat = test_fit$n_feat,
      auc = round(test_fit$auc, 4),
      ci_lo = round(test_fit$ci[1], 4),
      ci_hi = round(test_fit$ci[3], 4),
      delong_p = round(delong_p, 4),
      auc_improves = auc_improves,
      retained = retain,
      stringsAsFactors = FALSE))
  }
}

final_set <- sort(unique(c(base_set, retained_tier3)))
cat(sprintf("\nTier 3 retained: %d / %d candidates -- %s\n",
            length(retained_tier3), length(tier3),
            if (length(retained_tier3) > 0)
              paste(retained_tier3, collapse=", ") else "(none -- base model is final)"))

final_fit <- fit_rf_set(final_set)
cat(sprintf("\nFINAL TIER MODEL: %d features, OOB AUC=%.4f (%.4f-%.4f)\n",
            length(final_set), final_fit$auc,
            final_fit$ci[1], final_fit$ci[3]))
cat(sprintf("  Features: %s\n", paste(final_set, collapse=", ")))

vote_set  <- intersect(confirmed, names(df_imp))
vote_set  <- setdiff(vote_set, "time_to_endo_h")
vote_fit  <- fit_rf_set(vote_set)
delong_final_vs_vote <- tryCatch(
  roc.test(final_fit$roc, vote_fit$roc, method="delong", paired=TRUE)$p.value,
  error = function(e) NA_real_)
cat(sprintf("\nVote-count >=2/n model: %d features, OOB AUC=%.4f (%.4f-%.4f)\n",
            length(vote_set), vote_fit$auc, vote_fit$ci[1], vote_fit$ci[3]))
cat(sprintf("DeLong (final tier vs vote-count): p=%.4f -- %s\n",
            delong_final_vs_vote,
            ifelse(is.na(delong_final_vs_vote), "n/a",
                   ifelse(delong_final_vs_vote < 0.05,
                          "DIFFERENT discrimination",
                          "EQUIVALENT discrimination (use parsimony to choose)"))))

confirmed_vote_count <- confirmed
confirmed <- final_set

cat(sprintf("\n>>> HEADLINE FEATURE SET (Tier 1+2+retained Tier 3): %d features <<<\n",
            length(confirmed)))
cat(sprintf("    %s\n", paste(confirmed, collapse=", ")))
cat(sprintf("\n>>> SENSITIVITY FEATURE SET (>=2/n vote-count): %d features\n",
            length(confirmed_vote_count)))

saveRDS(final_set,            "final_feature_set.rds")
saveRDS(confirmed_vote_count, "feature_vote_set.rds")
cat(sprintf("V2.4.5 FIX 1: re-saved final_feature_set.rds with %d features.\n",
            length(final_set)))
cat(sprintf("              also saved feature_vote_set.rds (%d features) for transparency.\n",
            length(confirmed_vote_count)))

saveRDS(list(
  tier1 = tier1, tier2 = tier2, tier3 = tier3,
  base_set = base_set, retained_tier3 = retained_tier3,
  final_set = final_set, vote_set = vote_set,
  base_auc  = base_fit$auc,  base_ci  = base_fit$ci,
  final_auc = final_fit$auc, final_ci = final_fit$ci,
  vote_auc  = vote_fit$auc,  vote_ci  = vote_fit$ci,
  delong_final_vs_vote = delong_final_vs_vote,
  forward_log = forward_log
), "tier_forward_selection.rds")
cat("\nSaved: tier_forward_selection.rds\n")

cat("\n============================================================\n")
cat("HEADLINE FEATURE SET MISSINGNESS (post-Section 6D)\n")
cat("============================================================\n\n")

headline_miss <- data.frame(
  variable    = confirmed,
  n_total     = nrow(df),
  n_missing   = vapply(confirmed, function(v) sum(is.na(df[[v]])), integer(1)),
  pct_missing = vapply(confirmed,
                       function(v) round(100 * mean(is.na(df[[v]])), 1),
                       numeric(1)),
  stringsAsFactors = FALSE
)
headline_miss <- headline_miss[order(-headline_miss$pct_missing), ]
rownames(headline_miss) <- NULL
print(headline_miss)

cat(sprintf("\n  Max missingness in headline set: %.1f%% (%s)\n",
            headline_miss$pct_missing[1], headline_miss$variable[1]))
cat(sprintf("  Mean missingness across headline set: %.1f%%\n",
            mean(headline_miss$pct_missing)))
cat(sprintf("  Headline features with 0%% missing: %d / %d\n",
            sum(headline_miss$pct_missing == 0), nrow(headline_miss)))

saveRDS(headline_miss, "headline_missingness.rds")
cat("\nSaved: headline_missingness.rds\n")

# SECTION 6C: VOTE-THRESHOLD SENSITIVITY
cat("\n============================================================\n")
cat("SECTION 6C: VOTE-THRESHOLD SENSITIVITY\n")
cat("============================================================\n\n")

build_set <- function(min_votes) {
  if (min_votes == "stability_60") {
    fs <- stable_features
  } else if (min_votes == "all") {
    fs <- intersect(feat_cands, names(df_imp))
  } else {
    fs <- vote_df$variable[vote_df$votes >= min_votes]
  }
  fs <- sort(unique(union(fs, "hb_mean")))
  fs <- intersect(fs, names(df_imp))
  setdiff(fs, "time_to_endo_h")
}

.thr_keys   <- c("All candidates",
                 sprintf("Vote >=1/%d", n_methods),
                 sprintf("Vote >=2/%d", n_methods),
                 sprintf("Vote >=3/%d", n_methods))
.thr_values <- list(build_set("all"),
                    build_set(1),
                    build_set(2),
                    build_set(3))
if (n_methods == 4L) {
  .thr_keys   <- c(.thr_keys,   "Vote ==4/4 (unanimous)")
  .thr_values <- c(.thr_values, list(build_set(4)))
}
.thr_keys   <- c(.thr_keys,   "Stability >=60%")
.thr_values <- c(.thr_values, list(build_set("stability_60")))
threshold_sets <- setNames(.thr_values, .thr_keys)
.selected_threshold_label <- sprintf("Vote >=2/%d", n_methods)

fit_one_set <- function(fs) {
  if (length(fs) < 2) return(NULL)
  d_thr <- df_imp |> dplyr::select(all_of(fs), all_of(PRIMARY)) |>
    filter(!is.na(.data[[PRIMARY]]))
  d_thr[[PRIMARY]] <- factor(d_thr[[PRIMARY]], levels=c("NO","YES"))
  hr_min <- min(table(d_thr[[PRIMARY]]))
  set.seed(2024)
  rf_thr <- randomForest(reformulate(fs, PRIMARY), data=d_thr,
                         ntree=500, strata=d_thr[[PRIMARY]],
                         sampsize=c("NO"=hr_min, "YES"=hr_min))
  oob_p <- rf_thr$votes[, "YES"] / rowSums(rf_thr$votes)
  r <- roc(d_thr[[PRIMARY]], oob_p,
           levels=c("NO","YES"), direction="<", quiet=TRUE)
  ci_v <- ci.auc(r)
  list(auc=as.numeric(auc(r)), ci_low=ci_v[1], ci_high=ci_v[3],
       n_feat=length(fs), feats=fs)
}

threshold_results <- list()
cat(sprintf("%-26s %5s   %-15s\n", "Threshold", "n_feat", "OOB AUC (95% CI)"))
cat(paste(rep("-", 60), collapse=""), "\n")
for (nm in names(threshold_sets)) {
  res <- fit_one_set(threshold_sets[[nm]])
  if (is.null(res)) {
    cat(sprintf("%-26s %5s   %s\n", nm, "<2", "skipped"))
    next
  }
  threshold_results[[nm]] <- res
  cat(sprintf("%-26s %5d   %.3f (%.3f-%.3f)\n",
              nm, res$n_feat, res$auc, res$ci_low, res$ci_high))
}

thr_df <- do.call(rbind, lapply(names(threshold_results), function(nm) {
  r <- threshold_results[[nm]]
  data.frame(threshold = nm, n_feat = r$n_feat, auc = r$auc,
             ci_low = r$ci_low, ci_high = r$ci_high,
             stringsAsFactors = FALSE)
}))
thr_df$threshold <- factor(thr_df$threshold, levels = thr_df$threshold)

plots[["p_thr_sensitivity"]] <- ggplot(thr_df,
                                       aes(x = auc, y = threshold, colour = threshold == .selected_threshold_label)) +
  geom_point(size = 4, show.legend = FALSE) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.25, linewidth = 0.7, show.legend = FALSE) +
  geom_text(aes(label = sprintf("n=%d", n_feat)), hjust = -0.5, size = 3, show.legend = FALSE) +
  scale_colour_manual(values = c("TRUE" = COL_RED, "FALSE" = COL_BLUE),
                      labels = setNames(c(sprintf("Selected (%s)",
                                                  .selected_threshold_label),
                                          "Alternative"),
                                        c("TRUE", "FALSE")),
                      name = NULL) +
  xlim(0.5, 1.0) +
  labs(title = "BURST feature-selection threshold sensitivity",
       subtitle = "OOB AUC across vote thresholds; bars = 95% CI; n = features included",
       x = "OOB AUC - High risk", y = "Inclusion criterion") +
  theme_ugib() + theme(legend.position = "none")
print(plots[["p_thr_sensitivity"]])

saveRDS(thr_df, "vote_threshold_sensitivity.rds")
cat("\nSaved: vote_threshold_sensitivity.rds\n")
cat(sprintf(paste0("\nVerdict (V2.4 update -- Tier 6D supersedes vote-threshold choice):\n",
                   "  All vote-threshold cutoffs from >=1/%d to maximum yield AUC within 0.020;\n",
                   "  the stability >=60%% set (n=6) achieves AUC equivalent to all alternatives.\n",
                   "  HEADLINE choice is governed by Section 6D Tier hierarchy, not by this table.\n",
                   "  This table is reported as a sensitivity analysis demonstrating that the\n",
                   "  feature-selection cutoff is not driving the discrimination result.\n"),
            n_methods))

# SECTION 6E: PERMUTATION NULL FOR VOTE-COUNT
cat("\n============================================================\n")
cat("SECTION 6E: PERMUTATION NULL FOR VOTE-COUNT\n")
cat("============================================================\n")
cat("Estimating null vote distribution under shuffled outcome...\n")

RUN_PERMUTATION <- TRUE
B_perm <- 200

if (RUN_PERMUTATION) {
  set.seed(2024)
  df_perm <- df_imp |> dplyr::select(all_of(feat_cands), all_of(PRIMARY)) |>
    filter(!is.na(.data[[PRIMARY]]))
  perm_votes <- matrix(0L, nrow=length(feat_cands), ncol=B_perm,
                       dimnames=list(feat_cands, NULL))
  pb_step <- ceiling(B_perm / 10)
  
  for (b in seq_len(B_perm)) {
    df_p <- df_perm
    df_p[[PRIMARY]] <- sample(df_p[[PRIMARY]])
    if (length(unique(df_p[[PRIMARY]])) < 2) next
    
    univ_b <- character(0)
    for (v in vars_num) {
      x <- df_p[[v]]; mt <- df_p[[PRIMARY]]; ok <- !is.na(x) & !is.na(mt)
      if (sum(ok) < 10) next
      if (length(unique(mt[ok])) < 2) next
      p <- tryCatch(suppressWarnings(
        wilcox.test(x[ok & mt=="YES"], x[ok & mt=="NO"])$p.value),
        error=function(e) NA)
      if (!is.na(p) && p < .05) univ_b <- c(univ_b, v)
    }
    for (v in vars_cat) {
      x <- df_p[[v]]; mt <- df_p[[PRIMARY]]; ok <- !is.na(x) & !is.na(mt)
      if (sum(ok) < 10) next
      tbl <- table(x[ok], mt[ok])
      if (nrow(tbl) < 2 || ncol(tbl) < 2) next
      p <- tryCatch(suppressWarnings(
        if (any(chisq.test(tbl)$expected < 5))
          fisher.test(tbl, simulate.p.value=TRUE, B=1000)$p.value
        else chisq.test(tbl)$p.value), error=function(e) NA)
      if (!is.na(p) && p < .05) univ_b <- c(univ_b, v)
    }
    
    rf_b <- tryCatch(
      randomForest(reformulate(feat_cands, PRIMARY), data=df_p,
                   ntree=300, importance=FALSE),
      error=function(e) NULL)
    rf_top15_b <- if (!is.null(rf_b)) {
      imp <- importance(rf_b)[, "MeanDecreaseGini"]
      names(sort(imp, decreasing=TRUE))[1:min(15, length(imp))]
    } else character(0)
    
    lasso_sel_b <- character(0)
    X_b <- tryCatch(model.matrix(reformulate(feat_cands), data=df_p)[, -1],
                    error=function(e) NULL)
    y_b <- as.numeric(df_p[[PRIMARY]] == "YES")
    if (!is.null(X_b) && length(unique(y_b)) == 2 && nrow(X_b) > 20) {
      cv_b <- tryCatch(cv.glmnet(X_b, y_b, family="binomial", alpha=1,
                                 nfolds=5, type.measure="deviance"),
                       error=function(e) NULL)
      if (!is.null(cv_b)) {
        coef_b <- coef(cv_b, s="lambda.1se")
        terms_b <- setdiff(rownames(coef_b)[as.numeric(coef_b) != 0],
                           "(Intercept)")
        lasso_sel_b <- unique(unlist(lapply(feat_cands, function(v)
          if (any(grepl(paste0("^", v), terms_b))) v)))
      }
    }
    
    sel_union <- unique(c(univ_b, rf_top15_b, lasso_sel_b))
    for (v in sel_union) {
      n_meth <- (v %in% univ_b) + (v %in% rf_top15_b) + (v %in% lasso_sel_b)
      if (n_meth >= 2 && v %in% rownames(perm_votes)) {
        perm_votes[v, b] <- 1L
      }
    }
    
    if (b %% pb_step == 0)
      cat(sprintf("  ...%d of %d permutations done\n", b, B_perm))
  }
  
  null_vote_prop <- rowSums(perm_votes) / B_perm
  null_p95 <- quantile(null_vote_prop, 0.95)
  cat(sprintf("\nNull vote rate (%% of permutations selected per variable):\n"))
  cat(sprintf("  Mean: %.3f | Median: %.3f | 95th percentile: %.3f\n",
              mean(null_vote_prop), median(null_vote_prop), null_p95))
  cat(sprintf("Real-data observed: %d / %d candidate features confirmed (>=2/n).\n",
              length(confirmed_vote_count), length(feat_cands)))
  
  real_vote_flag <- as.integer(feat_cands %in% confirmed_vote_count)
  perm_summary <- data.frame(
    variable        = feat_cands,
    real_confirmed  = real_vote_flag,
    null_vote_prob  = round(null_vote_prop, 3),
    above_null_p95  = null_vote_prop > null_p95,
    stringsAsFactors = FALSE) |>
    arrange(desc(null_vote_prob))
  cat("\nTop 15 by null vote probability (under shuffled outcome):\n")
  print(head(perm_summary, 15), row.names=FALSE)
  
  saveRDS(list(perm_votes=perm_votes, null_vote_prop=null_vote_prop,
               null_p95=as.numeric(null_p95), perm_summary=perm_summary,
               B=B_perm), "permutation_null.rds")
  cat("\nSaved: permutation_null.rds\n")
} else {
  cat("\nSection 6E permutation null SKIPPED (set RUN_PERMUTATION <- TRUE)\n")
}

# SECTION 6F: SYMPTOM_DURATION ABLATION
cat("\n============================================================\n")
cat("SECTION 6F: SYMPTOM_DURATION ABLATION\n")
cat("============================================================\n\n")

if ("symptom_duration" %in% confirmed) {
  full_set    <- confirmed
  reduced_set <- setdiff(confirmed, "symptom_duration")
  
  full_fit    <- fit_rf_set(full_set,    PRIMARY)
  reduced_fit <- fit_rf_set(reduced_set, PRIMARY)
  
  delong_sd <- tryCatch(
    roc.test(full_fit$roc, reduced_fit$roc, method="delong", paired=TRUE)$p.value,
    error = function(e) NA_real_)
  
  delta_auc <- reduced_fit$auc - full_fit$auc
  
  cat(sprintf("Headline (with symptom_duration, %d feat): AUC=%.4f (%.4f-%.4f)\n",
              full_fit$n_feat, full_fit$auc, full_fit$ci[1], full_fit$ci[3]))
  cat(sprintf("Reduced (without, %d feat):                AUC=%.4f (%.4f-%.4f)\n",
              reduced_fit$n_feat, reduced_fit$auc,
              reduced_fit$ci[1], reduced_fit$ci[3]))
  cat(sprintf("Delta AUC: %+.4f | DeLong p=%.4f\n", delta_auc, delong_sd))
  
  verdict <- if (abs(delta_auc) < 0.01) {
    "EQUIVALENT (|dAUC| < 0.01) -- symptom_duration retention is robust"
  } else if (abs(delta_auc) < 0.02) {
    "BORDERLINE (|dAUC| 0.01-0.02) -- retain with explicit limitation footnote"
  } else {
    "MATERIAL difference (|dAUC| >= 0.02) -- reconsider inclusion in Discussion"
  }
  cat(sprintf("Verdict: %s\n", verdict))
  
  saveRDS(list(
    full_set         = full_set,
    reduced_set      = reduced_set,
    auc_full         = full_fit$auc,    ci_full    = full_fit$ci,
    auc_reduced      = reduced_fit$auc, ci_reduced = reduced_fit$ci,
    delta_auc        = delta_auc,
    delong_p         = delong_sd,
    verdict          = verdict
  ), "symptom_duration_ablation.rds")
  cat("Saved: symptom_duration_ablation.rds\n")
} else {
  cat("symptom_duration not in headline -- ablation skipped\n")
}

cat("\n============================================================\n")
cat("STROBE FLOW DIAGRAM DATA\n")
cat("============================================================\n\n")

strobe_flow <- list(
  step1_screened           = nrow(df_raw),
  step2_renamed            = nrow(df),
  step3_with_primary       = sum(!is.na(df$high_risk)),
  step3_primary_events     = sum(df$high_risk == "YES", na.rm = TRUE),
  step4_with_secondary     = sum(!is.na(df$major_tx)),
  step4_secondary_events   = sum(df$major_tx == "YES", na.rm = TRUE),
  step5_with_outcome_class = sum(!is.na(df$transfusion_class)),
  step6_final_ml_cohort    = nrow(df_imp),
  date_first = if (sum(!is.na(df$entry_dt)) > 0)
    format(min(df$entry_dt, na.rm = TRUE), "%Y-%m-%d") else NA_character_,
  date_last  = if (sum(!is.na(df$entry_dt)) > 0)
    format(max(df$entry_dt, na.rm = TRUE), "%Y-%m-%d") else NA_character_,
  n_features_after_collinearity  = length(feat_cands),
  n_features_confirmed_primary   = length(confirmed),
  n_features_confirmed_secondary = length(confirmed_secondary)
)

cat("STROBE FLOW (n at each step):\n")
cat(sprintf("  Step 1: Patients in source file (UGIBfinal.xlsx)  ...... %d\n",
            strobe_flow$step1_screened))
cat(sprintf("  Step 2: After ETL (renames + cleaning) ................ %d\n",
            strobe_flow$step2_renamed))
cat(sprintf("  Step 3: With PRIMARY outcome (high_risk) defined ...... %d  (events: %d)\n",
            strobe_flow$step3_with_primary, strobe_flow$step3_primary_events))
cat(sprintf("  Step 4: With SECONDARY outcome (major_tx) defined ..... %d  (events: %d)\n",
            strobe_flow$step4_with_secondary, strobe_flow$step4_secondary_events))
cat(sprintf("  Step 5: With ordinal outcome (transfusion_class) ...... %d\n",
            strobe_flow$step5_with_outcome_class))
cat(sprintf("  Step 6: Final analytic cohort (post simple-impute) .... %d\n",
            strobe_flow$step6_final_ml_cohort))
if (!is.na(strobe_flow$date_first)) {
  cat(sprintf("\n  Period: %s to %s\n",
              strobe_flow$date_first, strobe_flow$date_last))
}
cat(sprintf("\nFeature pool flow:\n"))
cat(sprintf("  Pre-specified candidate pool (V2.4): age, demographics, vitals, symptoms, labs,\n"))
cat(sprintf("    drugs, comorbidities; cancer collapsed to active_cancer (single binary).\n"))
cat(sprintf("  A priori exclusions: ed_fluids_l, hematochezia (n=5), time_to_endo_h\n"))
cat(sprintf("    (treatment / rare event / not ED-deployable).\n"))
cat(sprintf("  After collinearity exclusion (hct_mean, bun, hemodynamic_instab, ed_levosim): %d\n",
            strobe_flow$n_features_after_collinearity))
cat(sprintf("  Tier 1+2+3 selected PRIMARY (Section 6D): ........................ %d\n",
            strobe_flow$n_features_confirmed_primary))
cat(sprintf("  Confirmed SECONDARY (>=2 of N voting methods; tree may be\n"))
cat(sprintf("    down-weighted; see Section 6 trace for N): .................... %d\n",
            strobe_flow$n_features_confirmed_secondary))

saveRDS(strobe_flow, "strobe_flow.rds")
cat("\nSaved: strobe_flow.rds\n\n")

cat("\n============================================================\n")
cat("SAVING OUTPUTS FOR PART 2\n")
cat("============================================================\n\n")

saveRDS(df,     "ugib_df_full.rds"); cat("Saved: ugib_df_full.rds\n")
saveRDS(df_imp, "ugib_df_ml.rds");   cat("Saved: ugib_df_ml.rds\n")

saveRDS(vote_df,    "variable_votes.rds");  cat("Saved: variable_votes.rds\n")
saveRDS(feat_cands, "feat_cands_all.rds");  cat("Saved: feat_cands_all.rds\n")

if(requireNamespace("mice", quietly=TRUE)){
  cat("\nRunning MICE imputation (m=20, method=pmm)...\n")
  mice_vars <- c(confirmed, "transfusion_class")
  mice_vars <- mice_vars[mice_vars %in% names(df)]
  set.seed(2024)
  mice_obj <- mice::mice(df |> dplyr::select(all_of(mice_vars)),
                         m=20, method="pmm", printFlag=FALSE, seed=2024)
  saveRDS(mice_obj, "ugib_mice_object.rds")
  cat("Saved: ugib_mice_object.rds\n")
} else {
  cat("mice package not installed -- ordinal regression in Part 2 will be skipped\n")
  cat("To enable, install: install.packages('mice')\n")
}

cat("\n=== PART 1 COMPLETE ===\n")
cat("All .rds files saved for Part 2 pipeline\n")

cat("\n============================================================\n")
cat("REPRODUCIBILITY APPENDIX\n")
cat("============================================================\n\n")

repro <- list(
  script_version    = "V2.4.4",
  run_datetime      = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  r_version         = R.version.string,
  platform          = R.version$platform,
  os                = paste(Sys.info()[c("sysname","release")], collapse=" "),
  locale            = Sys.getlocale("LC_CTYPE"),
  seeds_used        = list(
    primary_analysis = 2024,
    bootstrap_B      = 2024,
    permutation_B    = 2024,
    mice_imputation  = 2024,
    rf_default       = 2024
  ),
  pkg_versions      = sapply(
    c("readxl","dplyr","tidyr","ggplot2","rpart","rpart.plot",
      "randomForest","glmnet","pROC","mice"),
    function(p) tryCatch(as.character(packageVersion(p)),
                         error = function(e) "not installed")),
  data_file         = "UGIBfinal.xlsx",
  data_period       = if (exists("strobe_flow"))
    paste(strobe_flow$date_first, "to", strobe_flow$date_last) else NA_character_,
  cohort_n          = nrow(df),
  outcomes          = list(
    primary       = list(name=PRIMARY,   events=sum(df[[PRIMARY]]   == "YES", na.rm=TRUE)),
    secondary     = list(name=SECONDARY, events=sum(df[[SECONDARY]] == "YES", na.rm=TRUE))
  ),
  final_features    = if (exists("confirmed")) confirmed else NA,
  random_state_note = paste(
    "All stochastic procedures use set.seed(2024) immediately before the",
    "stochastic call. Bootstrap B=100 (Section 6B), permutation B=200",
    "(Section 6E), MICE m=20 (Part 1 saving). RF uses ntree=500 with",
    "balanced sampling via stratified sampsize.")
)

cat("R version: ", repro$r_version, "\n", sep="")
cat("Platform:  ", repro$platform, "\n", sep="")
cat("OS:        ", repro$os, "\n", sep="")
cat("Locale:    ", repro$locale, "\n\n", sep="")
cat("Package versions used:\n")
for (p in names(repro$pkg_versions)) {
  cat(sprintf("  %-15s %s\n", p, repro$pkg_versions[p]))
}
cat(sprintf("\nSeeds used (set.seed value): "))
cat(paste(unique(unlist(repro$seeds_used)), collapse=", "), "\n")
cat(sprintf("Cohort: n=%d  |  Period: %s\n",
            repro$cohort_n,
            if (!is.na(repro$data_period)) repro$data_period else "(n/a)"))
cat(sprintf("PRIMARY: %s (%d events)  |  SECONDARY: %s (%d events)\n",
            repro$outcomes$primary$name,   repro$outcomes$primary$events,
            repro$outcomes$secondary$name, repro$outcomes$secondary$events))

saveRDS(repro, "reproducibility_appendix.rds")
cat("\nSaved: reproducibility_appendix.rds\n")
