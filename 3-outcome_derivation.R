#----------------------------------------------------------#
# Predicting outcomes after AKI using multistate models
# Code for outcome derivation
# Roemer J. Janse - Last updated on 2026-01-22
#----------------------------------------------------------#

# 0. Set-up ----
# Load packages
pacman::p_load("dplyr",          # Data wrangling
               "magrittr",       # Efficient pipelines
               "tidyr",          # Data cleaning
               "stringr",        # Working with strings
               "purrr",          # Functional programming
               "here",           # Local paths
               "ggplot2",        # Visual representation of outcomes
               "conflicted"      # Resolve function conflicts
)

# Resolve function conflicts
conflicts_prefer(dplyr::filter) # Between stats and dplyr

# Set data path, based on OS
if(.Platform[["OS.type"]] == "unix"){
  path <- "/Users/rjanse5/Networkshares/lab/lab_research/RES-Folder-UPOD/NOSTRADAMUS_SALTRO/E_ResearchData/2_ResearchData/CLEANED_for_Multistate_outcomes_project/20012026/dataframes/"
} else path <- "L:/lab_research/RES-Folder-UPOD/NOSTRADAMUS_SALTRO/E_ResearchData/2_ResearchData/CLEANED_for_Multistate_outcomes_project/20012026/dataframes/"

# Load functions
walk(list.files(here("funs")), \(x) source(paste0(here("funs"), "/", x)))

# Load spine data
load(paste0(path, "covariate_derivation.Rdata"))

# 1. Narrow (kidney) model ----
## 1.1. Recurrent AKI ----
# Load AKIs
load(paste0(path, "hosps.Rdata"))

# Load diagnosis codes
load(paste0(path, "diagnoses_procedures.Rdata"))

# Add recurrent AKIs to data
dat_spine %<>%
  # Join AKIs
  left_join(dat_hosp %>%
              # Keep only AKIs
              filter(!is.na(aki_episode_incl365d)) %>%
              # Cannot be the first episode
              filter(aki_episode_incl365d > 1) %>%
              # Arrange for grouping
              arrange(id) %>%
              # Group per inidvidual per AKI
              group_by(id) %>%
              # Keep one row per AKI
              slice(1L) %>% 
              # Remove grouping structure
              ungroup() %>%
              # Keep only relevant variables
              select(id, date, stage) %>%
              # Set date to AKI date in Date format (instead of POSIXct)
              mutate(lab_aki_dt = as.Date(date)) %>%
              # Rename variables
              rename(lab_aki_stage = stage),
            # Join by ID
            "id") %>%
  # Keep only AKIs after index date
  mutate(lab_aki_dt = if_else(lab_aki_dt <= discharge_dt, NA, lab_aki_dt)) %>%
  # Set first recurrent AKI on top
  arrange(id, lab_aki_dt) %>%
  # Keep only a single row per individual
  distinct(id, .keep_all = TRUE) %>% 
  # Add diagnoses of AKI
  left_join(dat_proc %>%
              # Only AKI diagnosis
              filter(code %in% c("aki_icd10_dbc_date",
                                 "aki_icd10_diagnosis_date",
                                 "aki_diagn_descr_date")),
            # By ID
            "id") %>%
  # Keep only AKI codes after index date
  mutate(code_dt = if_else(code_dt <= discharge_dt + 30, NA, code_dt)) %>%
  # Set first recurrent AKI on top
  arrange(id, code_dt) %>%
  # Keep only a single row per individual
  distinct(id, .keep_all = TRUE) %>%
  # Final AKI variables
  mutate(# Date of AKI
         aki_dt = pmin(code_dt, lab_aki_dt, na.rm = TRUE),
         # AKI indicator
         aki = if_else(!is.na(aki_dt), 1, 0),
         # Diagnosis type for AKI
         aki_type = case_when(aki == 1 & aki_dt == lab_aki_dt ~ "lab",
                              aki == 1 & aki_dt == code_dt ~ "icd",
                              aki == 0 ~ NA),
         # AKI stage for lab
         aki_stage = if_else(aki_type == "lab", lab_aki_stage, NA),
         # Time to AKI
         aki_tte = as.numeric(aki_dt - discharge_dt)) %>%
  # Drop redundant variables
  select(-c(lab_aki_dt, lab_aki_stage, code_dt, code))
  
## 1.2. CKD ----
# Load lab
load(paste0(path, "creats.Rdata"))

# Keep only relevant creatinine measurements
dat_scr_cohort <- dat_scr %>%
  # Keep only individuals in dat_spine
  filter(id %in% dat_spine[["id"]]) %>%
  # Drop any creatinine measurements taken during hospitalisation in which an AKI was present
  left_join(dat_hosp %>%
              # Keep only AKIs
              filter(!is.na(aki_episode_incl365d)) %>%
              # Keep only relevant variables
              select(id, admission_dt, discharge_dt),
            # Join by ID
            "id",
            # Both rows contain multiples of each individual
            relationship = "many-to-many") %>%
  # Drop creatinine measurements during AKI hospitalisations
  filter(lab_dt <= admission_dt | lab_dt >= discharge_dt) %>%
  # Keep only unique creatinine measurements
  distinct(id, creat, lab_dt) %>%
  # Add age and sex
  left_join(dat_spine %>%
              # Relevant variables
              select(id, age, female, discharge_dt),
            # By ID
            "id") %>%
  # Convert creatinine to eGFR
  mutate(egfr = ckd_epi(creat, female, age)) %>%
  # Keep only creatinines after index date
  filter(lab_dt >= discharge_dt)

# Fit linear model for each individual to get time to CKD
ttc_li <- sapply(group_split(dat_scr_cohort, id), \(x) linear_interpolation(x, threshold = 60))

# Get CKD based on two measurements with 90 days in-between
ttc_tm <- sapply(group_split(dat_scr_cohort, id), \(x) two_measurements(x, threshold = 60))

# Get IDs for all individuals and add ttc
dat_ttc <- distinct(dat_scr_cohort, id) %>%
  # Add TTC
  mutate(ttc_li = ttc_li,
         ttc_tm = ttc_tm,
         ckd_tte = pmin(ttc_li, ttc_tm, na.rm = TRUE),
         ckd_type = case_when(ckd_tte == ttc_li ~ "li", 
                              ckd_tte == ttc_tm ~ "tm",
                              .default = NA)) %>%
  # Drop irrelevant variables
  select(-c(ttc_li, ttc_tm))

# Add to data
dat_spine %<>%
  # Add ttc data
  left_join(dat_ttc, "id") %>%
  # Create date and indicator from time to CKD
  mutate(ckd_dt = as.Date(discharge_dt + ckd_tte),
         ckd = if_else(is.na(ckd_tte), 0, 1)) 

## 1.3. Kidney failure ----
# Fit linear model for each individual to get time to kidney failure
ttk_li <- sapply(group_split(dat_scr_cohort, id), \(x) linear_interpolation(x, threshold = 15))

# Get kidney failure based on two measurements with 90 days in-between
ttk_tm <- sapply(group_split(dat_scr_cohort, id), \(x) two_measurements(x, threshold = 15))

# Get IDs for all individuals and add ttc
dat_ttk <- distinct(dat_scr_cohort, id) %>%
  # Add TTC
  mutate(ttk_li = ttk_li,
         ttk_tm = ttk_tm,
         kf_tte = pmin(ttk_li, ttk_tm, na.rm = TRUE),
         kf_type = case_when(kf_tte == ttk_li ~ "li",
                             kf_tte == ttk_tm ~ "tm",
                             .default = NA))

# Add diagnoses for chronic dialysis or kidney transplantation
dat_ttk %<>% 
  # Add discharge date for each individual
  full_join(dat_spine %>%
              # Keep only relevant variables
              select(id, discharge_dt),
            # Join by ID
            "id") %>%
  # Join diagnoses and procedures data
  left_join(dat_proc %>%
              # Only KTx and chronic dialysis
              filter(code %in% c("surgery_ktx_procedure",
                                 "ktx_icd10_dbc_date",
                                 "ktx_icd10_diagnosis_date",
                                 "ktx_diagn_descr_date",
                                 "ktx_followup_procedure",
                                 "ktx_workup_procedure",
                                 "chronic_dialysis_procedure",
                                 "chronic_dialysis_icd10_dbc_date",
                                 "chronic_dialysis_icd10_diagnosis_date",
                                 "dialysis_shunt_procedure",
                                 "chronic_dialysis_diagn_descr_date")),
            # Join by ID
            "id") %>%
  # Keep only codes after discharge
  mutate(code_dt = if_else(code_dt <= discharge_dt, NA, code_dt)) %>%
  # Calculate date of kidney failure according to creatinine
  mutate(creat_kf_dt = discharge_dt + ceiling(kf_tte)) %>%
  # Arrange per individual with earliest outcomes first
  arrange(id, creat_kf_dt, code_dt) %>%
  # Group per individual
  group_by(id) %>%
  # Keep first row per individual
  slice(1L) %>%
  # Drop grouping structure
  ungroup() %>%
  # Calculate kidney failure variables
  mutate(# Kidney failure date
         kf_dt = pmin(code_dt, creat_kf_dt, na.rm = TRUE),
         # Kidney failure indicator
         kf = if_else(!is.na(kf_dt), 1, 0),
         # Kidney failure type
         kf_type = case_when(kf == 1 & kf_dt == creat_kf_dt ~ kf_type,
                             kf == 1 & kf_dt == code_dt ~ "icd",
                             kf == 0 ~ NA),
         # Time to kidney fialure
         kf_tte = as.numeric(kf_dt - discharge_dt)) %>%
  # Drop redundant variables
  select(-c(discharge_dt, code, code_dt, creat_kf_dt, ttk_li, ttk_tm))

# Add to dat_spine
dat_spine %<>% left_join(dat_ttk, "id")

## 1.4 Acute coronary heart disease ----
# Add diagnoses for acute coronary heart disease
dat_spine %<>%
  # Join diagnoses and procedures data
  left_join(dat_proc %>%
              # Only acute CHD events
              filter(code %in% c("chd_other_icd10_dbc_date",
                                 "chd_other_icd10_diagnosis_date",
                                 "chd_diagn_descr_date",
                                 "chd_icd10_dbc_date",
                                 "chd_icd10_diagnosis_date",
                                 "surgery_cardiothoracic_cabg_procedure",
                                 "cabg_diagn_descr_date",
                                 "angiogram_pci_procedure")),
            # Join by ID
            "id") %>%
  # Keep only codes after discharge and create indicator
  mutate(# Only dates after discharge
         chd_dt = if_else(code_dt <= discharge_dt, NA, code_dt),
         # Create indicator
         chd = if_else(is.na(chd_dt), 0, 1)) %>%
  # Arrange for grouping, putting first diagnosis on top
  arrange(id, code_dt) %>%
  # Group per individual
  group_by(id) %>%
  # Keep one row per individual
  slice(1L) %>%
  # Remove grouping structure
  ungroup() %>%
  # Drop redundant variables
  select(-c(code_dt, code))

## 1.5 Acute cerebrovascular disease ----
# Add diagnoses for acute cerebrovascular disease
dat_spine %<>%
  # Join diagnoses and procedures data
  left_join(dat_proc %>%
              # Only cerebrovascular disease
              filter(code %in% c("cvd_icd10_diagnosis_date", 
                                 "cvd_diagn_descr_date",
                                 "cvd_icd10_dbc_date",
                                 "comorb_cerebrovascular_disease_procedure")),
            # Join by ID
            "id") %>%
  # Keep only codes after discharge and create indicator
  mutate(# Only dates after discharge
    cvd_dt = if_else(code_dt <= discharge_dt, NA, code_dt),
    # Create indicator
    cvd = if_else(is.na(cvd_dt), 0, 1)) %>%
  # Arrange for grouping, putting first diagnosis on top
  arrange(id, code_dt) %>%
  # Group per individual
  group_by(id) %>%
  # Keep one row per individual
  slice(1L) %>%
  # Remove grouping structure
  ungroup() %>%
  # Drop redundant variables
  select(-c(code_dt, code))

## 1.6 Heart failure hospitalisation
# Load data
load(paste0(str_replace(path, "dataframes/", ""), "cleaned_DV_HF_hospitalization_dates.Rdata"))

# Extract heart failure hospitalisations
dat_spine %<>%
  # Join heart failure hospitalisation data
  left_join(cleaned_dataset %>%
              # Clean data
              set_colnames(tolower(colnames(.))) %>%
              # Join discharge date
              left_join(dat_spine %>%
                          # Keep only ID and discharge date
                          select(id, discharge_dt),
                        # Join by ID
                        "id") %>%
              # Set POSIXct to Date
              mutate(hfhosp_dt = as.Date(hf_hospitalization),
                     .keep = "unused") %>%
              # Keep only dates after discharge
              mutate(hfhosp_dt = if_else(hfhosp_dt <= discharge_dt, NA, hfhosp_dt)) %>%
              # Arrange per person on time
              arrange(id, hfhosp_dt) %>%
              # Group per individual
              group_by(id) %>%
              # Keep one row per individual
              slice(1L) %>%
              # Remove grouping structure
              ungroup() %>%
              # Drop discharge date
              select(-discharge_dt),
            # Join by ID
            "id")

## 1.7. Death ----
load(paste0(path, "death.Rdata"))

# Join to spine data
dat_spine %<>% left_join(dat_death, "id")

# 2. Finalise data ----
## 2.1. Final computations ----
# Finalise spine data
dat_spine %<>%
  # Adjust dates based on censoring, administrative censoring (2025-09-17) and 10-year final follow-up
  mutate(# Administrative censoring
         adm_cens_dt = as.Date("2025-09-17"),
         # 10-year max. follow-up
         end_fu_dt = discharge_dt + 10 * 365.25,
         # Final date for each individual
         final_dt = pmin(death_dt, adm_cens_dt, end_fu_dt, censor_dt, na.rm = TRUE),
         # If final date is discharge, add 1 day and censor there
         # We cannot exclude these individuals as lack of follow-up is future information;
         # So we just give them minimal follow-up
         final_dt = if_else(final_dt == discharge_dt, final_dt + 1, final_dt),
         # Total follow-up in years
         total_fu = as.numeric(final_dt - discharge_dt) / 365.25,
         # Death variables
         death = if_else(final_dt == death_dt, 1, 0, missing = 0),
         death_dt = if_else(death == 1, death_dt, final_dt),
         death_tte = as.numeric(death_dt - discharge_dt),
         # AKI variables
         aki = if_else(aki_dt == pmin(aki_dt, final_dt, na.rm = TRUE), 1, 0, missing = 0),
         aki_dt = if_else(aki == 1, aki_dt, final_dt),
         aki_tte = as.numeric(aki_dt - discharge_dt),
         aki_type = if_else(aki == 1, aki_type, NA),
         # CKD variables
         ckd = if_else(ckd_dt == pmin(ckd_dt, final_dt, na.rm = TRUE), 1, 0, missing = 0),
         ckd_dt = if_else(ckd == 1, ckd_dt, final_dt),
         ckd_tte = as.numeric(ckd_dt - discharge_dt),
         ckd_type = if_else(ckd == 1, ckd_type, NA),
         # Kidney failure variables
         kf = if_else(kf_dt == pmin(kf_dt, final_dt, na.rm = TRUE), 1, 0, missing = 0),
         kf_dt = if_else(kf == 1, kf_dt, final_dt),
         kf_tte = as.numeric(kf_dt - discharge_dt),
         kf_type = if_else(kf == 1, kf_type, NA),
         # Acute CHD variables
         chd = if_else(chd_dt == pmin(chd_dt, final_dt, na.rm = TRUE), 1, 0, missing = 0),
         chd_dt = if_else(chd == 1, chd_dt, final_dt),
         chd_tte = as.numeric(chd_dt - discharge_dt),
         # Acute CVD variables
         cvd = if_else(cvd_dt == pmin(cvd_dt, final_dt, na.rm = TRUE), 1, 0, missing = 0),
         cvd_dt = if_else(cvd == 1, cvd_dt, final_dt),
         cvd_tte = as.numeric(cvd_dt - discharge_dt),
         # Heart failure hospitalisation variables
         hfhosp = if_else(hfhosp_dt == pmin(hfhosp_dt, final_dt, na.rm = TRUE), 1, 0, missing = 0),
         hfhosp_dt = if_else(hfhosp == 1, hfhosp_dt, final_dt),
         hfhosp_tte = as.numeric(hfhosp_dt - discharge_dt),
         # MACCE
         macce = pmax(chd, cvd, hfhosp),
         macce_dt = if_else(macce == 1, pmin(chd_dt, cvd_dt, hfhosp_dt, na.rm = TRUE), final_dt),
            # If diagnosis is at the same moment, put MACCE 1 day later. This is the case for 0 CKD, 6 MACCE, and 15 AKI patients
         macce_dt = if_else(macce == 1 & 
                              (
                                (macce_dt == ckd_dt & ckd == 1) | 
                                (macce_dt == kf_dt & kf == 1) |
                                (macce_dt == aki_dt & aki == 1)
                                ), macce_dt + 1, macce_dt),
         macce_tte = as.numeric(macce_dt - discharge_dt),
         # MACCE with AKI
         macce_aki = if_else(macce == 1 & aki == 1, 1, 0),
         macce_aki_dt = if_else(macce_aki == 1, pmax(macce_dt, aki_dt), final_dt),
         macce_aki_tte = as.numeric(macce_aki_dt - discharge_dt),
         # MACCE with kidney failure
         macce_kf = if_else(macce == 1 & kf == 1, 1, 0),
         macce_kf_dt = if_else(macce_kf == 1, pmax(macce_dt, kf_dt), final_dt),
         macce_kf_tte = as.numeric(macce_kf_dt - discharge_dt),
         # MACCE with CKD
         macce_ckd = if_else(macce == 1 & ckd == 1, 1, 0),
         macce_ckd_dt = if_else(macce_ckd == 1, pmax(macce_dt, ckd_dt), final_dt),
         macce_ckd_tte = as.numeric(macce_ckd_dt - discharge_dt))

## 2.2. Save data ----
save(dat_spine,
     file = paste0(path, "outcome_derivation.Rdata"))
