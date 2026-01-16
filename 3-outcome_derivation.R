#----------------------------------------------------------#
# Predicting outcomes after AKI using multistate models
# Code for covariate derivation
# Roemer J. Janse - Last updated on 2025-12-22
#----------------------------------------------------------#

# 0. Set-up ----
# Load packages
pacman::p_load("dplyr",          # Data wrangling
               "magrittr",       # Efficient pipelines
               "tidyr",          # Data cleaning
               "stringr",        # Working with strings
               "purrr",          # Functional programming
               "here",           # Local paths
               "conflicted"      # Resolve function conflicts
)

# Resolve function conflicts
conflicts_prefer(dplyr::filter) # Between stats and dplyr

# Set data path
path <- "L:/lab_research/RES-Folder-UPOD/NOSTRADAMUS_SALTRO/E_ResearchData/2_ResearchData/CLEANED_for_Multistate_outcomes_project/17092025/dataframes/"

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
              filter(code == "ICD10_AKI"),
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
  mutate(egfr = ckd_epi(creat, female, age))

# Fit linear model for each individual to get time to CKD
ttc <- sapply(group_split(dat_scr_cohort, id), \(x) linear_interpolation(x, threshold = 60))

# Get IDs for all individuals and add ttc
dat_ttc <- distinct(dat_scr_cohort, id) %>%
  # Add TTC
  mutate(ckd_tte = ttc)

# Add to data
dat_spine %<>%
  # Add ttc data
  left_join(dat_ttc, "id") %>%
  # Create date and indicator from time to CKD
  mutate(ckd_dt = as.Date(discharge_dt + ckd_tte),
         ckd = if_else(is.na(ckd_tte), 0, 1)) 

## 1.3. Kidney failure ----
# Fit linear modeal for each individual to get time to kidney failure
ttk <- sapply(group_split(dat_scr_cohort, id), \(x) linear_interpolation(x, threshold = 15))

# Get IDs for all individuals and add ttc
dat_ttk <- distinct(dat_scr_cohort, id) %>%
  # Add TTC
  mutate(kf_tte = ttk)

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
              filter(code %in% c("Procedure_KTx_workup",
                                 "Procedure_surgery_KTx",
                                 "Procedure_KTx_followup",
                                 "Procedure_dialysis_shunt", 
                                 "Procedure_Chronic_dialysis",
                                 "ICD10_KTx", 
                                 "ICD10_Chronic_dialysis")),
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
         kf_type = case_when(kf == 1 & kf_dt == creat_kf_dt ~ "lab",
                             kf == 1 & kf_dt == code_dt ~ "icd",
                             kf == 0 ~ NA),
         # Time to kidney fialure
         kf_tte = as.numeric(kf_dt - discharge_dt)) %>%
  # Drop redundant variables
  select(-c(discharge_dt, code, code_dt, creat_kf_dt))

# Add to dat_spine
dat_spine %<>% left_join(dat_ttk, "id")

## 1.4. Save interim data ----
save(dat_spine, file = paste0(path, "narrow_outcome_derivation.Rdata"))

# 2. Broad (general) model ----
## 2.1. Kidney disease ----
# Determine kidney disease based on the narrow model
dat_spine %<>%
  # Create kidney disease variables
  mutate(kidney = pmax(aki, ckd, kf),
         kidney_dt = pmin(aki_dt,  ckd_dt, kf_dt, na.rm = TRUE),
         kidney_tte = pmin(aki_tte, ckd_tte, kf_tte, na.rm = TRUE),
         kidney_component = case_when(kidney_dt == aki_dt ~ "aki",
                                      kidney_dt == ckd_dt ~ "ckd",
                                      kidney_dt == kf_dt ~ "kf",
                                      .default = NA))

## 2.2. Oncological disease ----
# Determine oncological disease from diagnoses
test <- dat_spine %>%
  # Join diagnoses
  left_join(dat_proc %>%
