#----------------------------------------------------------#
# Predicting outcomes after AKI using multistate models
# Code for cohort derivation
# Roemer J. Janse - Last updated on 2026-02-24
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
path <- "L:/lab_research/RES-Folder-UPOD/NOSTRADAMUS_SALTRO/E_ResearchData/2_ResearchData/CLEANED_for_Multistate_outcomes_project/20012026/"

# Load functions
walk(list.files(here("funs")), \(x) source(paste0(here("funs"), "/", x)))

# Inclusion criteria
# - 18 years and older
# - episode of AKI in overlap with a hospitalisation
# - hospitalisation discharge from 2012-01-01 onwards
# - kidney function >= 60 mL/min/1.73m2
# - no maintenance dialysis
# - no prior kidney transplantation

# 1. Load and clean data ----
# Load hospitalisation data
load(paste0(path, "cleaned_DV_lab_SCr_AKI.Rdata"))

# Clean hospitalisation data
dat_hosp <- cleaned_dataset %>%
  # Set all column names to lower case
  set_colnames(tolower(colnames(.))) %>%
  # Arrange date based on hour
  arrange(date) %>%
  # Set POSIXct dates to Date format
  mutate(across(c("admission_date", "discharge_date"), as.Date)) %>% 
  # Rename variables
  rename(admission_dt = admission_date,
         discharge_dt = discharge_date,
         dialysis = any_dialysis,
         stage = aki_stage_incl365d) 

# Save clean hospitalisation data
save(dat_hosp,
     file = paste0(path, "dataframes/hosps.Rdata"))

# Load procedure data
load(paste0(path, "cleaned_DV_ICD10codes_procedures_freetext_comorbidities.Rdata"))

# Clean procedure and diagnosis data
dat_proc <- cleaned_dataset %>%
  # Set all column names to lower case
  set_colnames(tolower(colnames(.))) %>%
  # Set POSIXct date to Date object
  mutate(code_dt = as.Date(date),
         .keep = "unused") %>%
  # Set code to lowercase
  mutate(code = tolower(code))

# Save cleaned data
save(dat_proc,
     file = paste0(path, "dataframes/diagnoses_procedures.Rdata"))

# Load laboratory data
load(paste0(path, "cleaned_DV_lab_and_measurements.Rdata"))

# Clean laboratory data
dat_lab <- cleaned_dataset %>%
  # Set all column names to lower case
  set_colnames(tolower(colnames(.))) %>%
  # Set POSIXct date to Date object
  mutate(lab_dt = as.Date(date),
         .keep = "unused")

# Save cleaned data
save(dat_lab,
     file = paste0(path, "dataframes/laboratory.Rdata"))

# Load creatinine data
load(paste0(path, "cleaned_DV_lab_SCr_AKI.Rdata"))

# Clean creatinine data
dat_scr <- cleaned_dataset %>%
  # Set all column names to lower case
  set_colnames(tolower(colnames(.))) %>%
  # Keep only relevant columns
  select(id, date, creatinine_bl_umol_l, during_hospitalisation, specialty_hos) %>%
  # Set POSIXct date to Date object
  mutate(lab_dt = as.Date(date),
         .keep = "unused") %>%
  # Rename creatinine
  rename(creat = creatinine_bl_umol_l,
         ih = during_hospitalisation,
         dep = specialty_hos)

# Save cleaned data
save(dat_scr,
     file = paste0(path, "dataframes/creats.Rdata"))

# Load death data
load(paste0(path, "cleaned_DV_Death_or_last_contact_dates.Rdata"))

# Clean data
dat_death <- cleaned_dataset %>%
  # Set all column names to lower
  set_colnames(tolower(colnames(.))) %>%
  # Remove remaining grouping structure
  ungroup() %>%
  # Change dates from POSIXct to Date format
  mutate(across(death_date:last_contact_date, as.Date)) %>%
  # Change names
  rename(death_dt = death_date,
         censor_dt = last_contact_date)

# Save data
save(dat_death,
     file = paste0(path, "death.Rdata"))

# Load data frames
load(paste0(path, "dataframes/hosps.Rdata"))
load(paste0(path, "dataframes/diagnoses_procedures.Rdata"))
load(paste0(path, "dataframes/laboratory.Rdata"))
load(paste0(path, "dataframes/creats.Rdata"))
load(paste0(path, "dataframes/death.Rdata"))

# 2. Cohort spine ----
# We create the cohort spine based on AKI hospitalisation
# Number of individuals
n_distinct(dat_hosp[["id"]]) # n = 567,526

# Clean data 
# Save twice for later inclusion criterion
dat_spine <- (dat_ongoing_aki <- dat_hosp %>%
  # Drop non-AKI's
  filter(!is.na(aki_episode_incl365d)) %>%
  # Arrange per person per AKI episode
  arrange(id, aki_episode_incl365d) %>%
  # Group per person per AKI epsiode
  group_by(id, aki_episode_incl365d) %>%
  # Information about AKI stage is stored in the first row, so we want to keep that row; 
  # the overlapping hospitalisation is taken from rows further down
  fill(admission_dt, discharge_dt, .direction = "down") %>%
  # Keep only first AKI episode for each individual
  filter(aki_episode_incl365d == 1)) %>%
  # Keep first row per individual
  slice(1L) %>%
  # Remove grouping structure
  ungroup() %>%
  # Keep only relevant variables
  select(id, sex, dob, admission_dt, discharge_dt, stage)
         
# Number of individuals
n_distinct(dat_spine[["id"]]) # n = 51,058

# Exclude individuals with AKI that was never seen in-hospital 
dat_spine %<>% filter(!is.na(admission_dt) & !is.infinite(admission_dt))

# Number of individuals
n_distinct(dat_spine[["id"]]) # n = 37,199

# Keep hospitalisations only if discharge was from 2012 onwards
dat_spine %<>%
  # Drop earlier hospitalisations
  filter(discharge_dt >= as.Date("2012-01-01")) %>%
  # Create additional variables
  mutate(# Recode sex to female
         female = if_else(sex == "female", 1, 0),
         # Length of stay
         los = as.numeric(discharge_dt - admission_dt),
         # Age at discharge
         age = round(as.numeric(discharge_dt - dob) / 365.25)) %>%
  # Drop sex and reorder columns
  select(id, dob, age, female, admission_dt, discharge_dt, stage)

# Number of individuals
n_distinct(dat_spine[["id"]]) # n = 12,628

# Determine pre-AKI eGFR
dat_pae <- dat_scr %>%
  # Drop in-hospital creats
  filter(ih == 0) %>%
  # Join admission date
  left_join(dat_spine %>%
              # Keep only relevant variables
              select(id, admission_dt, female, dob),
            # Join by ID
            "id") %>%
  # Keep only creatinines in relevant time frame
  filter(lab_dt >= admission_dt - 365 & lab_dt <= admission_dt - 7) %>%
  # Calculate eGFR
  mutate(# Age at creatinine measurement
         age = round(as.numeric(lab_dt - dob) / 365.25),
         # eGFR
         pre_aki_egfr = ckd_epi(creat, female, age)) %>%
  # Arrange for grouping
  arrange(id) %>%
  # Group per individual
  group_by(id) %>%
  # Get median eGFR per individual
  summarise(pre_aki_egfr = median(pre_aki_egfr))

# Join to spine data
dat_spine %<>% left_join(dat_pae, "id")

# Determine max stage (separately as sometimes, two AKIs occur during hospitalisation and we want
# to take the max stage of any AKI during hospitalisation
# Data frame with all AKIs per individual
dat_max_stage <- dat_spine %>%
  # Keep only relevant information
  select(id, discharge_dt) %>%
  # Join all AKIs
  left_join(dat_hosp %>%
              # Keep only AKIs
              filter(aki_episode_incl365d >= 1) %>%
              # Keep only relevant information
              select(id, date, aki_episode_incl365d, stage),
            # Join by ID
            "id") %>%
  # Keep only AKIs prior to discharge
  filter(as.Date(date) <= discharge_dt) %>%
  # Arrange for grouping
  arrange(id) %>%
  # Group per individual
  group_by(id) %>%
  # Get max stage
  summarise(stage = max(stage))

# Add max stage to spine
dat_spine %<>%
  # Drop previous stage indicator
  select(-stage) %>%
  # Add stage indicator
  left_join(dat_max_stage, "id")

# Individuals with previous AKI based on ICD-10 codes
vec_paki <- dat_spine %>%
  # Join diagnostic codes
  left_join(dat_proc, "id") %>%
  # Keep only hypertension codes prior to admission
  filter(code %in% c("aki_icd10_dbc_date",
                     "aki_icd10_diagnosis_date",
                     "aki_diagn_descr_date") &
         code_dt <= admission_dt) %>%
  # Keep only one ID per individual
  distinct(id) %>%
  # Reduce to vector
  extract2("id")

# Exclude individuals with a previous AKI based on ICD-10
dat_spine %<>% filter(!(id %in% vec_paki))

# Number of individuals
n_distinct(dat_spine[["id"]]) # n = 12,534

# 3. Apply other inclusion criteria ----
# Age of 18 years and older
dat_spine %<>% filter(age >= 18)

# Number of individuals
n_distinct(dat_spine[["id"]]) # n = 10,165

## No previous maintenance dialysis
# Select individuals with maintenance dialysis
vec_main_dial <- dat_proc %>%
  # Keep only maintenance dialysis 
  filter(code %in% c("chronic_dialysis_procedure", 
                     "chronic_dialysis_icd10_dbc_date",
                     "chronic_dialysis_icd10_diagnosis_date",
                     "dialysis_shunt_procedure",
                     "chronic_dialysis_diagn_descr_date")) %>%
  # Add AKI admission date
  left_join(dat_spine %>%
              # Keep only ID and admission date
              select(id, discharge_dt),
            "id") %>%
  # Keep only observations with dialysis date prior to admission date
  filter(code_dt <= discharge_dt) %>%
  # Keep only unique individuals
  distinct(id) %>%
  # Keep only vector of IDs
  extract2("id")

# Remove individuals with prior maintenance dialysis
dat_spine %<>% filter(!(id %in% vec_main_dial))

# Number of individuals
n_distinct(dat_spine[["id"]]) # n = 9,643

## No previous kidney transplantation
# Select individuals with kidney transplantation
vec_ktx <- dat_proc %>%
  # Keep only maintenance dialysis 
  filter(code %in% c("surgery_ktx_procedure", 
                     "ktx_icd10_dbc_date",
                     "ktx_icd10_diagnosis_date",
                     "ktx_diagn_descr_date",
                     "ktx_followup_procedure",
                     "ktx_workup_procedure")) %>%
  # Add AKI admission date
  left_join(dat_spine %>%
              # Keep only ID and admission date
              select(id, discharge_dt),
            "id") %>%
  # Keep only observations with dialysis date prior to admission date
  filter(code_dt <= discharge_dt) %>%
  # Keep only unique individuals
  distinct(id) %>%
  # Keep only vector of IDs
  extract2("id")

# Remove individuals with prior kidney transplantation 
dat_spine %<>% filter(!(id %in% vec_ktx))

# Number of individuals
n_distinct(dat_spine[["id"]]) # n = 9,251

## Ongoing AKI
# Determine ongoing AKIs
vec_ongoing_aki <- dat_ongoing_aki %>%
  # Arrange for grouping
  arrange(id) %>%
  # Group per individual
  group_by(id) %>%
  # Get date of last AKI observation
  mutate(last_aki_measurement = as.Date(last(date))) %>%
  # Join index date of each individual
  left_join(dat_spine %>% 
              # Only relevant variables
              select(id, discharge_dt),
            # Join by ID
            "id") %>%
  # Keep individuals that had their last AKI measurement no later than discharge from their index hospitalisation
  filter(last_aki_measurement > discharge_dt.y) %>%
  # Keep only unique individuals
  distinct(id) %>%
  # Keep only vector of IDs
  extract2("id")

# Remove individuals whose AKI did not resolve at hospital discharge
dat_spine %<>% filter(!(id %in% vec_ongoing_aki))

# Number of individuals
n_distinct(dat_spine[["id"]]) # n = 7,665

## Kidney function at end of hospitalisation of 60 mL/min/1.73m2 or higher
dat_discharge_egfr <- dat_scr %>%
  # Keep only values for individuals in spine data
  filter(id %in% dat_spine[["id"]]) %>%
  # Add admission date discharge date, age, and sex for each individual
  left_join(dat_spine %>%
              # Keep specific variables
              select(id, age, female, admission_dt, discharge_dt), 
            "id") %>%
  # Drop all measurements that or not within one week before discharge
  filter(lab_dt >= discharge_dt - 7 & lab_dt <= discharge_dt) %>%
  # Calculate new variables
  mutate(# eGFr based on CKD-EPI 2009
         egfr = ckd_epi(creat, female, age),
         # Absolute time difference from discharge
         days_since_egfr = as.numeric(abs(discharge_dt - lab_dt))) %>%
  # Arrange data to put most recent measurement first
  arrange(id, days_since_egfr) %>%
  # Group per individual to reduce data to a single row per individual
  group_by(id) %>%
  # Keep first row per individual
  slice(1L) %>%
  # Remove grouping structure
  ungroup() %>%
  # Drop unnecessary columns
  select(-c(age:discharge_dt, lab_dt))

# Join kidney function to spine data
dat_spine %<>%
  # Join data
  left_join(dat_discharge_egfr, "id") %>%
  # Drop individuals with missing eGFR 
  filter(!is.na(egfr))

# Number of individuals
n_distinct(dat_spine[["id"]]) # n = 7,138

# Join kidney function to spine data
dat_spine %<>% filter(egfr >= 60)

# Number of individuals
n_distinct(dat_spine[["id"]]) # n = 4,874

## Censor individuals who died at hospital discharge
dat_spine %<>% 
  # Join death data
  left_join(dat_death %>% 
              # Keep only relevant variables
              select(id, death_dt),
            # Join by ID
            "id") %>%
  # Drop individuals who died at discharge
  filter(is.na(death_dt) | death_dt > discharge_dt) %>%
  # Drop death date
  select(-death_dt)

# Number of individuals
n_distinct(dat_spine[["id"]]) # n = 4,356

## Censor individuals who were censored prior to hospital discharge
dat_spine %<>% 
  # Join death data
  left_join(dat_death %>% 
              # Keep only relevant variables
              select(id, censor_dt),
            # Join by ID
            "id") %>%
  # Drop individuals who died at discharge
  filter(is.na(censor_dt) | censor_dt > discharge_dt) %>%
  # Drop death date
  select(-censor_dt)

# Number of individuals
n_distinct(dat_spine[["id"]]) # n = 4,023

# Remove vectors
rm(vec_ktx, vec_main_dial, vec_ongoing_aki)

# 4. Finalise data ----
# Clean up data and calculate last variables
dat_spine <- dat_spine %>%
  # Calculate last variables
  mutate(# Relative eGFR change (negative sign equals increase)
         egfr_change = (pre_aki_egfr - egfr) / pre_aki_egfr,
         # Length of hospital stay
         los = as.numeric(discharge_dt - admission_dt),
         # Missing pre-AKI-eGFR
         missing_pre_aki_egfr = if_else(is.na(pre_aki_egfr), 1, 0))

# Save data
save(dat_spine,
     file = paste0(path, "dataframes/cohort_derivation.Rdata"))
