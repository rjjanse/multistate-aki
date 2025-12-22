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
load(paste0(path, "cohort_derivation.Rdata"))

# 1. Prepare data ----
# Load diagnoses data
load(paste0(path, "diagnoses_procedures.Rdata"))

# Load laboratory data
load(paste0(path, "laboratory.Rdata"))

# Load medication data
load(paste0(str_replace(path, "dataframes/", ""), "cleaned_DV_medication_start_stop_dates.Rdata"))

# Clean medication data
dat_med <- cleaned_dataset %>%
  # Set column names to lower case
  set_colnames(tolower(colnames(.))) %>%
  # Clean variables
  mutate(# Dates from POSIXct format to Date format
         across(start_date:stop_date, as.Date),
         # Medication type to lower case values only
         drug = tolower(medication_type)) %>%
  # Drop medication type column
  select(-medication_type) %>%
  # Rename date variables
  rename(start_dt = start_date,
         stop_dt = stop_date)

# Save medication data
save(dat_med,
     file = paste0(path, "medication.Rdata"))

# 2. Comorbidities ----
## 2.1. Diabetes mellitus ----
# ICD-10 codes E10-E14
# Any glucose lowering agent
# HbA1c mmol/mol >= 48 or HbA1c >= 6.5%
# Glucose >= 7 as fasting or >= 11.1 as non-fasting

# Individuals with diagnostic codes
vec_diag <- dat_spine %>%
  # Join diagnostic codes
  left_join(dat_proc, "id") %>%
  # Keep only diabetes diagnoses prior to discharge date
  filter(code == "ICD10_DM" &
         code_dt <= discharge_dt) %>%
  # Keep only one ID per individual
  distinct(id) %>%
  # Reduce to vector
  extract2("id")

# Individuals with medication
vec_med <- dat_spine %>%
  # Join diagnostic codes
  left_join(dat_med, "id") %>%
  # Keep only glucose lowering drugs that were an active prescription in the last year prior to discharge
  filter(drug == "gluc_lowering_med" &
         stop_dt >= discharge_dt - 365 &
         stop_dt <= discharge_dt) %>%
  # Keep only one ID per individual
  distinct(id) %>%
  # Reduce to vector
  extract2("id")

# HbA1c measurements
vec_hba1c <- dat_spine %>%
  # Join lab data
  left_join(dat_lab, "id") %>%
  # Keep only HbA1c >= 6.5% within three months prior to discharge date
  filter(code == "HbA1c_bl_procent_lab" &
         lab_dt >= discharge_dt - 90 &
         lab_dt <= discharge_dt &
         result >= 6.5) %>%
  # Keep only one ID per individual
  distinct(id) %>%
  # Reduce to vector
  extract2("id")

# Glucose fasting measurements
vec_gluc_f <- dat_spine %>%
  # Join lab data
  left_join(dat_lab, "id") %>%
  # Keep only HbA1c >= 6.5% within three months prior to discharge date
  filter(code == "Glucose_f_bl_mmol_L_lab" &
         lab_dt >= discharge_dt - 90 &
         lab_dt <= discharge_dt &
         result >= 7) %>%
  # Keep only one ID per individual
  distinct(id) %>%
  # Reduce to vector
  extract2("id")

# Glucose non-fasting measurements
vec_gluc_nf <- dat_spine %>%
  # Join lab data
  left_join(dat_lab, "id") %>%
  # Keep only HbA1c >= 6.5% within three months prior to discharge date
  filter(code == "Glucose_nf_bl_mmol_L_lab" &
           lab_dt >= discharge_dt - 90 &
           lab_dt <= discharge_dt &
           result >= 11.1) %>%
  # Keep only one ID per individual
  distinct(id) %>%
  # Reduce to vector
  extract2("id")

# Final vector for individuals with diabetes
vec_diab <- unique(c(vec_diag, vec_med, vec_hba1c, vec_gluc_f, vec_gluc_nf))

# Drop other vectors
rm(vec_diag, vec_med, vec_hba1c, vec_gluc_f, vec_gluc_nf)

## 2.2. Heart failure ----
# ICD10-codes I11.0, I13.0, I13.2, I50
# Procedures: left ventricular assist device

# Diagnostic codes & procedures
vec_hf <- dat_spine %>%
  # Join diagnostic codes
  left_join(dat_proc, "id") %>%
  # Keep only diabetes diagnoses prior to discharge date
  filter(code %in% c("ICD10_HF", "Procedure_surgery_cardiothoracic_artificial_heart") &
         code_dt <= discharge_dt) %>%
  # Keep only one ID per individual
  distinct(id) %>%
  # Reduce to vector
  extract2("id")

## 2.3. Hypertension ----
# ICD10 codes I10-I15
# Any antihypertensive agent
# Two consecutive measurements, at least 90 days apart, with either SBP >= 140 mmHg or DBP >= 90 mmHg

# Diagnostic codes
vec_diag <- dat_spine %>%
  # Join diagnostic codes
  left_join(dat_proc, "id") %>%
  # Keep only hypertension codes in the year prior to discharge
  filter(code == "ICD10_HT" &
         code_dt >= discharge_dt - 365 &
         code_dt <= discharge_dt) %>%
  # Keep only one ID per individual
  distinct(id) %>%
  # Reduce to vector
  extract2("id")

# Any antihypertensive agent
vec_med <- dat_spine %>%
  # Join diagnostic codes
  left_join(dat_med, "id") %>%
  # Keep only glucose lowering drugs that were an active prescription in the last year prior to discharge
  filter(drug == "any_antihypertensives_med" &
         stop_dt >= discharge_dt - 365) %>%
  # Keep only one ID per individual
  distinct(id) %>%
  # Reduce to vector
  extract2("id")

# Final vector for individuals with hypertension
vec_hyp <- unique(c(vec_diag, vec_med))

# Remove other vectors
rm(vec_diag, vec_med)

## 2.4. Coronary artery/heart disease ----
# ICD-10 Codes: I20; I21; I22; I23; I24; I25
# Procedures: Coronary artery bypass grafting or percutaneous coronary Intervention

# Individuals with CAD
vec_cad <- dat_spine %>%
  # Join diagnostic codes
  left_join(dat_proc, "id") %>%
  # Keep only hypertension codes prior to discharge
  filter(code %in% c("ICD10_CHD", 
                     "Procedure_surgery_cardiothoracic_CABG",
                     "Procedure_angiogram_PCI") &
         code_dt <= discharge_dt) %>%
  # Keep only one ID per individual
  distinct(id) %>%
  # Reduce to vector
  extract2("id")

## 2.5. Cerebrovascular disease ----
# ICD-10 Codes: I60-64; I69.3; I69.8; I69.4
# Procedures: Carotid endarterectomy or carotid artery stenting

# Individuals with CBVD
vec_cbvd <- dat_spine %>%
  # Join diagnostic codes
  left_join(dat_proc, "id") %>%
  # Keep only hypertension codes prior to discharge
  filter(code %in% c("ICD10_CVD"
                     ) &
           code_dt <= discharge_dt) %>%
  # Keep only one ID per individual
  distinct(id) %>%
  # Reduce to vector
  extract2("id")

## 2.6. COPD ----
# ICD-10 Codes: J43; J44

# Individuals with COPD
vec_copd <- dat_spine %>%
  # Join diagnostic codes
  left_join(dat_proc, "id") %>%
  # Keep only hypertension codes prior to discharge
  filter(code == "ICD10_COPD" &
         code_dt <= discharge_dt) %>%
  # Keep only one ID per individual
  distinct(id) %>%
  # Reduce to vector
  extract2("id")

## 2.7. Liver cirrhosis ----
# ICD-10 Codes: K74.3; K74.4; K74.5; K74.6; K70.3; K70.4; K76.6; I85.0; I85.9; I98.2

# Individuals with liver cirrhosis
vec_liver <- dat_spine %>%
  # Join diagnostic codes
  left_join(dat_proc, "id") %>%
  # Keep only hypertension codes prior to discharge
  filter(code == "ICD10_Liverdisease" &
         code_dt <= discharge_dt) %>%
  # Keep only one ID per individual
  distinct(id) %>%
  # Reduce to vector
  extract2("id")

## 2.8. Malignancy ----
# ICD-10 Codes: C00-97
# Procedures: surgical removal of a tumor, the administration of cytostatic drugs, bone marrow transplantation

# Individuals with a malignancy
vec_mal <- dat_spine %>%
  # Join diagnostic codes
  left_join(dat_proc, "id") %>%
  # Keep only hypertension codes prior to discharge
  filter(code %in% c("ICD10_Malignancy",
                     "Procedure_BMT") &
         code_dt <= discharge_dt) %>%
  # Keep only one ID per individual
  distinct(id) %>%
  # Reduce to vector
  extract2("id")

## 2.9. Arrhythmias/conduction disorders ----
# ICD-10 Codes: I44-49; Z99.4; Z95.0; T82.1; Z45.0
# Procedures: ICD/pacemaker implantation, cardiac defibrillation, and pulmonary vein isolation (ablation)

# Individuals with arrhythmias
vec_arr <- dat_spine %>%
  # Join diagnostic codes
  left_join(dat_proc, "id") %>%
  # Keep only hypertension codes prior to discharge
  filter(code %in% c("ICD10_Arrhythmia",
                     "Procedure_arrhythmia") &
           code_dt <= discharge_dt) %>%
  # Keep only one ID per individual
  distinct(id) %>%
  # Reduce to vector
  extract2("id")

## 2.10. Peripheral artery disease ----
# ICD-10 Codes: I70.2
# Procedures: revascularization or bypass of peripheral vessels

# Individuals with PAD
vec_pad <- dat_spine %>%
  # Join diagnostic codes
  left_join(dat_proc, "id") %>%
  # Keep only hypertension codes prior to discharge
  filter(code %in% c("ICD10_PAD") &
           code_dt <= discharge_dt) %>%
  # Keep only one ID per individual
  distinct(id) %>%
  # Reduce to vector
  extract2("id")

## 2.11. Autoimmune disease ----
# ICD-10 Codes: M05-14; M30-36; D50-89; K75.4; E06.3; G35; E27.1; M45; M35.3

# Individuals with autoimmune disease
vec_aid <- dat_spine %>%
  # Join diagnostic codes
  left_join(dat_proc, "id") %>%
  # Keep only hypertension codes prior to discharge
  filter(code == "ICD10_Autoimmunedisease" &
         code_dt <= discharge_dt) %>%
  # Keep only one ID per individual
  distinct(id) %>%
  # Reduce to vector
  extract2("id")

## 2.12. Polycystic kidney disease ----
# ICD-10 Codes: Q60-63; Q85.1

# Individuals with PCKD
vec_pckd <- dat_spine %>%
  # Join diagnostic codes
  left_join(dat_proc, "id") %>%
  # Keep only hypertension codes prior to discharge
  filter(code == "ICD10_Hereditary_congenital_and_polycystic_kidneydisease" &
         code_dt <= discharge_dt) %>%
  # Keep only one ID per individual
  distinct(id) %>%
  # Reduce to vector
  extract2("id")

## 2.13. Heart/lung transplantation ----
# ICD-10 Codes: Z94.1; Z94.2; T86.2; T86.3
# Procedures: heart or lung transplantation

# Individuals with Heart/lung Tx
vec_hltx <- dat_spine %>%
  # Join diagnostic codes
  left_join(dat_proc, "id") %>%
  # Keep only hypertension codes prior to discharge
  filter(code %in% c("ICD10_Heart_or_LungTx", 
                     "Procedure_surgery_cardiothoracic_Heart_or_LungTx") &
         code_dt <= discharge_dt) %>%
  # Keep only one ID per individual
  distinct(id) %>%
  # Reduce to vector
  extract2("id")

## 2.14. Previous AKI ----
# Individuals with previous AKI based on ICD-10 codes
vec_paki <- dat_spine %>%
  # Join diagnostic codes
  left_join(dat_proc, "id") %>%
  # Keep only hypertension codes prior to admission
  filter(code == "ICD10_AKI" &
         code_dt <= admission_dt) %>%
  # Keep only one ID per individual
  distinct(id) %>%
  # Reduce to vector
  extract2("id")

# 3. Prescriptions ----
## 3.1. Glucose-lowering agents ----
# Derive vector
vec_gla <- dat_spine %>%
  # Join diagnostic codes
  left_join(dat_med, "id") %>%
  # Keep only glucose lowering drugs that were an active prescription in the last year prior to discharge
  filter(drug == "gluc_lowering_med" &
         stop_dt >= discharge_dt - 365) %>%
  # Keep only one ID per individual
  distinct(id) %>%
  # Reduce to vector
  extract2("id")

## 3.2. Antihypertensives ----
# Derive vector
vec_aht <- dat_spine %>%
  # Join diagnostic codes
  left_join(dat_med, "id") %>%
  # Keep only glucose lowering drugs that were an active prescription in the last year prior to discharge
  filter(drug == "any_antihypertensives_med" &
         stop_dt >= discharge_dt - 365) %>%
  # Keep only one ID per individual
  distinct(id) %>%
  # Reduce to vector
  extract2("id")

## 3.3. Lipod-lowering agents ----
# Derive vector
vec_lla <- dat_spine %>%
  # Join diagnostic codes
  left_join(dat_med, "id") %>%
  # Keep only glucose lowering drugs that were an active prescription in the last year prior to discharge
  filter(drug == "statin_med" &
         stop_dt >= discharge_dt - 365) %>%
  # Keep only one ID per individual
  distinct(id) %>%
  # Reduce to vector
  extract2("id")

## 3.4. Antithrombotic agents ----
# Derive vector
vec_ata <- dat_spine %>%
  # Join diagnostic codes
  left_join(dat_med, "id") %>%
  # Keep only glucose lowering drugs that were an active prescription in the last year prior to discharge
  filter(drug == "antithrom_med" &
         stop_dt >= discharge_dt - 365) %>%
  # Keep only one ID per individual
  distinct(id) %>%
  # Reduce to vector
  extract2("id")

## 3.5. Immunosuppressors ----
# Derive vector
vec_ims <- dat_spine %>%
  # Join diagnostic codes
  left_join(dat_med, "id") %>%
  # Keep only glucose lowering drugs that were an active prescription in the last year prior to discharge
  filter(drug %in% c("imsup_med",
                     "imsup_nephotoxic_med") &
         stop_dt >= discharge_dt - 365) %>%
  # Keep only one ID per individual
  distinct(id) %>%
  # Reduce to vector
  extract2("id")

## 3.6. Any nephrotoxic drug ----
# Derive vector
vec_ntd <- dat_spine %>%
  # Join diagnostic codes
  left_join(dat_med, "id") %>%
  # Keep only glucose lowering drugs that were an active prescription in the last year prior to discharge
  filter(drug %in% c("antivir_nephtoxic_med",
                     "any_nephtoxic_med",
                     "imsup_nephotoxic_med") &
           stop_dt >= discharge_dt - 365) %>%
  # Keep only one ID per individual
  distinct(id) %>%
  # Reduce to vector
  extract2("id")

# 4. Hospitalisation characteristics ----
## 4.1. Length of stay ----
# Calculate LOS immediately in spine data
dat_spine %<>% mutate(los = as.numeric(discharge_dt - admission_dt))

## 4.2. ICU admission ----

## 4.3. Primary reason for admission ----

## 4.4. Cardiothoracic surgery ----
# Individuals with cardiothoracic surgery
vec_cts <- dat_spine %>%
  # Join diagnostic codes
  left_join(dat_proc, "id") %>%
  # Keep only surgery codes during admission
  filter(str_detect(code, "surgery_cardiothoracic") &
         code_dt >= admission_dt &
         code_dt <= discharge_dt) %>%
  # Keep only one ID per individual
  distinct(id) %>%
  # Reduce to vector
  extract2("id")

## 4.5. Non-cardiothoracic surgery ----
# Individuals with non-cardiothoracic surgery
vec_ncts <- dat_spine %>%
  # Join diagnostic codes
  left_join(dat_proc, "id") %>%
  # Keep only surgery codes during admission
  filter(str_detect(code, "surgery_(?!cardiothoracic)") &
         code_dt <= discharge_dt) %>%
  # Keep only one ID per individual
  distinct(id) %>%
  # Reduce to vector
  extract2("id")

## 4.6. Resuscitation ----
# Individuals with resuscitation
vec_resc <- dat_spine %>%
  # Join diagnostic codes
  left_join(dat_proc, "id") %>%
  # Keep only codes during admission
  filter(code %in% c("ICD10_Resuscitation", 
                     "Procedure_Resuscitation") &
         code_dt >= admission_dt &
         code_dt <= discharge_dt) %>%
  # Keep only one ID per individual
  distinct(id) %>%
  # Reduce to vector
  extract2("id")

# 5. Combine all data ----
# Vector with names of all data vectors
vec_data <- ls()[str_detect(ls(), "vec_")]

# Add vectors to data
for(i in vec_data){
  # Add indicator to data
  dat_spine %<>% 
    # Create indicator if ID is present in vector with IDs
    mutate(ind = if_else(id %in% get(i), 1, 0)) %>%
    # Rename indicator to vector
    set_colnames(c(colnames(.)[1:length(colnames(.)) - 1], i))
}

# Rename and re-order columns
dat_spine %<>%
  # Renaming
  rename(med_aht = vec_aht,       # Anti-hypertensives    
         com_aid = vec_aid,       # Autoimmune disease    
         com_arr = vec_arr,       # Arrhythmias
         med_ata = vec_ata,       # Anti-thrombotic agents   
         com_cad = vec_cad,       # Cororary artery disease    
         com_cbvd = vec_cbvd,     # Cerebrovascular disease      
         com_copd = vec_copd,     # COPD      
         ihe_cts = vec_cts,       # Cardiothoracic surgery (in-hospital event) 
         com_diab = vec_diab,     # Diabetes      
         med_gla = vec_gla,       # Glucose-lowering agents    
         com_hf = vec_hf,         # Heart failure  
         com_hltx = vec_hltx,     # Heart/lung Tx      
         com_hyp = vec_hyp,       # Hypertension    
         med_ims = vec_ims,       # Immunosuppressive medication    
         com_liver = vec_liver,   # Liver disease        
         med_lla = vec_lla,       # Lipid-lowering agents    
         com_mal = vec_mal,       # Malignancy    
         ihe_ncts = vec_ncts,     # Non-cardiothoracic surgery (in-hospital event)      
         med_ntd = vec_ntd,       # Nephrotoxic drugs    
         com_pad = vec_pad,       # Peripheral artery disease    
         com_paki = vec_paki,     # Previous AKI      
         com_pckd = vec_pckd,     # Polycystic kidney disease      
         ihe_resc = vec_resc      # Resuscitation (in-hospital event)    
  ) %>%# 
  # Change order
  relocate(com_arr, com_aid, com_cad, com_cbvd, com_copd, com_diab, com_hf, com_hltx, com_hyp, com_liver, com_mal, com_pad,
           com_paki, com_pckd, med_aht, med_ata, med_gla, med_ims, med_lla, med_ntd, ihe_cts, ihe_ncts, ihe_resc,
           .after = los)

# Save data
save(dat_spine,
     file = paste0(path, "covariate_derivation.Rdata"))
