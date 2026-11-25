#----------------------------------------------------------#
# Predicting outcomes after AKI using multistate models
# Code for model development
# Roemer J. Janse - Last updated on 2026-02-12
#----------------------------------------------------------#

# 0. Set-up ----
# Load packages
pacman::p_load("dplyr",          # Data wrangling
               "magrittr",       # Efficient pipelines
               "tidyr",          # Data cleaning
               "stringr",        # Working with strings
               "purrr",          # Functional programming
               "here",           # Local paths
               "survival",       # Survival models
               "mstate",         # Multistate models
               "conflicted"      # Resolve function conflicts
)

# Resolve function conflicts
conflicts_prefer(dplyr::filter) # Between stats and dplyr

# Load functions
walk(list.files(here("funs")), \(x) source(paste0(here("funs"), "/", x)))

# Set data path
path <- psp(unix = "/volumes/Data/lab/lab_research/RES-Folder-UPOD/NOSTRADAMUS_SALTRO/E_ResearchData/2_ResearchData/CLEANED_for_Multistate_outcomes_project/20012026/dataframes/",
            windows = "L:/lab_research/RES-Folder-UPOD/NOSTRADAMUS_SALTRO/E_ResearchData/2_ResearchData/CLEANED_for_Multistate_outcomes_project/20012026/dataframes/")

## Load spine dataframes
loadr("outcome_derivation_em.Rdata")   # Endpoint model
loadr("outcome_derivation_tm.Rdata")   # Trajectory models

# 1. Endpoint model ----
## 1.1. Prepare data ----
### 1.1.1. Transition matrix ----
# Set-up transition matrix 
mat_trans <- transMat(x = list(c(2:3, 5),     # From baseline
                               4:5,           # From KF
                               4:5,           # From MACCE
                               5,             # From MACCE + KF
                               c()),          # From death           
                      names =  c("Baseline",  # State 1
                                 "KF",        # State 2
                                 "MACCE",     # State 3
                                 "MACCE_KF",  # State 4
                                 "Death"))    # State 5

### 1.1.2. Data transformation ----
# Get vector of covariates that are to be modelled
vec_covars <- setdiff(colnames(select(dat_spine_em, age, female, stage, los:ihe_resc)), 
                      c("com_aid", "med_aht", "med_ata", "med_gla", "com_hltx", "com_liver", "com_pad", "com_pckd", "ihe_resc"))

# Save vec_covars for calibration (script 5)
saver(vec_covars, "vec_covars.Rdata")

# Get vector of TTE variables
vec_tte <- c(NA, "kf_tte", "macce_tte", "macce_kf_tte", "death_tte")

# Get vector of status indicators
vec_stat <- c(NA, "kf", "macce", "macce_kf", "death")

# Transform data
dat_mstate <- msprep(time = vec_tte,
                     status = vec_stat,
                     data = dat_spine_em,
                     trans = mat_trans,
                     keep = vec_covars,
                     id = "id") %>%
  # For individuals with Tstop at Tstart, add one day to Tstop
  mutate(Tstop = if_else(Tstop == Tstart, Tstop + 1, Tstop),
         # Also recalculate time
         time = Tstop - Tstart) 

# Mutate drops transition matrix attribute, here we add it back
attr(dat_mstate, "trans") <- mat_trans

# Expand covariates
dat_mstate <- expand.covs(dat_mstate, vec_covars, append = TRUE)

# Event count
events(dat_mstate)[["Frequencies"]]

# Save mstate data (for calibration in script 5)
saver(dat_mstate, "dat_mstate_em.Rdata")

## 1.2. Fit model ----
# Model right hand side
model_rhs <- paste0(colnames(dat_mstate)[str_detect(colnames(dat_mstate), "\\.\\d")], collapse = " + ")

# Fit model on mstate data
fit <- coxph(as.formula(paste0("Surv(Tstart, Tstop, status) ~ strata(trans) + ", model_rhs)),
             data = dat_mstate, 
             method = "breslow",
             iter.max = 1e4)

# Save model fit
save(fit, file = paste0(path, "mstate_fit_em.Rdata"))

## 1.3. Get model predictions ----
# Load model
loadr("mstate_fit_em.Rdata")

# Get predictions from model
dat_prds <- prd_mstate(dat_spine_em, vec_covars, mat_trans)

# Save predictions
saver(dat_prds, "predictions_em.Rdata")

# Fix predictions
dat_prds_fixed <- prd_fix(dat_prds)

# Save fixed predictions
saver(dat_prds_fixed, "predictions_fixed_em.Rdata")

# 2. Transition model without AKI ----
## 2.1. Prepare data ----
### 2.1.1. Transition matrix ----
# Set-up transition matrix 
mat_trans <- transMat(x = list(c(2, 4, 7),    # From baseline
                               c(3, 5, 7),    # From CKD
                               c(6, 7),       # From KF
                               c(5, 7),       # From MACCE
                               6:7,           # From MACCE + CKD
                               7,             # From MACCE + KF
                               c()),          # From death           
                      names =  c("Baseline",  # State 1
                                 "CKD",       # State 2
                                 "KF",        # State 3
                                 "MACCE",     # State 4
                                 "MACCE_CKD", # State 5
                                 "MACCE_KF",  # State 6
                                 "Death"))    # State 7

### 2.1.2. Data transformation ----
# Get vector of covariates that are to be modelled
vec_covars <- setdiff(colnames(select(dat_spine_tm, age, female, stage, los:ihe_resc)), 
                      c("com_aid", "med_aht", "med_ata", "med_gla", "com_hltx", "com_liver", "com_pad", "com_pckd", "ihe_resc"))

# Get vector of TTE variables
vec_tte <- c(NA, "ckd_tte", "kf_tte", "macce_tte", "macce_ckd_tte", "macce_kf_tte", "death_tte")

# Get vector of status indicators
vec_stat <- c(NA, "ckd", "kf", "macce", "macce_ckd", "macce_kf", "death")

# Transform data
dat_mstate <- msprep(time = vec_tte,
                     status = vec_stat,
                     data = dat_spine_tm,
                     trans = mat_trans,
                     keep = vec_covars,
                     id = "id") %>%
  # For individuals with Tstop at Tstart, add one day to Tstop
  mutate(Tstop = if_else(Tstop == Tstart, Tstop + 1, Tstop),
         # Also recalculate time
         time = Tstop - Tstart) 

# Mutate drops transition matrix attribute, here we add it back
attr(dat_mstate, "trans") <- mat_trans

# Expand covariates
dat_mstate <- expand.covs(dat_mstate, vec_covars, append = TRUE)

# Event count
events(dat_mstate)[["Frequencies"]]

# Save mstate data (for calibration in script 5)
saver(dat_mstate, "dat_mstate_tm.Rdata")

## 2.2. Fit model ----
# Model right hand side
model_rhs <- paste0(colnames(dat_mstate)[str_detect(colnames(dat_mstate), "\\.\\d")], collapse = " + ")

# Fit model on mstate data
fit <- coxph(as.formula(paste0("Surv(Tstart, Tstop, status) ~ strata(trans) + ", model_rhs)),
             data = dat_mstate, 
             method = "breslow",
             iter.max = 100)

# Save model fit
save(fit, file = paste0(path, "mstate_fit_tm.Rdata"))

## 2.3. Get model predictions ----
# Load model
loadr("mstate_fit_tm.Rdata")

# Get predictions from model
dat_prds <- prd_mstate(dat_spine_tm, vec_covars, mat_trans)

# Save predictions
saver(dat_prds, "predictions_tm.Rdata")

# Fix predictions
dat_prds_fixed <- prd_fix(dat_prds)
  
# Save fixed predictions
saver(dat_prds_fixed, "predictions_fixed_tm.Rdata")
