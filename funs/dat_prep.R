# Function to prepare data into format for prd_mstate function
dat_prep <- function(
    age,                  # Age in years
    female,               # 1 if female, 0 if male
    aki_stage,            # AKI Stage (1-3 or 4 if dialysis was required)
    los,                  # Length of stay in hospital (days)
    # All below arguments are 1 if yes or 0 if no
    arrhythmia,           # History of arrhythmia
    cad,                  # History of coornary artery disease
    cbvd,                 # History of cerebrovascular disease
    copd,                 # History of chronic obstructive pulmonary disease
    diabetes,             # History of diabetes mellitus
    heartfailure,         # History of heart failure
    hypertension,         # History of hypertension
    malignancy,           # History of malignancy
    immunosuppresives,    # Use of immunosuppresive medication
    lipidloweringagents,  # Use of lipid lowering agents
    nephrotoxicdrugs,     # Use of nehprotoxic medication
    cardiothoracicsurg,   # Cardiothoracic surgery during index hospitalisation
    noncardiothoracicsurg # Non-cardiothoracic surgery during index hospitalisation
){
  # Create data frame
  dat_predictors <- data.frame(id = 1:length(age),
                               age = age,
                               female = female,
                               stage = aki_stage,
                               los = los,
                               com_arr = arrhythmia,
                               com_cad = cad,
                               com_cbvd = cbvd,
                               com_copd = copd,
                               com_diab = diabetes,
                               com_hf = heartfailure,
                               com_hyp = hypertension,
                               com_mal = malignancy,
                               med_ims = immunosuppresives,
                               med_lla = lipidloweringagents,
                               med_ntd = nephrotoxicdrugs,
                               ihe_cts = cardiothoracicsurg,
                               ihe_ncts = noncardiothoracicsurg)
  
  # Return data 
  return(dat_predictors)
}
