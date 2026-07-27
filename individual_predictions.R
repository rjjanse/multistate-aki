# Create predictor object
vec_predictors <- c("age", "female", "stage", "los", "com_arr", "com_cad", 
                    "com_cbvd", "com_copd", "com_diab", "com_hf", "com_hyp", 
                    "com_mal", "med_ims", "med_lla", "med_ntd", "ihe_cts", 
                    "ihe_ncts")

# Create transition matrix (for endpoint model)
mat_trans <- structure(c(NA, NA, NA, NA, NA, 1L, NA, NA, NA, 
                         NA, 2L, NA, NA, NA, NA, NA, 4L, 6L, 
                         NA, NA, 3L, 5L, 7L, 8L, NA), 
                       dim = c(5L, 5L), 
                       dimnames = 
                         list(from = c("Baseline", "KF", "MACCE", 
                                       "MACCE_KF", "Death"), 
                              to = c("Baseline", "KF", "MACCE", 
                                     "MACCE_KF", "Death")))

# Source prd_mstate function
source("https://github.com/rjjanse/multistate-aki/raw/refs/heads/main/funs/prd_mstate.R")

# Source dat_prep function
source("https://github.com/rjjanse/multistate-aki/raw/refs/heads/main/funs/dat_prep.R")