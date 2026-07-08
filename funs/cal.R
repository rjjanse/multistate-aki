# Function to get calibration plot data and calibration metrics
cal <- function(dat_mstate, 
                dat_raw, 
                dat_prd, 
                timepoint, 
                j = 1, # From state
                s = 0, # From time
                method = "mlr", 
                smoother = "loess", 
                w_covs){
  # Subset dat_prd to timepoint of interest
  dat_prd_subset <- filter(dat_prd, time == timepoint)
  
  # Determine individuals with invalid predictions (set to missing in script 4)
  vec_invalid_ids <- dat_prd_subset %>%
    # Keep only individuals with any missing prediction
    filter(!(id %in% na.omit(.)[["id"]])) %>%
    # Extract IDs
    extract2("id")
  
  # Subset dat_raw to timepoint of interest and drop individuals with an invalid prediction
  dat_raw_subset <- filter(dat_raw, time == timepoint & !(id %in% vec_invalid_ids))
  
  # Subset dat_prd_subset further based on invalid IDs 
  dat_prd_subset %<>% filter(!(id %in% vec_invalid_ids)) %>%
    # Change to data frame
    # A tibble will throw: Error in Ops.data.frame(tp_pred_mlr[, 2:ncol(tp_pred_mlr)], tp_pred_mlr[,  : 
    # ‘/’ only defined for equally-sized data frames
    as.data.frame() %>%
    # Keep only columns with predictions
    select(contains("pstate"))
  
  # Subset dat_mstate based on invalid IDs
  dat_mstate_subset <- filter(dat_mstate, !(id %in% vec_invalid_ids))
  
  # Calculate calibration
  lst_cal_dat <- calib_msm(
    data_ms = dat_mstate_subset,
    data_raw = dat_raw_subset,
    j = j, # From state j
    s = s, # From time s
    t = timepoint, # At time t
    tp_pred = dat_prd_subset,
    calib_type = method,
    curve_type = smoother,
    w_covs = w_covs,
    assess_moderate = TRUE,   # Calibration plots
    assess_mean = TRUE        # Calibration-in-the-large
  )
  
  # Extract calibration-in-the-large
  vec_citl <- lst_cal_dat[["mean"]]
  
  # Extract data for calibration plots
  dat_cal <- lst_cal_dat[["plotdata"]]
  
  # Number of individuals dropped due to invalid predictions
  n_dropped <- n_distinct(dat_prd[["id"]]) - n_distinct(dat_raw_subset[["id"]])
  
  # Create output object
  output <- list(
    citl = vec_citl,
    dat_plot = dat_cal,
    n_dropped = n_dropped
  )
  
  # Return output
  return(output)
}
