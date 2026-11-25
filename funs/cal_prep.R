# Function to prepare data for calibration plots
cal_prep <- function(dat_outcome, dat_prd, timepoints, outcomes){
  # Drop censoring from outcome vector (calib_msm already takes it into account)
  outcomes <- setdiff(outcomes, "censored")
  
  # Add tte variables to outcome vector
  outcomes_tte <- c(outcomes, paste0(outcomes, "_tte"))
  
  # Prepare outcome data with censoring indicator
  dat_tmp <- dat_outcome %>%
    # Create tte for censoring
    mutate(dtcens = as.numeric(final_dt - discharge_dt)) %>%
    # Pivot to long format 
    # This gives a warning because tte column will be NA, but we fix this below
    pivot_longer(all_of(outcomes_tte),
                 names_sep = "_(?=tte|$)",
                 names_to = c("event", "tte"),
                 values_to = "occurred") %>%
    # Set NA in tte column to indicate occurrence
    mutate(tte = if_else(is.na(tte), "occurred", tte)) %>%
    # Pivot tte back to wide format
    pivot_wider(names_from = "tte",
                values_from = "occurred") %>%
    # Suppress warning that occurrs due to pivot_longer (see above annotation)
    suppressWarnings()
  
  # Add timepoints to create raw data
  dat_raw <- expand.grid(id = unique(dat_tmp[["id"]]),
                         time = timepoints) %>%
    # Join dat_tmp completely
    left_join(dat_tmp, 
              "id",
              relationship = "many-to-many") %>%
    # Arrange for grouping
    arrange(id, time) %>%
    # Group per individual per timepoint
    group_by(id, time) %>%
    # Create indicators for calibmsm package
    mutate(
      # First, check if each outcome occurred at the timepoint of interest
      occurred = if_else(tte > time, 0, occurred),
      # Then, reset tte
      tte = if_else(occurred == 1, tte, pmin(time, dtcens)),
      # Now determine if someone was censored
      dtcens_s = if_else(sum(occurred) == 0, 1, 0)) %>%
    # Keep one row per individual
    # It does not matter which one as we only need the covariates for weighting and dtcens & dtcens_s, not event & occurred
    slice(1L) %>%
    # Remove grouping structure
    ungroup()
  
  # Determine timepoints closest before timepoint of interest
  tps <- sapply(timepoints, \(x) dat_prd[["time"]][which.max(if_else(dat_prd[["time"]] - x <= 0, dat_prd[["time"]] - x, NA))])
  
  # Subset predictions data to relevant timepoints
  dat_tp_pred <- filter(dat_prd, time %in% tps) 
  
  # Vector of current times
  vec_times <- sort(unique(dat_tp_pred[["time"]]))
  
  # Change available times to intended timepoints
  for(i in 1:length(timepoints)) dat_tp_pred %<>% mutate(time = if_else(time == vec_times[[i]], timepoints[[i]], time))
  
  # Create output object
  output <- list(
    dat_raw = dat_raw,
    dat_prd = dat_tp_pred
  )
  
  # Return output
  return(output)
}
