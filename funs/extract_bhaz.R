# Function to extract baseline hazards at relevant timepoints from model fit object
extract_bhaz <- function(.fit, timepoints){
  # Get baseline hazards
  dat_bhaz <- basehaz(.fit) 
  
  # Further prepare data per stratum (as the closest available timepoint differs per transition
  dat_bhaz <- bind_rows(lapply(unique(dat_bhaz[["strata"]]), \(x){
    # Subset dat_bhaz to stratum
    dat_bhaz_strat <- filter(dat_bhaz, strata == x)
    
    # Get timepoints
    tps <- sapply(timepoints, \(x) dat_bhaz_strat[["time"]][which.max(if_else(dat_bhaz_strat[["time"]] - x <= 0, dat_bhaz_strat[["time"]] - x, NA))])
    
    # Keep only timepoints of interest
    dat_bhaz_strat %<>% 
      # Filter timepoints of interest
      filter(time %in% tps) %>%
      # Clean up transition strata indicator
      mutate(strata = as.numeric(str_replace(strata, "trans=", ""))) 
    
    # Change available times to intended timepoints
    for(i in 1:length(timepoints)) dat_bhaz_strat %<>% mutate(time = if_else(time == tps[[i]], timepoints[[i]], time))
      
    # Return timepoints
    return(x = dat_bhaz_strat)
  })) %>%
    # Pivot data to wide format
    pivot_wider(id_cols = time,
                names_from = strata,
                values_from = hazard)
  
  # Return data
  return(dat_bhaz)
}
