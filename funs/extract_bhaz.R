# Function to extract baseline hazards at relevant timepoints from model fit object
extract_bhaz <- function(.fit, timepoints){
  # Get baseline hazards
  dat_bhaz <- basehaz(.fit) 
  
  # Determine timepoints closest before timepoint of interest
  tps <- sapply(timepoints, \(x) dat_bhaz[["time"]][which.max(if_else(dat_bhaz[["time"]] - x <= 0, dat_bhaz[["time"]] - x, NA))])
  
  # Further prepare data
  dat_bhaz %<>% 
    # Filter timepoints of interest
    filter(time %in% tps) %>%
    # Clean up transition strata indicator
    

  # Change available times to intended timepoints
  for(i in 1:length(timepoints)) dat_bhaz %<>% mutate(time = if_else(time == tps[[i]], timepoints[[i]], time))
  
  
}