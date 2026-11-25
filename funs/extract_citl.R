# Function to extract and format CITL data after cal()
extract_citl <- function(.list, x, outcomes){
  # Prepare data
  dat_citl <- .list[[x]][["citl"]] %>%
    # Set to tibble
    as_tibble() %>%
    # Change column name
    rename(citl = 1) %>%
    # Add indicators
    mutate(
      # State indicator
      state = names(.list[[x]][["citl"]]),
      # Timepoint indicator
      time = timepoints[[x]]
    )
  # All states
  states <- sort(unique(dat_citl[["state"]]))
  
  # Rename states
  for(i in 1:length(outcomes)) dat_citl %<>% mutate(state = if_else(state == states[[i]], outcomes[[i]], state))
  
  # Return data
  return(dat_citl)
}
