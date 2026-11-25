# Function to calculate, format, and rescale PDI
pdi <- function(.data, outcomes, timepoint){
  # Calculate PDI
  tab_pdi <- pdiest(
    # Use prepped PDI data
    .data %>%
      # Subset data to timepoint of interest
      filter(time == timepoint) %>%
      # Drop time variable
      select(-time)
  ) %>%
    # Change to data frame
    as.data.frame() %>%
    # Rename PDI column
    rename(pdi = 1) %>%
    # Add additional columns
    mutate(
      # Indicator for outcome
      outcome = c("overall", outcomes),
      # Indicator for timeopint
      time = timepoint) %>%
    # Relocate indicator
    relocate(outcome, time, .before = pdi) %>%
    # Set to tibble
    as_tibble()
    
  # Return data
  return(tab_pdi)
}
