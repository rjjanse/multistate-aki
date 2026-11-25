# Function to reshape prediction data into format for PDI calculation
pdi_prep <- function(.data, timepoints, dat_outcome, outcomes){
  # Get columns
  vec_cols <- setdiff(colnames(.data), c("time", "id"))
  
  # Add tte variables to outcome vector
  outcomes_tte <- c(outcomes, paste0(outcomes, "_tte"))
  
  # Determine timepoints closest before timepoint of interest
  tps <- sapply(timepoints, \(x) .data[["time"]][which.max(if_else(.data[["time"]] - x <= 0, .data[["time"]] - x, NA))])
  
  # Add outcomes that occurred
  dat_pdi <- .data %>%
    # Keep only timepoint of interest
    filter(time %in% tps) %>%
    # Join outcome
    left_join(dat_outcome %>%
                # Determine first event time
                mutate(first_event_tte = apply(select(., all_of(setdiff(paste0(outcomes, "_tte"), "censored_tte"))), 1, min)) %>%
                # Select relevant columns
                select(id, all_of(outcomes_tte), first_event_tte) %>%
                # Pivot to long format 
                # This gives a warning because tte column will be NA, but we fix this below
                pivot_longer(all_of(outcomes_tte),
                             names_sep = "_(?=tte|$)",
                             names_to = c("event", "tte"),
                             values_to = "occurred") %>%
                # Set NA in tte column to indicate occurrence
                mutate(tte = if_else(is.na(tte), "occurred", tte)) %>%
                # Pivot tte back to wide format
                pivot_wider(id_cols = c(id, event, first_event_tte),
                            names_from = "tte",
                            values_from = "occurred"),
              # Join on ID
              "id",
              # Individuals occur multiple times in both data frames
              relationship = "many-to-many") %>%
    # Keep only events and censoring
    filter((event != "censored" & occurred == 1 & tte <= time) |   # Occurence of event
           (event == "censored" & time <= first_event_tte)) %>%    # First event has not occurred yet
    # Outcome placeholder
    mutate(outcome = NA) %>%
    # Suppress warning that occurrs due to pivot_longer (see above annotation)
    suppressWarnings()
  
  # Change outcomes (important: outcome vector should be in same order as pstate)
  for(i in 1:length(outcomes)) dat_pdi %<>% mutate(outcome = if_else(event == outcomes[[i]], as.character(i), outcome))

  # Change available times to intended timepoints
  for(i in 1:length(timepoints)) dat_pdi %<>% mutate(time = if_else(time == tps[[i]], timepoints[[i]], time))
  
  # Finalise PDI data
  dat_pdi %<>%
    # Drop first event tte
    select(-first_event_tte) %>%
    # Set outcome to numeric
    mutate(outcome = as.numeric(outcome)) %>%
    # Rename pstate columns and drop pstate1 (this is the censoring state)
    set_colnames(c("time", paste0("p", 1:length(outcomes)), "id", "event", "occurred", "tte", "outcome")) %>%
    # Keep only relevant columns
    select(time, outcome, all_of(str_replace(vec_cols, "state", "")[1:length(vec_cols)])) %>%
    # Set to data frame (instead of tibble)
    as.data.frame()
  
  # Return PDI data
  return(dat_pdi)
}
