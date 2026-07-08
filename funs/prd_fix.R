# Function to remove invalid predictions from multistate model
prd_fix <- function(.data){
  # Columns
  vec_cols <- setdiff(colnames(.data), c("time", "id"))
  
  # Create new data frame
  dat_prds_fixed <- .data
  
  ## Fix invalid predictions
  for(i in vec_cols){
    # Fix per column
    dat_prds_fixed %<>%
      # Rename target column
      set_colnames(str_replace(colnames(.), i, "target_pstate")) %>%
      # Arrange per person per time
      arrange(id, time) %>%
      # Group per individual
      group_by(id) %>%
      # Fix predictions
      mutate(# Get row number
        row = row_number(),
        # Register row number if prediction is invalid
        invalid_index = if_else(target_pstate < 0 | target_pstate > 1, row, Inf),
        # Get first occurrence
        invalid_index = min(invalid_index, na.rm = TRUE),
        # Set predictions to missing
        target_pstate = if_else(!is.na(invalid_index) & row >= invalid_index, NA, target_pstate)) %>%
      # Remove grouping structure
      ungroup() %>%
      # Name target pstate back
      set_colnames(str_replace(colnames(.), "target_pstate", i))
  }
  
  # Drop irrelevant variables
  dat_prds_fixed %<>% select(-c(row, invalid_index))
  
  # Return data
  return(dat_prds_fixed)
}
