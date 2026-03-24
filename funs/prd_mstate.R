# Create function to calculate multistate predicted risks per individual 
prd_mstate <- function(.data, predictors, transition_matrix){
  # Select all predictor columns
  dat_new <- select(.data, id, all_of(predictors))
  
  # Maximum transition
  max_trans <- max(mat_trans, na.rm = TRUE)
  
  # Create new column names
  new_cols <- sort(do.call("c", lapply(1:max_trans, function(x) paste0(predictors, ".", x))))
  
  # Duplicate new data for each transition
  dat_new <- do.call("rbind", replicate(max_trans, dat_new, simplify = FALSE)) %>%
    # Arrange on individual
    arrange(id) %>%
    # Group on individuals
    group_by(id) %>%
    # Compute strata and trans
    mutate(strata = 1:max_trans,
           trans = 1:max_trans) %>%
    # Ungroup again
    ungroup()
  
  # Add new (empty) columns
  for(i in new_cols) dat_new[[i]] <- NA
  
  # Set new column values to covariate value if stratum, else 0
  dat_new %<>%
    # Set column values
    # If the number in the column name equals the stratum, get that value from the original columns, else set to 0
    mutate(across(all_of(new_cols), function(x) x = ifelse(as.numeric(str_extract(cur_column(), "(?<=\\.)\\d{1,2}")) == trans, 
                                                           dat_new[[str_replace(cur_column(), "\\.\\d{1,2}", "")]], 0))) %>%
    # Remove unnecessary variables
    select(-all_of(predictors))
  
  # Get unique studynrs
  ids <- unique(dat_new[["id"]])
  
  # Total studynrs
  tot <- length(ids)
  
  # Get individual predictions
  all_probs <- do.call("rbind", lapply(1:tot, function(x){
    # Get ID
    id_iteration = ids[[x]]
    
    # Get relevant data
    dat_tmp <- filter(dat_new, id == id_iteration)
    
    # Prepare individual data
    dat_fit <- msfit(fit, newdata = dat_tmp, trans = mat_trans)
    
    # Get individual probabilities for each timepoint
    probs <- probtrans(dat_fit, 
                       predt = 0, 
                       method = "aalen", 
                       direction = "forward", 
                       variance = FALSE)[[1]] %>%
      # Add studynr
      mutate(id = id_iteration)
    
    # Iteration counter
    cat("\rIteration ", x, "/",tot, "       ")
    
    # Return probs
    return(probs)
  }))
  
  # Return all probabilities
  return(all_probs)
}
