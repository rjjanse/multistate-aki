# Function to extarct and tabulate coefficients from mstate model fit object
extract_coefs <- function(.fit){
  # Extract coefficients
  dat_fit <- .fit[["coefficients"]] %>%
    # Set to tibble
    as_tibble() %>%
    # Rename column
    rename(coef = 1) %>%
    # Add indicators
    mutate(
      # Coefficient name
      var = names(.fit[["coefficients"]]),
      # Transition indicator
      trans = str_extract(var, "(?<=\\.)\\d+"),
      # Remove transition indicator from var
      var = str_replace(var, "\\.\\d+$", "")) %>%
    # Pivot to wide format
    pivot_wider(id_cols = var,
                names_from = trans,
                values_from = coef)
  
  # Return data
  return(dat_fit)
}
