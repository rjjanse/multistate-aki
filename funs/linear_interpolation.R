# Function for linear interpolation of CKD
linear_interpolation <- function(.data, threshold = 15, visual_representation = FALSE){
  # Set date to time based on index date
  dat_tmp <- mutate(.data, time = as.numeric(lab_dt - discharge_dt))
  
  # Get max time
  max_time <- max(dat_tmp[["time"]])
  
  # Fit linear model on data
  fit <- lm(egfr ~ time,
            data = dat_tmp)
  
  # Get parameters
  y <- threshold
  a <- fit[["coefficients"]][["(Intercept)"]]
  b <- fit[["coefficients"]][["time"]]
  
  # Calculate time to threshold
  ttt <- (y - a) / b
  
  # Set to NA if it is beyond max time (we do not extrapolate) or before index date
  ttt <- if_else(ttt > max_time | ttt <= 0, NA, ttt)
  
  # Create output
  output <- ttt
  
  # Create visual representation if requested (requires ggplot2)
  if(visual_representation){
    # Create plot
    p <- ggplot(dat_tmp, 
           aes(x = time,
               y = egfr)) +
      # Geometries
      geom_vline(xintercept = 0,
                 linetype = "dashed",
                 colour = "black",
                 alpha = 0.3) +
      geom_vline(xintercept = ttt,
                 linetype = "longdash",
                 colour = "black",
                 alpha = 0.3) +
      geom_point() +
      geom_smooth(method = "lm",
                  colour = "darkorange",
                  formula = "y ~ x") +
      geom_hline(yintercept = threshold,
                 colour = "darkred") +
      # Aesthethics
      theme_minimal()
    
    # Create output list
    output <- list(ttt, p)
  }
  
  # Return output
  return(output)
}
