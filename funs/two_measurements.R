# Function to assess kidney disease based on two subequent measurements below a threshold
two_measurements <- function(.data, threshold = 15, visual_representation = FALSE){
  # Set date to time based on index date
  dat_tmp <- mutate(.data, time = as.numeric(lab_dt - discharge_dt))
  
  # Set empty time_stop
  time_stop <- NA
  
  ## Determine first sustained drop of 90 days
  # For loop for all rows in data
  for(i in 1:nrow(dat_tmp)){
    # If eGFR is above threshold, immediately go to the next iteration
    if(dat_tmp[i, "egfr"] > threshold) next
    
    # If eGFR is below threshold, store start time
    time_start <- dat_tmp[i, "time"][[1]]
    
    # Set iterator for next loop, starts at row after current
    j <- i + 1
    
    # Set stopping condition for for-loop
    stop <- FALSE
    
    # Iterate over all rows after the current one
    while(j <= nrow(dat_tmp)){
      # If eGFR is above threshold, immediately stop checking
      # Loop will continue to next iteration of i
      if(dat_tmp[j, "egfr"] > threshold) break
      
      # If time period is long enough, we can diagnose CKD
      if(dat_tmp[j, "time"] - time_start >= 90){
        # Store stoping time
        time_stop <- dat_tmp[j, "time"][[1]]
        
        # Add stopping indicator for upper for-loop
        stop <- TRUE
        
        # Stop while-loop
        break
      }
      
      # If the while-loop was not stopped, continue to next row
      j <- j + 1
    }
    
    # Stop for-loop if diagnosis was possible
    if(stop) break
  }
  
  # Store time stop in output (this is the diagnosis time)
  output <- time_stop
  
  # Visual representation if requested
  if(visual_representation){
    # Create indicator for diagnosis
    dat_tmp %<>% mutate(# Indicator for measurements that contributed to diagnosis
                        diag = if_else(time >= time_start & time <= time_stop, 1, 0,
                                       missing = 0),
                        # Set diag to factor
                        diag = factor(diag, levels = 0:1))
    
    # Create plot
    p <- ggplot(dat_tmp,
                aes(x = time,
                    y = egfr,
                    colour = diag,
                    shape = diag)) +
      # Geometries
      geom_point() +
      geom_vline(xintercept = time_stop,
                 linetype = "longdash",
                 colour = "black",
                 alpha = 0.3) +
      geom_hline(yintercept = threshold,
                 colour = "darkred") +
      # Colour scale
      scale_colour_manual(values = c("darkgreen", "darkblue")) +
      # Aesthetics
      theme_minimal()
    
    # Add plot to output
    output <- list(output, p)
  }
  
  # Return output
  return(output)
}
