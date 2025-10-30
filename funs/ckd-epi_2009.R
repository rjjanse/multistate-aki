# Function to calculate eGFR based on CKD-EPI 2009 creatinine formula
ckd_epi <- function(scr,
                    female,
                    age,
                    black = NULL,
                    scr_unit = "umol/l"){
  # If black is NULL, ethnicity correction factor (ECF) is not taken into account in the calculation
  if(is.null(black)) ecf <- 1 else ecf <- if_else(black == 1, 1.159, 1)
  
  # Define kappa
  k <- if_else(female == 1, 0.7, 0.9)
  
  # Define alpha
  a <- if_else(female == 1, -0.329, -0.411)
  
  # Convert serum creatinine to mg/dL if necessary
  if(scr_unit == "umol/l") cr <- scr / 88.42 else cr <- scr
  
  # Calculate eGFR
  egfr <- 141 * pmin(cr / k, 1) ^ a * pmax(cr / k, 1) ^ -1.209 * 0.993 ^ age * if_else(female == 1, 1.018, 1) * ecf
  
  # Return eGFR
  return(egfr)
}
