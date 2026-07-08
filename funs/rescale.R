# Linear rescaling function
rescale <- function(value, xmin, xmax = 1, ymin = 0.5, ymax = 1) return(ymin + (ymax - ymin) / (xmax - xmin) * (value - xmin))
