# Function to quickly save object to prespecified path
saver <- function(object, file, prefix_path = path) save(object, file = paste0(prefix_path, file))
