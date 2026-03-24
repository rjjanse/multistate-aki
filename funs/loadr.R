# Shortcut function to load with a pre-specified path, named path
loadr <- function(file, prefix_path = path) load(paste0(prefix_path, file))
