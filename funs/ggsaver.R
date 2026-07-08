# Function to quickly save figure to prespecified figure path
ggsaver <- function(plot, file, width, height, dpi = 1200, device = c("png", "pdf"), prefix_path = path_fig){
  # Save plot for each device
  for(i in device) ggsave(paste0(file, ".", i), plot, path = path_fig, width= width, height = height, dpi = dpi)
}
