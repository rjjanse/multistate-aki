# Function to draw calibration plots
plot_cal <- function(.list, 
                     outcomes){
  # Prepare plotting data
  dat_plot <- bind_rows(lapply(1:length(.list), \(x) .list[[x]] %<>% mutate(state = names(.list)[[x]]))) 
  
  # All states
  states <- sort(unique(dat_plot[["state"]]))
  
  # Nr. of states
  n_states <- length(states)
  
  # Rename states
  for(i in 1:length(outcomes)) dat_plot %<>% mutate(state = if_else(state == states[[i]], outcomes[[i]], state))

  # Define outcome labels
  outcomes_labs <- recode_values(outcomes,
                                 "censored" ~ "No event",
                                 c("kf", "macce", "ckd") ~ toupper(outcomes),
                                 "aki" ~ "Recurent AKI",
                                 "death" ~ "Death",
                                 c("macce_aki", "macce_ckd", "macce_kf") ~ toupper(str_replace(outcomes, "_", " + ")))
  
  # Set states to factor
  dat_plot %<>% mutate(state = factor(state,
                                      levels = outcomes,
                                      labels = outcomes_labs))
  
  # Palette
  palette <- c("#4B6BC0",  # 1: Censored
               "#648FFF",  # 2: CKD
               "#785EF0",  # 3: KF
               "#DC267F",  # 4: MACCE
               "#ED4440",  # 5: MACCE + CKD
               "#FE6100",  # 6: MACCE + KF
               "#FFB000"   # 7: Death
               )
  
  # Selected colours for different nr. of categories
  if(n_states == 5) pal <- palette[c(1, 3, 4, 6, 7)] 
  if(n_states == 7) pal <- palette
  
  # Create plot
  p <- ggplot(dat_plot,
              aes(x = pred,
                  y = obs,
                  colour = state,
                  fill = state)) +
    # Geometries
    geom_abline(colour = "#000000", alpha = 0.3) +
    geom_point(colour = "#000000", alpha = 0.3) +
    geom_smooth(method = "gam",
                formula = y ~ s(x, bs = "cs")) +
    # Scaling
    scale_x_continuous(breaks = seq(0, 1, 0.2),
                       name = "Predicted risk") +
    scale_y_continuous(breaks = seq(0, 1, 0.2),
                       name = "Observed risk") +
    scale_fill_manual(values = pal,
                      aesthetics = c("colour", "fill")) +
    # Transformations
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    facet_wrap(vars(state),
               axes = "all") +
    # Aesthetics
    theme(plot.background = element_rect(fill = "white", 
                                         colour = "white"),
          panel.background = element_rect(fill = "white", 
                                          colour = NA),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "#000000"),
          strip.text = element_text(face = "bold", 
                                    hjust = 0.5),
          strip.background = element_blank(),
          legend.position = "none")
  
  # Return plot
  return(p)
} 

