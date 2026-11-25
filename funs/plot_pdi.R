# Function to plot PDI
plot_pdi <- function(.data, outcomes, c1 = "#785EF0", c2 = "#FE6100"){
  # Calculate lower bound of PDI
  lb <- 1 / (n_distinct(.data[["outcome"]]) - 1)
  
  # Define outcome labels
  outcomes_labs <- case_match(outcomes,
                              "censored" ~ "No event",
                              c("kf", "macce", "ckd") ~ toupper(outcomes),
                              "aki" ~ "Recurent AKI",
                              "death" ~ "Death",
                              c("macce_aki", "macce_ckd", "macce_kf") ~ toupper(str_replace(outcomes, "_", " + ")))
  
  # Prettify data for plotting
  dat_plot <- mutate(.data,
                     # Set outcome to factor based on order of outcomes vector
                     outcome = factor(outcome, 
                                      levels = c("overall", outcomes),
                                      labels = c("Overall", outcomes_labs)),
                     # Set time to years
                     time = paste0("Year ", round(time / 365.25)),
                     # Set time to factor
                     time = factor(time,
                                   levels = paste0("Year ", 1:10)))
    
  # Create plot for PDI
  p <- ggplot(dat_plot,
              aes(x = pdi,
                  y = outcome)) +
    # Geometries
    geom_vline(xintercept = lb,
               linetype = "dashed",
               colour = "black",
               alpha = 0.3) +
    geom_point(alpha = 0.75,
               colour = c1,
               shape = "circle") +
    # Scales
    scale_x_continuous(name = "PDI",
                       limits = c(lb, 1),
                       breaks = round(seq(lb, 1, length.out = 6), digits = 2)) +
    scale_y_discrete(limits = rev) +
    # Transformations
    facet_wrap(vars(time),
               axes = "all",
               axis.labels = "margins") +
    # Aesthetics
    theme(plot.background = element_rect(colour = "white",
                                         fill = "white"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(),
          axis.title.y = element_blank(),
          strip.text = element_text(face = "bold"),
          strip.background = element_blank())

  # Return plot
  return(p)
}
