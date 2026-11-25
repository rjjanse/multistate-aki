# Function to plot risks from multistate model for individual patient
plot_rsk <- function(.data, 
                     outcomes, 
                     timepoints,
                     type = c("area", "line"),
                     show_censored = FALSE,
                     verbose_legend = FALSE){
  # Match type
  type <- match.arg(type)
  
  
  # Define verbose legend
  outcomes_verbose <- c("censored", "ckd", "kf", "macce", "macce_ckd", "macce_kf", "death")
  
  # Number of outcomes
  if(verbose_legend) n_outcomes <-  length(outcomes_verbose) else n_outcomes <- length(outcomes)
                                                                        
  # Palette
  palette <- c("#4B6BC0",  # 1: Censored
               "#648FFF",  # 2: CKD
               "#785EF0",  # 3: KF
               "#DC267F",  # 4: MACCE
               "#ED4440",  # 5: MACCE + CKD
               "#FE6100",  # 6: MACCE + KF
               "#FFB000"   # 7: Death
  )
  
  ## Subset palette
  # If censored should be included
  if(show_censored){
    # Based on number of outcomes
    if(n_outcomes == 5) pal <- palette[c(1, 3, 4, 6, 7)] 
    if(n_outcomes == 7) pal <- palette
  }
  
  # If censored should not be included
  else { 
    # Based on number of states
    if(n_outcomes == 5) pal <- palette[c(3, 4, 6, 7)] 
    if(n_outcomes == 7) pal <- palette[2:7]
  }
  
  # Prepare plotting data
  dat_plt <- .data %>%
    # Set pstate colnames to outcomes (outcomes and pstate should be in the same order)
    set_colnames(c("time", outcomes, "id")) %>%
    # Pivot data to long format
    pivot_longer(cols = all_of(outcomes),
                 names_to = "outcome",
                 values_to = "prob")
  
  # Drop censored if requested 
  if(!show_censored){
    # From data
    dat_plt %<>% filter(outcome != "censored")
    
    # From outcome vectors
    outcomes <- setdiff(outcomes, "censored")
    outcomes_verbose <- setdiff(outcomes_verbose, "censored")
  }
  
  # Define outcome labels
  outcomes_labs <- recode_values(outcomes,
                                 "censored" ~ "No event",
                                 c("kf", "macce", "ckd") ~ toupper(outcomes),
                                 "aki" ~ "Recurent AKI",
                                 "death" ~ "Death",
                                 c("macce_aki", "macce_ckd", "macce_kf") ~ toupper(str_replace(outcomes, "_", " + ")))
  
  # Verbose outcome labels
  outcomes_labs_verbose <- recode_values(outcomes_verbose,
                                         "censored" ~ "No event",
                                         c("kf", "macce", "ckd") ~ toupper(outcomes_verbose),
                                         "aki" ~ "Recurent AKI",
                                         "death" ~ "Death",
                                         c("macce_aki", "macce_ckd", "macce_kf") ~ toupper(str_replace(outcomes_verbose, "_", " + ")))
  
  # Set outcome to factor for verbose legend
  if(verbose_legend) dat_plt %<>% mutate(outcome = factor(outcome, 
                                                          levels = outcomes_verbose,
                                                          labels = outcomes_labs_verbose))
  
  # Set outcome to factor for normal legend
  else dat_plt %<>% mutate(outcome = factor(outcome, 
                                            levels = outcomes,
                                            labels = outcomes_labs))
  
  # Draw area plot
  p <- ggplot(dat_plt,
              aes(x = time,
                  y = prob,
                  fill = outcome)) +
    # Geometries
    geom_area(alpha = 0.8,
              colour = if_else(show_censored, "black", NA),
              show.legend = TRUE) +
    # Scaling
    scale_fill_manual(values = pal,
                      drop = FALSE) +
    scale_x_continuous(name = "Time (years)",
                       breaks = c(0, timepoints),
                       labels = 0:10,
                       expand = c(0, 0)) +
    scale_y_continuous(name = "Probability of event",
                       breaks = seq(0, 1, 0.2),
                       labels = paste0(seq(0, 100, 20), "%"),
                       expand = c(0, 0)) +
    # Transformations
    coord_cartesian(ylim = c(0, 1)) +
    # Aesthetics
    theme(plot.background = element_rect(fill = "white",
                                         colour = NA),
          panel.background = element_blank(),
          panel.border = element_rect(),
          panel.grid.major = element_line(colour = "lightgrey"),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          legend.position = "bottom")
  
  # Overwrite with lineplot if requested
  if(type == "line"){
    # Draw plot
    p <- ggplot(dat_plt,
                aes(x = time,
                    y = prob,
                    colour = outcome,
                    linetype = outcome)) +
      # Geometries
      geom_line() +
      # Scaling
      scale_colour_manual(values = pal) +
      scale_x_continuous(name = "Time (years)",
                         breaks = c(0, timepoints),
                         labels = 0:10,
                         expand = c(0, 0)) +
      scale_y_continuous(name = "Probability of event",
                         breaks = seq(0, 1, 0.2),
                         labels = paste0(seq(0, 100, 20), "%"),
                         expand = c(0, 0)) +
      # Transformations
      coord_cartesian(ylim = c(0, 1)) +
      # Aesthetics
      theme(plot.background = element_rect(fill = "white",
                                           colour = NA),
            panel.background = element_blank(),
            panel.border = element_rect(),
            panel.grid.major = element_line(colour = "lightgrey"),
            panel.grid.minor = element_blank(),
            legend.title = element_blank(),
            legend.position = "bottom")
  }
  
  # Return plot
  return(p)
}
