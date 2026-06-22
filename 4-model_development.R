#----------------------------------------------------------#
# Predicting outcomes after AKI using multistate models
# Code for model development
# Roemer J. Janse - Last updated on 2026-02-12
#----------------------------------------------------------#

# 0. Set-up ----
# Load packages
pacman::p_load("dplyr",          # Data wrangling
               "magrittr",       # Efficient pipelines
               "tidyr",          # Data cleaning
               "stringr",        # Working with strings
               "purrr",          # Functional programming
               "here",           # Local paths
               "survival",       # Survival models
               "mstate",         # Multistate models
               "conflicted"      # Resolve function conflicts
)

# Resolve function conflicts
conflicts_prefer(dplyr::filter) # Between stats and dplyr

# Load functions
walk(list.files(here("funs")), \(x) source(paste0(here("funs"), "/", x)))

# Set data path
path <- psp(unix = "/Users/rjanse5/Networkshares/Data/lab_research/RES-Folder-UPOD/NOSTRADAMUS_SALTRO/E_ResearchData/2_ResearchData/CLEANED_for_Multistate_outcomes_project/20012026/dataframes/",
            windows = "L:/lab_research/RES-Folder-UPOD/NOSTRADAMUS_SALTRO/E_ResearchData/2_ResearchData/CLEANED_for_Multistate_outcomes_project/20012026/dataframes/")

## Load spine dataframes
loadr("outcome_derivation_em.Rdata")   # Endpoint model
loadr("outcome_derivation_tm.Rdata")   # Trajectory models

# 1. Endpoint model ----
## 1.1. Prepare data ----
### 1.1.1. Transition matrix ----
# Set-up transition matrix 
mat_trans <- transMat(x = list(c(2:3, 5),     # From baseline
                               4:5,           # From KF
                               4:5,           # From MACCE
                               5,             # From MACCE + KF
                               c()),          # From death           
                      names =  c("Baseline",  # State 1
                                 "KF",        # State 2
                                 "MACCE",     # State 3
                                 "MACCE_KF",  # State 4
                                 "Death"))    # State 5

### 1.1.2. Data transformation ----
# Get vector of covariates that are to be modelled
vec_covars <- setdiff(colnames(select(dat_spine_em, age, female, stage, los:ihe_resc)), 
                      c("com_aid", "med_aht", "med_ata", "med_gla", "com_hltx", "com_liver", "com_pad", "com_pckd", "ihe_resc"))

# Get vector of TTE variables
vec_tte <- c(NA, "kf_tte", "macce_tte", "macce_kf_tte", "death_tte")

# Get vector of status indicators
vec_stat <- c(NA, "kf", "macce", "macce_kf", "death")

# Transform data
dat_mstate <- msprep(time = vec_tte,
                     status = vec_stat,
                     data = dat_spine_em,
                     trans = mat_trans,
                     keep = vec_covars,
                     id = "id") %>%
  # For individuals with Tstop at Tstart, add one day to Tstop
  mutate(Tstop = if_else(Tstop == Tstart, Tstop + 1, Tstop),
         # Also recalculate time
         time = Tstop - Tstart) 

# Mutate drops transition matrix attribute, here we add it back
attr(dat_mstate, "trans") <- mat_trans

# Expand covariates
dat_mstate <- expand.covs(dat_mstate, vec_covars, append = TRUE)

# Event count
events(dat_mstate)[["Frequencies"]]

## 1.2. Fit model ----
# Model right hand side
model_rhs <- paste0(colnames(dat_mstate)[str_detect(colnames(dat_mstate), "\\.\\d")], collapse = " + ")

# Fit model on mstate data
fit <- coxph(as.formula(paste0("Surv(Tstart, Tstop, status) ~ strata(trans) + ", model_rhs)),
             data = dat_mstate, 
             method = "breslow",
             iter.max = 1e4)

# Save model fit
save(fit, file = paste0(path, "mstate_fit_em.Rdata"))

## 1.3. Get model predictions ----
# Load model
loadr("mstate_fit_em.Rdata")

# Get predictions from model
dat_prds <- prd_mstate(dat_spine_em, vec_covars, mat_trans)

# Save predictions
save(dat_prds, file = paste0(path, "predictions_em.Rdata"))

# 2. Transition model without AKI ----
## 2.1. Prepare data ----
### 2.1.1. Transition matrix ----
# Set-up transition matrix 
mat_trans <- transMat(x = list(c(2, 4, 7),    # From baseline
                               c(3, 5, 7),    # From CKD
                               c(6, 7),       # From KF
                               c(5, 7),       # From MACCE
                               6:7,           # From MACCE + CKD
                               7,             # From MACCE + KF
                               c()),          # From death           
                      names =  c("Baseline",  # State 1
                                 "CKD",       # State 2
                                 "KF",        # State 3
                                 "MACCE",     # State 4
                                 "MACCE_CKD", # State 5
                                 "MACCE_KF",  # State 6
                                 "Death"))    # State 7

### 2.1.2. Data transformation ----
# Get vector of covariates that are to be modelled
vec_covars <- setdiff(colnames(select(dat_spine_tm, age, female, stage, los:ihe_resc)), 
                      c("com_aid", "med_aht", "med_ata", "med_gla", "com_hltx", "com_liver", "com_pad", "com_pckd", "ihe_resc"))

# Get vector of TTE variables
vec_tte <- c(NA, "ckd_tte", "kf_tte", "macce_tte", "macce_ckd_tte", "macce_kf_tte", "death_tte")

# Get vector of status indicators
vec_stat <- c(NA, "ckd", "kf", "macce", "macce_ckd", "macce_kf", "death")

# Transform data
dat_mstate <- msprep(time = vec_tte,
                     status = vec_stat,
                     data = dat_spine_tm,
                     trans = mat_trans,
                     keep = vec_covars,
                     id = "id") %>%
  # For individuals with Tstop at Tstart, add one day to Tstop
  mutate(Tstop = if_else(Tstop == Tstart, Tstop + 1, Tstop),
         # Also recalculate time
         time = Tstop - Tstart) 

# Mutate drops transition matrix attribute, here we add it back
attr(dat_mstate, "trans") <- mat_trans

# Expand covariates
dat_mstate <- expand.covs(dat_mstate, vec_covars, append = TRUE)

# Event count
events(dat_mstate)[["Frequencies"]]

## 2.2. Fit model ----
# Model right hand side
model_rhs <- paste0(colnames(dat_mstate)[str_detect(colnames(dat_mstate), "\\.\\d")], collapse = " + ")

# Fit model on mstate data
fit <- coxph(as.formula(paste0("Surv(Tstart, Tstop, status) ~ strata(trans) + ", model_rhs)),
             data = dat_mstate, 
             method = "breslow",
             iter.max = 100)

# Save model fit
save(fit, file = paste0(path, "mstate_fit_tm.Rdata"))

## 2.3. Get model predictions ----
# Load model
loadr("mstate_fit_tm.Rdata")

# Get predictions from model
dat_prds <- prd_mstate(dat_spine_tm, vec_covars, mat_trans)

# Save predictions
saver(dat_prds, "predictions_tm.Rdata")






















library(calibmsm)
loadr("predictions_em.Rdata")

# Columns
vec_cols <- setdiff(colnames(dat_prds), c("time", "id"))

# Create new data frame
dat_prds_fixed <- dat_prds

## Fix invalid predictions
for(i in vec_cols){
  # Fix per column
  dat_prds_fixed %<>%
    # Rename target column
    set_colnames(str_replace(colnames(.), i, "target_pstate")) %>%
    # Arrange per person per time
    arrange(id, time) %>%
    # Group per individual
    group_by(id) %>%
    # Fix predictions
    mutate(# Get row number
           row = row_number(),
           # Register row number if prediction is invalid
           invalid_index = if_else(target_pstate < 0 | target_pstate > 1, row, NA),
           # Get first occurrence
           invalid_index = min(invalid_index, na.rm = TRUE),
           # Set predictions to missing
           target_pstate = if_else(!is.na(invalid_index) & row >= invalid_index, NA, target_pstate)) %>%
    # Remove grouping structure
    ungroup() %>%
    # Name target pstate back
    set_colnames(str_replace(colnames(.), "target_pstate", i))
}


# Change data for PDI calculation
dat_pdi <- dat_prds_fixed %>%
  # Join outcome
  left_join(dat_spine_em %>%
              # Select relevant columns
              select(id, kf, macce, macce_kf, death) %>%
              # Pivot to long format 
              pivot_longer(c(kf, macce, macce_kf, death),
                           values_to = "occurred",
                           names_to = "event") %>%
              # Keep only events
              filter(occurred == 1),
            # Join on ID
            "id",
            # Individuals occur multiple times in both data frames
            relationship = "many-to-many") %>%
  # Change outcomes
  mutate(outcome = case_match(event,
                              "kf" ~ 1,
                              "macce" ~ 2,
                              "macce_kf" ~ 3,
                              "death" ~ 4)) %>%
  # Rename columns with predictions
  rename(p1 = pstate2,   # Outcome 1 (KF) is stored in pstate 2
         p2 = pstate3,
         p3 = pstate4,
         p4 = pstate5) %>%
  # Keep only individuals with an event
  filter(!is.na(outcome)) %>%
  # Arrange for grouping
  arrange(id, time) %>%
  # Group per individual
  group_by(id) %>%
  # Keep last observation per individual (10 year prediction)
  # Use filter to keep multiple rows for individuals with multiple outcomes
  filter(time == last(time)) %>%
  # Remove grouping structure
  ungroup() %>%
  # Keep only relevant columns
  select(outcome, all_of(str_replace(vec_cols, "state", "")[1:length(vec_cols)-1])) %>%
  # Set to data frame (instead of tibble)
  as.data.frame()

    
pdis <- pdiest(dat_pdi)

# Rescale function
rescale <- function(value, xmin = 1/6, xmax = 1, ymin = 0.5, ymax = 1) return(ymin + (ymax - ymin) / (xmax - xmin) * (value - xmin))

rescale(pdis) %>% as.data.frame() %>% mutate(pdi = pdis[, 1], outcome = c("overall", "kf", "macce", "macce_kf", "death"))

# Calibration
pacman::p_load("calibmsm")
      
# Prepare data raw
dat_raw <- dat_spine_em %>%
  mutate(cens = 1 - pmax(kf, macce, macce_kf, death, na.rm = TRUE)) %>%
  pivot_longer(c(kf, macce, macce_kf, death, cens),
               names_to = "occurred",
               values_to = "event") %>%
  mutate(event_dt = case_match(occurred,
                               "kf" ~ kf_dt,
                               "macce" ~ macce_dt,
                               "macce_kf" ~ macce_kf_dt,
                               "death" ~ death_dt,
                               "cens" ~ final_dt)) %>%
  # Keep only events
  filter(event == 1) %>%
  mutate(dtcens_s = if_else(occurred == "cens", 1, 0),
         dtcens = as.numeric(event_dt - discharge_dt)) %>%
         #occurred = factor(occurred,
         #                  levels = c("macce_kf", "macce", "kf", "death", "cens"))) %>% # order of least to most frequent
  arrange(id, event_dt) %>%
  group_by(id) %>% 
  slice(1L) %>%
  ungroup()

tp_prds <- filter(dat_prds, time == 365.25 * 10) %>% select(vec_cols)

# Calibration plots                                  
p_cal <- calib_msm(data_ms = dat_mstate, 
                   data_raw = dat_raw, 
                   j = 1, 
                   s = 0, 
                   t = 3652,
                   tp_pred = tp_prds,
                   calib_type = "mlr",
                   curve_type = "loess",
                   w_covs = vec_covars)

plot(p_cal)


dat_plot <- bind_rows(p_cal[["plotdata"]]) %>%
  mutate(state = c(rep("Censored", 2201), 
                   rep("Kidney failure", 2201), 
                   rep("MACCE", 2201),
                   rep("MACCE + Kidney Failure", 2201),
                   rep("Death", 2201))) 


# Create plot
p <- ggplot(dat_plot, aes(x = pred, y = obs, colour = state, fill = state)) +
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
  scale_fill_manual(values = c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#FFB000"),
                    aesthetics = c("colour", "fill")) +
  # Transformations
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  facet_wrap(vars(state),
             axes = "all") +
  # Aesthetics
  theme(plot.background = element_rect(fill = "transparent", colour = "transparent"),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "#000000"),
        strip.text = element_text(face = "bold", hjust = 0.5),
        strip.background = element_blank(),
        legend.position = "none")

p

ggsave("calibration_em.png",
       plot = p,
       path = "/Users/rjanse5/Library/CloudStorage/OneDrive-UMCUtrecht/Documenten/Projects/202508 - Multistate AKI - Denise/Figures/",
       width = 9,
       height = 6,
       dpi = 1200)














## For abstract
library(patchwork)
library(ggplot2)

# Function to plot
# Individual predictions
fun_p <- function(.data, state){
  # Get states that need to be dropped (starting state & before)
  drop <- paste0("pstate", 1:state)

  # Derive plotting data
  dat_tmp <- pivot_longer(.data,
                          pstate1:pstate5) %>%
    filter(!(name %in% drop)) %>%
    mutate(outcome = factor(name,
                            levels = c("pstate2", "pstate3", "pstate4", "pstate5"),
                            labels = c("Kidney failure",
                                       "MACCE",
                                       "MACCE + Kidney failure",
                                       "Death")))
  
  # First break
  first_break <- if_else(min(dat_tmp[["time"]]) > 0, min(dat_tmp[["time"]]), 0)
  
  # Breaks
  breaks <- (0:10) * 365.25
  
  # Truncated breaks
  t_breaks <- unique(c(0, first_break, breaks[breaks > first_break]))
  
  # Labels
  t_labels <- round(t_breaks / 365.25, 1)
  
  # Create plot
  p <- ggplot(dat_tmp,
              aes(x = time,
                  y = value,
                  fill = outcome,
                  colour = outcome)) +
    # Geometries
    geom_area(colour = NA, 
              alpha = 0.8,
              show.legend = TRUE) +
    # Scales
    scale_x_continuous(breaks = t_breaks,
                       labels = t_labels,
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2),
                       expand = c(0, 0),
                       labels = paste0(seq(0, 100, 20), "%")) +
    scale_fill_manual(values = c("#785EF0",
                                 "#DC267F",
                                 "#FE6100",
                                 "#FFB000"),
                      drop = FALSE) +
    # Labels & guides
    xlab("Time (years)") + 
    ylab("Probability of outcome") +
    guides(fill = guide_legend(title = NULL,
                               ncol = 1,
                               position = "right")) +
    # Transformations
    coord_cartesian(xlim = c(0, 365.25 * 10),
                    ylim = c(0, 1)) +
    # Aesthetics
    theme(plot.background = element_rect(colour = "transparent", fill = "white"),
          plot.title = element_text(hjust = 0.5,
                                    face = "italic"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.line.x.top = element_line(colour = "black"),
          axis.line.y.right = element_line(colour = "black"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "darkgrey"))
  
  return(p)
}

# Get three nice individuals
ids_abstract <- c(11917530550, # KF & MACCE, no death
                  555729182,    # KF & Death, no MACCE
                  225297902)     # Only death

# Keep only relevant data
dat_abstract <- dat_spine_em %>%
  # Filter IDs
  filter(id %in% ids_abstract)

# Get predictions from model from baseline
dat_prds <- prd_mstate(dat_abstract, vec_covars, mat_trans)

# Individual data frames
dat_abstract_id1 <- filter(dat_abstract, id == ids_abstract[[1]])
dat_abstract_id2 <- filter(dat_abstract, id == ids_abstract[[2]])
dat_abstract_id3 <- filter(dat_abstract, id == ids_abstract[[3]])

### Individual 1
## Get predictions for first individual from moment of kidney failure
# Moment of kidney failure (state 2)
t_kf <- extract2(dat_abstract_id1, "kf_tte")

# Moment of MACCE_KF (state 4)
t_macce <- extract2(dat_abstract_id1, "macce_kf_tte")

# Predictions for this individual from kidney failure onwards
dat_prds_id1_kf <- prd_mstate(dat_abstract_id1, vec_covars, mat_trans, t = t_kf, state = 2)

# Predictions for this individual from MACCE onwards
dat_prds_id1_macce <- prd_mstate(dat_abstract_id1, vec_covars, mat_trans, t = t_macce, state = 4)

# Store data
dat_id1 <- list(
  "1" = filter(dat_prds, id == ids_abstract[[1]]),
  "2" = dat_prds_id1_kf,
  "4" = dat_prds_id1_macce
)

# Make plots
p_id1_baseline <- fun_p(dat_id1[[1]], 1) + ggtitle("From baseline")
p_id1_kf <- fun_p(dat_id1[[2]], 2) + ggtitle("After kidney failure")
p_id1_macce_kf <- fun_p(dat_id1[[3]], 4) + ggtitle("After MACCE + kidney failure")

### Individual 2
## Get predictions for first individual from moment of kidney failure
# Moment of kidney failure (state 2)
t_kf <- extract2(dat_abstract_id2, "kf_tte")

# Predictions for this individual from kidney failure onwards
dat_prds_id2_kf <- prd_mstate(dat_abstract_id2, vec_covars, mat_trans, t = t_kf, state = 2) 

# Store data
dat_id2 <- list(
  "1" = filter(dat_prds, id == ids_abstract[[2]]),
  "2" = dat_prds_id2_kf
)

# Store plots
p_id2_baseline <- fun_p(dat_id2[[1]], 1) + ggtitle("From baseline")
p_id2_kf <- fun_p(dat_id2[[2]], 2) + ggtitle("After kidney failure")

# Add all plots together
p <- wrap_plots(p_id1_baseline, p_id1_kf, p_id1_macce_kf,
                p_id2_baseline, p_id2_kf, plot_spacer(), 
                fun_p(filter(dat_prds, id == ids_abstract[[3]]), state = 0) + ggtitle("From baseline"), plot_spacer(), plot_spacer(),
                ncol = 3,
                nrow = 3,
                byrow = FALSE,
                guides = "collect") &
  # Aesthetics
  theme(legend.position = "left",
        plot.background = element_rect(colour = "transparent",
                                       fill = "transparent"),
        panel.background = element_rect(colour = NA,
                                        fill = "white"),
        legend.background = element_rect(colour = "black",
                                         fill = "white"))

# Save figure
ggsave("preds_em.png", p, path = "/Users/rjanse5/OneDrive - UMC Utrecht/Documenten/Projects/202508 - Multistate AKI - Denise/figures/",
       dpi = 1200, width = 15, height = 12)




# EPPerMed figure
# Add all plots together
p_eppermed <- wrap_plots(p_id1_baseline, p_id1_kf, p_id1_macce_kf,
                nrow = 1,
                guides = "collect") &
  guides(fill = guide_legend(nrow = 1,
                             title = NULL)) &
  # Aesthetics
  theme(legend.position = "bottom")

# Save figure
ggsave("preds_eppermed.png", 
       p_eppermed, 
       path = "/Users/rjanse5/OneDrive - UMC Utrecht/Documenten/Grants/2026/CARMEN2026/full proposal/",
       dpi = 1200, width = 10, height = 4)

