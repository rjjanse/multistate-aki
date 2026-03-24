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

# Set data path
if(.Platform[["OS.type"]] == "unix"){
  path <- "/Users/rjanse5/Networkshares/lab/lab_research/RES-Folder-UPOD/NOSTRADAMUS_SALTRO/E_ResearchData/2_ResearchData/CLEANED_for_Multistate_outcomes_project/20012026/dataframes/"
} else path <- "L:/lab_research/RES-Folder-UPOD/NOSTRADAMUS_SALTRO/E_ResearchData/2_ResearchData/CLEANED_for_Multistate_outcomes_project/20012026/dataframes/"

# Load functions
walk(list.files(here("funs")), \(x) source(paste0(here("funs"), "/", x)))

# Load spine data
load(paste0(path, "outcome_derivation.Rdata"))

# 1. Prepare data ----
## 1.1. Transition matrix ----
# Set-up transition matrix (we are only interested in baseline transitions
mat_trans <- transMat(x = list(c(2:5, 9), c(3:4, 6, 9), c(4, 7, 9), 8:9, 6:9, 9, 9, 9, c()),
                      names = c("Baseline", "AKI", "CKD", "KF", "MACCE", "MACCE_AKI", "MACCE_CKD", "MACCE_KF", "Death"))

# Transitions
# 1 = Baseline to AKI
# 2 = Baseline to CKD
# 3 = Baseline to KF
# 4 = Baseline to MACCE
# 5 = Baseline to death
# 6 = AKI to CKD
# 7 = AKI to KF
# 8 = AKI to MACCE + AKI
# 9 = AKI to death
# 10 = CKD to KF
# 11 = CKD to MACCE + CKD
# 12 = CKD to death
# 13 = KF to MACCE + KF
# 14 = KF to death
# 15 = MACCE to MACCE + AKI
# 16 = MACCE to MACCE + CKD
# 17 = MACCE to MACCE + KF
# 18 = MACCE to death
# 19 = MACCE + AKI to death
# 20 = MACCE + CKD to death
# 21 = MACCE + KF to death

## 1.2. Data transformation ----
# Get vector of covariates that are to be modelled
vec_covars <- setdiff(colnames(select(dat_spine, age, female, stage, los:ihe_resc)), 
                      c("com_aid", "med_aht", "med_ata", "med_gla", "com_hltx", "com_liver", "com_pad", "com_pckd", "ihe_resc"))

# Get vector of TTE variables
vec_tte <- c(NA, "aki_tte", "ckd_tte", "kf_tte", "macce_tte", "macce_aki_tte", "macce_ckd_tte", "macce_kf_tte", "death_tte")

# Get vector of status indicators
vec_stat <- c(NA, "aki", "ckd", "kf", "macce", "death")

# Transform data
dat_mstate <- msprep(time = vec_tte,
                     status = vec_stat,
                     data = dat_spine,
                     trans = mat_trans,
                     keep = vec_covars,
                     id = "id") %>%
  # For individuals with Tstop at Tstart, add one day to Tstop
  mutate(Tstop = if_else(Tstop == Tstart, Tstop + 1, Tstop),
         # Also recalculate time
         time = Tstop - Tstart) %>%
  # We get warnings about simultaneous state transitions: these are fine as they choose the smallest receiving
  # state and these are already ordered by importance
  suppressWarnings()

# Mutate drops transition matrix attribute, here we add it back
attr(dat_mstate, "trans") <- mat_trans

# Expand covariates
dat_mstate <- expand.covs(dat_mstate, vec_covars, append = TRUE)

# Event count
events(dat_mstate)[["Frequencies"]]

# 2. Fit model ----
# Model right hand side
model_rhs <- paste0(colnames(dat_mstate)[str_detect(colnames(dat_mstate), "\\.\\d")], collapse = " + ")

# Fit model on mstate data
fit <- coxph(as.formula(paste0("Surv(Tstart, Tstop, status) ~ strata(trans) + ", model_rhs)),
             data = dat_mstate, 
             method = "breslow",
             iter.max = 100)

# Save model fit
save(fit, file = paste0(path, "mstate_fit.Rdata"))

# 3. Get model predictions ----
# Get predictions from model
dat_prds <- prd_mstate(dat_spine, vec_covars, mat_trans)

# Save predictions
save(dat_prds, file = paste0(path, "predictions.Rdata"))











## Metrics for abstract
# Individuals
n_distinct(dat_spine$id)

# Median age
summary(dat_spine$age) # 63 [52-72]

# Proportion female
proportions(table(dat_spine$female)) # 41.4% = 58.6% male

# Median follow-up
summary(dat_spine$total_fu) # 2.2 [0.6-5.4]

## Prepare data for predictions
# Get predictions from model
dat_prds <- mstate_pred(dat_spine, vec_covars, mat_trans)

# Save predictions
save(dat_prds, file = paste0(path, "predictions.Rdata"))









# Calculate PDI
# Change data for PDI calculation
dat_pdi <- dat_prds %>%
  # Join outcome
  left_join(dat_spine %>%
              # Calculate first event
              mutate(event_dt = pmin(aki_dt, ckd_dt, kf_dt, death_dt, chd_dt, cvd_dt, na.rm = TRUE),
                     outcome = case_when(event_dt == aki_dt & aki == 1 ~ 1,
                                       event_dt == ckd_dt & ckd == 1 ~ 2,
                                       event_dt == kf_dt & kf == 1 ~ 3,
                                       event_dt == chd_dt &chd == 1 ~ 4,
                                       event_dt == cvd_dt & cvd == 1~ 5,
                                       event_dt == death_dt & death == 1 ~ 6,
                                       .default = NA)) %>%
              select(id, outcome),
            "id") %>%
  # Rename prediction columns
  rename(p1 = pstate2,   # Here, p1 corresponds to outcome 1 which is pstate2
         p2 = pstate3,
         p3 = pstate4,
         p4 = pstate5,
         p5 = pstate6,
         p6 = pstate7) %>%
  # Keep only predictions for outcomes (not censored)
  filter(!is.na(outcome))  %>% 
  # Keep last observation per individual
  arrange(id, desc(time)) %>%
  group_by(id) %>%
  slice(1L) %>%
  ungroup() %>%
  # Keep only relevant columns
  select(outcome, p1:p6)%>%
  # Explicitly set to data frame for pdiest function
  as.data.frame()

pdis <- pdiest(dat_pdi)

# Rescale function
rescale <- function(value, xmin = 1/6, xmax = 1, ymin = 0.5, ymax = 1) return(ymin + (ymax - ymin) / (xmax - xmin) * (value - xmin))

rescale(pdis) %>% as.data.frame() %>% mutate(outcome = c("overall", "aki", "ckd", "kf", "chd", "cvd", "death"))

# Calibration
pacman::p_load("calibmsm")
      
# Prepare data raw
dat_raw <- dat_spine %>%
  mutate(event_dt = pmin(aki_dt, ckd_dt, kf_dt, chd_dt, cvd_dt, death_dt, final_dt, na.rm = TRUE),
         event = pmax(aki, ckd, kf, chd, cvd, death, na.rm = TRUE),
         dtcens_s = if_else(event == 1, 0, 1),
         dtcens = as.numeric(event_dt - discharge_dt))
       
tp_prds <- filter(dat_prds, time == 365.25 * 10) %>% select(pstate1:pstate7)

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

fun_p <- function(.data){
  dat_tmp <- pivot_longer(.data,
                          pstate1:pstate7) %>%
    filter(name != "pstate1") %>%
    mutate(outcome = factor(name,
                            labels = c("Recurrent AKI",
                                       "CKD",
                                       "Kidney failure",
                                       "Acute coronary heart disease",
                                       "Acute cerebrovascular disease",
                                       "Death")))
  
  p <- ggplot(dat_tmp,
              aes(x = time,
                  y = value,
                  fill = outcome,
                  colour = outcome)) +
    geom_area(colour = NA, alpha = 0.8) +
    scale_x_continuous(breaks = (0:10) * 365.25,
                       labels = 0:10,
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2),
                       expand = c(0, 0),
                       labels = paste0(seq(0, 100, 20), "%")) +
    scale_fill_manual(values = c("#648FFF",
                                 "#785EF0",
                                 "#DC267F",
                                 "#FE6100",
                                 "#FFB000",
                                 "#FFD885")) +
    xlab("Time (years)") + ylab("Probability of outcome") +
    coord_cartesian(ylim = c(0, 1)) +
    theme(plot.background = element_rect(colour = "transparent", fill = "white"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.line.x.top = element_line(colour = "black"),
          axis.line.y.right = element_line(colour = "black"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "darkgrey"),
          legend.position = "bottom",
          legend.title = element_blank())
  
  return(p)
}

ids <- unique(dat_spine$id)

p1 <- fun_p(filter(dat_prds, id == ids[[1]]))
p2 <- fun_p(filter(dat_prds, id == ids[[2]]))  

fun_p(filter(dat_prds, id == ids[[3000]]))
#541 acute coronary disease
# 500 kidney failure

library(patchwork)

final <- wrap_plots(p1, p2, nrow = 1, axis_titles = "collect", guides = "collect") + plot_annotation(tag_levels = "A") & 
  theme(legend.position = "bottom",
        legend.title = element_blank())

final

ggsave("fig.png", final, path = "C:/Users/rjanse5/OneDrive - UMC Utrecht/Documenten/Projects/202508 - Multistate AKI - Denise/figures/",
       dpi = 1200, width = 8, height = 5)
