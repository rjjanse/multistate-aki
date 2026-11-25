## Functions to check risk plots for specific individuals, which helps to find examples
# Endpoint model
foo_em <- function(x) plot_rsk(filter(dat_prds_fixed_em, id == unique(dat_prds_fixed_em[["id"]])[[x]]), outcomes_em, timepoints)

# Endpoint + CKD model
foo_tm <- function(x) plot_rsk(filter(dat_prds_fixed_tm, id == unique(dat_prds_fixed_em[["id"]])[[x]]), outcomes_tm, timepoints)
