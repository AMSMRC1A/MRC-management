# working directory should be "MRC-management/code/R"

library(deSolve)
library(tidyverse)
library(pracma)
# libraries for parallelization
library(doParallel)
library(foreach)
# update this if you want to be less aggressive with your core use
registerDoParallel(detectCores() - 2) 

source("implementation/optimal_control_functions.R")

# implement single optimal control analysis ------------------------------------
test <- oc_optim(model = "cholera")

# plot state variables and controls over time in each patch
test$trajectories %>%
  melt(c("time")) %>%
  mutate(
    patch = substr(variable, 2, 2),
    variable = substr(variable, 1, 1)
  ) %>%
  mutate(
    variable = factor(variable, levels = c("S", "I", "R", "W", "v", "u"))
  ) %>%
  ggplot(aes(x = time, y = value, color = patch)) +
  geom_line(size = 1) +
  facet_wrap(vars(variable), scales = "free") +
  #scale_color_manual(values = patch_colors2) +
  scale_linetype_manual(values = c("dotted", "solid")) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank()
  )


# implement across multiple optimal control parameters -------------------------
# define parameters to change
test_params <- data.frame(
  control_type = c("unique", "uniform")
)

# iterate over multiple
# EH Note: is there a more efficient way to do this?
test <- foreach(
  i = 1:nrow(test_params),
  .packages = c("deSolve", "tidyverse", "pracma"),
  # explicitly give 'foreach' the functions and data it needs
  .export = c(
    # optimal control functions
    "oc_optim", "run_no_optim","norm_oc",
    "setup_model", "set_old_variables","define_interp_fns", 
    "update_optimal_solution", "param_changer")
              # EH: KYLE CAN YOU MAKE SURE WE DON'T NEED THESE NOW?
              # "adj","calc_opt_v", "calc_opt_u",
              # "calc_j","eval_j_integrand",
              # "chol","params", "oc_params")
) %dopar% {
  oc_optim(
    model = "cholera",
    change_params = test_params[i, ]
  )
}


test <- setDT(change_params)[, oc_optim(model = "cholera", 
                                        change_params = .SD), by = .I]

