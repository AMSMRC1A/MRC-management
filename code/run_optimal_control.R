# working directory should be "MRC-management/code/R"

library(deSolve)
library(tidyverse)
library(pracma)
#library(reshape2)
library(data.table)
# libraries for parallelization
library(doParallel)
library(foreach)
# use all of your CPU cores except 2 for processing 
registerDoParallel(detectCores() - 2) # decrease for less aggressive core use

source("implementation/optimal_control_functions.R")


# implement across multiple optimal control parameters -------------------------
# define baseline parameters to change
test_params <- data.frame(
  control_type = c("unique", "uniform"), 
  m1 = 5E-4,
  m2 = 5E-4
)


# iterate over each parameter set in test_params
# EH Note: is there a more efficient way to do this?
# KD: maybe. You could create a full parameter set data table, then use mutate
#     to add the appropriate columns from oc_optim. But I think mutate doesn't
#     play well with multi-output functions. Maybe one of the apply functions?
#     Brandon might know more
test2 <- foreach(
  i = 1:nrow(test_params),
  .packages = c("deSolve", "tidyverse", "pracma")
) %dopar% {
  oc_optim(
    model = "ebola",
    change_params = test_params[i, 1:3]
  )
}
test_params$test_case <- 1:nrow(test_params)


# reformat for plotting
states <- lapply(
  1:nrow(test_params),
  function(i) { # browser();
    return(data.frame(
      test_case = i,
      reshape2::melt(test2[[i]]$trajectories, "time")
    ))
  }
)
states <- as.data.frame(do.call(rbind, states))
states <- left_join(test_params, states)
# reformat for easy plotting
states <- states %>%
  select(-test_case) %>%
  mutate(
    patch = substr(variable, 2, 2),
    variable = substr(variable, 1, 1)
  ) %>%
  mutate(
    plot_var = ifelse(control_type == "unique", paste("patch", patch), "uniform")
  )

# plot state variables and controls over time in each patch
states %>%
  mutate(
    variable = factor(variable, levels = c("S", "E", "I", "R", "H", "D", "v", "u"))
  ) %>%
  ggplot(aes(x = time, y = value, color = patch, linetype = control_type)) +
  geom_line(size = 1) +
  facet_wrap(vars(variable), scales = "free") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank()
  )


## TESTING WITH DAT.TABLE, NOT WORKING RIGHT NOW
# KD: This isn't working for me. I had to add the "data.table" library to get 
#     the setDT function. "Error: Cannot find symbol change_params"
# test3 <- setDT(change_params)[, oc_optim(model = "cholera", 
#                                          change_params = .SD), by = .I]

