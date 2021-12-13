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
test_params <- expand.grid(
  control_type = c("unique", "uniform"),
  beta_I1 = c(2.64E-6,2.64E-5)
)

# iterate over multiple
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
    model = "cholera",
    change_params = test_params[i, ]
  )
}
test_params$test_case <- 1:nrow(test_params)

# KD: the following makes the output from the above "nice"
#     could/should this be built into the oc_optim function?

states <- lapply(1:nrow(test_params), function(i) {
  return(data.frame(test_case = i, test2[[i]]$trajectories))
})
states <- as.data.frame(do.call(rbind, states))
states <- left_join(test_params, states)
j_vals <- lapply(1:nrow(test_params), function(i) {
  return(data.frame(test_case = i, test2[[i]][["j"]]))
})
j_vals <- as.data.frame(do.call(rbind, j_vals))
j_vals <- left_join(test_params, j_vals)
j_vals$j <- apply(j_vals[, c("j_case1", "j_case2", "j_vacc1", "j_vacc2")], 1, sum)


# KD: This isn't working for me. I had to add the "data.table" library to get 
#     the setDT function. "Error: Cannot find symbol change_params"
test3 <- setDT(change_params)[, oc_optim(model = "cholera", 
                                        change_params = .SD), by = .I]

