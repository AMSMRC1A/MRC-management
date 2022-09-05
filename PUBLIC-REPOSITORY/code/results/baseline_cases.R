# required packages
library(ggplot2)
library(deSolve)
library(reshape2)
library(tidyverse)
library(cowplot)
library(pracma)
library(RColorBrewer)
library(scales)

# to paralellize
library(doParallel)
library(foreach)
registerDoParallel(detectCores() - 2) # update this if you want to use more cores

# load optimal control files
source("code/optimal_control_functions.R")

#### RUN BASELINE OPTIMAL CONTROL CASES ----------------------------------------
chol_mod_details <- setup_model("cholera")
ebol_mod_details <- setup_model("ebola")

#### setup parameter sets to test
test_params <- expand.grid(
  control_type = c("unique", "uniform"),
  model = c("cholera", "ebola")
)
## add no control scenario
# add min/max control bounds to test_params
for (i in c("v1_max", "v2_max", "u1_max", "u2_max")) {
  test_params[i] <- ifelse(test_params$model == "cholera",
                           unlist(chol_mod_details$params[i]),
                           unlist(ebol_mod_details$params[i])
  )
}
# only use no movement scenario for now
test_params <- bind_rows(
  test_params,
  expand.grid(
    control_type = "uniform", # name uniform for now so code will run, fix this later
    model = c("cholera", "ebola"),
    v1_max = 0, v2_max = 0,
    u1_max = 0, u2_max = 0
  )
)
test_params <- as.data.frame(test_params)


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
    model = test_params[i, 2],
    change_params = test_params[i, c(1, 3:5)]
  )
}

# change test_params control_type to none where necessary -- can remove when update code (see comment above)
test_params$control_type <- ifelse(test_params$v1_max == 0,
                                   "none", as.character(test_params$control_type)
)

test_params$test_case <- 1:nrow(test_params)

# print number of iterations for each run
sapply(test2, function(i){i$n_iterations})



#### MANIPULATE OUTPUT FOR PLOTTING --------------------------------------------
## states
states <- lapply(
  1:nrow(test_params),
  function(i) { # browser();
    return(data.frame(
      test_case = i,
      melt(test2[[i]]$trajectories, "time")
    ))
  }
)
states <- as.data.frame(do.call(rbind, states))
# initial conditions
ICs <- lapply(
  1:nrow(test_params),
  function(i) { # browser();
    return(data.frame(
      test_case = i,
      melt(as.data.frame(test2[[i]]$uncontrolled), "time")
    ))
  }
)
ICs <- as.data.frame(do.call(rbind, ICs))
ICs = ICs %>%
  group_by(test_case) %>%
  mutate(time = time - max(time))
states <- bind_rows(states, ICs)
states <- left_join(test_params, states)
# reformat for easy plotting
states <- states %>%
  select(-test_case) %>%
  mutate(
    patch = substr(variable, 2, 2),
    variable = substr(variable, 1, 1)
  ) %>%
  mutate(
    plot_var = ifelse(control_type == "unique", paste("patch", patch), control_type)
  )
rm(ICs)


## relative costs
j_vals <- lapply(
  1:nrow(test_params),
  function(i) {
    j <- test2[[i]][["j"]]
    j <- c(j, j_tot = sum(j[substr(names(j), 1, 1) == "j"]))
    return(data.frame(test_case = i, melt(j)))
  }
)
j_vals <- as.data.frame(do.call(rbind, j_vals))
j_vals <- left_join(test_params, j_vals, by = "test_case") %>%
  rename(variable = L1)



