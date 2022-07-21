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
source("implementation/optimal_control_functions.R")

## create multi-panel plot with control dynamics for a variety of cost params

test_params <- expand.grid(
  control_type = c("unique", "uniform"),
  Cv1 = c(0.01,0.1),
  Cv2 = c(0.01,0.1),
  Cu1 = c(0.1,1), 
  Cu2 = c(0.1,1)
)
# filter out multiple changes
test_params$rm_flag = with(test_params,ifelse(Cv2/Cv1 == 10 & Cu2/Cu1 == 10, 1, 0))
test_params$rm_flag = with(test_params,ifelse(Cv1/Cv2 == 10 & Cu1/Cu2 == 10, 1, rm_flag))
test_params$rm_flag = with(test_params,ifelse(Cv2/Cv1 == 10 & Cu1/Cu2 == 10, 1, rm_flag))
test_params$rm_flag = with(test_params,ifelse(Cv1/Cv2 == 10 & Cu2/Cu1 == 10, 1, rm_flag))
test_params$rm_flag = with(test_params, ifelse(Cv1 + Cv2 == 0.2, 1, rm_flag))
test_params$rm_flag = with(test_params, ifelse(Cu1 + Cu2 == 2, 1, rm_flag))
test_params <- test_params %>%
  filter(rm_flag == 0) %>%
  select(-rm_flag)


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
    change_params = test_params[i, 1:5]
  )
}
test_params$test_case <- 1:nrow(test_params)

# print number of iterations for each run
sapply(test2, function(i){i$n_iterations})

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
# initial conditions
ICs <- lapply(
  1:nrow(test_params),
  function(i) { # browser();
    return(data.frame(
      test_case = i,
      reshape2::melt(as.data.frame(test2[[i]]$uncontrolled), "time")
    ))
  }
)
ICs <- as.data.frame(do.call(rbind, ICs))
ICs$time = ICs$time - max(ICs$time)
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
    plot_var = ifelse(control_type == "unique", paste("patch", patch), "uniform")
  )
rm(ICs)


# plot state variables over time in each patch
states_plot <- states %>%
  mutate(
    variable = factor(variable, levels = c("S", "E", "I", "R", "H", "D", "v", "u"))
  ) %>%
  filter(variable %in% c("I"), time >= 0) %>%
  ggplot(aes(x = time, y = value, color = patch, linetype = control_type)) +
  geom_line(size = 1) +
  facet_grid(rows = vars(paste0(variable,patch)),
             cols = vars(paste0("Cv1: ", Cv1, ", Cv2: ", Cv2, 
                                "\nCu1: ", Cu1, ", Cu2: ",Cu2)),
             scales = "free") +
  scale_linetype_manual(values = c("dotted", "solid"))+
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank()
  )

# plot controls over time in each patch
p1 = states %>%
  mutate(
    variable = factor(variable, levels = c("S", "E", "I", "R", "H", "D", "v", "u"))
  ) %>%
  filter(variable %in% c("u"), time >= 0, Cv1 == .01, Cv2 == .01) %>%
  ggplot(aes(x = time, y = value, color = patch, linetype = control_type)) +
  geom_line(size = 1) +
  facet_grid(cols = vars(Cu1),
             rows = vars(Cu2),
             labeller = label_both,
             scales = "free") +
  labs(y = "hosptialization", 
       subtitle = "Ebola: change Cu") +
  scale_linetype_manual(values = c("dotted", "solid"))+
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank()
  )

p1a = states %>%
  mutate(
    variable = factor(variable, levels = c("S", "E", "I", "R", "H", "D", "v", "u"))
  ) %>%
  filter(variable %in% c("v"), time >= 0, Cv1 == .01, Cv2 == .01) %>%
  ggplot(aes(x = time, y = value, color = patch, linetype = control_type)) +
  geom_line(size = 1) +
  facet_grid(cols = vars(Cu1),
             rows = vars(Cu2),
             labeller = label_both,
             scales = "free") +
  labs(y = "vaccination", 
       subtitle = "Ebola: change Cu") +
  scale_linetype_manual(values = c("dotted", "solid"))+
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank()
  )


# plot controls over time in each patch
p2 =  states %>%
  mutate(
    variable = factor(variable, levels = c("S", "E", "I", "R", "H", "D", "v", "u"))
  ) %>%
  filter(variable %in% c("v"), time >= 0, Cu1 == 0.1, Cu2 == 0.1) %>%
  ggplot(aes(x = time, y = value, color = patch, linetype = control_type)) +
  geom_line(size = 1) +
  facet_grid(cols = vars(Cv1),
             rows = vars(Cv2),
             labeller = label_both,
             scales = "free") +
  labs(y = "vaccination", 
       subtitle = "Ebola: change Cv") +
  scale_linetype_manual(values = c("dotted", "solid"))+
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank()
  )
p2a =  states %>%
  mutate(
    variable = factor(variable, levels = c("S", "E", "I", "R", "H", "D", "v", "u"))
  ) %>%
  filter(variable %in% c("u"), time >= 0, Cu1 == 0.1, Cu2 == 0.1) %>%
  ggplot(aes(x = time, y = value, color = patch, linetype = control_type)) +
  geom_line(size = 1) +
  facet_grid(cols = vars(Cv1),
             rows = vars(Cv2),
             labeller = label_both,
             scales = "free") +
  labs(y = "hospitalization", 
       subtitle = "Ebola: change Cv") +
  scale_linetype_manual(values = c("dotted", "solid"))+
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank()
  )

ebola_plots <- plot_grid(plot_grid(p2, p2a, ncol = 1), plot_grid(p1a, p1, ncol = 1), nrow = 1)



######## SECOND EXAMPLE ########
test_params <- expand.grid(
  control_type = c("unique", "uniform"),
  u1_max = 1:5/10
)
test_params$u2_max = test_params$u1_max


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

# print number of iterations for each run
sapply(test2, function(i){i$n_iterations})

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
# initial conditions
ICs <- lapply(
  1:nrow(test_params),
  function(i) { # browser();
    return(data.frame(
      test_case = i,
      reshape2::melt(as.data.frame(test2[[i]]$uncontrolled), "time")
    ))
  }
)
ICs <- as.data.frame(do.call(rbind, ICs))
ICs$time = ICs$time - max(ICs$time)
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
    plot_var = ifelse(control_type == "unique", paste("patch", patch), "uniform")
  )
rm(ICs)


# plot state variables over time in each patch
states_plot2 <- states %>%
  mutate(
    variable = factor(variable, levels = c("S", "E", "I", "R", "H", "D", "v", "u"))
  ) %>%
  filter(variable %in% c("I"), time >= 0) %>%
  ggplot(aes(x = time, y = value, color = patch, linetype = control_type)) +
  geom_line(size = 1) +
  facet_grid(rows = vars(paste0(variable,patch)),
             cols = vars(u1_max),
             scales = "free") +
  scale_linetype_manual(values = c("dotted", "solid"))+
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank()
  )

######## RERUN FOR CHOLERA ##########

## create multi-panel plot with control dynamics for a variety of cost params

test_params <- expand.grid(
  control_type = c("unique", "uniform"),
  C1 = c(0.125,1.25),
  C2 = c(0.125,1.25),
  D1 = c(0.0125,0.125), 
  D2 = c(0.0125,0.125)
)
# filter out multiple changes
test_params$rm_flag = with(test_params,ifelse(C2/C1 == 10 & D2/D1 == 10, 1, 0))
test_params$rm_flag = with(test_params,ifelse(C1/C2 == 10 & D1/D2 == 10, 1, rm_flag))
test_params$rm_flag = with(test_params,ifelse(C2/C1 == 10 & D1/D2 == 10, 1, rm_flag))
test_params$rm_flag = with(test_params,ifelse(C1/C2 == 10 & D2/D1 == 10, 1, rm_flag))
test_params$rm_flag = with(test_params, ifelse(D1 + D2 == 0.25, 1, rm_flag))
test_params$rm_flag = with(test_params, ifelse(C1 + C2 == 2.5, 1, rm_flag))
test_params <- test_params %>%
  filter(rm_flag == 0) %>%
  select(-rm_flag)


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
    model = "cholera",
    change_params = test_params[i, 1:5]
  )
}
test_params$test_case <- 1:nrow(test_params)

# print number of iterations for each run
sapply(test2, function(i){i$n_iterations})

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
# initial conditions
ICs <- lapply(
  1:nrow(test_params),
  function(i) { # browser();
    return(data.frame(
      test_case = i,
      reshape2::melt(as.data.frame(test2[[i]]$uncontrolled), "time")
    ))
  }
)
ICs <- as.data.frame(do.call(rbind, ICs))
ICs$time = ICs$time - max(ICs$time)
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
    plot_var = ifelse(control_type == "unique", paste("patch", patch), "uniform")
  )
rm(ICs)


p1 = states %>%
  mutate(
    variable = factor(variable, levels = c("S", "E", "I", "R", "H", "D", "v", "u"))
  ) %>%
  filter(variable %in% c("u"), time >= 0, C1 == 0.125, C2 == 0.125) %>%
  ggplot(aes(x = time, y = value, color = patch, linetype = control_type)) +
  geom_line(size = 1) +
  facet_grid(cols = vars(D1),
             rows = vars(D2),
             labeller = label_both,
             scales = "free") +
  labs(y = "sanitation", 
       subtitle = "Cholera: change Cu") +
  scale_linetype_manual(values = c("dotted", "solid"))+
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank()
  )
p1a = states %>%
  mutate(
    variable = factor(variable, levels = c("S", "E", "I", "R", "H", "D", "v", "u"))
  ) %>%
  filter(variable %in% c("v"), time >= 0, C1 == 0.125, C2 == 0.125) %>%
  ggplot(aes(x = time, y = value, color = patch, linetype = control_type)) +
  geom_line(size = 1) +
  facet_grid(cols = vars(D1),
             rows = vars(D2),
             labeller = label_both,
             scales = "free") +
  labs(y = "vaccination", 
       subtitle = "Cholera: change Cu") +
  scale_linetype_manual(values = c("dotted", "solid"))+
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank()
  )


# plot controls over time in each patch
p2 =  states %>%
  mutate(
    variable = factor(variable, levels = c("S", "E", "I", "R", "H", "D", "v", "u"))
  ) %>%
  filter(variable %in% c("v"), time >= 0, D1 == 0.0125, D2 == 0.0125) %>%
  ggplot(aes(x = time, y = value, color = patch, linetype = control_type)) +
  geom_line(size = 1) +
  facet_grid(cols = vars(C1),
             rows = vars(C2),
             labeller = label_both,
             scales = "free") +
  labs(y = "vaccination", 
       subtitle = "Cholera: change Cv") +
  scale_linetype_manual(values = c("dotted", "solid"))+
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank()
  )
p2a =  states %>%
  mutate(
    variable = factor(variable, levels = c("S", "E", "I", "R", "H", "D", "v", "u"))
  ) %>%
  filter(variable %in% c("u"), time >= 0, D1 == 0.0125, D2 == 0.0125) %>%
  ggplot(aes(x = time, y = value, color = patch, linetype = control_type)) +
  geom_line(size = 1) +
  facet_grid(cols = vars(C1),
             rows = vars(C2),
             labeller = label_both,
             scales = "free") +
  labs(y = "sanitation", 
       subtitle = "Cholera: change Cv") +
  scale_linetype_manual(values = c("dotted", "solid"))+
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank()
  )

cholera_plots <- plot_grid(plot_grid(p2, p2a, ncol = 1), plot_grid(p1a, p1, ncol = 1), nrow = 1)

