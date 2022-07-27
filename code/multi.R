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

# load model details to global environment
chol_mod_details <- setup_model("cholera")
ebol_mod_details <- setup_model("ebola")


#### FUNCTIONS -----------------------------------------------------------------
sens_analysis_setup <- function(change_params, 
                                multiplier, 
                                id_col = NA){
  # create data.frame with parameters to change (start with baseline)
  test_params <- data.frame(change_params)
  # expand to have one row per parameter to change
  test_params <- test_params[rep(1, times = length(change_params) + 1),]
  # change each parameter by x-fold once
  for(i in 1:length(change_params)){
    test_params[i,i] = multiplier*test_params[i,i]
  }
  # create id column if desired
  if(!any(is.na(id_col))){
    test_params$id = id_col
  }
  # add rows for unique and uniform
  test_params <- bind_rows(
    test_params %>% mutate(control_type = "unique"),
    test_params %>% mutate(control_type = "uniform")
  )
}

### SETUP PLOTTING -------------------------------------------------------------
## colors
patch_colors <- c("#1b9e77", "#d95f02")
# cost_colors6 <- c(rev(brewer.pal(3, "Reds")), rev(brewer.pal(3, "Blues")))
control_type_labs <- c("uniform", "non-uniform") # << change name of unique here
names(control_type_labs) <- c("uniform", "unique")
# labels
ebola_control_labs <- c("Vaccination effort", "Hospitalization effort")
names(ebola_control_labs) <- c("v", "u")

#### EBOLA ---------------------------------------------------------------------
# setup test parameters
# test_params <- sens_analysis_setup(change_params = ebol_mod_details$params[c("Cv1", "Cv2", "Cu1", "Cu2")], 
#                                    multiplier = 10, 
#                                    id_col = c("Cv1", "Cv2", "Cu1", "Cu2", "base"))

test_params <- sens_analysis_setup(change_params = chol_mod_details$params[c("C1", "C2", "D1", "D2")],
                                   multiplier = 10,
                                   id_col = c("C1", "C2", "D1", "D2", "base"))

# iterate over each parameter set in test_params
test2 <- foreach(
  i = 1:nrow(test_params),
  .packages = c("deSolve", "tidyverse", "pracma")
) %dopar% {
  oc_optim(
    model = "cholera",
    change_params = test_params[i, 1:6]
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
states <- left_join(test_params, states)


# reformat for easy plotting
states <- states %>%
  #select(-test_case) %>%
  mutate(
    patch = substr(variable, 2, 2),
    variable = substr(variable, 1, 1)
  ) %>%
  mutate(
    plot_var = ifelse(control_type == "unique", paste("patch", patch), "uniform")
  ) %>%
  mutate(control_type = ifelse(control_type == "unique", "non-uniform", "uniform"))
  



######## Ebola plots ########
ebol_titles <- c("(a) increase cost of vaccination in patch 1", 
                 "(b) increase cost of vaccination in patch 2",
                 "(c) increase cost of hospitalization in patch 1", 
                 "(d) increase cost of hospitalization in patch 2")
names(ebol_titles) <- c("Cv1", "Cv2", "Cu1", "Cu2")

# plot controls over time in each patch
p<- list()
for(i in 1:(nrow(test_params)/2 - 1)){
  p[[i]] <- states %>%
    filter(test_case %in% c(i, nrow(test_params)/2 + i), 
           variable %in% c("v", "u")) %>%
    # mutate(
    #   variable = factor(variable, levels = c("S", "E", "I", "R", "H", "D", "v", "u"))
    # ) %>%
    ggplot(aes(x = time, y = value, color = patch, linetype = control_type)) +
    geom_line(size = 1) +
    facet_grid(rows = vars(variable),
               labeller = labeller(variable = ebola_control_labs),
               switch = "y",
               scales = "free") +
    labs(color = "Patch:", 
         linetype = "Control type:",
         subtitle = ebol_titles[i]
         ) +
    scale_color_manual(values = patch_colors) +
    scale_linetype_manual(values = c("dotted", "solid")) +
    theme_minimal_grid(12) +
    theme(
      axis.title.y = element_blank(),
      legend.position = "bottom",
      panel.spacing = unit(c(0), "lines"),
      plot.margin = unit(c(0.1,0.1,0.1,2), "lines"),
      strip.background = element_blank(),
      strip.placement = "outside" 
    ) +
  panel_border()
}


leg <- get_legend(p[[1]])
p <- lapply(p, function(i){i + theme(legend.position = "none")})

plot_grid(plotlist = p, ncol = 4)


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

######## Cholera plots ########

## Create Cholera plots
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
  labs(y = "Sanitation", 
       subtitle = "",
       linetype = "Control type:",
       color = "Patch:") +
  scale_linetype_manual(values = c("dotted", "solid"))+
  theme_bw(12) +
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
  labs(y = "Vaccination", 
       subtitle = "Cholera: Effect of changing the patch-specific cost of sanitation") +
  scale_linetype_manual(values = c("dotted", "solid"))+
  theme_bw(12) +
  theme(
    legend.position = "none",
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
  labs(y = "Vaccination",
       subtitle = "Cholera: Effect of changing the patch-specific cost of vaccination") +
  scale_linetype_manual(values = c("dotted", "solid"))+
  theme_bw(12) +
  theme(
    legend.position = "none",
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
  labs(y = "Sanitation", 
       subtitle = "",
       linetype = "Control type:",
       color = "Patch:") +
  scale_linetype_manual(values = c("dotted", "solid"))+
  theme_bw(12) +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank()
  )

cholera_plots_Cv <- plot_grid(p2, p2a, ncol = 1)

cholera_plots_Cu <- plot_grid(p1a, p1, ncol = 1)

cholera_plots <- plot_grid(plot_grid(p2, p2a, ncol = 1), plot_grid(p1a, p1, ncol = 1), nrow = 1)

