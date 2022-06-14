# required packages
library(ggplot2)
library(deSolve)
library(reshape2)
library(tidyverse)
library(cowplot)
library(pracma)
library(RColorBrewer)

# to paralellize
library(doParallel)
library(foreach)
registerDoParallel(detectCores() - 2) # update this if you want to use more cores

# load optimal control files
source("implementation/optimal_control_functions.R")

#### OUTLINE -------------------------------------------------------------------
# group comments
# show infectious classes, show no control case somewhere
# relative change plot: include only cases, deaths, total cost
# Figure out why ebola no movement scenario has cases in patch 2

# plan
# Figures 1/2: model diagrams
# Figure 3: cholera infectious classes no control, unique, uniform (one panel per patch)
#           cholera control over time unique/uniform (one panel per patch) - both controls 
#           (check not worth showing single control)
# Figure 4: repeat for ebola
# Figure 5: relative change (cases, deaths, total cost)  (one panel per model)
# Figure 6: relative change + movement (no movement, symmetric movement, asymmetric movement)
#           play around with faceting
#           decide names of costs
#
# supplementary
# all classes  no control 
# single control cases
#
# theme
# theme_minimal (cowplot)
# colors: colorbrewer - patches
# linetype: no control, unique, uniform

#### RUN OC --------------------------------------------------------------------
#ebola_m_test <- c(0, 0.005, .05)
test_params <- expand.grid(
  control_type = c("unique", "uniform"), 
  m1 = c(0, .05),
  m2 = c(0, .05), 
  model = c("cholera", "ebola")
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
    model = test_params[i,4],
    change_params = test_params[i, 1:3]
  )
}
test_params$test_case <- 1:nrow(test_params)

# KD: the following makes the output from the above "nice"
#     could/should this be built into the oc_optim function?
states <- lapply(1:nrow(test_params), function(i) {
  return(data.frame(test_case = i, test2[[i]]$trajectories[,c("v1", "v2", "u1", "u2")]))
})
states <- as.data.frame(do.call(rbind, states))
states <- left_join(test_params, states)
j_vals <- lapply(1:nrow(test_params), function(i) {
  return(data.frame(test_case = i, test2[[i]][["j"]]))
})
j_vals <- as.data.frame(do.call(rbind, j_vals))
j_vals <- left_join(test_params, j_vals)
j_vals$j <- apply(j_vals[, c("j_case1", "j_case2", "j_vacc1", "j_vacc2")], 1, sum)

### PLOTTING -------------------------------------------------------------------
patch_colors2 <- c("#DE2D26", "#3182BD")
cost_colors6 <- c(rev(brewer.pal(3, "Reds")), rev(brewer.pal(3, "Blues")))

# time series
states %>%
  group_by(test_case) %>%
  mutate(time = 1:n()) %>%
  ungroup() %>%
  select(-test_case) %>%
  reshape2::melt(c("control_type", "m1", "m2", "model", "time")) %>%
  mutate(
    patch = substr(variable, 2, 2),
    variable = substr(variable, 1, 1)
  ) %>%
  mutate(
    plot_var = ifelse(control_type == "unique", paste("patch", patch), "uniform")#,
    #variable = factor(variable, levels = level_codes)
  ) %>%
  ggplot(aes(x = time, y = value, color = patch, linetype = control_type)) +
  geom_line(size = 1) +
  facet_grid(cols = vars(paste(m1, m2)), rows = vars(paste(model, variable)), scales = "free") +
  scale_color_manual(values = patch_colors2) +
  scale_linetype_manual(values = c("dotted", "solid")) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank()
    )

# relative costs
j_vals %>%
  unique() %>%
  select(-test_case) %>%
  melt(c("control_type","m1", "m2", "model")) %>% 
  dcast(variable + m1 + m2 + model ~ control_type) %>%
  mutate(rel_change = (unique - uniform)/uniform) %>%
  ggplot(aes(x = variable, y = rel_change)) +
  geom_col() +
  facet_grid(cols = vars(paste(m1, m2)), rows = vars(paste(model)), scales = "free") +
  labs(y = "% change", x = "") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


