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

#### RUN OC --------------------------------------------------------------------
chol_mod_details <- setup_model("cholera")
ebol_mod_details <- setup_model("ebola")

# ebola_m_test <- c(0, 0.005, .05)
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

# KD: the following makes the output from the above "nice"
#     could/should this be built into the oc_optim function?
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

### SETUP PLOTTING -------------------------------------------------------------
## colors
patch_colors <- c("#1b9e77", "#d95f02")
# cost_colors6 <- c(rev(brewer.pal(3, "Reds")), rev(brewer.pal(3, "Blues")))
control_type_labs <- c("uniform", "non-uniform") # << change name of unique here
names(control_type_labs) <- c("uniform", "unique")

#### TIME SERIES PLOT FUNCTIONS ------------------------------------------------

# plot_df should have columns for time, value, patch, control_type
plot_timeseries <- function(plot_df, patch_colors, facet_type, facet_labs,
                            lty_lab, y_lab, leg_pos) {
  # define linetype and size by control type
  ltys = c("longdash", "solid", "dotted")
  names(ltys) = c("none", "uniform", "unique")
  szs = c(0.6, 0.9, 0.9)
  names(szs) = c("none", "uniform", "unique")
  # filter to only those in plot_df
  ltys = ltys[names(ltys) %in% unique(plot_df$control_type)]
  szs = szs[names(szs) %in% unique(plot_df$control_type)]
  # plot
  p <- ggplot(
    data = plot_df,
    aes(
      x = time, y = value,
      color = patch, linetype = control_type, size = control_type
    )
  ) +
    geom_line() +
    scale_color_manual(name = "Patch:", values = patch_colors) +
    scale_linetype_manual(values = ltys, labels = lty_lab) +
    scale_size_manual(values = szs, labels = lty_lab) +
    scale_x_continuous( expand = expansion(mult = c(0,0.05))) +
    scale_y_continuous( expand = expansion(mult = c(0,0.05))) +
    labs(
      x = "Days since start of control",
      y = y_lab,
      linetype = "Control type:",
      size = "Control type:"
    ) +
    theme_minimal_grid() +
    theme(
      legend.position = leg_pos,
      legend.key.width = unit(1, "cm"),
      panel.grid.minor = element_blank(),
      strip.placement = "outside"
    ) +
    panel_border()
  if (facet_type == "states") {
    p <- p +
      facet_wrap(vars(variable),
        labeller = labeller(variable = facet_labs),
        nrow = 1,
        scales = "free", strip.position = "top"
      )
  } else if (facet_type == "Patch:") {
    p <- p +
      facet_wrap(vars(patch),
        labeller = labeller(patch = facet_labs),
        nrow = 2,
        scales = "free",
        strip.position = "top"
      )
  }
  return(p)
}

create_multipanel_ts_plot <- function(model_name, states, patch_colors,
                                      I_labs, control_labs, control_type_labs) {
  # plot infectious class
  p_Istates <- states %>%
    filter(
      model == model_name,
      variable == "I"
    ) %>%
    plot_timeseries(
      patch_colors = patch_colors,
      facet_labs = I_labs,
      facet_type = "Patch:",
      lty_lab = control_type_labs,
      y_lab = "",
      leg_pos = "none"
    )
  # plot controls
  p_controls <- states %>%
    filter(
      model == model_name,
      variable %in% c("v", "u")
    ) %>%
    plot_timeseries(
      patch_colors = patch_colors,
      facet_labs = control_labs,
      facet_type = "states",
      lty_lab = control_type_labs,
      y_lab = "",
      leg_pos = "bottom"
    )
  # put together
  p <- plot_grid(p_Istates, p_controls,
    rel_widths = c(0.3, 0.7),
    align = "h",
    axis = "b"
  )
  return(p)
}

#### FIGURE 3: Ebola infection trajectories + control effort -------------------
# repeat for ebola
# define labels
ebola_control_labs <- c("Vaccination effort", "Hospitalization effort")
names(ebola_control_labs) <- c("v", "u")
# ebola states
ebola_I_labs <- c("Infectives in Patch 1", "Infectives in Patch 2")
names(ebola_I_labs) <- 1:2

# plot figure
fig3 <- create_multipanel_ts_plot(
  model_name = "ebola",
  states = states %>% filter(# do not show before control
                             time >= 0,
                             # exclude no control
                             v1_max != 0),
  patch_colors = patch_colors,
  control_type_labs = control_type_labs,
  I_labs = ebola_I_labs,
  control_labs = ebola_control_labs
)
ggsave("../results/figures/Ebola_trajectories_control.pdf",
       plot = fig3,
       width = 6, height = 3, scale = 2)

## values for text:
# number vaccinated
# number hospitalized/sanitized
# number deaths
# j values
temp_df_unique <- states %>%
  filter(time >= 0) %>% 
  filter(model == "ebola") %>%
  filter(control_type == "unique") %>%
  select(-c(model, control_type,v1_max, v2_max, u1_max, u2_max, plot_var)) %>%
  mutate(compartment = paste0(variable,patch)) %>%
  select(-c(patch, variable)) %>%
  relocate(value, .after = compartment) %>% 
  unique() %>%
  pivot_wider(names_from = compartment)
  
j_vals_unique <- calc_j_ebola(temp_df_unique$time,
                              select(temp_df_unique,`S1`:`u2`),
                              ebol_mod_details$params)

temp_df_uniform <- states %>%
  filter(time >= 0) %>% 
  filter(model == "ebola") %>%
  filter(control_type == "uniform") %>%
  select(-c(model, control_type,v1_max, v2_max, u1_max, u2_max, plot_var)) %>%
  mutate(compartment = paste0(variable,patch)) %>%
  select(-c(patch, variable)) %>%
  relocate(value, .after = compartment) %>% 
  unique() %>%
  pivot_wider(names_from = compartment)

j_vals_uniform <- calc_j_ebola(temp_df_uniform$time,
                              select(temp_df_uniform,`S1`:`u2`),
                              ebol_mod_details$params)
  

# days hosp switches which higher
test <- states %>% 
  filter(# do not show before control
    time >= 0,
    # exclude no control
    v1_max != 0, 
    model == "ebola", 
    variable == "u") %>%
  dcast(time ~ paste0(control_type,patch), value.var = "value") %>%
  mutate(flag = ifelse(unique1 < uniform1, 1, 0)) %>%
  filter(flag ==1) %>%
  pull(time) %>%
  min()




#### FIGURE 5: Cholera infection trajectories + control effort -----------------
# define labels
chol_control_labs <- c("Vaccination effort", "Sanitation effort")
names(chol_control_labs) <- c("v", "u")
# cholera states
chol_I_labs <- c("Infectives in Patch 1", "Infectives in Patch 2")
names(chol_I_labs) <- 1:2

# plot figure
fig5 <- create_multipanel_ts_plot(
  model_name = "cholera",
  states = states %>%
    filter(# do not show before control
           time >= 0,
           # exclude no control
           v1_max != 0),
  patch_colors = patch_colors,
  control_type_labs = control_type_labs,
  I_labs = chol_I_labs,
  control_labs = chol_control_labs
)
ggsave("../results/figures/Cholera_trajectories_control.pdf",
       plot = fig5,
       width = 6, height = 3, scale = 2)

# Get data for tables
temp_df_unique <- states %>%
  filter(time >= 0) %>% 
  filter(model == "cholera") %>%
  filter(control_type == "unique") %>%
  select(-c(model, control_type,v1_max, v2_max, u1_max, u2_max, plot_var)) %>%
  mutate(compartment = paste0(variable,patch)) %>%
  select(-c(patch, variable)) %>%
  relocate(value, .after = compartment) %>% 
  unique() %>%
  pivot_wider(names_from = compartment)

j_vals_unique <- calc_j_cholera(temp_df_unique$time,
                              select(temp_df_unique,`S1`:`u2`),
                              chol_mod_details$params)

temp_df_uniform <- states %>%
  filter(time >= 0) %>% 
  filter(model == "cholera") %>%
  filter(control_type == "uniform") %>%
  select(-c(model, control_type,v1_max, v2_max, u1_max, u2_max, plot_var)) %>%
  mutate(compartment = paste0(variable,patch)) %>%
  select(-c(patch, variable)) %>%
  relocate(value, .after = compartment) %>% 
  unique() %>%
  pivot_wider(names_from = compartment)

j_vals_uniform <- calc_j_cholera(temp_df_uniform$time,
                               select(temp_df_uniform,`S1`:`u2`),
                               chol_mod_details$params)


#### FIGURE A1: Cholera compartments time series -------------------------------
chol_all_states_labs <- c("Susceptible", "Infected", "Recovered", "Water")
names(chol_all_states_labs) <- c("S", "I", "R", "W")

figA1 <- states %>%
  filter(
    model == "cholera",
    variable %in% c("S", "I", "R", "W"),
    #m1 == 0,
    #m2 == 0,
    control_type == "none"
  ) %>%
  mutate(variable = factor(variable, levels = c("S", "I", "R", "W", "u", "v"))) %>%
  plot_timeseries(
    patch_colors = patch_colors,
    facet_type = "states",
    facet_labs = chol_all_states_labs,
    lty_lab = control_type_labs,
    y_lab = "",
    leg_pos = "bottom"
  ) +
  #  We don't need to show a legend for linetype since there's only one type of control.
  guides(linetype = "none") + # KD: This command isn't working.
  scale_linetype_manual(values = "solid", labels = control_type_labs)

ggsave("../results/figures/Appendix_cholera_no control.pdf",
       plot = figA1,
       width = 6, height = 2, scale = 2)

#### COST PLOT FUNCTIONS -------------------------------------------------------

# relative costs
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

patch_labs <- c("Patch 1", "Patch 2", "Total")
names(patch_labs) <- c("1", "2", "t")


#### FIGURE 4: Ebola relative costs ---------------------------------------------
var_labs <- c("Vaccination", "Hospitalization", "Cases", "Total cost")
names(var_labs) <- c("vacc", "sani", "case", "to")

fig4 <- j_vals %>%
  select(-test_case) %>%
  dcast(variable + model ~ control_type) %>%
  mutate(
    type = ifelse(substr(variable, 1, 1) == "j", "cost",
                  ifelse(substr(variable, 1, 1) == "e", "epi", "res")
    ),
    patch = substr(variable, nchar(variable), nchar(variable)),
    variable_short = substr(variable, unlist(gregexpr("_", variable)) + 1, nchar(variable) - 1),
    rel_change = (unique / uniform)-1 # this treats "unique" as the before and "uniform" as the after
  ) %>%
  filter(
    model == "ebola",
    !(variable %in% c(paste0("j_", c("case1", "case2", "vacc1", "vacc2", "sani1", "sani2")), 
                      paste0("epi_", c("hosp1", "hosp2", "death1", "death2"))))
  ) %>%
  mutate(variable_short = factor(variable_short, levels = c("vacc", "sani", "case", "to")),
         text_pos = ifelse(rel_change < 0, 1, 0)) %>%
  ggplot(aes(x = variable_short, y = rel_change)) +
  geom_col(position = "dodge", color = "black") +
  geom_text(aes(label = paste0(round(rel_change*100,1), "%"), vjust = text_pos)) +
  geom_hline(yintercept = 0) +
  facet_grid(
    cols = vars(patch),
    labeller = labeller(patch = patch_labs),
    scales = "free",
    space = "free_x"
  ) +
  labs(
    x = "",
    y = "",
    title = "Ebola: Percent change from uniform to non-uniform policy"
  ) +
  # scale_fill_brewer(palette = "Greys", ) +
  scale_x_discrete(labels = var_labs) +
  scale_y_continuous(labels = percent) +
  theme_minimal(18) +
  theme(
    legend.position = "bottom",
    panel.border = element_rect(color = "lightgrey", fill = NA),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(0, "cm")
  )

ggsave("../results/figures/Ebola_relative_costs.pdf",
       plot = fig4,
       width = 6, height = 3, scale = 2)


## values for text
j_vals %>%
  select(-test_case) %>%
  dcast(variable + model ~ control_type) %>%
  filter(model == "ebola", 
         substr(variable, 1,3) == "epi")

#### FIGURE 6: Cholera relative costs ------------------------------------------
var_labs <- c("Vaccination", "Sanitation", "Cases", "Total cost")
names(var_labs) <- c("vacc", "sani", "case", "to")
fig6 <- j_vals %>%
  select(-test_case) %>%
  dcast(variable + model ~ control_type) %>%
  mutate(
    type = ifelse(substr(variable, 1, 1) == "j", "cost",
      ifelse(substr(variable, 1, 1) == "e", "epi", "res")
    ),
    patch = substr(variable, nchar(variable), nchar(variable)),
    variable_short = substr(variable, unlist(gregexpr("_", variable)) + 1, nchar(variable) - 1),
    rel_change = (unique/uniform)-1 # this treats "unique" as the before and "uniform" as the after
  ) %>%
  filter(
    model == "cholera",
    !(variable %in% paste0("j_", c("case1", "case2", "vacc1", "vacc2", "sani1", "sani2")))
  ) %>%
  mutate(variable_short = factor(variable_short, levels = c("vacc", "sani", "case", "to")),
         text_pos = ifelse(rel_change < 0, 1, 0)) %>%
  ggplot(aes(x = variable_short, y = rel_change)) +
  geom_col(position = "dodge", color = "black") +
  geom_text(aes(label = paste0(round(rel_change*100,1), "%"), vjust = text_pos)) +
  geom_hline(yintercept = 0) +
  facet_grid(
    cols = vars(patch),
    labeller = labeller(patch = patch_labs),
    scales = "free",
    space = "free_x"
  ) +
  labs(
    x = "",
    y = "",
    title = "Cholera: Percent change from uniform to non-uniform policy"
  ) +
  # scale_fill_brewer(palette = "Greys", ) +
  scale_x_discrete(labels = var_labs) +
  scale_y_continuous(labels = percent) +
  theme_minimal(18) +
  theme(
    legend.position = "bottom",
    panel.border = element_rect(color = "lightgrey", fill = NA),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(0, "cm")
  )
ggsave("../results/figures/Cholera_relative_costs.pdf",
       plot = fig6,
       width = 6, height = 3, scale = 2)



#### FIGURE 7: Ebola - Effect of varying patch cost of Ebola control------------
## Build test set for movement parameters
## create multi-panel plot with control dynamics for a variety of cost params

## colors
patch_colors <- c("#1b9e77", "#d95f02")
# cost_colors6 <- c(rev(brewer.pal(3, "Reds")), rev(brewer.pal(3, "Blues")))
control_type_labs <- c("uniform", "non-uniform") # << change name of unique here
names(control_type_labs) <- c("uniform", "unique")
# labels
ebola_control_labs <- c("Vaccination effort", "Hospitalization effort")
names(ebola_control_labs) <- c("v", "u")
cholera_control_labs <- c("Vaccination effort", "Sanitation effort")
names(cholera_control_labs) <- c("v", "u")

# setup test parameters
test_params <- sens_analysis_setup(change_params = ebol_mod_details$params[c("Cv1", "Cv2", "Cu1", "Cu2")],
                                   multiplier = 10,
                                   id_col = c("Cv1", "Cv2", "Cu1", "Cu2", "base"))
test_params <- test_params[-c(2,7),1:6]

# iterate over each parameter set in test_params
test2 <- foreach(
  i = 1:nrow(test_params),
  .packages = c("deSolve", "tidyverse", "pracma")
) %dopar% {
  oc_optim(
    model = "ebola",
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


ebol_titles <- c("(a) increase cost of vaccination in patch 1", 
                 #"(b) increase cost of vaccination in patch 2",
                 "(c) increase cost of hospitalization in patch 1", 
                 "(d) increase cost of hospitalization in patch 2")
names(ebol_titles) <- c("Cv1",
                        #"Cv2",
                        "Cu1",
                        "Cu2")

# plot controls over time in each patch
p<- list()
for(i in 1:(nrow(test_params)/2 - 1)){
  p[[i]] <- states %>%
    filter(test_case %in% c(i, nrow(test_params)/2 + i), 
           variable %in% c("v", "u")) %>%
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

Ebola_vary_cost_plots <- plot_grid(plotlist = p, ncol = 4)
Ebola_vary_cost_plots

ggsave("../results/figures/Ebola_vary_cost.pdf",
       Ebola_vary_cost_plots,
       width = 6, height = 3, scale = 2)

#### FIGURE 8: Cholera - Effect of varying patch cost of control----------------
# setup test parameters
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
  mutate(
    patch = substr(variable, 2, 2),
    variable = substr(variable, 1, 1)
  ) %>%
  mutate(
    plot_var = ifelse(control_type == "unique", paste("patch", patch), "uniform")
  ) %>%
  mutate(control_type = ifelse(control_type == "unique", "non-uniform", "uniform"))

chol_titles <- c("(a) increase cost of vaccination in patch 1", 
                 "(b) increase cost of vaccination in patch 2",
                 "(c) increase cost of sanitation in patch 1", 
                 "(d) increase cost of sanitation in patch 2")
names(chol_titles) <- c("Cv1", "Cv2", "Cu1", "Cu2")

# plot controls over time in each patch
p<- list()
for(i in 1:(nrow(test_params)/2 - 1)){
  p[[i]] <- states %>%
    filter(test_case %in% c(i, nrow(test_params)/2 + i), 
           variable %in% c("v", "u")) %>%
    ggplot(aes(x = time, y = value, color = patch, linetype = control_type)) +
    geom_line(size = 1) +
    facet_grid(rows = vars(variable),
               labeller = labeller(variable = cholera_control_labs),
               switch = "y",
               scales = "free") +
    labs(color = "Patch:", 
         linetype = "Control type:",
         subtitle = chol_titles[i]
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


Cholera_vary_cost_plots <- plot_grid(plotlist = p, ncol = 4)

ggsave("../results/figures/Cholera_vary_cost.pdf",
       Cholera_vary_cost_plots,
       width = 6, height = 3, scale = 2)