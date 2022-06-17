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

#### RUN OC --------------------------------------------------------------------
chol_mod_details <- setup_model("cholera")
ebol_mod_details <- setup_model("ebola")

#ebola_m_test <- c(0, 0.005, .05)
test_params <- expand.grid(
  control_type = c("unique", "uniform"), 
  m1 = c(0, .005),
  m2 = c(0, .005), 
  model = c("cholera", "ebola")
)

## add no control scenario
# add min/max control bounds to test_params
test_params$v1_max <- ifelse(test_params$model == "cholera", 
                             chol_mod_details$params["v1_max"], 
                             ebol_mod_details$params["v1_max"])
test_params$v2_max <- ifelse(test_params$model == "cholera", 
                             chol_mod_details$params["v2_max"], 
                             ebol_mod_details$params["v2_max"])
test_params$u1_max <- ifelse(test_params$model == "cholera", 
                             chol_mod_details$params["u1_max"], 
                             ebol_mod_details$params["u1_max"])
test_params$u2_max <- ifelse(test_params$model == "cholera", 
                             chol_mod_details$params["u2_max"], 
                             ebol_mod_details$params["u2_max"])
# only use no movement scenario for now 
test_params <- bind_rows(test_params, 
                         data.frame(control_type = "none", 
                         ))


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
states <- lapply(1:nrow(test_params), 
                 function(i){#browser();
                   return(data.frame(test_case = i, 
                                     melt(test2[[i]]$trajectories, "time")))
                 })
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

### SETUP PLOTTING -------------------------------------------------------------
## colors
patch_colors <- c("#1b9e77", "#d95f02")
#cost_colors6 <- c(rev(brewer.pal(3, "Reds")), rev(brewer.pal(3, "Blues")))


#### TIME SERIES PLOT FUNCTIONS ------------------------------------------------

# plot_df should have columns for time, value, patch, control_type
plot_timeseries <- function(plot_df, patch_colors, facet_type, facet_labs,
                            y_lab, leg_pos){
  p <-   ggplot(data = plot_df, 
                aes(x = time, y = value, 
                    color = patch, linetype = control_type)) +
    geom_line(size = 1) +
    scale_color_manual(values = patch_colors) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    labs(x = "time since start of control (days)", 
         y = y_lab, 
         linetype = "control type") +
    theme_minimal() +
    theme(
      legend.position = leg_pos, 
      legend.key.width = unit(1, "cm"),
      panel.grid.minor = element_blank(),
      strip.placement = "outside"
    )
  if(facet_type == "states"){
    p <- p + 
      facet_wrap(vars(variable), 
               labeller = labeller(variable = facet_labs),
               nrow = 1,
               scales = "free", strip.position = "left")
  }
  else if(facet_type == "patch"){
    p <- p + 
      facet_wrap(vars(patch), 
                 labeller = labeller(patch = facet_labs), 
                 nrow = 2,
                 strip.position = "left")
  }
  return(p)
}

create_multipanel_ts_plot <- function(model_name, states, patch_cholors,
                                 I_labs, control_labs){
  # plot infectious class
  p_Istates <- states %>%
    filter(model == model_name, 
           m1 == 0, 
           m2 == 0,
           variable == "I") %>%
    plot_timeseries(patch_colors = patch_colors, 
                    facet_labs = I_labs, 
                    facet_type = "patch",
                    y_lab = "", 
                    leg_pos = "none")
  # plot controls
  p_controls <- states %>%
    filter(model == model_name, 
           m1 == 0, 
           m2 == 0,
           variable %in% c("v", "u")) %>%
    plot_timeseries(patch_colors = patch_colors, 
                    facet_labs = control_labs, 
                    facet_type = "states",
                    y_lab = "", 
                    leg_pos = "bottom")
  # put together
  p <- plot_grid(p_Istates, p_controls, 
                    rel_widths = c(0.3,0.7), 
                    align = "h", 
                    axis = "b")
  return(p)
}

#### FIGURE 3 ------------------------------------------------------------------

fig3 <- create_multipanel_ts_plot(model_name = "cholera", 
                                  states = states, 
                                  patch_cholors = patch_cholors, 
                                  I_labs = chol_I_labs, 
                                  control_labs = chol_control_labs)
#ggsave("results/figures/Figure3.pdf", width = 6, height = 3, scale = 1.5)


#### FIGURE 4 ------------------------------------------------------------------
#repeat for ebola

#### COST PLOT FUNCTIONS -------------------------------------------------------


# relative costs
j_vals <- lapply(1:nrow(test_params), 
                 function(i) {j = test2[[i]][["j"]];
                 j = c(j, j_tot = sum(j[substr(names(j),1,1) == "j"]));
  return(data.frame(test_case = i, melt(j)))
})
j_vals <- as.data.frame(do.call(rbind, j_vals))
j_vals <- left_join(test_params, j_vals, by = "test_case") %>%
  rename(variable = L1)

patch_labs <- c("patch 1", "patch 2", "total")
names(patch_labs) <- c('1',"2","t")

var_labs <- c("vaccination", "sanitation", "cases", "total cost")
names(var_labs) <- c("vacc", "sani", "case", 'to')

j_vals %>%
  select(-test_case) %>%
  dcast(variable + m1 + m2 + model ~ control_type) %>%
  mutate(type = ifelse(substr(variable, 1,1) == "j", "cost", 
                       ifelse(substr(variable, 1,1)== "e", "epi", "res")),
         patch =substr(variable, nchar(variable), nchar(variable)),
         variable_short = substr(variable, unlist(gregexpr("_", variable)) + 1, nchar(variable) - 1),
         rel_change = log(uniform/unique)) %>%
  filter(m1 == 0, 
         m2 == 0, 
         model == "cholera",
         !(variable %in% paste0("j_", c("case1", "case2", "vacc1", "vacc2", "sani1", "sani2")))) %>%
  mutate(variable_short = factor(variable_short, levels = c("vacc", "sani", "case", "to"))) %>%
  ggplot(aes(x = variable_short, y = rel_change)) +
  geom_col(position = "dodge", color = "black") +
  geom_hline(yintercept = 0) +
  facet_grid(cols = vars(patch), 
             labeller = labeller(patch = patch_labs),
             scales = "free",
             space = "free_x") +
  labs( x = "", 
        y = "% change (log scale)\nlog(unique/uniform)", 
        title = "Cholera model") +
  #scale_fill_brewer(palette = "Greys", ) + 
  scale_x_discrete(labels = var_labs) +
  scale_y_continuous(labels = percent) +
  theme_minimal() +
  theme(legend.position = "bottom", 
        panel.border = element_rect(color = "lightgrey", fill = NA), 
        panel.grid.major.x = element_blank(), 
        panel.spacing = unit(0, "cm"))
#ggsave("results/figures/Figure5.pdf", width = 6, height = 3, scale = 1.5)
  


