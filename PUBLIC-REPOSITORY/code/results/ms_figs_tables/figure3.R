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

# load baseline cases
source("code/results/baseline_cases.R")
# load plotting utilities
source("code/results/plotting_utils.R")

#### FIGURE 3: Cholera infection trajectories + control effort -----------------
# plot figure
fig <- create_multipanel_ts_plot(
  model_name = "cholera",
  states = states %>%
    filter(# do not show before control
      time >= 0,
      # exclude no control
      v1_max != 0) %>% 
    mutate(patch_control = paste(patch, control_type, sep = "-")),
  patch_colors = ebola_cols,
  control_type_labs = ebola_gov_labs,
  I_labs = chol_I_labs,
  control_labs = chol_control_labs,
  patch_ltys = ebola_ltys
)

ggsave("figures/Cholera_trajectories_control.pdf",
       plot = fig,
       width = 6, height = 3, scale = 2)

