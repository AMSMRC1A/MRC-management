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


#### FIGURE A1: Cholera compartments time series -------------------------------
chol_all_states_labs <- c("Susceptible", "Infected", "Recovered", "Water")
names(chol_all_states_labs) <- c("S", "I", "R", "W")
ltys = c("solid", "solid")
names(ltys) = c("1-none", "2-none")
cols = patch_colors
names(cols) = c("1-none", "2-none")

figA1 <- states %>%
  filter(
    model == "cholera",
    variable %in% c("S", "I", "R", "W"),
    #m1 == 0,
    #m2 == 0,
    control_type == "none"
  ) %>%
  mutate(variable = factor(variable, levels = c("S", "I", "R", "W", "u", "v")), 
         patch_control = paste(patch, control_type, sep = "-")) %>%
  plot_timeseries(
    patch_colors = cols,
    facet_type = "states",
    facet_labs = chol_all_states_labs,
    lty_lab = control_type_labs,
    y_lab = "",
    leg_pos = "bottom", 
    patch_ltys = ltys
  )

# fix legend
figA1 <- figA1 + 
  geom_vline(xintercept = 0, linetype = "dotted")+ 
  scale_color_manual(values = patch_colors, 
                     labels = c("Patch 1", "Patch 2"), 
                     name = "") +
  scale_linetype_manual(values = c("solid", "solid"), 
                        labels = c("Patch 1", "Patch 2"), 
                        name = "")

ggsave("figures/appendix_cholera_no control.pdf",
       plot = figA1,
       width = 6, height = 2, scale = 2)

