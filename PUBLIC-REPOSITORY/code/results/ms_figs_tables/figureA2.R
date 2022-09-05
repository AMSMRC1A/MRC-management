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

#### FIGURE A2: Ebola compartments time series -------------------------------
ebol_all_states_labs <- c("Susceptible",
                          "Exposed",
                          "Infected",
                          "Hospitalized",
                          "Dead",
                          "Recovered")
names(ebol_all_states_labs) <- c("S", "E", "I", "H", "D", "R")

figA2 <- states %>%
  filter(
    model == "ebola",
    variable %in% c("S", "E", "I", "H", "D", "R"),
    #m1 == 0,
    #m2 == 0,
    control_type == "none"
  ) %>%
  mutate(variable = factor(variable, levels = c("S", "E", "I", "H", "D", "R", "u", "v")),
         patch_control = paste(patch, control_type, sep = "-")) %>%
  plot_timeseries(
    patch_colors = cols,
    facet_type = "states",
    facet_labs = ebol_all_states_labs,
    lty_lab = control_type_labs,
    y_lab = "",
    leg_pos = "bottom", 
    patch_ltys = ltys
  )

# fix legend, add line to indicate start of control
figA2 <- figA2 + 
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_manual(values = patch_colors, 
                     labels = c("Patch 1", "Patch 2"), 
                     name = "") +
  scale_linetype_manual(values = c("solid", "solid"), 
                        labels = c("Patch 1", "Patch 2"), 
                        name = "")



ggsave("figures/appendix_ebola_no control.pdf",
       plot = figA2,
       width = 6, height = 1.5, scale = 2.3)
