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

#### FIGURE 6: Ebola relative costs --------------------------------------------
patch_labs <- c("Patch 1", "Patch 2", "Total")
names(patch_labs) <- c("1", "2", "t")

var_labs <- c("Vaccinations", "Hospitalizations", "Cases", "Total cost")
names(var_labs) <- c("vacc", "hosp", "case", "to")

fig <- j_vals %>%
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
    variable %in% c(paste0("j_", c("tot")),
                    paste0("epi_", c("case1", "case2", "hosp1", "hosp2")),
                    paste0("res_", c("vacc1", "vacc2")))
  ) %>%
  mutate(variable_short = factor(variable_short, levels = c("vacc", "hosp", "case", "to")),
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

ggsave("figures/Ebola_relative_costs.pdf",
       plot = fig,
       width = 6, height = 3, scale = 2)

