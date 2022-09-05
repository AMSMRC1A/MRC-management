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

# load cost change cases
source("code/results/vary_cost_params_ebola.R")
# load plotting utilities
source("code/results/plotting_utils.R")

#### FIGURE A3: Ebola controls with changing costs -----------------------------
ebol_titles <- c("(a) increase cost of vaccination in Patch 1", 
                 #"(b) increase cost of vaccination in patch 2",
                 "(b) increase cost of hospitalization in Patch 1", 
                 "(c) increase cost of hospitalization in Patch 2")
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
    mutate(patch_control = paste(patch, control_type, sep = "-")) %>%
    ggplot(aes(x = time, y = value, color = patch_control, linetype = patch_control)) +
    geom_line() +
    facet_grid(rows = vars(variable),
               labeller = labeller(variable = ebola_control_labs),
               switch = "y",
               scales = "free") +
    labs(color = "Patch:", 
         linetype = "Control type:",
         subtitle = ebol_titles[i]
    ) +
    scale_color_manual(values = ebola_cols, labels = ebola_gov_labs, name = "Control:") +
    scale_linetype_manual(values = ebola_ltys, labels = ebola_gov_labs, name = "Control:") +
    theme_minimal_grid(12) +
    theme(
      axis.title.y = element_blank(),
      legend.position = "bottom",
      legend.key.width = unit(2,"lines"),
      legend.title = element_blank(),
      panel.spacing = unit(c(0), "lines"),
      plot.margin = unit(c(0.1,0.1,0.1,2), "lines"),
      strip.background = element_blank(),
      strip.placement = "outside" 
    ) +
    panel_border()
}

# fix legend
# make dummy plot with the correct legend
new_lty <-  c("solid", "23", "solid", "23")
names(new_lty) <- c("1-unique", "1-uniform",  "2-unique","2-uniform")
leg <- states  %>% 
  filter(variable == "I" & control_type != "none") %>%
  mutate(patch_control = paste(patch, control_type, sep = "-")) %>% 
  ggplot(aes(x = time, y = value, color = patch_control, linetype = patch_control)) +
  geom_line() +
  scale_color_manual(values = ebola_cols, labels = ebola_gov_labs, name = "Control:") +
  scale_linetype_manual(values = new_lty, labels = ebola_gov_labs, name = "Control:") +
  theme_bw() +
  theme(legend.position = "bottom")
leg <- get_legend(leg)


p <- lapply(p, function(i){i + theme(legend.position = "none")})
p[[2]] <- p[[2]] + theme(strip.text.y = element_blank())
p[[3]] <- p[[3]] + theme(strip.text.y = element_blank())
#p <- lapply(p[[2:3]], function(i){i + theme(axis.title.y = element_blank())})
bl <- ggplot() + theme_nothing()

Ebola_vary_cost_plots <- plot_grid(plot_grid(plotlist = p, ncol = 3),
                                   leg,
                                   ncol = 1, 
                                   rel_heights = c(0.95,0.05))
Ebola_vary_cost_plots

ggsave("figures/Ebola_vary_cost.pdf",
       Ebola_vary_cost_plots,
       width = 6, height = 3, scale = 2)