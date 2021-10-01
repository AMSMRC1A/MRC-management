# libraries
library(deSolve)
library(reshape2)
library(tidyverse)
library(cowplot)
library(pracma)
library(RColorBrewer)

# to paralellize
library(doParallel)
library(foreach)
registerDoParallel(4) # update this 4 if you want to use more cores


# load files--------------------------------------------------------------------
source("CholeraSIRW_ODE.R")
source("Cholera_params.R")
source("Cholera_OCfunc.R")
source("Cholera_adjoints.R")
source("Cholera_analysisfunc.R")

# initial guesses - controls----------------------------------------------------
guess_v1 = rep(0,length(times))
guess_v2 = rep(0, length(times))

# setup baseline optimal control parameters-------------------------------------
tol = 0.01 # tolerance parameter for optimization
oc_params <- c(b1 = 1, b2 = 1, # cost of cases
               C1 = 0.125, C2 = 0.125,  # cost of vaccinations
               epsilon1  = 200000, epsilon2 = 10000) # non-linearity

# Experiment 1: vary movement and control type----------------------------------
# define parameters to sweep across
test_params <- expand.grid(m1 = c(0,0.025), # max movement rate 5% / day
                           m2 = c(0,0.025), # max movement rate 5% / day
                           control_type = c("unique", "uniform", "max", "none"))

# calculate OC
# begin system clock 
# to keep track of run-time, guide decisions about code optimization
start_time <- Sys.time()

# run OC across multiple parameters
exper <- test_mult_params(test_params = test_params, 
                          return_type = c("X", "j", "v"), 
                          base_params = c(params, oc_params), 
                          guess_v1 = guess_v1, guess_v2 = guess_v2, 
                          IC = IC, bounds = bounds, times = times, tol = tol)

# calculate run-time and print it
end_time <- Sys.time()
end_time - start_time

## reformat output
exper <- reformat_mult_params_output(output = exper, test_params = test_params)



#### PLOT RESULTS: EH to put into functions ####

# change to long for plotting---------------------------------------------------
mult_oc_params <- melt(exper$v %>% select(-test_case), 
                       id = c("m1", "m2", "control_type","time"))
mult_oc_params$scenario = with(mult_oc_params,
                               ifelse(control_type == "unique",
                                      paste0(control_type,": ", variable),
                                      as.character(control_type)))

# plot 1: vaccination strategies------------------------------------------------
p1 = ggplot(data = mult_oc_params %>% filter(control_type %in% c("unique", "uniform"))) + 
  geom_line(aes(x = time, y = value, linetype = control_type, color = scenario), size = 1) +
  guides(linetype = "none") +
  facet_grid(rows = vars(m2), cols = vars(m1), labeller = label_both)+
  labs(y = "optimal vaccination rate") +
  scale_color_manual(values = c("black", "red", "blue")) +
  theme_bw()+
  theme(legend.position = "bottom")


# compute relative changes in cost
j_vals <- exper$j
j_vals <- j_vals %>%
  # add column representing cost relative to "no control"
  mutate(rel_j = j/j[control_type=="none"])

# separately compute the change in cost due to having separate controls
percent_change_unique_to_uniform <- j_vals %>%
  filter(control_type %in% c("unique","uniform")) %>%
  mutate(j_change = 100*(j-j[control_type=="uniform"])/j[control_type=="uniform"])

# plot 2: relative costs of strategies------------------------------------------
p2 = ggplot(data = filter(j_vals,control_type != "none")) + 
  geom_bar(aes(x = paste0("m1: ", m1,", m2: ", m2), 
               y = 1-rel_j, fill=as.factor(control_type)), 
           size = 3, position = "dodge", stat='identity') + 
  #geom_path(aes(x = as.factor(control_type), y = j, group = paste0("m1: ", m1,", m2: ", m2)), size = 1) +
  labs(x = "control strategy", y = "cost relative to no control") +
  theme_bw()+
  theme(legend.position = "bottom",
        legend.title = element_blank())

j_vals_long = melt(j_vals %>% select(-j, -rel_j), c("m1", "m2", "control_type","test_case"))

# plot 3: direct cost comparisons-----------------------------------------------
p3 = ggplot(data = j_vals_long, )+
  geom_bar(aes(x = as.factor(control_type), y = value, fill=variable), 
           size = 3, position = "stack", stat='identity') + 
  facet_grid(rows = vars(m2), cols = vars(m1), labeller = label_both)+
  labs(x = "control strategy", y = "absolute cost") +
  scale_fill_brewer(palette = "Set2") +
  theme_bw()+
  theme(legend.position = "bottom",
        legend.title = element_blank())

# put plots together and save---------------------------------------------------
plot_grid(p1, p3, p2, nrow = 1)
ggsave("figures/vary_control_strategies.pdf", width = 14, height = 6)

# plot X: infections time series, unique controls and none----------------------
states_long <- melt(exper$X %>% select(-test_case), c("m1", "m2", "control_type", "time"))
states_long$state <- substr(states_long$variable, 1,1)
states_long$patch <- substr(states_long$variable, 2,2)
ggplot(data = states_long %>% filter(state == "I", 
                                     control_type %in% c("none", "unique"),
                                     m1 == 0.05,
                                     m2 == 0), 
       aes(x = time, y = value, color = patch, linetype = control_type)) +
  geom_line() + 
  labs(y = "infections") + 
  scale_color_manual(values =c("red", "blue")) +
  theme_bw() + 
  theme(legend.position = "bottom")


