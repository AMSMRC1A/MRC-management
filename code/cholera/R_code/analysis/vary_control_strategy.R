# libraries
library(deSolve)
library(reshape2)
library(tidyverse)
library(cowplot)
library(pracma)

# to paralellize
library(doParallel)
library(foreach)
registerDoParallel(4) # update this 4 if you want to use more cores


# load files
source("CholeraSIRW_ODE.R")
source("Cholera_params.R")
source("Cholera_OCfunc.R")
source("Cholera_adjoints.R")

# initial guesses - controls
guess_v1 = rep(0,length(times))
guess_v2 = rep(0, length(times))

# setup optimal control parameters
delta = 0.01
oc_params <- c(b1 = 1, b2 = 1,
               C1 = 0.125, C2 = 0.125, 
               epsilon1  = 10000, epsilon2 = 10000)

# define parameters to sweep across
test_params <- expand.grid(m1 = c(0,0.05), 
                           m2 = c(0,0.05),
                           control_type = c("unique", "equiv", "max", "none"))
test_params$test_case <- 1:nrow(test_params)

# calculate OC
start_time <- Sys.time()
vary_params <- foreach (i=1:nrow(test_params), .packages = c("deSolve","tidyverse", "pracma")) %dopar% { 
  apply_oc(change_params = test_params[i,],
           guess_v1 = guess_v1, guess_v2 = guess_v2, 
           init_x = IC, bounds = bounds,
           ode_fn = chol, adj_fn = adj,
           times = times, params = c(params, oc_params), delta = delta, control_type = test_params[i,"control_type"])
}
vary_params <- do.call(rbind, vary_params)
end_time <- Sys.time()
end_time - start_time

# reformat output
j_vals <- lapply(1:length(vary_params), function(i){return(data.frame(test_case = i, vary_params[[i]][["j"]]))})
j_vals <- as.data.frame(do.call(rbind, j_vals))
j_vals <- left_join(test_params,j_vals)
j_vals$j = apply(j_vals[,5:8],1,sum)
mult_oc_params <- lapply(1:length(vary_params), function(i){return(data.frame(test_case = i, vary_params[[i]][["ts"]]))})
mult_oc_params <- as.data.frame(do.call(rbind, mult_oc_params))
mult_oc_params <- left_join(test_params,mult_oc_params)


# change to long for plotting
mult_oc_params <- melt(mult_oc_params %>% select(-test_case), 
                       id = c("m1", "m2", "control_type","time"))
mult_oc_params$scenario = with(mult_oc_params,ifelse(control_type == "unique", paste0(control_type,": ", variable), as.character(control_type)))

p1 = ggplot(data = mult_oc_params %>% filter(control_type %in% c("unique", "equiv"))) + 
  geom_line(aes(x = time, y = value, linetype = control_type, color = scenario), size = 1) +
  guides(linetype = FALSE) +
  facet_grid(rows = vars(m2), cols = vars(m1), labeller = label_both)+
  labs(y = "optimal vaccination rate") +
  scale_color_manual(values = c("black", "red", "blue")) +
  theme_bw()+
  theme(legend.position = "bottom")

# compute relative changes in cost
j_vals <- j_vals %>%
  # add column representing cost relative to "no control"
  mutate(rel_j = j/j[control_type=="none"])

  # separately compute the change in cost due to having separate controls
percent_change_unique_to_equiv <- j_vals %>%
  filter(control_type %in% c("unique","equiv")) %>%
  mutate(j_change = 100*(j-j[control_type=="equiv"])/j[control_type=="equiv"])

p2 = ggplot(data = filter(j_vals,control_type != "none")) + 
  geom_bar(aes(x = paste0("m1: ", m1,", m2: ", m2), y = 1-rel_j, fill=as.factor(control_type)), size = 3, position = "dodge", stat='identity') + 
  #geom_path(aes(x = as.factor(control_type), y = j, group = paste0("m1: ", m1,", m2: ", m2)), size = 1) +
  labs(x = "control strategy", y = "cost relative to no control") +
  theme_bw()+
  theme(legend.position = "bottom",
        legend.title = element_blank())

j_vals_long = melt(j_vals %>% select(-j, -rel_j), c("m1", "m2", "control_type","test_case"))

p3 = ggplot(data = j_vals_long, )+
  geom_bar(aes(x = , y = 1-rel_j, fill=as.factor(control_type)), size = 3, position = "dodge", stat='identity') + 
  facet_grid(cols = vars(paste0("m1: ", m1,", m2: ", m2)))
  labs(x = "control strategy", y = "absolute cost") +
  theme_bw()+
  theme(legend.position = "bottom",
        legend.title = element_blank())

plot_grid(p1, p2)
ggsave("figures/vary_control_strategies.pdf", width = 14, height = 6)
