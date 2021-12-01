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
registerDoParallel(detectCores() - 2) # update this if you want to be less aggressive with your core use


# load files--------------------------------------------------------------------
source("CholeraSIRW_ODE.R")
source("Cholera_params.R")
source("Cholera_OCfunc.R")
source("Cholera_adjoints.R")
source("Cholera_analysisfunc.R")


# Helper function to combine parameter sets we're testing-----------------------
test_param_builder <- function(change_params, all_params) {
  change_names <- names(select(change_params, -control_type, -experiment_type))
  test_params <- full_join(select(as.data.frame(all_params), -all_of(change_names)),
    as.data.frame(change_params),
    by = character()
  )
}

# initial guesses - controls----------------------------------------------------
guess_v1 <- rep(0, length(times))
guess_v2 <- rep(0, length(times))
guess_u1 <- rep(0, length(times))
guess_u2 <- rep(0, length(times))

# setup baseline optimal control parameters-------------------------------------
tol = 0.01 # tolerance parameter for optimization
oc_params <- list(b1 = 1, b2 = 1, # cost of cases
               C1 = 0.125, C2 = 0.125,  # cost of vaccinations
               epsilon1  = 10000, epsilon2 = 10000, # non-linearity for vacc
               D1 = 0.125, D2 = 0.125, # cost of sanitation 
               eta1 = 100, eta2 = 100 ) # non-linearity for sanitation

# Experiment 0: vary control type only------------------------------------------
# define parameters to sweep across
test_params <- expand.grid(
  m1 = 0, 
  m2 = 0,
  control_type = c("unique", "uniform")
)

# calculate OC
# begin system clock
# to keep track of run-time, guide decisions about code optimization
start_time <- Sys.time()

# run OC across multiple parameters
exper0 <- test_mult_params(
  test_params = test_params,
  base_params = c(params, oc_params),
  guess_v1 = guess_v1, guess_v2 = guess_v2,
  guess_u1 = guess_u1, guess_u2 = guess_u2,
  IC = IC, bounds = bounds, times = times, tol = tol
)

# calculate run-time and print it
end_time <- Sys.time()
end_time - start_time

exper0$states %>% 
  melt(c("m1", "m2", "control_type", "test_case", "time")) %>%
  mutate(patch = substr(variable, 2,2), 
         variable = substr(variable, 1,1)) %>%
  filter(control_type %in% c("uniform", "unique")) %>%
  mutate(plot_var = ifelse(control_type == "unique", paste("patch", patch), "uniform")) %>%
  ggplot(aes(x = time, y = value, color = plot_var)) + 
  geom_point(size = 1)+
  facet_wrap(~ variable, scale = "free_y") +
  scale_color_manual(values = c("red", "blue", "black")) +
  theme_bw()

# Experiment 1: vary movement and control type----------------------------------
# define parameters to sweep across
test_params <- expand.grid(
  m1 = c(0, 0.025), # max movement rate 5% / day
  m2 = c(0, 0.025), # max movement rate 5% / day
  control_type = c("unique", "uniform", "max", "none"),
  experiment_type = "movement"
)

test_params <- right_join(select(as.data.frame(all_params), -m1, -m2), as.data.frame(test_params), by = character())

# calculate OC
# begin system clock
# to keep track of run-time, guide decisions about code optimization
start_time <- Sys.time()

# run OC across multiple parameters
exper1 <- test_mult_params(
  test_params = test_params,
  base_params = all_params,
  guess_v1 = guess_v1, guess_v2 = guess_v2,
  guess_u1 = guess_u1, guess_u2 = guess_u2,
  IC = IC, bounds = bounds, times = times, tol = tol
)

# calculate run-time and print it
end_time <- Sys.time()
end_time - start_time

# Experiment 2: vary only disease dynamics parameters----------------------------------
# assume gamma (recovery) and delta (mortality) are characteristics of pathogen
# and therefore fixed
# define parameters to sweep across (+/- 10% from baseline parameter)
transmission_params <- expand.grid(
  beta_I1 = params$beta_I1 * c(0.9, 1, 1.1),
  beta_I2 = params$beta_I2 * c(0.9, 1, 1.1),
  beta_W1 = params$beta_W1 * c(0.9, 1, 1.1),
  beta_W2 = params$beta_W2 * c(0.9, 1, 1.1),
  control_type = c("unique", "uniform"),
  experiment_type = "transmission"
)
test_params <- test_param_builder(transmission_params, all_params)

# calculate OC
# begin system clock
# to keep track of run-time, guide decisions about code optimization
start_time <- Sys.time()

# run OC across multiple parameters
exper2 <- test_mult_params(
  test_params = test_params,
  base_params = c(params, oc_params),
  guess_v1 = guess_v1, guess_v2 = guess_v2,
  guess_u1 = guess_u1, guess_u2 = guess_u2,
  IC = IC, bounds = bounds, times = times, tol = tol
)

# calculate run-time and print it
end_time <- Sys.time()
end_time - start_time



ggplot(data = exper2$j_vals) +
  geom_line(
    data = data.frame(x = with(exper2$j_vals, c(min(c(j_case1, j_case2)), max(c(j_case1, j_case2))))),
    aes(x = x, y = x)
  ) +
  geom_point(aes(x = j_case1, y = j_case2, color = as.factor(paste(beta_I1, beta_I2, beta_W1, beta_W2)))) +
  geom_point(aes(x = j_vacc1, y = j_vacc2, color = as.factor(paste(beta_I1, beta_I2, beta_W1, beta_W2))), shape = 2) +
  facet_wrap(vars(control_type)) +
  theme_bw() +
  theme(legend.position = "none")


# Experiment 3: vary many parameters, organized by type-------------------------
# set up basic parameters with no contorl
base_params <- expand.grid(
  control_type = c("unique", "uniform", "none"),
  experiment_type = "control"
)
test_params <- test_param_builder(base_params, all_params)


# assume recovery and mortality are characteristics of pathogen (therefore fixed)
# define movement parameters to sweep across
movement_params <- expand.grid(
  m1 = c(0, 0.01, 0.005),
  m2 = c(0, 0.01, 0.005),
  control_type = c("unique", "uniform", "none"),
  experiment_type = "movement"
)
movement_params <- test_param_builder(movement_params, all_params)
test_params <- rbind(test_params, movement_params)
# define disease parameters to sweep across (+/- 10% from baseline parameter)
transmission_params <- expand.grid(
  beta_W1 = params$beta_W1 * c(0.9, 1, 1.1),
  beta_W2 = params$beta_W2 * c(0.9, 1, 1.1),
  control_type = c("unique", "uniform", "none"),
  experiment_type = "transmission"
)
transmission_params <- test_param_builder(transmission_params, all_params)
test_params <- rbind(test_params, transmission_params)
# define environmental parameters to sweep across (+/- 10% from baseline parameter)
env_params <- expand.grid(
  rho1 = params$rho1 * c(0.9, 1, 1.1),
  rho2 = params$rho2 * c(0.9, 1, 1.1),
  control_type = c("unique", "uniform", "none"),
  experiment_type = "environmental"
)
env_params <- test_param_builder(env_params, all_params)
test_params <- rbind(test_params, env_params)
# define cost parameters to sweep across (+/- 10% from baseline parameter)
cost_params <- expand.grid(
  b2 = oc_params$b2 * c(0.9, 1, 1.1),
  C1 = oc_params$C1 * c(0.9, 1, 1.1),
  C2 = oc_params$C2 * c(0.9, 1, 1.1),
  epsilon1 = oc_params$epsilon1 * c(0.9, 1, 1.1),
  epsilon2 = oc_params$epsilon2 * c(0.9, 1, 1.1),
  control_type = c("unique", "uniform", "none"),
  experiment_type = "cost"
)
cost_params <- test_param_builder(cost_params, all_params)
test_params <- rbind(test_params, cost_params)

# calculate OC
# begin system clock
# to keep track of run-time, guide decisions about code optimization
start_time <- Sys.time()
# run OC across multiple parameters
exper3 <- test_mult_params(
  test_params = test_params,
  base_params = c(params, oc_params),
  guess_v1 = guess_v1, guess_v2 = guess_v2,
  IC = IC, bounds = bounds, times = times, tol = tol
)
# calculate run-time and print it
end_time <- Sys.time()
end_time - start_time

## reformat output
j_vals <- exper3$j_vals

# plot cost outcomes (unique/uniform) based on the type of parameter change
max_j_unif <- max(j_vals %>% filter(control_type == "uniform") %>% pull(j))
max_j_uni <- max(j_vals %>% filter(control_type == "unique") %>% pull(j))
j_vals %>%
  select(-test_case) %>%
  reshape2::melt(c(
    "m1", "m2", "beta_W1", "beta_W2", "rho1", "rho2",
    "b2", "C1", "C2", "epsilon1", "epsilon2",
    "control_type", "experiment_type"
  )) %>%
  reshape2::dcast(m1 + m2 + beta_W1 + beta_W2 + rho1 + rho2 +
    b2 + C1 + C2 + epsilon1 + epsilon2 +
    variable + experiment_type ~ control_type, value.var = "value") %>%
  mutate(test_case = 1:n()) %>%
  filter(variable %in% c("j_case1", "j_case2", "j_vacc1", "j_vacc2", "j")) %>%
  ggplot() +
  geom_abline() +
  geom_abline(slope = 1.1, color = "grey") +
  geom_abline(slope = 0.9, color = "grey") +
  geom_abline(slope = 1.05, color = "grey") +
  geom_abline(slope = 0.95, color = "grey") +
  geom_point(aes(x = unique, y = uniform, color = as.factor(experiment_type)), size = 2, alpha = 0.5) +
  geom_label(
    data = data.frame(
      y = c(max_j_unif, max_j_unif, max_j_uni * 0.9, max_j_uni * 0.95),
      x = c(max_j_unif * 0.9, max_j_unif * .95, max_j_uni, max_j_uni),
      text = c("+10%", "+5%", "-10%", "-5%"),
      variable = c("j", "j")
    ),
    aes(x = x, y = y, label = text), color = "grey", label.size = NA
  ) +
  facet_wrap(vars(variable), scales = "free") +
  scale_size_continuous(range = c(0.5, 1.5)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid = element_blank()
  )

#### PLOT RESULTS: EH to put into functions ####

# change to long for plotting---------------------------------------------------
v_states <- select(exper3$states, m1, m2, control_type, time, v1, v2) # will streamline this later

mult_oc_params <- melt(v_states,
  id = c("m1", "m2", "control_type", "time")
)
mult_oc_params$scenario <- with(
  mult_oc_params,
  ifelse(control_type == "unique",
    paste0(control_type, ": ", variable),
    as.character(control_type)
  )
)

# plot 1: vaccination strategies------------------------------------------------
p1 <- mult_oc_params %>%
  ggplot() +
  geom_line(aes(x = time, y = value, linetype = control_type, color = scenario), size = 1) +
  guides(linetype = "none") +
  facet_grid(rows = vars(m2), cols = vars(m1), labeller = label_both) +
  labs(y = "optimal vaccination rate") +
  scale_color_manual(values = c("black", "red", "blue", "orange")) +
  theme_bw() +
  theme(legend.position = "bottom")


# compute relative changes in cost
j_vals <- exper3$j_vals
j_vals <- j_vals %>%
  # add column representing cost relative to "no control"
  mutate(rel_j = j / j[control_type == "none"])

# separately compute the change in cost due to having separate controls
percent_change_unique_to_uniform <- j_vals %>%
  filter(control_type %in% c("unique", "uniform")) %>%
  mutate(j_change = 100 * (j - j[control_type == "uniform"]) / j[control_type == "uniform"])

# plot 2: relative costs of strategies------------------------------------------
p2 <- ggplot(data = filter(j_vals, control_type != "none")) +
  geom_bar(aes(
    x = paste0("m1: ", m1, ", m2: ", m2),
    y = 1 - rel_j, fill = as.factor(control_type)
  ),
  size = 3, position = "dodge", stat = "identity"
  ) +
  # geom_path(aes(x = as.factor(control_type), y = j, group = paste0("m1: ", m1,", m2: ", m2)), size = 1) +
  labs(x = "control strategy", y = "cost relative to no control") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )

j_vals_long <- melt(j_vals %>% select(-j, -rel_j), c("m1", "m2", "control_type", "test_case"))

# plot 3: direct cost comparisons-----------------------------------------------
p3 <- ggplot(data = j_vals_long, ) +
  geom_bar(aes(x = as.factor(control_type), y = value, fill = variable),
    size = 3, position = "stack", stat = "identity"
  ) +
  facet_grid(rows = vars(m2), cols = vars(m1), labeller = label_both) +
  labs(x = "control strategy", y = "absolute cost") +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )

# put plots together and save---------------------------------------------------
#plot_grid(p1, p3, p2, nrow = 1)
#ggsave("figures/vary_control_strategies.pdf", width = 14, height = 6)

# plot X: infections time series, unique controls and none----------------------
states_long <- melt(exper3$states %>% select(-test_case), c("m1", "m2", "control_type", "time"))
states_long$state <- substr(states_long$variable, 1, 1)
states_long$patch <- substr(states_long$variable, 2, 2)
ggplot(
  data = states_long %>% filter(
    state == "I",
    control_type %in% c("none", "unique"),
    m1 == 0.05,
    m2 == 0
  ),
  aes(x = time, y = value, color = patch, linetype = control_type)
) +
  geom_line() +
  labs(y = "infections") +
  scale_color_manual(values = c("red", "blue")) +
  theme_bw() +
  theme(legend.position = "bottom")
