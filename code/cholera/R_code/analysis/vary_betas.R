# libraries
library(deSolve)
library(reshape2)
library(tidyverse)
library(cowplot)

# load files
source("CholeraSIRW_ODE.R")
source("Cholera_params.R")
source("Cholera_OCfunc.R")
source("Cholera_runODE.R")
source("Cholera_adjoints.R")

# initial guesses - controls
guess_v1 = rep(0,length(times))
guess_v2 = rep(0, length(times))

#### test multiple oc parameter values ####

## testing multiple betas

betas <- expand.grid(beta_I1 = c(0,2.64E-5),
                     beta_I2 = c(0,2.64E-5),
                     beta_W1 = 1.21E-4*c(0.1,1,10),
                     beta_W2 = 1.21E-4*c(0.1,1,10))
betas$scenario <- 1

# combine into single data.frame to run all combinations
test_params <- bind_rows(betas)
test_params$counter <- 1:nrow(test_params)

# calculate trajectories

# replicate oc analysis across multiple parameters
get.trajectories <- function(change_params,init_x,times,ode_fn,params,v1_interp,v2_interp) {
  if("counter" %in% names(change_params)){print(as.numeric(change_params["counter"]))}
  # update parameters
  new_params <- params 
  p_loc <- match(names(change_params), names(new_params))
  new_params[p_loc[!is.na(p_loc)]] = change_params[!is.na(p_loc)]
  x <- as.data.frame(ode(y = init_x, times = times, func = ode_fn, parms = new_params, 
           v1_interp = v1_interp, v2_interp = v2_interp))
  
  ret <- data.frame(t(change_params)) %>%
    bind_cols(time=times, I1=x$I1, I2=x$I2)
#    bind_cols(time=times, S1=x$S1, I1=x$I1, R1=x$R1, W1=x$W1,
#              S2=x$S2, I2=x$I2, R2=x$R2, W2=x$W2)
  return(ret)
  
}

start_time <- Sys.time()
mult_trajectories <- apply(test_params, 1, get.trajectories,
                            init_x = IC,
                            times = times,
                            ode_fn = chol,
                            params = c(params),
                            v1_interp = guess_v1,
                            v2_interp = guess_v2)

mult_trajectories <- do.call(rbind, mult_trajectories)
end_time <- Sys.time()
end_time - start_time

# change to long for plotting
mult_trajectories <- melt(mult_trajectories %>% select(-counter), 
                       id = c("beta_I1", "beta_I2", "beta_W1", "beta_W2","scenario","time"))


# plot outcomes and save in a single file
plts <- list()
plts <- ggplot(data = mult_trajectories) + 
    geom_line(aes(x = time, y = value, color = variable)) + 
    facet_grid(cols = vars(beta_W1, beta_W2), rows = vars(beta_I1, beta_I2), labeller = "label_both") +
    theme_half_open(12) +
    theme(legend.position = "bottom")

plts
