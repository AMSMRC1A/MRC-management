## create data.frame of parameters to test
# assume costs are equal in both patches
equal <- expand.grid(b1 = c(1,10), 
                     C1 = c(0.125, 0.625, 1.25), 
                     epsilon1 = c(1000, 100000)) # consider changing to 3E5
equal$b2 <- equal$b1
equal$C2 <- equal$C1
equal$epsilon2 <- equal$epsilon1
equal$scenario <- 1
# assume patch 2 has 5X cost of patch 1
p2greater <- expand.grid(b1 = c(1,10), 
                         C1 = c(0.125, 0.625, 1.25), 
                         epsilon1 = c(1000, 100000)) # consider changing to 3E5
p2greater$b2 <- 5*p2greater$b1
p2greater$C2 <- 5*p2greater$C1
p2greater$epsilon2 <- 5*p2greater$epsilon1
p2greater$scenario <- 2
# assume patch 1 has 5X cost of patch 2
p1greater <- expand.grid(b1 = c(1,10)*5, 
                         C1 = c(0.125, 0.625, 1.25)*5, 
                         epsilon1 = c(1000, 100000)*5) # consider changing to 3E5
p1greater$b2 <- 1/5*p1greater$b1
p1greater$C2 <- 1/5*p1greater$C1
p1greater$epsilon2 <- 1/5*p1greater$epsilon1
p1greater$scenario <- 3
# combine into single data.frame to run all combinations
test_params <- bind_rows(equal, p1greater, p2greater)
test_params$counter <- 1:nrow(test_params)

# calculate OC
start_time <- Sys.time()
mult_oc_params <- foreach (i=1:nrow(test_params), .packages = c("deSolve","tidyverse")) %dopar% { 
  apply_oc(change_params = test_params[i,],
           guess_v1 = guess_v1, guess_v2 = guess_v2, 
           init_x = IC, bounds = bounds,
           ode_fn = chol, adj_fn = adj,control_type = "unique",
           times = times, params = c(params, oc_params), delta = delta)
}
mult_oc_params <- do.call(rbind, mult_oc_params)
end_time <- Sys.time()
end_time - start_time

# change to long for plotting
mult_oc_params <- melt(mult_oc_params %>% select(-counter), 
                       id = c("b1", "b2", "C1", "C2", "epsilon1", "epsilon2","j","scenario","time"))
# rename scenarios to be meaningful
mult_oc_params$scenario <- sapply(mult_oc_params$scenario, function(i){
  switch(i, "Equal cost", "Greater cost in P2", "Greater cost in P1")})

# add j values to plot
j_vals <- mult_oc_params %>% 
  group_by(b1, b2, C1, C2, epsilon1, epsilon2, scenario) %>% 
  summarise(j = unique(j)) %>% 
  mutate(label = paste0("j = ", round(j)), 
         x = 200, 
         y = max(bounds)*9/10)

# plot outcomes and save in a single file
plts <- list()
for(i in unique(mult_oc_params$scenario)){
  p <- ggplot(data = mult_oc_params %>% filter(scenario == i)) + 
    geom_line(aes(x = time, y = value, color = variable)) + 
    geom_text(data = j_vals %>% filter(scenario == i), aes(x = x, y = y, label = label), hjust = 1, vjust = 0) +
    ggtitle(i) +
    facet_grid(cols = vars(b1,b2, epsilon1, epsilon2), rows = vars(C1, C2), labeller = "label_both") +
    theme_half_open(12) +
    theme(legend.position = "bottom")
  plts[[i]] <- p
}

pdf("vary_cost_parameters.pdf")
invisible(lapply(plts, print))
dev.off()