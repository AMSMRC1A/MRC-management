#### BASELINE CASES ------------------------------------------------------------
# load baseline cases
source("code/results/baseline_cases.R")

#### SECTION 3.1: CHOLERA BASELINE 
# initial conditions at start of control
test2[[1]]$uncontrolled[nrow(test2[[1]]$uncontrolled),]
# cases before control implemented
200000 - sum(test2[[2]]$uncontrolled[nrow(test2[[2]]$uncontrolled),c("S1", "S2")])
# deaths before control implemented
200000 - sum(test2[[1]]$uncontrolled[nrow(test2[[1]]$uncontrolled), c("S1", "S2", "I1", "I2", "R1", "R2")])

# % change vs. no control
j_vals %>%
  filter(variable %in% c("epi_case1", "epi_case2"), 
         model == "cholera") %>%
  group_by(control_type, model) %>%
  summarise(value = sum(value)) %>%
  group_by(model) %>%
  mutate(change = (value - value[control_type == "none"])/value[control_type == "none"])


#### TABLE 2: EPI/CONTROL OUTCOMES 
j_vals %>%
  filter(variable %in% c("epi_case1", "epi_case2", "epi_death1", "epi_death2", 
                         "res_vacc1", "res_vacc2", "res_sani1","res_sani2"), 
         control_type != "none") %>%
  mutate(patch = substr(variable, nchar(variable), nchar(variable)), 
         variable = substr(variable, 1, nchar(variable)-1)) %>%
  mutate(value = round(value)) %>%
  dcast(model + control_type + patch ~ variable, value.var = "value")
  

#### SECTION 3.2: EBOLA BASELINE
# initial conditions at start of control
round(test2[[3]]$uncontrolled[nrow(test2[[3]]$uncontrolled),])
# cases before control implemented
round(200000 - sum(test2[[3]]$uncontrolled[nrow(test2[[3]]$uncontrolled),c("S1", "S2")]))
# deaths before control implemented
round(200000 - sum(test2[[3]]$uncontrolled[nrow(test2[[3]]$uncontrolled),c("S1", "S2","E1", "E2", "I1", "I2", "H1", "H2", "R1", "R2")]))

# days hosp switches which higher
states %>% 
  filter(# do not show before control
    time >= 0,
    # exclude no control
    v1_max != 0, 
    model == "ebola", 
    variable == "u") %>%
  dcast(time ~ paste0(control_type,patch), value.var = "value") %>%
  mutate(flag = ifelse(unique1 < uniform1, 1, 0)) %>%
  filter(flag ==1) %>%
  pull(time) %>%
  min()


#### CHOLERA COST CHANGE -------------------------------------------------------
# load cost change cases
source("code/results/vary_cost_params_cholera.R")

# start of sanitation control in (b) of Figure 7
states %>%
  filter(control_type == "unique",
         id == "C2",
         variable == "u",
         patch == 2, 
         value > 0.01) %>%
  summarise(m = min(time))

# Figure 7 caption 
j_vals %>% 
  filter(variable %in% c("j_vacc1", "j_vacc2", "j_sani1", "j_sani2")) %>%
  group_by(id, control_type) %>%
  summarise(value = sum(value)) %>%
  mutate(value = round(value)) %>%
  dcast(id ~ control_type, value.var = "value")


#### EBOLA COST CHANGE -------------------------------------------------------
# load cost change cases
source("code/results/vary_cost_params_ebola.R")
  
# Figure A3 caption 
j_vals %>% 
  filter(variable %in% c("j_vacc1", "j_vacc2", "j_sani1", "j_sani2")) %>%
  group_by(id, control_type) %>%
  summarise(value = sum(value)) %>%
  mutate(value = round(value)) %>%
  dcast(id ~ control_type, value.var = "value")






## values for text:
# number vaccinated
# number hospitalized/sanitized
# number deaths
# j values
temp_df_unique <- states %>%
  filter(time >= 0) %>% 
  filter(model == "ebola") %>%
  filter(control_type == "unique") %>%
  select(-c(model, control_type,v1_max, v2_max, u1_max, u2_max, plot_var)) %>%
  mutate(compartment = paste0(variable,patch)) %>%
  select(-c(patch, variable)) %>%
  relocate(value, .after = compartment) %>% 
  unique() %>%
  pivot_wider(names_from = compartment)

j_vals_unique <- calc_j_ebola(temp_df_unique$time,
                              select(temp_df_unique,`S1`:`u2`),
                              ebol_mod_details$params)

temp_df_uniform <- states %>%
  filter(time >= 0) %>% 
  filter(model == "ebola") %>%
  filter(control_type == "uniform") %>%
  select(-c(model, control_type,v1_max, v2_max, u1_max, u2_max, plot_var)) %>%
  mutate(compartment = paste0(variable,patch)) %>%
  select(-c(patch, variable)) %>%
  relocate(value, .after = compartment) %>% 
  unique() %>%
  pivot_wider(names_from = compartment)

j_vals_uniform <- calc_j_ebola(temp_df_uniform$time,
                               select(temp_df_uniform,`S1`:`u2`),
                               ebol_mod_details$params)




## values for text
j_vals %>%
  select(-test_case) %>%
  dcast(variable + model ~ control_type) %>%
  filter(model == "ebola", 
         substr(variable, 1,3) == "epi")




# Get data for tables
temp_df_unique <- states %>%
  filter(time >= 0) %>% 
  filter(model == "cholera") %>%
  filter(control_type == "unique") %>%
  select(-c(model, control_type,v1_max, v2_max, u1_max, u2_max, plot_var)) %>%
  mutate(compartment = paste0(variable,patch)) %>%
  select(-c(patch, variable)) %>%
  relocate(value, .after = compartment) %>% 
  unique() %>%
  pivot_wider(names_from = compartment)

j_vals_unique <- calc_j_cholera(temp_df_unique$time,
                                select(temp_df_unique,`S1`:`u2`),
                                chol_mod_details$params)

temp_df_uniform <- states %>%
  filter(time >= 0) %>% 
  filter(model == "cholera") %>%
  filter(control_type == "uniform") %>%
  select(-c(model, control_type,v1_max, v2_max, u1_max, u2_max, plot_var)) %>%
  mutate(compartment = paste0(variable,patch)) %>%
  select(-c(patch, variable)) %>%
  relocate(value, .after = compartment) %>% 
  unique() %>%
  pivot_wider(names_from = compartment)

j_vals_uniform <- calc_j_cholera(temp_df_uniform$time,
                                 select(temp_df_uniform,`S1`:`u2`),
                                 chol_mod_details$params)







#### values for text
states %>%
  filter(control_type == "unique",
         id == "C2",
         variable == "u",
         patch == 2, 
         value > 0.01) %>%
  summarise(m = min(time))

# caption 
j_vals %>% 
  filter(variable %in% c("j_vacc1", "j_vacc2", "j_sani1", "j_sani2")) %>%
  group_by(id, control_type) %>%
  summarise(value = sum(value)) %>%
  mutate(value = round(value)) %>%
  dcast(id ~ control_type, value.var = "value")







#### values for text
# relative costs
j_vals <- lapply(
  1:nrow(test_params),
  function(i) {
    j <- test2[[i]][["j"]]
    j <- c(j, j_tot = sum(j[substr(names(j), 1, 1) == "j"]))
    return(data.frame(test_case = i, melt(j)))
  }
)
j_vals <- as.data.frame(do.call(rbind, j_vals))
j_vals <- left_join(test_params, j_vals, by = "test_case") %>%
  rename(variable = L1)

# caption 
j_vals %>% 
  filter(variable %in% c("j_vacc1", "j_vacc2", "j_sani1", "j_sani2")) %>%
  group_by(id, control_type) %>%
  summarise(value = sum(value)) %>%
  mutate(value = round(value)) %>%
  dcast(id ~ control_type, value.var = "value")
