require(deSolve)
require(tidyverse)
require(reshape2)
require(ggpubr)
require(RColorBrewer)
source("ebola_functions_mrc.R")
source("ebola_params_mrc.R")

##Look at changes in optimal vaccination rate and costs associated with disease
##and vaccination as movement rates change

#set movement rates
m_values<-seq(-5,-1,by=1)
m_grid<-expand.grid(m1=m_values,m2=m_values)  #Expands 2 (or more) vectors into a grid

#Set maximum control levels M1,M2
maxControl<-c(.1,.1)

#Set timeframe
times=seq(0,730,by=1)

#set intial conditions (initial outbreak) - 10 I in patch 1
initial=c(S1=parm$N1-10,E1=0,I1=10,H1=0,D1=0,R1=0,S2=parm$N2,E2=0,I2=0,H2=0,D2=0,R2=0)

v1 <- data.frame(times = times, v1 = rep(0,length(times)))
v1_interp <- approxfun(v1, rule = 2)
v1 = v1[,2]
v2 = data.frame(times = times, v2 = rep(0,length(times)))
v2_interp <- approxfun(v2, rule = 2)
v2 = v2[,2]
parm_temp<-c(parm,v1_interp=v1_interp,v2_interp=v2_interp)
##Simulate out uncontrolled outbreak------------------------
outbreak_sens<-lapply(1:dim(m_grid)[1],function(x){
  parm_temp$m1<-10^(m_grid[x,"m1"]) #overwrite variables of interest in the temp
  parm_temp$m2<-10^(m_grid[x,"m2"]) #parameter set
  temp<-ode(y=initial,
            times=times,
            func=ebola.ode,
            parms=parm_temp,
            method = "ode45")%>%
    data.frame()
  size<-list(SolX=temp,
             `Patch 1`=with(c(parm,temp),sum(b1*(betaI1*S1*I1 + betaD1*S1*D1))),
             `Patch 2`=with(c(parm,temp),sum(b2*(betaI2*S2*I2 + betaD2*S2*D2))))
  print(x)
  return(size)
})
#Extract outbreak size information
size_sens<-lapply(1:dim(m_grid)[1],function(x){
  temp<-data.frame("Patch 1"=outbreak_sens[[x]]$`Patch 1`,
                   "Patch 2"=outbreak_sens[[x]]$`Patch 2`)
  return(temp)
})%>%bind_rows()%>%
  bind_cols(m_grid)%>%
  mutate(Total=Patch.1+Patch.2)%>%
  pivot_longer(cols=c("Patch.1","Patch.2","Total"),
               names_to = "Location",values_to = "Outbreak_Size")
size_sens$Location<-size_sens$Location%>%
  factor(levels=c("Patch.1","Patch.2","Total"),
         labels=c("Patch 1","Patch 2","Total"))

##Solve Optimal Control Problem for each parameter combination in m_grid--------------
parm_temp<-parm
movement_sens<-lapply(1:dim(m_grid)[1],function(x){
  parm_temp$m1<-10^(m_grid[x,"m1"]) #overwrite variables of interest in the temp
  parm_temp$m2<-10^(m_grid[x,"m2"]) #parameter set
  temp<-ebola.optim(inits=initial,M=maxControl,params=parm_temp,
                    times=times,maxIter=100)
  return(temp)
})

#Extract Optimal Vaccination Rates for each parameter set------------------
movement_sens_vacc<-lapply(1:dim(m_grid)[1],function(x){
  temp<-data.frame(m1=m_grid[x,"m1"],m2=m_grid[x,"m2"],
                   movement_sens[[x]]$SolX)
  return(temp)
})%>%bind_rows()%>%
  select(m1,m2,time,v1,v2)%>%
  pivot_longer(cols=c("v1","v2"),names_to = "patch",values_to = "vacc_rate")%>%
  mutate(m1=paste0("10^",m1),m2=paste0("10^",m2))

#Renaming parameters so that they are more interpretable
movement_sens_vacc$m1<-movement_sens_vacc$m1%>%
  factor(levels=paste0("10^",seq(-5,1)),ordered = T)
movement_sens_vacc$m2<-movement_sens_vacc$m2%>%
  factor(levels=paste0("10^",seq(-5,1)),ordered = T)

#Extract cost values---------------------
movement_sens_costs<-lapply(1:dim(m_grid)[1],function(x){
  temp<-data.frame(m1=m_grid[x,"m1"],m2=m_grid[x,"m2"],
                   J11=movement_sens[[x]]$J11,
                   J12=movement_sens[[x]]$J12,
                   J21=movement_sens[[x]]$J21,
                   J22=movement_sens[[x]]$J22)
  return(temp)
})%>%bind_rows()%>%
  mutate(Total=J11+J12+J21+J22,Total_1=J11+J21,Total_2=J12+J22)%>%
  pivot_longer(cols=c("J11","J12","J21","J22","Total","Total_1","Total_2"),names_to = "Type",values_to = "Cost")%>%
  mutate(m1=10^m1,m2=10^m2)
  
#Rename costs to be interpretable
movement_sens_costs$Type<-movement_sens_costs$Type%>%
  factor(levels=c("J11","J12","J21","J22","Total_1","Total_2","Total"),
         labels=c("Patch 1 Disease Cost","Patch 2 Disease Cost",
                  "Patch 1 Vaccine Cost","Patch 2 Vaccine Cost",
                  "Patch 1 Total Cost","Patch 2 Total Cost",
                  "Total Cost"),
         ordered = T)
##Plot uncontrolled outbreak sizes-------------------------------------
outbreak_plot_1<-ggplot(filter(size_sens,Location!="Total"),
                             aes(x=10^m1,y=10^m2,fill=Outbreak_Size))+
  geom_tile()+
  xlab("Movement Out of Patch 1")+
  ylab("Movement Out of Patch 2")+
  scale_fill_distiller(name="Total Infections",palette="YlOrRd",direction = 1)+
  scale_x_log10()+
  scale_y_log10()+
  facet_wrap(~Location,ncol = 2)+
  theme_bw()
outbreak_plot_2<-ggplot(filter(size_sens,Location=="Total"),
                        aes(x=10^m1,y=10^m2,fill=Outbreak_Size))+
  geom_tile()+
  xlab("Movement Out of Patch 1")+
  ylab("Movement Out of Patch 2")+
  scale_fill_distiller(name="Total Infections",palette="YlOrRd",direction = 1)+
  scale_x_log10()+
  scale_y_log10()+
  theme_bw()

(outbreak_plot<-ggarrange(outbreak_plot_1,outbreak_plot_2,
                          ncol=1,nrow=2))
rm(outbreak_plot_1,outbreak_plot_2)
##Plot Cost Values - Done this way so that each figure has its color range defined------------------
##by its range of costs. Otherwise, the discrepancy in magnitude makes differences in 
##figures hard to see
movement_cost_plot_1<-ggplot(filter(movement_sens_costs,Type=="Patch 1 Disease Cost"|Type=="Patch 2 Disease Cost"),
                              aes(x=m1,y=m2,fill=Cost))+
  geom_tile()+
  xlab("Movement Out of Patch 1")+
  ylab("Movement Out of Patch 2")+
  scale_fill_distiller(trans="log10",palette="YlOrRd",direction = 1)+
  scale_x_log10()+
  scale_y_log10()+
  facet_wrap(~Type,ncol = 2)+
  theme_bw()
movement_cost_plot_2<-ggplot(filter(movement_sens_costs,Type=="Patch 1 Vaccine Cost"|Type=="Patch 2 Vaccine Cost"),
                             aes(x=m1,y=m2,fill=Cost))+
  geom_tile()+
  xlab("Movement Out of Patch 1")+
  ylab("Movement Out of Patch 2")+
  scale_x_log10()+
  scale_y_log10()+
  scale_fill_distiller(trans="log10",palette="YlOrRd",direction = 1)+
  facet_wrap(~Type,ncol = 2)+
  theme_bw()
movement_cost_plot_3<-ggplot(filter(movement_sens_costs,Type=="Patch 1 Total Cost"|Type=="Patch 2 Total Cost"),
                             aes(x=m1,y=m2,fill=Cost))+
  geom_tile()+
  xlab("Movement Out of Patch 1")+
  ylab("Movement Out of Patch 2")+
  scale_x_log10()+
  scale_y_log10()+
  scale_fill_distiller(trans="log10",palette="YlOrRd",direction = 1)+
  facet_wrap(~Type,ncol = 2)+
  theme_bw()
movement_cost_plot_4<-ggplot(filter(movement_sens_costs,Type=="Total Cost"),
                             aes(x=m1,y=m2,fill=Cost))+
  geom_tile()+
  xlab("Movement Out of Patch 1")+
  ylab("Movement Out of Patch 2")+
  scale_x_log10()+
  scale_y_log10()+
  scale_fill_distiller(trans="log10",palette="YlOrRd",direction = 1)+
  facet_wrap(~Type,ncol = 2)+
  theme_bw()
(movement_cost_plot<-ggarrange(movement_cost_plot_1,movement_cost_plot_2,
                              movement_cost_plot_3,movement_cost_plot_4,
                              ncol=2,nrow=2))

rm(movement_cost_plot_1,movement_cost_plot_2,
   movement_cost_plot_3,movement_cost_plot_4)
##Plot optimal vaccination rates for different parameter combinations-----------------
(movement_sens_vacc_plot<-ggplot(movement_sens_vacc,aes(x=time,y=vacc_rate,color=patch))+
    geom_line()+
    xlab("Days Since Outbreak Began")+
    ylab("Optimal Vaccination Rate")+
    scale_color_discrete(name="Patch",breaks=c("v1","v2"),labels=c("1","2"))+
    facet_grid(m1~m2,labeller = "label_both")+
    theme_bw())

##Select specific simulations to plot--------------------------
#get array of simulations that match our choices
getThese<-c(which(m_grid$m1==-1&m_grid$m2==-1),
            which(m_grid$m1==-5&m_grid$m2==-3),
            which(m_grid$m1==-3&m_grid$m2==-5),
            which(m_grid$m1==-1&m_grid$m2==-5),
            which(m_grid$m1==-5&m_grid$m2==-1))
#Get particular simulations, add in information on parameter values, 
#and join as onedataframe
#OC Simulations
simulations_oc<-lapply(getThese,function(x){
  temp<-movement_sens[[x]]$SolX%>%
    data.frame(m1=m_grid[x,"m1"],m2=m_grid[x,"m2"])
  return(temp)
})%>%
  bind_rows()%>%
  mutate(simulation="OC")

#No Control Simulations
simulations_nc<-lapply(getThese,function(x){
  temp<-outbreak_sens[[x]]$SolX%>%
    data.frame(m1=m_grid[x,"m1"],m2=m_grid[x,"m2"])
  return(temp)
})%>%
  bind_rows()%>%
  mutate(simulation="NC",v1=0,v2=0)

simulations<-bind_rows(simulations_nc,
                       simulations_oc)
rm(simulations_oc,simulations_nc)

simulations_for_plots<-simulations%>%
  select("time","I1","I2","H1","H2","D1","D2","v1","v2","m1","m2","simulation")%>%
  pivot_longer(cols=c("I1","I2","H1","H2","D1","D2","v1","v2"),
               names_to="State",
               values_to = "Value")%>%
  mutate(Class=substr(State,1,1),Patch=paste("Patch",substr(State,2,2),simulation),
         Situation=paste0("m1=10^",m1,", m2=10^",m2))
simulations_for_plots$Situation<-factor(simulations_for_plots$Situation)
simulations_for_plots$Class<-simulations_for_plots$Class%>%
  factor(levels=c("I","H","D","v"),
         labels=c("Infectious","Hospitalized","Infectious Dead","Vaccination Rate"))
#plot the selected simulations
simulation_plot_1<-ggplot(filter(simulations_for_plots,
                                 Class!="Vaccination Rate"),
                          aes(x=time,y=Value,color=Class))+
  geom_line(size=1.1)+
  scale_y_continuous(breaks=c(0,20000,40000,60000),limits = c(0,60000))+
  ylab("Number of Individuals")+
  xlab("Day")+
  facet_grid(Patch~Situation,scales = "free_y")+
  theme_bw()

simulation_plot_2<-ggplot(filter(simulations_for_plots,
                                 Class=="Vaccination Rate",simulation=="OC"),
                          aes(x=time,y=Value,color=Patch))+
  geom_line(size=1.1)+
  facet_grid(1~Situation,scales = "free_y")+
  theme_bw()+
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank())
(simulation_plot<-ggarrange(simulation_plot_1,
                            simulation_plot_2,
                            ncol=1,nrow=2))


##Save Figures-----------------------------------------------
ggsave(file="Movement_cost.tiff",plot=movement_cost_plot)
ggsave(file="Movement_vacc.tiff",plot=movement_sens_vacc_plot)
ggsave(file="Oubreak_size.tiff",plot=outbreak_plot)
