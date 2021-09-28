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
m_values<-seq(-5,1,by=1)
m_grid<-expand.grid(m1=m_values,m2=m_values)  #Expands 2 (or more) vectors into a grid

#set intial conditions (initial outbreak) - 10 I in patch 1
initial=c(S1=parm$N1-10,E1=0,I1=10,H1=0,D1=0,R1=0,S2=parm$N2,E2=0,I2=0,H2=0,D2=0,R2=0)

parm_temp<-parm

##Solve Optimal Control Problem for each parameter combination in m_grid--------------
movement_sens<-lapply(1:dim(m_grid)[1],function(x){
  parm_temp$m1<-10^(m_grid[x,"m1"]) #overwrite variables of interest in the temp
  parm_temp$m2<-10^(m_grid[x,"m2"]) #parameter set
  temp<-ebola.optim(inits=initial,M=c(.001,.001),params=parm_temp,times=seq(0,730,by=1))
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

##Plot optimal vaccination rates for different parameter combinations-----------------
(movement_sens_vacc_plot<-ggplot(movement_sens_vacc,aes(x=time,y=vacc_rate,color=patch))+
    geom_line()+
    xlab("Days Since Outbreak Began")+
    ylab("Optimal Vaccination Rate")+
    scale_color_discrete(name="Patch",breaks=c("v1","v2"),labels=c("1","2"))+
    facet_grid(m1~m2,labeller = "label_both")+
    theme_bw())
ggsave(file="Movement_cost.tiff",plot=movement_cost_plot)
ggsave(file="Movement_vacc.tiff",plot=movement_sens_vacc_plot)
