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
m_values<-seq(-5,-3,by=1)
m_grid<-expand.grid(m1=m_values,m2=m_values)  #Expands 2 (or more) vectors into a grid
m_grid<-m_grid%>%filter(abs(m1-m2)<2)
#Set maximum control levels M1,M2
maxControl<-c(.015,.015)

#Set timeframe
times<-seq(0,365,by=1)
preTimes<-seq(0,100,by=1)

#set intial conditions (initial outbreak) - 10 I in patch 1
initial=c(S1=parm$N1-1,E1=0,I1=1,H1=0,D1=0,R1=0,S2=parm$N2,E2=0,I2=0,H2=0,D2=0,R2=0)

v1 <- data.frame(times = times, v1 = rep(0,length(times)))
v1_interp <- approxfun(v1, rule = 2)
v1 = v1[,2]
v2 = data.frame(times = times, v2 = rep(0,length(times)))
v2_interp <- approxfun(v2, rule = 2)
v2 = v2[,2]
parm_temp<-c(parm,v1_interp=v1_interp,v2_interp=v2_interp)
##Simulate out uncontrolled outbreak------------------------
uncontrolled<-lapply(1:dim(m_grid)[1],function(x){
  parm_temp$m1<-10^(m_grid[x,"m1"]) #overwrite variables of interest in the temp
  parm_temp$m2<-10^(m_grid[x,"m2"]) #parameter set
  prePeriod<-ode(y=initial,
                 times=preTimes,
                 func=ebola.ode,
                 parms=parm_temp,
                 method = "ode45")
  temp<-ode(y=prePeriod[dim(prePeriod)[1],2:dim(prePeriod)[2]],
            times=times,
            func=ebola.ode,
            parms=parm_temp,
            method = "ode45")%>%
    data.frame()
  size<-list(SolX=temp,
             `Patch 1`=with(c(parm,temp),sum(b1*(betaI1*S1*I1 + betaD1*S1*D1)*(times[2]-times[1]))),
             `Patch 2`=with(c(parm,temp),sum(b2*(betaI2*S2*I2 + betaD2*S2*D2)*(times[2]-times[1]))))
  print(x)
  return(size)
})
#Extract outbreak size information
uncontrolled_cost<-lapply(1:dim(m_grid)[1],function(x){
  temp<-data.frame("Patch 1"=uncontrolled[[x]]$`Patch 1`,
                   "Patch 2"=uncontrolled[[x]]$`Patch 2`)
  return(temp)
})%>%bind_rows()%>%
  bind_cols(m_grid)%>%
  mutate(Total=Patch.1+Patch.2)%>%
  pivot_longer(cols=c("Patch.1","Patch.2","Total"),
               names_to = "Type",values_to = "Cost")%>%
  mutate(m1=10^m1,m2=10^m2,control_type="Uncontrolled")
uncontrolled_cost$Type<-uncontrolled_cost$Type%>%
  factor(levels=c("Patch.1","Patch.2","Total"),
         labels=c("Patch 1 Disease Cost","Patch 2 Disease Cost","Total Cost"))

##Solve Optimal Control Problem for unique controls--------------
parm_temp<-parm
unique_control<-lapply(1:dim(m_grid)[1],function(x){
  parm_temp$m1<-10^(m_grid[x,"m1"]) #overwrite variables of interest in the temp
  parm_temp$m2<-10^(m_grid[x,"m2"]) #parameter set
  prePeriod<-ode(y=initial,
            times=preTimes,
            func=ebola.ode,
            parms=parm_temp,
            method = "ode45")
  temp<-ebola.optim(inits=prePeriod[dim(prePeriod)[1],2:dim(prePeriod)[2]],
                    M=maxControl,params=parm_temp,
                    times=times,maxIter = 100,strictConv = F)
  return(temp)
})

#Extract optimal Vaccination for unique plans------------------
unique_control_vacc<-lapply(1:dim(m_grid)[1],function(x){
  temp<-data.frame(m1=m_grid[x,"m1"],m2=m_grid[x,"m2"],
                   unique_control[[x]]$SolX)
  return(temp)
})%>%bind_rows()%>%
  select(m1,m2,time,v1,v2)%>%
  pivot_longer(cols=c("v1","v2"),names_to = "patch",values_to = "vacc_rate")%>%
  mutate(m1=10^m1,m2=10^m2,control_type="Unique")


#Extract cost values---------------------
unique_control_costs<-lapply(1:dim(m_grid)[1],function(x){
  temp<-data.frame(m1=m_grid[x,"m1"],m2=m_grid[x,"m2"],
                   J11=unique_control[[x]]$J11,
                   J12=unique_control[[x]]$J12,
                   J21=unique_control[[x]]$J21,
                   J22=unique_control[[x]]$J22)
  return(temp)
})%>%bind_rows()%>%
  mutate(Total=J11+J12+J21+J22,Total_1=J11+J21,Total_2=J12+J22)%>%
  pivot_longer(cols=c("J11","J12","J21","J22","Total","Total_1","Total_2"),names_to = "Type",values_to = "Cost")%>%
  mutate(m1=10^m1,m2=10^m2,control_type="Unique")

#Rename costs to be interpretable
unique_control_costs$Type<-unique_control_costs$Type%>%
  factor(levels=c("J11","J12","J21","J22","Total_1","Total_2","Total"),
         labels=c("Patch 1 Disease Cost","Patch 2 Disease Cost",
                  "Patch 1 Vaccine Cost","Patch 2 Vaccine Cost",
                  "Patch 1 Total Cost","Patch 2 Total Cost",
                  "Total Cost"))

##Solve Optimal Control Problem for same controls--------------
parm_temp<-parm
same_control<-lapply(1:dim(m_grid)[1],function(x){
  parm_temp$m1<-10^(m_grid[x,"m1"]) #overwrite variables of interest in the temp
  parm_temp$m2<-10^(m_grid[x,"m2"]) #parameter set
  prePeriod<-ode(y=initial,
                 times=preTimes,
                 func=ebola.ode,
                 parms=parm_temp,
                 method = "ode45")
  temp<-ebola.optim(inits=prePeriod[dim(prePeriod)[1],2:dim(prePeriod)[2]],
                    M=maxControl,params=parm_temp,
                    times=times,maxIter = 100,strictConv = F,sameRate = T)
  return(temp)
})


#Extract optimal Vaccination for same plans------------------
same_control_vacc<-lapply(1:dim(m_grid)[1],function(x){
  temp<-data.frame(m1=m_grid[x,"m1"],m2=m_grid[x,"m2"],
                   same_control[[x]]$SolX)
  return(temp)
})%>%bind_rows()%>%
  select(m1,m2,time,v1,v2)%>%
  pivot_longer(cols=c("v1","v2"),names_to = "patch",values_to = "vacc_rate")%>%
  mutate(m1=10^m1,m2=10^m2,control_type="Same")

#Extract cost values---------------------
same_control_costs<-lapply(1:dim(m_grid)[1],function(x){
  temp<-data.frame(m1=m_grid[x,"m1"],m2=m_grid[x,"m2"],
                   J11=same_control[[x]]$J11,
                   J12=same_control[[x]]$J12,
                   J21=same_control[[x]]$J21,
                   J22=same_control[[x]]$J22)
  return(temp)
})%>%bind_rows()%>%
  mutate(Total=J11+J12+J21+J22,Total_1=J11+J21,Total_2=J12+J22)%>%
  pivot_longer(cols=c("J11","J12","J21","J22","Total","Total_1","Total_2"),names_to = "Type",values_to = "Cost")%>%
  mutate(m1=10^m1,m2=10^m2,control_type="Same")

#Rename costs to be interpretable
same_control_costs$Type<-same_control_costs$Type%>%
  factor(levels=c("J11","J12","J21","J22","Total_1","Total_2","Total"),
         labels=c("Patch 1 Disease Cost","Patch 2 Disease Cost",
                  "Patch 1 Vaccine Cost","Patch 2 Vaccine Cost",
                  "Patch 1 Total Cost","Patch 2 Total Cost",
                  "Total Cost"))

##Combine Data tables------------------------------
vacc_rates<-bind_rows(same_control_vacc,unique_control_vacc)
costs<-bind_rows(uncontrolled_cost,same_control_costs,unique_control_costs)
costs$control_type<-costs$control_type%>%
  factor(levels=c("Uncontrolled","Same","Unique"))
vacc_rates$m1<-vacc_rates$m1%>%
  factor(ordered = T)
vacc_rates$m2<-vacc_rates$m2%>%
  factor(ordered = T)
costs$m1<-costs$m1%>%
  factor(ordered = T)
costs$m2<-costs$m2%>%
  factor(ordered = T)


#Plot figure of vaccination rates
rates_plot<-ggplot(vacc_rates,aes(x=time,y=vacc_rate,color=patch,linetype=control_type))+
  geom_line(size=1.25)+
  facet_grid(m1~m2,labeller = "label_both")+
  scale_linetype_discrete(name="Control")+
  scale_color_discrete(name="Patch",labels=c("1","2"))+
  ylab("Vaccination Rate")+
  xlab("Day")+
  theme_bw()+
  guides(linetype=guide_legend(nrow=2),
         color=guide_legend(nrow=2))+
  theme(legend.position = "bottom",
        legend.text=element_text(size=8),
        legend.title = element_text(size=10),
        axis.title = element_text(size=10),
        axis.text = element_text(size=8),
        strip.text = element_text(size=8))

#Plot Cost Figure
costs_plot<-ggplot(filter(costs,!Type%in%c("Total Cost","Patch 1 Total Cost","Patch 2 Total Cost")),
       aes(x=control_type,y=Cost,fill=Type))+
  geom_col()+
  scale_x_discrete(breaks=c("Uncontrolled","Same","Unique"),labels=c("None","Same","Unique"))+
  facet_grid(m1~m2,labeller = "label_both")+
  ylab("Cost")+
  xlab(NULL)+
  scale_fill_discrete(name=NULL,guide=guide_legend(nrow=2))+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.text=element_text(size=8),
        legend.title = element_text(size=10),
        axis.title = element_text(size=10),
        axis.text.y = element_text(size=8),
        axis.text.x=element_text(size=8,angle=45,hjust=1),
        strip.text = element_text(size=8))

#Plot Relative Costs
costs_relative<-costs%>%
  filter(Type%in%c("Total Cost"))%>%
  pivot_wider(names_from = control_type,values_from = Cost)%>%
  mutate(Unique_rel=Unique/Uncontrolled,Same_rel=Same/Uncontrolled,.keep="unused")%>%
  pivot_longer(cols=c("Unique_rel","Same_rel"),names_to = "Cost_Type",values_to = "Cost")%>%
  mutate(movement=paste("m1:",m1,"\n","m2:",m2))

rel_costs_plot<-ggplot(costs_relative,aes(x=movement,y=Cost,color=Cost_Type,fill=Cost_Type))+
  geom_col(position="dodge")+
  scale_color_discrete(name="Control",breaks=c("Same_rel","Unique_rel"),labels=c("Same","Unique"))+
  scale_fill_discrete(name="Control",breaks=c("Same_rel","Unique_rel"),labels=c("Same","Unique"))+
  xlab(NULL)+
  ylab("Cost Relative to No Control")+
  theme_bw()+
  theme(legend.position = "bottom",
      legend.text=element_text(size=8),
      legend.title = element_text(size=10),
      axis.title = element_text(size=10),
      axis.text.y = element_text(size=8),
      axis.text.x = element_text(size=8,angle=45,hjust=1),
      strip.text = element_text(size=8))

full_figure<-ggarrange(rates_plot,costs_plot,rel_costs_plot,nrow=1)
full_figure
ggsave("Example_Costs.tiff",plot=full_figure)
