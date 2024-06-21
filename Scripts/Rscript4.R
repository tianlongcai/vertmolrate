rm(list=ls())
gc()
#define work direction
workdir <- "/Users/tianlong/Desktop/VertMolRate"
setwd(workdir)

library(tidyverse)
library(MCMCglmm)
library(MuMIn)
library(patchwork)

sisters <- read.csv("./DataFiles/sisters_family/sisters.csv", row.names = 1)

groups <- c(unique(sisters$Group), unique(sisters$ThermoMode))

ols_list <- vector("list", length = length(groups) * 2)
names(ols_list) <- paste(rep(groups, each = 2), c("dN", "dS"), sep = ".")


for(group in groups){
  if(group %in% unique(sisters$ThermoMode)){
    subdata <- sisters%>% 
      filter(ThermoMode == group)
  }else{
    subdata <- sisters %>% 
      filter(Group == group)
  }
  
  fit_dn <- lm(dN~diff.sp+0, subdata)
  fit_ds <- lm(dS~diff.sp+0, subdata)
  
  ols_list[[paste(group, "dS", sep = ".")]] <- fit_ds
  ols_list[[paste(group, "dN", sep = ".")]] <- fit_dn
}

ols_out <- NULL
for(i in 1:length(groups)){
  group <- groups[i]
  fit.ds.out<- summary(ols_list[[paste(group, "dS", sep = ".")]])$coefficients[,4]
  fit.dn.out<- summary(ols_list[[paste(group, "dN", sep = ".")]])$coefficients[,4]
  
  fit.out <- data.frame(MolRate=c("dS", "dN"), p =c(fit.ds.out, fit.dn.out))%>%
    mutate(Group=group, ThermoMode=ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"))%>%
    select(ThermoMode, Group, MolRate, p)
  ols_out <- rbind(ols_out, fit.out)
}

ols_out

f1 <- sisters%>%
  filter(Group=="Fishes")%>%
  ggplot(aes(x=diff.sp, y=dS, colour=Group))+
  geom_point(size=1)+
  geom_smooth(method="lm", formula = y~x+0)+
  geom_errorbar(aes(ymin=dS-dS.SD, ymax=dS+dS.SD), width=.2, linewidth=0.2)+
  scale_color_manual(values = c("#2a6aaf"))+
  scale_x_continuous(limits = c(0,6), breaks = seq(0,6,2))+
  scale_y_continuous(limits = c(-2,2), breaks = seq(-2,2,1))+
  labs(x="Diff. in ln(Clade Size)", y="Diff. in ln(dS)", title="Fishes")+
  annotate("text", x=4,y=2, label="p = 0.94", size=2.5)+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f2 <- sisters%>%
  filter(Group=="Amphibians")%>%
  ggplot(aes(x=diff.sp, y=dS, colour=Group))+
  geom_point(size=1)+
  geom_smooth(method="lm", formula = y~x+0)+
  geom_errorbar(aes(ymin=dS-dS.SD, ymax=dS+dS.SD), width=.2, linewidth=0.2)+
  scale_color_manual(values = c("#2a6aaf"))+
  scale_x_continuous(limits = c(0,6), breaks = seq(0,6,2))+
  scale_y_continuous(limits = c(-1,0.5), breaks = seq(-1,0.5,0.5))+
  labs(x="Diff. in ln(Clade Size)", y="", title="Amphibians")+
  annotate("text", x=4,y=0.5, label="p = 0.42", size=2.5)+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))
f3 <- sisters%>%
  filter(Group=="Reptiles")%>%
  ggplot(aes(x=diff.sp, y=dS, colour=Group))+
  geom_point(size=1)+
  geom_smooth(method="lm", formula = y~x+0)+
  geom_errorbar(aes(ymin=dS-dS.SD, ymax=dS+dS.SD), width=.2, linewidth=0.2)+
  scale_color_manual(values = c("#2a6aaf"))+
  scale_x_continuous(limits = c(0,4), breaks = seq(0,4,1))+
  scale_y_continuous(limits = c(-1.5,1), breaks = seq(-1.5,1,0.5))+
  labs(x="Diff. in ln(Clade Size)", y="", title="Reptiles")+
  annotate("text", x=3,y=1, label="p = 0.08", size=2.5)+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f4 <- sisters%>%
  filter(Group=="Mammals")%>%
  ggplot(aes(x=diff.sp, y=dS, colour=Group))+
  geom_point(size=1)+
  geom_smooth(method="lm", formula = y~x+0)+
  geom_errorbar(aes(ymin=dS-dS.SD, ymax=dS+dS.SD), width=.2, linewidth=0.2)+
  scale_color_manual(values = c("#d3292f"))+
  scale_x_continuous(limits = c(0,4), breaks = seq(0,4,1))+
  scale_y_continuous(limits = c(-1,1), breaks = seq(-1,1,0.5))+
  labs(x="Diff. in ln(Clade Size)", y="", title="Mammals")+
  annotate("text", x=3,y=1, label="p = 0.99", size=2.5)+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))
f5 <- sisters%>%
  filter(Group=="Birds")%>%
  ggplot(aes(x=diff.sp, y=dS, colour=Group))+
  geom_point(size=1)+
  geom_smooth(method="lm", formula = y~x+0)+
  geom_errorbar(aes(ymin=dS-dS.SD, ymax=dS+dS.SD), width=.2, linewidth=0.2)+
  scale_color_manual(values = c("#d3292f"))+
  scale_x_continuous(limits = c(0,6), breaks = seq(0,6,2))+
  scale_y_continuous(limits = c(-0.8,1.2), breaks = seq(-0.8,1.2,0.4))+
  labs(x="Diff. in ln(Clade Size)", y="", title="Birds")+
  annotate("text", x=4,y=1.2, label="p = 0.043", size=2.5)+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f6 <- sisters%>%
  filter(Group=="Fishes")%>%
  ggplot(aes(x=diff.sp, y=dN, colour=Group))+
  geom_point(size=1)+
  geom_smooth(method="lm", formula = y~x+0)+
  geom_errorbar(aes(ymin=dN-dN.SD, ymax=dN+dN.SD), width=.2, linewidth=0.2)+
  scale_color_manual(values = c("#2a6aaf"))+
  scale_x_continuous(limits = c(0,6), breaks = seq(0,6,2))+
  scale_y_continuous(limits = c(-2,2), breaks = seq(-2,2,1))+
  labs(x="Diff. in ln(Clade Size)", y="Diff. in ln(dN)", title="Fishes")+
  annotate("text", x=4,y=2, label="p = 0.70", size=2.5)+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f7 <- sisters%>%
  filter(Group=="Amphibians")%>%
  ggplot(aes(x=diff.sp, y=dN, colour=Group))+
  geom_point(size=1)+
  geom_smooth(method="lm", formula = y~x+0)+
  geom_errorbar(aes(ymin=dN-dN.SD, ymax=dN+dN.SD), width=.2, linewidth=0.2)+
  scale_color_manual(values = c("#2a6aaf"))+
  scale_x_continuous(limits = c(0,6), breaks = seq(0,6,2))+
  scale_y_continuous(limits = c(-1,0.5), breaks = seq(-1,0.5,0.5))+
  labs(x="Diff. in ln(Clade Size)", y="", title="Amphibians")+
  annotate("text", x=4,y=0.5, label="p = 0.39", size=2.5)+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f8 <- sisters%>%
  filter(Group=="Reptiles")%>%
  ggplot(aes(x=diff.sp, y=dN, colour=Group))+
  geom_point(size=1)+
  geom_smooth(method="lm", formula = y~x+0)+
  geom_errorbar(aes(ymin=dN-dN.SD, ymax=dN+dN.SD), width=.2, linewidth=0.2)+
  scale_color_manual(values = c("#2a6aaf"))+
  scale_x_continuous(limits = c(0,4), breaks = seq(0,4,1))+
  scale_y_continuous(limits = c(-0.5,1), breaks = seq(-0.5,1,0.5))+
  labs(x="Diff. in ln(Clade Size)", y="", title="Reptiles")+
  annotate("text", x=3,y=1, label="p = 0.60", size=2.5)+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f9 <- sisters%>%
  filter(Group=="Mammals")%>%
  ggplot(aes(x=diff.sp, y=dN, colour=Group))+
  geom_point(size=1)+
  geom_smooth(method="lm", formula = y~x+0)+
  geom_errorbar(aes(ymin=dN-dN.SD, ymax=dN+dN.SD), width=.2, linewidth=0.2)+
  scale_color_manual(values = c("#d3292f"))+
  scale_x_continuous(limits = c(0,4), breaks = seq(0,4,1))+
  scale_y_continuous(limits = c(-1.5,1.5), breaks = seq(-1.5,1.5,0.5))+
  labs(x="Diff. in ln(Clade Size)", y="", title="Mammals")+
  annotate("text", x=3,y=1.5, label="p = 0.45", size=2.5)+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))
f10 <- sisters%>%
  filter(Group=="Birds")%>%
  ggplot(aes(x=diff.sp, y=dN, colour=Group))+
  geom_point(size=1)+
  geom_smooth(method="lm", formula = y~x+0)+
  geom_errorbar(aes(ymin=dN-dN.SD, ymax=dN+dN.SD), width=.2, linewidth=0.2)+
  scale_color_manual(values = c("#d3292f"))+
  scale_x_continuous(limits = c(0,6), breaks = seq(0,6,2))+
  scale_y_continuous(limits = c(-0.8,1.2), breaks = seq(-0.8,1.2,0.4))+
  labs(x="Diff. in ln(Clade Size)", y="", title="Birds")+
  annotate("text", x=4,y=1.2, label="p = 0.39", size=2.5)+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+
  plot_layout(ncol = 5, nrow = 2)

ggsave(filename="./Outputs/MainFigures/Fig4.pdf", height=3.8, width=8.27)


#################
rm(list=ls())
gc()
library(caper)
library(phytools)
workdir <- "/Users/tianlong/Desktop/VertMolRate"
source(paste0(workdir, "/SourceFunctions/source_functions.r"))

######
subrate <- read.csv(paste0(workdir, "/DataFiles/MolEvolRate/VertMolRate.csv"))
phy <- read.tree(paste0(workdir, "/DataFiles/trees/phy_all_sampled_vertebrates.tre"))
phylo <- force.ultrametric(phy, method="nnls")#nnls or extend

groups <- unique(subrate$Group)
results <-vector("list", length = length(groups) * 2)
names(results) <- paste(rep(groups, each = 2), c("dN", "dS"), sep = ".")
prior <- list(R=list(V=diag(1),nu=0.002), 
              G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*1000)))

for(group in groups){
  data1 <- subrate%>%filter(Group == group)%>%
    dplyr::select(Species, dN, dS, Div.Rate:NTR.DNA)%>%
    filter(!is.na(NTR))
  data2 <- subrate%>%filter(Group == group)%>%
    dplyr::select(Species, dN, dS, Div.Rate:NTR.DNA)%>%
    filter(!is.na(NTR.DNA))
  
  phy <- drop.tip(phylo, setdiff(phylo$tip.label, data1$Species))
  
  
  fit1 <- MCMCglmm(log(ClaDS.Rate)~log(dS), 
                   random=~Species, prior=prior,
                   ginverse = list(Species=inverseA(phy)$Ainv),
                   data=data1,verbose=TRUE, 
                   nitt=11000, burnin=1000, thin=10, family = c("gaussian"))
  fit2 <- MCMCglmm(log(ClaDS.Rate.DNA)~log(dS), 
                   random=~Species, prior=prior,
                   ginverse = list(Species=inverseA(phy)$Ainv),
                   data=data2,verbose=TRUE, 
                   nitt=11000, burnin=1000, thin=10, family = c("gaussian"))
  fit3 <- MCMCglmm(log(NTR)~log(dS), 
                   random=~Species, prior=prior,
                   ginverse = list(Species=inverseA(phy)$Ainv),
                   data=data1,verbose=TRUE, 
                   nitt=11000, burnin=1000, thin=10, family = c("gaussian"))
  fit4 <- MCMCglmm(log(NTR.DNA)~log(dS), 
                   random=~Species, prior=prior,
                   ginverse = list(Species=inverseA(phy)$Ainv),
                   data=data2,verbose=TRUE, 
                   nitt=11000, burnin=1000, thin=10, family = c("gaussian"))
  fit5 <- MCMCglmm(log(ClaDS.Rate)~log(dN), 
                   random=~Species, prior=prior,
                   ginverse = list(Species=inverseA(phy)$Ainv),
                   data=data1,verbose=TRUE, 
                   nitt=11000, burnin=1000, thin=10, family = c("gaussian"))
  fit6 <- MCMCglmm(log(ClaDS.Rate.DNA)~log(dN), 
                   random=~Species, prior=prior,
                   ginverse = list(Species=inverseA(phy)$Ainv),
                   data=data2,verbose=TRUE, 
                   nitt=11000, burnin=1000, thin=10, family = c("gaussian"))
  fit7 <- MCMCglmm(log(NTR)~log(dN), 
                   random=~Species, prior=prior,
                   ginverse = list(Species=inverseA(phy)$Ainv),
                   data=data1,verbose=TRUE, 
                   nitt=11000, burnin=1000, thin=10, family = c("gaussian"))
  fit8 <- MCMCglmm(log(NTR.DNA)~log(dN), 
                   random=~Species, prior=prior,
                   ginverse = list(Species=inverseA(phy)$Ainv),
                   data=data2,verbose=TRUE, 
                   nitt=11000, burnin=1000, thin=10, family = c("gaussian"))
  
  
  results[[paste(group, "dS", sep = ".")]]  <- data.frame(Group=group, ClaDS.Slope=summary(fit1)$solutions[2,1], ClaDS.P=summary(fit1)$solutions[2,5], 
                                                          ClaDS.DNA.Slope=summary(fit2)$solutions[2,1], ClaDS.DNA.P=summary(fit2)$solutions[2,5], 
                                                          Nnode.Slope=summary(fit3)$solutions[2,1], Nnode.P=summary(fit3)$solutions[2,5], 
                                                          Nnode.DNA.Slope=summary(fit4)$solutions[2,1], Nnode.DNA.P=summary(fit4)$solutions[2,5])
  
  results[[paste(group, "dN", sep = ".")]]  <- data.frame(Group=group, ClaDS.Slope=summary(fit5)$solutions[2,1], ClaDS.P=summary(fit5)$solutions[2,5], 
                                                          ClaDS.DNA.Slope=summary(fit6)$solutions[2,1], ClaDS.DNA.P=summary(fit6)$solutions[2,5], 
                                                          Nnode.Slope=summary(fit7)$solutions[2,1], Nnode.P=summary(fit7)$solutions[2,5], 
                                                          Nnode.DNA.Slope=summary(fit8)$solutions[2,1], Nnode.DNA.P=summary(fit8)$solutions[2,5])
  print(group)
}
results%>%bind_rows()%>%mutate(MolRate=rep(c("dN", "dS"), 5))

