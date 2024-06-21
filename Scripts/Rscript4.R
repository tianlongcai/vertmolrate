################################
# R script to examine relationships between diversification rates and molecular rates
# Author: Tianlong Cai
# Email: caitianlong@westlake.edu.cn

#####################
rm(list=ls())
gc()
#define work direction
workdir <- "/Users/tianlong/VertMolRate"
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




