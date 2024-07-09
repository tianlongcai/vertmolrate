################################
# An R script to predict molecular rates using multiple PGLMMs
# Author: Tianlong Cai
# Email: caitianlong@westlake.edu.cn

#####################
rm(list=ls())
gc()
#define work direction
workdir <- "/Users/tianlong/VertMolRate"
setwd(workdir)

library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(phytools)
library(MCMCglmm)
library(MuMIn)
library(geiger)
library(patchwork)
library(corrplot)


source(paste0(workdir, "/Scripts/source_functions.r"))

######
subrate <- read.csv(paste0(workdir, "/DataFiles/MolEvolRate/VertMolRate.csv"))%>%
  mutate(AnnualTemp=ifelse(ThermoMode=="Ectotherms"&AnnualTemp<0.01, 0.01, AnnualTemp))


row.names(subrate) <- subrate$Species

#input phylogeny
phy <- read.tree(paste0(workdir, "/DataFiles/trees/phy_all_sampled_vertebrates.tre"))
phylo <- force.ultrametric(phy, method="nnls")#nnls or extend


#get most recent common ancestor (node) of species pair
Ainv <- inverseA(phylo)$Ainv


####################################
#1. Correlations between molecular rates and life-history traits
groups <- c("Fishes","Amphibians","Reptiles","Mammals","Birds")

pdf("./Outputs/Supplementary/Extended Fig8x.pdf", width=9.27, height = 2.8)
layout(matrix(1:10, nrow = 2, ncol = 5, byrow = FALSE))

for(group in groups){
  mydata<-subrate%>%
    filter(Group==group)%>%
    mutate(dS=log(dS), dN=log(dN),BodyMass=log(BodyMass), Fecundity=log(Fecundity), 
           MaturityAge=log(MaturityAge), Longevity=log(Longevity))%>%
    select(dS, BodyMass, Fecundity, MaturityAge, Longevity)
  corr_matrix <- cor(mydata)
  corrplot(corr_matrix, method = "circle", 
           type="lower",
           col = colorRampPalette(brewer.pal(11, "Spectral"))(200),
           addCoef.col = "black", # 
           tl.col = "black", tl.srt = 25, # 
           number.cex = 0.4, cl.cex = 0.5, # 
           tl.cex =0.6,
           cex.main = 0.6,
           title = group,
           cl.length = 5,
           mar = c(1, 1, 1, 1)) # 
  
  mydata<-subrate%>%
    filter(Group==group)%>%
    mutate(dS=log(dS), dN=log(dN),BodyMass=log(BodyMass), Fecundity=log(Fecundity), 
           MaturityAge=log(MaturityAge), Longevity=log(Longevity))%>%
    select(dN, dS, BodyMass, Fecundity, MaturityAge, Longevity)
  corr_matrix <- cor(mydata)
  corrplot(corr_matrix, method = "circle", 
           type="lower",
           col = colorRampPalette(brewer.pal(11, "Spectral"))(200),
           addCoef.col = "black", # 
           tl.col = "black", tl.srt = 25, # 
           number.cex = 0.4, cl.cex = 0.5, # 
           tl.cex =0.6,
           cex.main = 0.6,
           title = group,
           cl.length = 5,
           mar = c(1, 1, 1, 1)) # 
  
}
dev.off()


######################
#2. Model selection
#PLMM with tips data
scale_subrate <- subrate %>%
  mutate(dS=scale(log(dS)), dN=scale(log(dN)),AnnualTemp=scale(AnnualTemp),
         BodyMass=scale(log(BodyMass)), Fecundity=scale(log(Fecundity)), 
         MaturityAge=scale(log(MaturityAge)), Longevity=scale(log(Longevity)))


#################################
prior <- list(R=list(V=diag(1),nu=0.002), 
               G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*1000)))


full.group.interact.ds<- MCMCglmm(dS~AnnualTemp*Group+BodyMass*Group+Fecundity*Group+MaturityAge*Group+Longevity*Group, 
                                  random=~Species, prior = prior, ginverse=list(Species=Ainv), 
                                  data=scale_subrate, verbose=TRUE,
                                  nitt=60000, burnin=10000, thin=25, family = c("gaussian"))

full.thermo.interact.ds<- MCMCglmm(dS~AnnualTemp*ThermoMode+BodyMass*ThermoMode+Fecundity*ThermoMode+MaturityAge*ThermoMode+Longevity*ThermoMode, 
                                   random=~Species, prior = prior, ginverse=list(Species=Ainv), 
                                   data=scale_subrate, verbose=TRUE,
                                   nitt=60000, burnin=10000, thin=25, family = c("gaussian"))

full.no.interact.ds <- MCMCglmm(dS~AnnualTemp+BodyMass+Fecundity+MaturityAge+Longevity, 
                                random=~Species, prior = prior, ginverse=list(Species=Ainv),
                                data=scale_subrate, verbose=TRUE,
                                nitt=60000, burnin=10000, thin=25, family = c("gaussian"))

fixed.only.group.interact.ds<- MCMCglmm(dS~AnnualTemp*Group+BodyMass*Group+Fecundity*Group+MaturityAge*Group+Longevity*Group,
                                        data=scale_subrate, verbose=TRUE,
                                        nitt=60000, burnin=10000, thin=25, family = c("gaussian"))

fixed.only.thermo.interact.ds<- MCMCglmm(dS~AnnualTemp*ThermoMode+BodyMass*ThermoMode+Fecundity*ThermoMode+MaturityAge*ThermoMode+Longevity*ThermoMode,
                                         data=scale_subrate, verbose=TRUE,
                                         nitt=60000, burnin=10000, thin=25, family = c("gaussian"))


fixed.only.no.interact.ds <- MCMCglmm(dS~AnnualTemp+BodyMass+Fecundity+MaturityAge+Longevity, 
                                      data=scale_subrate, verbose=TRUE, 
                                      nitt=60000, burnin=10000, thin=25, family = c("gaussian"))

phylo.only.ds<- MCMCglmm(dS~1,
                         random=~Species, prior = prior, ginverse=list(Species=Ainv),
                         data=scale_subrate, verbose=TRUE, 
                         nitt=60000, burnin=10000, thin=25, family = c("gaussian"))
  
  
full.group.interact.dn<- MCMCglmm(dN~dS*Group+AnnualTemp*Group+BodyMass*Group+Fecundity*Group+MaturityAge*Group+Longevity*Group, 
                                  random=~Species, prior = prior, ginverse=list(Species=Ainv), 
                                  data=scale_subrate, verbose=TRUE,
                                  nitt=60000, burnin=10000, thin=25, family = c("gaussian"))

full.thermo.interact.dn<- MCMCglmm(dN~dS*ThermoMode+AnnualTemp*ThermoMode+BodyMass*ThermoMode+Fecundity*ThermoMode+MaturityAge*ThermoMode+Longevity*ThermoMode, 
                                   random=~Species, prior = prior, ginverse=list(Species=Ainv), 
                                   data=scale_subrate, verbose=TRUE,
                                   nitt=60000, burnin=10000, thin=25, family = c("gaussian"))

full.no.interact.dn <- MCMCglmm(dN~dS+AnnualTemp+BodyMass+Fecundity+MaturityAge+Longevity, 
                                random=~Species, prior = prior, ginverse=list(Species=Ainv),
                                data=scale_subrate, verbose=TRUE,
                                nitt=60000, burnin=10000, thin=25, family = c("gaussian"))

fixed.only.group.interact.dn<- MCMCglmm(dN~dS*Group+AnnualTemp*Group+BodyMass*Group+Fecundity*Group+MaturityAge*Group+Longevity*Group,
                                        data=scale_subrate, verbose=TRUE,
                                        nitt=60000, burnin=10000, thin=25, family = c("gaussian"))

fixed.only.thermo.interact.dn<- MCMCglmm(dN~dS*ThermoMode+AnnualTemp*ThermoMode+BodyMass*ThermoMode+Fecundity*ThermoMode+MaturityAge*ThermoMode+Longevity*ThermoMode,
                                         data=scale_subrate, verbose=TRUE,
                                         nitt=60000, burnin=10000, thin=25, family = c("gaussian"))


fixed.only.no.interact.dn <- MCMCglmm(dN~dS+AnnualTemp+BodyMass+Fecundity+MaturityAge+Longevity, 
                                      data=scale_subrate, verbose=TRUE, 
                                      nitt=60000, burnin=10000, thin=25, family = c("gaussian"))

phylo.only.dn<- MCMCglmm(dN~1,
                         random=~Species, prior = prior, ginverse=list(Species=Ainv),
                         data=scale_subrate, verbose=TRUE, 
                         nitt=60000, burnin=10000, thin=25, family = c("gaussian"))


pglmm_model_sel <-list(full.group.interact.ds=full.group.interact.ds, 
                           full.thermo.interact.ds=full.thermo.interact.ds, 
                           full.no.interact.ds=full.no.interact.ds, 
                           fixed.only.group.interact.ds=fixed.only.group.interact.ds, 
                           fixed.only.thermo.interact.ds=fixed.only.thermo.interact.ds, 
                           fixed.only.no.interact.ds=fixed.only.no.interact.ds, 
                           phylo.only.ds=phylo.only.ds,
                           full.group.interact.dn=full.group.interact.dn,
                           full.thermo.interact.dn=full.thermo.interact.dn,
                           full.no.interact.dn=full.no.interact.dn,
                           fixed.only.group.interact.dn=fixed.only.group.interact.dn, 
                           fixed.only.thermo.interact.dn=fixed.only.thermo.interact.dn,
                           fixed.only.no.interact.dn=fixed.only.no.interact.dn,
                           phylo.only.dn=phylo.only.dn)
  

save(pglmm_model_sel, file="./Outputs/Data/Multiple_PGLMM_Model_Selection.rdata")
load(file="./Outputs/Data/Multiple_PGLMM_Model_Selection.rdata")


ds <- model.sel(full.group.interact.ds, full.thermo.interact.ds, full.no.interact.ds, fixed.only.group.interact.ds, 
                fixed.only.thermo.interact.ds, fixed.only.no.interact.ds, phylo.only.ds,  rank=DIC,
                extra = c(R2m=function(x) r.squaredMCMCglmm(x)$R2m,
                          R2m.sd=function(x) r.squaredMCMCglmm(x)$R2m.sd,
                          R2c=function(x) r.squaredMCMCglmm(x)$R2c,
                          R2c.sd=function(x) r.squaredMCMCglmm(x)$R2c.sd,
                          R2r=function(x) r.squaredMCMCglmm(x)$R2r,
                          R2r.sd=function(x) r.squaredMCMCglmm(x)$R2r.sd))
dn <- model.sel(full.group.interact.dn, full.thermo.interact.dn, full.no.interact.dn, fixed.only.group.interact.dn, 
          fixed.only.thermo.interact.dn, fixed.only.no.interact.dn, phylo.only.dn,  rank=DIC,
          extra = c(R2m=function(x) r.squaredMCMCglmm(x)$R2m,
                    R2m.sd=function(x) r.squaredMCMCglmm(x)$R2m.sd,
                    R2c=function(x) r.squaredMCMCglmm(x)$R2c,
                    R2c.sd=function(x) r.squaredMCMCglmm(x)$R2c.sd,
                    R2r=function(x) r.squaredMCMCglmm(x)$R2r,
                    R2r.sd=function(x) r.squaredMCMCglmm(x)$R2r.sd))


model_sel <- list(ds,dn)
save(model_sel, file="./Outputs/Data/model_sel.rdata")
load(file="./Outputs/Data/model_sel.rdata")
ds <- model_sel[[1]]
dn <- model_sel[[2]]
df1 <- data.frame(
  model = c("Fixed-only", "Fixed-Interact", "Best Model", "Phylo-only"),
  Mol.Rate = rep(c("dS", "dN"),each=4),
  R2c = c(ds["fixed.only.no.interact.ds","R2c"], 
          ds["fixed.only.group.interact.ds","R2c"], 
          ds["full.group.interact.ds","R2c"], 
          ds["phylo.only.ds", "R2c"],
          dn["fixed.only.no.interact.dn","R2c"], 
          dn["fixed.only.group.interact.dn","R2c"], 
          dn["full.group.interact.dn","R2c"],
          dn["phylo.only.dn", "R2c"]), 
  R2c.sd = c(ds["fixed.only.no.interact.ds","R2c.sd"], 
             ds["fixed.only.group.interact.ds","R2c.sd"], 
             ds["full.group.interact.ds","R2c.sd"], 
             ds["phylo.only.ds", "R2c.sd"],
             dn["fixed.only.no.interact.dn","R2c.sd"], 
             dn["fixed.only.group.interact.dn","R2c.sd"], 
             dn["full.group.interact.dn","R2c.sd"],
             dn["phylo.only.dn", "R2c.sd"]))%>%
  mutate(model=factor(model, levels = c("Fixed-only", "Fixed-Interact", "Best Model", "Phylo-only")))
  
df2 <- data.frame(R2=c("R2.Phy", "R2m", "R2.Phy", "R2m"),
                  Mean=c(ds[1,"R2r"], ds[1,"R2m"], dn[1,"R2r"], dn[1,"R2m"]),
                  SD=c(ds[1,"R2r.sd"], ds[1,"R2m.sd"], dn[1,"R2r.sd"], dn[1,"R2m.sd"]),
                  group=c("dS", "dS", "dN", "dN"))%>%
  mutate(group=factor(group, levels=c("dS","dN")))

  
#Fig.3a-b
f3a <- ggplot(df1, aes(model, R2c, group=Mol.Rate, colour = Mol.Rate))+
  geom_line(size=0.4)+
  geom_point(size=1.2)+
  geom_errorbar(aes(ymin = R2c-R2c.sd, ymax = R2c+R2c.sd), width=.3, size=0.4)+
  scale_color_manual(values=c("#93c47d", "#fdb96b"))+
  labs(x="", y=expression(R[c] ^  "2" ), title = "Model selection", tag="a")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0,1))+
  theme_classic()+
  theme(axis.title = element_text(size=8),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=9),
        legend.text = element_text(size=7),
        legend.position = c(0.20, 0.8),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key.height = unit(0.4, 'cm'),
        legend.key.width = unit(0.6, 'cm'),
        axis.text.x = element_text(size=8, angle = 20, vjust = 1, hjust=0.8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))


f3b<-df2%>%
  ggplot(aes(x=group, y=Mean, fill=R2)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar(aes(ymin = Mean-SD, ymax = Mean+SD), width=.3, position=position_dodge(0.9), size=0.2)+
  labs(x="", y=expression("R" ^ 2), title = "Best model", tag="b")+
  scale_fill_manual(labels=expression("R"["phy"]^2,"R"["fixed"]^2), values=c("#7570b3", "#e7298a"))+
  guides(y=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0,1))+
  theme_classic()+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        legend.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=9),
        legend.text.align = 0,
        legend.position = c(0.80, 0.8),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key.height = unit(0.2, 'cm'),
        legend.key.width = unit(0.2, 'cm'),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))



############################################################
#3.PGLMMs for each class, endotherm and ectotherms using all data
#################################################################
#3.1 Using phylogenetic signal as a random effect
pred1 <- c("AnnualTemp", "BodyMass", "Longevity","MaturityAge","Fecundity")
pred2 <- c("AnnualTemp", "BodyMass", "Longevity","MaturityAge","Fecundity","dS")

prior <- list(R=list(V=diag(1),nu=0.002), 
              G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*1000)))

groups <- c(unique(subrate$Group), unique(subrate$ThermoMode))

multi_pglmm_list <- vector("list", length = length(groups) * 2)
names(multi_pglmm_list) <- paste(rep(groups, each = 2), c("dN", "dS"), sep = ".")

# Prior specification for the model
for (group in groups) {
  # Choose data based on group
  if(group %in% unique(subrate$ThermoMode)){
    subdata <- subrate%>% 
      filter(ThermoMode == group)%>%
      mutate(dS=scale(log(dS)), dN=scale(log(dN)),AnnualTemp=scale(AnnualTemp),
             BodyMass=scale(log(BodyMass)), Fecundity=scale(log(Fecundity)), 
             MaturityAge=scale(log(MaturityAge)), Longevity=scale(log(Longevity)))
    
  }else{
    subdata <- subrate %>% 
      filter(Group == group)%>%
      mutate(dS=scale(log(dS)), dN=scale(log(dN)),AnnualTemp=scale(AnnualTemp),
             BodyMass=scale(log(BodyMass)), Fecundity=scale(log(Fecundity)), 
             MaturityAge=scale(log(MaturityAge)), Longevity=scale(log(Longevity)))
  }
  
  
  #Priors
  fit_ds <- MCMCglmm(formula(paste("dS~", paste(pred1, collapse="+"))), 
                     random=~Species, prior=prior,
                     ginverse = list(Species=inverseA(phylo)$Ainv), 
                     data=subdata,verbose=TRUE, 
                     nitt=60000, burnin=10000, thin=25, family = c("gaussian"))
  fit_dn <- MCMCglmm(formula(paste("dN~", paste(pred2, collapse="+"))), 
                     random=~Species, prior=prior,
                     ginverse = list(Species=inverseA(phylo)$Ainv),
                     data=subdata,verbose=TRUE, 
                     nitt=60000, burnin=10000, thin=25, family = c("gaussian"))
  
  # Store models in list
  multi_pglmm_list[[paste(group, "dS", sep = ".")]] <- fit_ds
  multi_pglmm_list[[paste(group, "dN", sep = ".")]] <- fit_dn
  
  print(group)
}

#save(multi_pglmm_list, file="./Outputs/Data/Multiple_PGLMM.rdata")
load(file="./Outputs/Data/Multiple_PGLMM.rdata")




summary_ds <-NULL
summary_dn <-NULL
for(group in groups){
  fit_ds <- multi_pglmm_list[[paste(group, "dS", sep = ".")]]
  fit_dn <- multi_pglmm_list[[paste(group, "dN", sep = ".")]]
  
  sum_ds <- summary(fit_ds)
  sum_dn <- summary(fit_dn)
  
  ds_fixed <- as.data.frame(sum_ds$solutions)
  ds_random <- as.data.frame(sum_ds$Gcovariances)%>%
    mutate(pMCMC=NA)
  ds_r <- r.squaredMCMCglmm(fit_ds)
  
  ds_fixed <- data.frame(Var=row.names(ds_fixed), ds_fixed)
  ds_random <- data.frame(Var=row.names(ds_random), ds_random)
  ds_r <- data.frame(Var=c("Rm2", "Rc2","Rr2"),rbind(as.numeric(ds_r[1:3]), as.numeric(ds_r[5:7]), as.numeric(ds_r[9:11])), NA,NA)
  colnames(ds_r) <- colnames(ds_random)
  
  summary_ds <- rbind(summary_ds, rbind(ds_fixed, ds_random,ds_r)%>%mutate(Group=group))
  
  dn_fixed <- as.data.frame(sum_dn$solutions)
  dn_random <- as.data.frame(sum_dn$Gcovariances)%>%
    mutate(pMCMC=NA)
  dn_r <- r.squaredMCMCglmm(fit_dn)
  
  dn_fixed <- data.frame(Var=row.names(dn_fixed), dn_fixed)
  dn_random <- data.frame(Var=row.names(dn_random), dn_random)
  dn_r <- data.frame(Var=c("Rm2", "Rc2","Rr2"),rbind(as.numeric(dn_r[1:3]), as.numeric(dn_r[5:7]), as.numeric(dn_r[9:11])), NA,NA)
  colnames(dn_r) <- colnames(dn_random)
  
  summary_dn <- rbind(summary_dn, rbind(dn_fixed, dn_random,dn_r)%>%mutate(Group=group))
}

summary_ds 
summary_dn 

#extract parameter values from PGLMMs
pglmm_out <- NULL
for(i in 1:length(groups)){
  group <- groups[i]
  fit.ds.out <- as.data.frame(summary(multi_pglmm_list[[paste(group, "dS", sep = ".")]])$solutions[-1,])
  fit.dn.out <- as.data.frame(summary(multi_pglmm_list[[paste(group, "dN", sep = ".")]])$solutions[-1,])
  fit.ds.out <- fit.ds.out%>%mutate(env=row.names(fit.ds.out), Subrate="dS")
  fit.dn.out <- fit.dn.out%>%mutate(env=row.names(fit.dn.out), Subrate="dN")
  
  fit.out <- rbind(fit.ds.out, fit.dn.out)%>%
    mutate(Group=group, ThermoMode=ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"),
           sig=ifelse(pMCMC<0.05, "1", "0"))
  pglmm_out <- rbind(pglmm_out, fit.out)
}



#Fig.3c-f
f3c <- pglmm_out%>%
  filter(Subrate=="dS")%>%
  filter(Group =="Ectotherms")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="dS", title = "Ectotherms", tag="c")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.45,0.21), breaks = seq(-0.4,0.2,0.2))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=9),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()


f3d <- pglmm_out%>%
  filter(Subrate=="dS")%>%
  filter(Group =="Endotherms")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="dS", title = "Endotherms", tag="d")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.45,0.21), breaks = seq(-0.4,0.2,0.2))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=9),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()



f3e <- pglmm_out%>%
  filter(Subrate=="dN")%>%
  filter(Group =="Ectotherms")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="dN", title = "Ectotherms", tag="e")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.21,0.4), breaks = seq(-0.2,0.4,0.2))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=9),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()


f3f <- pglmm_out%>%
  filter(Subrate=="dN")%>%
  filter(Group =="Endotherms")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="dN", title = "Endotherms", tag="f")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.21,0.4), breaks = seq(-0.2,0.4,0.2))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=9),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()





f3a+f3b+f3c+f3d+f3e+f3f
ggsave(filename="./Outputs/MainFigures/Fig3af.pdf", height=5, width=8.27)

#Fig.3gh
f3g1 <- pglmm_out%>%
  filter(Subrate=="dS")%>%
  filter(Group =="Fishes")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="dS", title = "Fishes")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.5,0.5), breaks = seq(-0.5,0.5,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f3g2 <- pglmm_out%>%
  filter(Subrate=="dS")%>%
  filter(Group =="Amphibians")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="", title = "Amphibians")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.5,0.5), breaks = seq(-0.5,0.5,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f3g3 <- pglmm_out%>%
  filter(Subrate=="dS")%>%
  filter(Group =="Reptiles")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="", title = "Reptiles")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.5,0.5), breaks = seq(-0.5,0.5,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f3g4 <- pglmm_out%>%
  filter(Subrate=="dS")%>%
  filter(Group =="Mammals")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="", title = "Mammals")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.5,0.5), breaks = seq(-0.5,0.5,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f3g5 <- pglmm_out%>%
  filter(Subrate=="dS")%>%
  filter(Group =="Birds")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="", title = "Birds")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.5,0.5), breaks = seq(-0.5,0.5,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f3h1 <- pglmm_out%>%
  filter(Subrate=="dN")%>%
  filter(Group =="Fishes")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="dN", title = "Fishes")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.5,0.75), breaks = seq(-0.5,0.75,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f3h2 <- pglmm_out%>%
  filter(Subrate=="dN")%>%
  filter(Group =="Amphibians")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="", title = "Amphibians")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.5,0.75), breaks = seq(-0.5,0.75,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f3h3 <- pglmm_out%>%
  filter(Subrate=="dN")%>%
  filter(Group =="Reptiles")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="", title = "Reptiles")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.5,0.75), breaks = seq(-0.5,0.75,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f3h4 <- pglmm_out%>%
  filter(Subrate=="dN")%>%
  filter(Group =="Mammals")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="", title = "Mammals")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.5,0.75), breaks = seq(-0.5,0.75,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f3h5 <- pglmm_out%>%
  filter(Subrate=="dN")%>%
  filter(Group =="Birds")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="", title = "Birds")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.5,0.75), breaks = seq(-0.5,0.75,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()


f3g1+f3g2+f3g3+f3g4+f3g5+f3h1+f3h2+f3h3+f3h4+f3h5+
  plot_layout(ncol = 5, nrow = 2)

ggsave(filename="./Outputs/MainFigures/Fig3gh.pdf", height=4, width=8.27)


#################################################
#3.2 Accounting for phylogenetic signal using the best-fit evolutionary model
load(file="./Outputs/Data/trait_models.rdata")

pred1 <- c("AnnualTemp", "BodyMass", "Longevity","MaturityAge","Fecundity")
pred2 <- c("AnnualTemp", "BodyMass", "Longevity","MaturityAge","Fecundity","dS")

prior <- list(R=list(V=diag(1),nu=0.002), 
              G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*1000)))

groups <- c(unique(subrate$Group), unique(subrate$ThermoMode))

multi_pglmm_list <- vector("list", length = length(groups) * 2)
names(multi_pglmm_list) <- paste(rep(groups, each = 2), c("dN", "dS"), sep = ".")

# Prior specification for the model
for (group in groups) {
  # Choose data based on group
  if(group %in% unique(subrate$ThermoMode)){
    subdata <- subrate%>% 
      filter(ThermoMode == group)%>%
      mutate(dS=scale(log(dS)), dN=scale(log(dN)),AnnualTemp=scale(AnnualTemp),
             BodyMass=scale(log(BodyMass)), Fecundity=scale(log(Fecundity)), 
             MaturityAge=scale(log(MaturityAge)), Longevity=scale(log(Longevity)))
    
  }else{
    subdata <- subrate %>% 
      filter(Group == group)%>%
      mutate(dS=scale(log(dS)), dN=scale(log(dN)),AnnualTemp=scale(AnnualTemp),
             BodyMass=scale(log(BodyMass)), Fecundity=scale(log(Fecundity)), 
             MaturityAge=scale(log(MaturityAge)), Longevity=scale(log(Longevity)))
  }
  
  phy <- drop.tip(phylo, setdiff(phylo$tip.label, subdata$Species))
  
  #The best traits evolution model of dS
  model.ds <- trait_models[[paste(group, "dS", sep = ".")]]
  
  #The best traits evolution model of dN
  model.dn <- trait_models[[paste(group, "dN", sep = ".")]]
  
  
  #rescale tree
  if(model.ds$best_model=="OU"){phy_ds <- phytools::rescale(phy, model=model.ds$best_model, alpha=model.ds$opt, sigsq=model.ds$sigsq)}
  if(model.ds$best_model=="BM"){phy_ds <- phytools::rescale(phy, model=model.ds$best_model, sigsq=model.ds$sigsq)}
  if(model.ds$best_model=="EB"){phy_ds <- phytools::rescale(phy, model=model.ds$best_model, a=model.ds$opt, sigsq=model.ds$sigsq)}
  
  #rescale tree
  if(model.dn$best_model=="OU"){phy_dn <- phytools::rescale(phy, model=model.dn$best_model, alpha=model.dn$opt, sigsq=model.dn$sigsq)}
  if(model.dn$best_model=="BM"){phy_dn <- phytools::rescale(phy, model=model.dn$best_model, sigsq=model.dn$sigsq)}
  if(model.dn$best_model=="EB"){phy_dn <- phytools::rescale(phy, model=model.dn$best_model, a=model.dn$opt, sigsq=model.dn$sigsq)}
  
  
  #modify edge.length in case of small values
  for (i in 1:length(phy_ds$edge.length)) {
    if (phy_ds$edge.length[i] < 1e-16) {
      phy_ds$edge.length[i] = 1e-16
    }
  }
  for (i in 1:length(phy_dn$edge.length)) {
    if (phy_dn$edge.length[i] < 1e-16) {
      phy_dn$edge.length[i] = 1e-16
    }
  }
  
  #Priors
  fit_ds <- MCMCglmm(formula(paste("dS~", paste(pred1, collapse="+"))), 
                     random=~Species, prior=prior,
                     ginverse = list(Species=inverseA(phy_ds)$Ainv), 
                     data=subdata,verbose=TRUE, 
                     nitt=60000, burnin=10000, thin=25, family = c("gaussian"))
  fit_dn <- MCMCglmm(formula(paste("dN~", paste(pred2, collapse="+"))), 
                     random=~Species, prior=prior,
                     ginverse = list(Species=inverseA(phy_dn)$Ainv),
                     data=subdata,verbose=TRUE, 
                     nitt=60000, burnin=10000, thin=25, family = c("gaussian"))
  
  # Store models in list
  multi_pglmm_list[[paste(group, "dS", sep = ".")]] <- fit_ds
  multi_pglmm_list[[paste(group, "dN", sep = ".")]] <- fit_dn
  
  print(group)
}

#save(multi_pglmm_list, file="./Outputs/Data/Multiple_PGLMM_best_traits_model.rdata")
load(file="./Outputs/Data/Multiple_PGLMM_best_traits_model.rdata")




summary_ds <-NULL
summary_dn <-NULL
for(group in groups){
  fit_ds <- multi_pglmm_list[[paste(group, "dS", sep = ".")]]
  fit_dn <- multi_pglmm_list[[paste(group, "dN", sep = ".")]]
  
  sum_ds <- summary(fit_ds)
  sum_dn <- summary(fit_dn)
  
  ds_fixed <- as.data.frame(sum_ds$solutions)
  ds_random <- as.data.frame(sum_ds$Gcovariances)%>%
    mutate(pMCMC=NA)
  ds_r <- r.squaredMCMCglmm(fit_ds)
  
  ds_fixed <- data.frame(Var=row.names(ds_fixed), ds_fixed)
  ds_random <- data.frame(Var=row.names(ds_random), ds_random)
  ds_r <- data.frame(Var=c("Rm2", "Rc2","Rr2"),rbind(as.numeric(ds_r[1:3]), as.numeric(ds_r[5:7]), as.numeric(ds_r[9:11])), NA,NA)
  colnames(ds_r) <- colnames(ds_random)
  
  summary_ds <- rbind(summary_ds, rbind(ds_fixed, ds_random,ds_r)%>%mutate(Group=group))
  
  dn_fixed <- as.data.frame(sum_dn$solutions)
  dn_random <- as.data.frame(sum_dn$Gcovariances)%>%
    mutate(pMCMC=NA)
  dn_r <- r.squaredMCMCglmm(fit_dn)
  
  dn_fixed <- data.frame(Var=row.names(dn_fixed), dn_fixed)
  dn_random <- data.frame(Var=row.names(dn_random), dn_random)
  dn_r <- data.frame(Var=c("Rm2", "Rc2","Rr2"),rbind(as.numeric(dn_r[1:3]), as.numeric(dn_r[5:7]), as.numeric(dn_r[9:11])), NA,NA)
  colnames(dn_r) <- colnames(dn_random)
  
  summary_dn <- rbind(summary_dn, rbind(dn_fixed, dn_random,dn_r)%>%mutate(Group=group))
}

summary_ds
summary_dn

#extract parameter values from PGLMMs
pglmm_out <- NULL
for(i in 1:length(groups)){
  group <- groups[i]
  fit.ds.out <- as.data.frame(summary(multi_pglmm_list[[paste(group, "dS", sep = ".")]])$solutions[-1,])
  fit.dn.out <- as.data.frame(summary(multi_pglmm_list[[paste(group, "dN", sep = ".")]])$solutions[-1,])
  fit.ds.out <- fit.ds.out%>%mutate(env=row.names(fit.ds.out), Subrate="dS")
  fit.dn.out <- fit.dn.out%>%mutate(env=row.names(fit.dn.out), Subrate="dN")
  
  fit.out <- rbind(fit.ds.out, fit.dn.out)%>%
    mutate(Group=group, ThermoMode=ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"),
           sig=ifelse(pMCMC<0.05, "1", "0"))
  pglmm_out <- rbind(pglmm_out, fit.out)
}



#Extended Fig.10
f3 <- pglmm_out%>%
  filter(Subrate=="dS")%>%
  filter(Group =="Ectotherms")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="dS", title = "Ectotherms")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.45,0.21), breaks = seq(-0.4,0.2,0.2))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=9),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()


f4 <- pglmm_out%>%
  filter(Subrate=="dS")%>%
  filter(Group =="Endotherms")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="", title = "Endotherms")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.45,0.21), breaks = seq(-0.4,0.2,0.2))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=9),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()



f5 <- pglmm_out%>%
  filter(Subrate=="dN")%>%
  filter(Group =="Ectotherms")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="dN", title = "Ectotherms")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.21,0.4), breaks = seq(-0.2,0.4,0.2))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=9),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()


f6 <- pglmm_out%>%
  filter(Subrate=="dN")%>%
  filter(Group =="Endotherms")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="", title = "Endotherms")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.21,0.4), breaks = seq(-0.2,0.4,0.2))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=9),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()






#
f7 <- pglmm_out%>%
  filter(Subrate=="dS")%>%
  filter(Group =="Fishes")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="dS", title = "Fishes")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.5,0.5), breaks = seq(-0.5,0.5,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f8 <- pglmm_out%>%
  filter(Subrate=="dS")%>%
  filter(Group =="Amphibians")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="", title = "Amphibians")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.5,0.5), breaks = seq(-0.5,0.5,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f9 <- pglmm_out%>%
  filter(Subrate=="dS")%>%
  filter(Group =="Reptiles")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="", title = "Reptiles")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.5,0.5), breaks = seq(-0.5,0.5,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f10 <- pglmm_out%>%
  filter(Subrate=="dS")%>%
  filter(Group =="Mammals")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="", title = "Mammals")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.5,0.5), breaks = seq(-0.5,0.5,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f11 <- pglmm_out%>%
  filter(Subrate=="dS")%>%
  filter(Group =="Birds")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="", title = "Birds")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.5,0.5), breaks = seq(-0.5,0.5,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f12 <- pglmm_out%>%
  filter(Subrate=="dN")%>%
  filter(Group =="Fishes")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="dN", title = "Fishes")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.5,0.75), breaks = seq(-0.5,0.75,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f13 <- pglmm_out%>%
  filter(Subrate=="dN")%>%
  filter(Group =="Amphibians")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="", title = "Amphibians")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.5,0.75), breaks = seq(-0.5,0.75,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f14 <- pglmm_out%>%
  filter(Subrate=="dN")%>%
  filter(Group =="Reptiles")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="", title = "Reptiles")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.5,0.75), breaks = seq(-0.5,0.75,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f15 <- pglmm_out%>%
  filter(Subrate=="dN")%>%
  filter(Group =="Mammals")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="", title = "Mammals")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.5,0.75), breaks = seq(-0.5,0.75,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f16 <- pglmm_out%>%
  filter(Subrate=="dN")%>%
  filter(Group =="Birds")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="", title = "Birds")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.5,0.75), breaks = seq(-0.5,0.75,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f3+f4+f5+f6+plot_layout(nrow=1)

ggsave(filename="./Outputs/Supplementary/Extended Fig10a.pdf", height=2.2, width=8.27)

f7+f8+f9+f10+f11+f12+f13+f14+f15+f16+
  plot_layout(ncol = 5, nrow = 2)

ggsave(filename="./Outputs/Supplementary/Extended Fig10b.pdf", height=4, width=8.27)



################################################
#Multiple PGLMM using randoms of phylogeny, habitat and lineage
groups <- c(unique(subrate$Group), unique(subrate$ThermoMode))

multi_pglmm_list <- vector("list", length = length(groups) * 2)
names(multi_pglmm_list) <- paste(rep(groups, each = 2), c("dN", "dS"), sep = ".")

prior <- list(R=list(V=diag(1),nu=0.002), 
              G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*1000),
                     G2=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*1000),
                     G3=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*1000)))

for (group in groups) {
  # Choose data based on group
  if(group %in% unique(subrate$ThermoMode)){
    subdata <- subrate%>% 
      filter(ThermoMode == group)%>%
      mutate(dS=scale(log(dS)), dN=scale(log(dN)),AnnualTemp=scale(AnnualTemp),
             BodyMass=scale(log(BodyMass)), Fecundity=scale(log(Fecundity)), 
             MaturityAge=scale(log(MaturityAge)), Longevity=scale(log(Longevity)))
    
  }else{
    subdata <- subrate %>% 
      filter(Group == group)%>%
      mutate(dS=scale(log(dS)), dN=scale(log(dN)),AnnualTemp=scale(AnnualTemp),
             BodyMass=scale(log(BodyMass)), Fecundity=scale(log(Fecundity)), 
             MaturityAge=scale(log(MaturityAge)), Longevity=scale(log(Longevity)))
  }
  
  
  fit_ds <- MCMCglmm(formula(paste("dS~", paste(pred1, collapse="+"))), 
                     random=~Species+Habitat+Clade, prior=prior,
                     ginverse = list(Species=inverseA(phylo)$Ainv), 
                     data=subdata,verbose=TRUE, 
                     nitt=60000, burnin=10000, thin=25, family = c("gaussian"))
  fit_dn <- MCMCglmm(formula(paste("dN~", paste(pred2, collapse="+"))), 
                     random=~Species+Habitat+Clade, prior=prior,
                     ginverse = list(Species=inverseA(phylo)$Ainv),
                     data=subdata,verbose=TRUE, 
                     nitt=60000, burnin=10000, thin=25, family = c("gaussian"))
  
  # Store models in list
  multi_pglmm_list[[paste(group, "dS", sep = ".")]] <- fit_ds
  multi_pglmm_list[[paste(group, "dN", sep = ".")]] <- fit_dn
  
  print(group)
}

#save(multi_pglmm_list, file="./Outputs/Data/Multiple_PGLMM_three_random_effects.rdata")
load(file="./Outputs/Data/Multiple_PGLMM_three_random_effects.rdata")



for(group in groups){
  fit_ds <- multi_pglmm_list[[paste(group, "dS", sep = ".")]]
  fit_dn <- multi_pglmm_list[[paste(group, "dN", sep = ".")]]
  
  sum_ds <- summary(fit_ds)
  sum_dn <- summary(fit_dn)
  
  ds_fixed <- as.data.frame(sum_ds$solutions)
  ds_random <- as.data.frame(sum_ds$Gcovariances)%>%
    mutate(pMCMC=NA)
  ds_r <- r.squaredMCMCglmm(fit_ds)
  
  ds_fixed <- data.frame(Var=row.names(ds_fixed), ds_fixed)
  ds_random <- data.frame(Var=row.names(ds_random), ds_random)
  ds_r <- data.frame(Var=c("Rm2", "Rc2","Rr2"),rbind(as.numeric(ds_r[1:3]), as.numeric(ds_r[5:7]), as.numeric(ds_r[9:11])), NA,NA)
  colnames(ds_r) <- colnames(ds_random)
  
  summary_ds <- rbind(summary_ds, rbind(ds_fixed, ds_random,ds_r)%>%mutate(Group=group))
  
  dn_fixed <- as.data.frame(sum_dn$solutions)
  dn_random <- as.data.frame(sum_dn$Gcovariances)%>%
    mutate(pMCMC=NA)
  dn_r <- r.squaredMCMCglmm(fit_dn)
  
  dn_fixed <- data.frame(Var=row.names(dn_fixed), dn_fixed)
  dn_random <- data.frame(Var=row.names(dn_random), dn_random)
  dn_r <- data.frame(Var=c("Rm2", "Rc2","Rr2"),rbind(as.numeric(dn_r[1:3]), as.numeric(dn_r[5:7]), as.numeric(dn_r[9:11])), NA,NA)
  colnames(dn_r) <- colnames(dn_random)
  
  summary_dn <- rbind(summary_dn, rbind(dn_fixed, dn_random,dn_r)%>%mutate(Group=group))
}

summary_ds 
summary_dn 

#extract parameter values from PGLMMs
pglmm_out <- NULL
for(i in 1:length(groups)){
  group <- groups[i]
  fit.ds.out <- as.data.frame(summary(multi_pglmm_list[[paste(group, "dS", sep = ".")]])$solutions[-1,])
  fit.dn.out <- as.data.frame(summary(multi_pglmm_list[[paste(group, "dN", sep = ".")]])$solutions[-1,])
  fit.ds.out <- fit.ds.out%>%mutate(env=row.names(fit.ds.out), Subrate="dS")
  fit.dn.out <- fit.dn.out%>%mutate(env=row.names(fit.dn.out), Subrate="dN")
  
  fit.out <- rbind(fit.ds.out, fit.dn.out)%>%
    mutate(Group=group, ThermoMode=ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"),
           sig=ifelse(pMCMC<0.05, "1", "0"))
  pglmm_out <- rbind(pglmm_out, fit.out)
}

################
f3 <- pglmm_out%>%
  filter(Subrate=="dS")%>%
  filter(Group =="Ectotherms")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="dS", title = "Ectotherms")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.45,0.21), breaks = seq(-0.4,0.2,0.2))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=9),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()


f4 <- pglmm_out%>%
  filter(Subrate=="dS")%>%
  filter(Group =="Endotherms")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="", title = "Endotherms")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.45,0.21), breaks = seq(-0.4,0.2,0.2))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=9),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()



f5 <- pglmm_out%>%
  filter(Subrate=="dN")%>%
  filter(Group =="Ectotherms")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="dN", title = "Ectotherms")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.21,0.4), breaks = seq(-0.2,0.4,0.2))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=9),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()


f6 <- pglmm_out%>%
  filter(Subrate=="dN")%>%
  filter(Group =="Endotherms")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="", title = "Endotherms")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.21,0.4), breaks = seq(-0.2,0.4,0.2))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=9),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()




f3+f4+f5+f6+
  plot_layout(nrow=1)

ggsave(filename="./Outputs/Supplementary/Extended Fig11a.pdf", height=2.2, width=8.27)


#
f7 <- pglmm_out%>%
  filter(Subrate=="dS")%>%
  filter(Group =="Fishes")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="dS", title = "Fishes")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.5,0.5), breaks = seq(-0.5,0.5,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f8 <- pglmm_out%>%
  filter(Subrate=="dS")%>%
  filter(Group =="Amphibians")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="", title = "Amphibians")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.5,0.5), breaks = seq(-0.5,0.5,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f9 <- pglmm_out%>%
  filter(Subrate=="dS")%>%
  filter(Group =="Reptiles")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="", title = "Reptiles")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.5,0.5), breaks = seq(-0.5,0.5,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f10 <- pglmm_out%>%
  filter(Subrate=="dS")%>%
  filter(Group =="Mammals")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="", title = "Mammals")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.5,0.5), breaks = seq(-0.5,0.5,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f11 <- pglmm_out%>%
  filter(Subrate=="dS")%>%
  filter(Group =="Birds")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="", title = "Birds")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.5,0.5), breaks = seq(-0.5,0.5,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f12 <- pglmm_out%>%
  filter(Subrate=="dN")%>%
  filter(Group =="Fishes")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="dN", title = "Fishes")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.5,0.75), breaks = seq(-0.5,0.75,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f13 <- pglmm_out%>%
  filter(Subrate=="dN")%>%
  filter(Group =="Amphibians")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="", title = "Amphibians")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.5,0.75), breaks = seq(-0.5,0.75,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f14 <- pglmm_out%>%
  filter(Subrate=="dN")%>%
  filter(Group =="Reptiles")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="", title = "Reptiles")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.5,0.75), breaks = seq(-0.5,0.75,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f15 <- pglmm_out%>%
  filter(Subrate=="dN")%>%
  filter(Group =="Mammals")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="", title = "Mammals")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.5,0.75), breaks = seq(-0.5,0.75,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f16 <- pglmm_out%>%
  filter(Subrate=="dN")%>%
  filter(Group =="Birds")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="", title = "Birds")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.5,0.75), breaks = seq(-0.5,0.75,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()


f7+f8+f9+f10+f11+f12+f13+f14+f15+f16+
  plot_layout(ncol = 5, nrow = 2)

ggsave(filename="./Outputs/Supplementary/Extended Fig11b.pdf", height=4, width=8.27)


##############
#Multiple PGLMM using data without missing traits
groups <- c(unique(subrate$Group), unique(subrate$ThermoMode))

multi_pglmm_list <- vector("list", length = length(groups) * 2)
names(multi_pglmm_list) <- paste(rep(groups, each = 2), c("dN", "dS"), sep = ".")


prior <- list(R=list(V=diag(1),nu=0.002), 
              G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*1000)))

for (group in groups) {
  # Choose data based on group
  if(group %in% unique(subrate$ThermoMode)){
    subdata <- subrate%>% 
      filter(ThermoMode == group)%>%
      filter(!is.na(Ref.Fecundity), !is.na(Ref.Longevity),!is.na(Ref.MaturityAge))%>%
      mutate(dS=scale(log(dS)), dN=scale(log(dN)),AnnualTemp=scale(AnnualTemp),
             BodyMass=scale(log(BodyMass)), Fecundity=scale(log(Fecundity)), 
             MaturityAge=scale(log(MaturityAge)), Longevity=scale(log(Longevity)))
    
  }else{
    subdata <- subrate %>% 
      filter(Group == group)%>%
      filter(!is.na(Ref.Fecundity), !is.na(Ref.Longevity),!is.na(Ref.MaturityAge))%>%
      mutate(dS=scale(log(dS)), dN=scale(log(dN)),AnnualTemp=scale(AnnualTemp),
             BodyMass=scale(log(BodyMass)), Fecundity=scale(log(Fecundity)), 
             MaturityAge=scale(log(MaturityAge)), Longevity=scale(log(Longevity)))
  }
  
  
  fit_ds <- MCMCglmm(formula(paste("dS~", paste(pred1, collapse="+"))), 
                     random=~Species, prior=prior,
                     ginverse = list(Species=inverseA(phylo)$Ainv), 
                     data=subdata,verbose=TRUE, 
                     nitt=60000, burnin=10000, thin=25, family = c("gaussian"))
  fit_dn <- MCMCglmm(formula(paste("dN~", paste(pred2, collapse="+"))), 
                     random=~Species, prior=prior,
                     ginverse = list(Species=inverseA(phylo)$Ainv),
                     data=subdata,verbose=TRUE, 
                     nitt=60000, burnin=10000, thin=25, family = c("gaussian"))
  
  # Store models in list
  multi_pglmm_list[[paste(group, "dS", sep = ".")]] <- fit_ds
  multi_pglmm_list[[paste(group, "dN", sep = ".")]] <- fit_dn
  
  print(group)
}

#save(multi_pglmm_list, file="./Outputs/Data/Multiple_PGLMM_without_missing_traits.rdata")
load(file="./Outputs/Data/Multiple_PGLMM_without_missing_traits.rdata")

summary_ds <-NULL
summary_dn <-NULL
for(group in groups){
  fit_ds <- multi_pglmm_list[[paste(group, "dS", sep = ".")]]
  fit_dn <- multi_pglmm_list[[paste(group, "dN", sep = ".")]]
  
  sum_ds <- summary(fit_ds)
  sum_dn <- summary(fit_dn)
  
  ds_fixed <- as.data.frame(sum_ds$solutions)
  ds_random <- as.data.frame(sum_ds$Gcovariances)%>%
    mutate(pMCMC=NA)
  ds_r <- r.squaredMCMCglmm(fit_ds)
  
  ds_fixed <- data.frame(Var=row.names(ds_fixed), ds_fixed)
  ds_random <- data.frame(Var=row.names(ds_random), ds_random)
  ds_r <- data.frame(Var=c("Rm2", "Rc2","Rr2"),rbind(as.numeric(ds_r[1:3]), as.numeric(ds_r[5:7]), as.numeric(ds_r[9:11])), NA,NA)
  colnames(ds_r) <- colnames(ds_random)
  
  summary_ds <- rbind(summary_ds, rbind(ds_fixed, ds_random,ds_r)%>%mutate(Group=group))
  
  dn_fixed <- as.data.frame(sum_dn$solutions)
  dn_random <- as.data.frame(sum_dn$Gcovariances)%>%
    mutate(pMCMC=NA)
  dn_r <- r.squaredMCMCglmm(fit_dn)
  
  dn_fixed <- data.frame(Var=row.names(dn_fixed), dn_fixed)
  dn_random <- data.frame(Var=row.names(dn_random), dn_random)
  dn_r <- data.frame(Var=c("Rm2", "Rc2","Rr2"),rbind(as.numeric(dn_r[1:3]), as.numeric(dn_r[5:7]), as.numeric(dn_r[9:11])), NA,NA)
  colnames(dn_r) <- colnames(dn_random)
  
  summary_dn <- rbind(summary_dn, rbind(dn_fixed, dn_random,dn_r)%>%mutate(Group=group))
}

summary_ds 
summary_dn 

#extract parameter values from PGLMMs
pglmm_out <- NULL
for(i in 1:length(groups)){
  group <- groups[i]
  fit.ds.out <- as.data.frame(summary(multi_pglmm_list[[paste(group, "dS", sep = ".")]])$solutions[-1,])
  fit.dn.out <- as.data.frame(summary(multi_pglmm_list[[paste(group, "dN", sep = ".")]])$solutions[-1,])
  fit.ds.out <- fit.ds.out%>%mutate(env=row.names(fit.ds.out), Subrate="dS")
  fit.dn.out <- fit.dn.out%>%mutate(env=row.names(fit.dn.out), Subrate="dN")
  
  fit.out <- rbind(fit.ds.out, fit.dn.out)%>%
    mutate(Group=group, ThermoMode=ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"),
           sig=ifelse(pMCMC<0.05, "1", "0"))
  pglmm_out <- rbind(pglmm_out, fit.out)
}

################
f3 <- pglmm_out%>%
  filter(Subrate=="dS")%>%
  filter(Group =="Ectotherms")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="dS", title = "Ectotherms")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.6,0.4), breaks = seq(-0.6,0.3,0.3))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=9),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

#Extended Fig.12
f4 <- pglmm_out%>%
  filter(Subrate=="dS")%>%
  filter(Group =="Endotherms")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="", title = "Endotherms")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.6,0.3), breaks = seq(-0.6,0.3,0.3))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=9),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()



f5 <- pglmm_out%>%
  filter(Subrate=="dN")%>%
  filter(Group =="Ectotherms")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="dN", title = "Ectotherms")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.6,0.6), breaks = seq(-0.6,0.6,0.3))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=9),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()


f6 <- pglmm_out%>%
  filter(Subrate=="dN")%>%
  filter(Group =="Endotherms")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="", title = "Endotherms")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.6,0.6), breaks = seq(-0.6,0.6,0.3))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=9),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()




f3+f4+f5+f6+
  plot_layout(nrow=1)

ggsave(filename="./Outputs/Supplementary/Extended Fig12a.pdf", height=2.2, width=8.27)


#
f7 <- pglmm_out%>%
  filter(Subrate=="dS")%>%
  filter(Group =="Fishes")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="dS", title = "Fishes")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.8,0.8), breaks = seq(-0.8,0.8,0.4))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f8 <- pglmm_out%>%
  filter(Subrate=="dS")%>%
  filter(Group =="Amphibians")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="", title = "Amphibians")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.8,0.8), breaks = seq(-0.8,0.8,0.4))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f9 <- pglmm_out%>%
  filter(Subrate=="dS")%>%
  filter(Group =="Reptiles")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="", title = "Reptiles")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.8,0.8), breaks = seq(-0.8,0.8,0.4))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f10 <- pglmm_out%>%
  filter(Subrate=="dS")%>%
  filter(Group =="Mammals")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="", title = "Mammals")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.8,0.8), breaks = seq(-0.8,0.8,0.4))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f11 <- pglmm_out%>%
  filter(Subrate=="dS")%>%
  filter(Group =="Birds")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="", title = "Birds")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.8,0.8), breaks = seq(-0.8,0.8,0.4))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f12 <- pglmm_out%>%
  filter(Subrate=="dN")%>%
  filter(Group =="Fishes")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="dN", title = "Fishes")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.8,0.8), breaks = seq(-0.8,0.8,0.4))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f13 <- pglmm_out%>%
  filter(Subrate=="dN")%>%
  filter(Group =="Amphibians")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="", title = "Amphibians")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.8,0.8), breaks = seq(-0.8,0.8,0.4))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f14 <- pglmm_out%>%
  filter(Subrate=="dN")%>%
  filter(Group =="Reptiles")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="", title = "Reptiles")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.8,0.8), breaks = seq(-0.8,0.8,0.4))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f15 <- pglmm_out%>%
  filter(Subrate=="dN")%>%
  filter(Group =="Mammals")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="", title = "Mammals")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.8,0.8), breaks = seq(-0.8,0.8,0.4))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f16 <- pglmm_out%>%
  filter(Subrate=="dN")%>%
  filter(Group =="Birds")%>%
  ggplot(aes(x=env, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="", title = "Birds")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.8,0.8), breaks = seq(-0.8,0.8,0.4))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()


f7+f8+f9+f10+f11+f12+f13+f14+f15+f16+
  plot_layout(ncol = 5, nrow = 2)

ggsave(filename="./Outputs/Supplementary/Extended Fig12b.pdf", height=4, width=8.27)





