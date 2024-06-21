################################
# An R script for analyzing latitudinal gradient in molecular rates at species and assemblage levels
# Author: Tianlong Cai
# Email: caitianlong@westlake.edu.cn

#####################
rm(list=ls())
gc()
#define work direction
workdir <- "/Users/Tianlong/VertMolRate"
setwd(workdir)

### Load packages
library(ape)
library(tidyverse)
library(phytools)
library(MCMCglmm)
library(geiger)
library(RColorBrewer)
library(agricolae)
library(nlme)
# Load source functions
source(paste0(workdir, "/Scripts/source_functions.r"))

###########################################################
# Load input data for molecular evolution rates
subrate <- read.csv(paste0(workdir, "/DataFiles/MolEvolRate/VertMolRate.csv"))

# Load input phylogeny
phylo <- read.tree(paste0(workdir, "/DataFiles/trees/phy_all_sampled_vertebrates.tre"))

# Check whether the phylogeny is ultrametric
is.ultrametric(phylo)

# Force the tree to be ultrametric
phylo <- force.ultrametric(phylo, method = "nnls")  # Use nnls or extend

################################################################################
#1. Analyzing latitudinal gradients in molecular rates across species using
# a Phylogenetic Generalized Linear Mixed Model (PGLMM) that accounts for 
# phylogenetic relatedness as a random effect.

#groups
groups <- c(unique(subrate$Group), unique(subrate$ThermoMode), "Non-long migrants", "Non migrants")

pglmm_list <- vector("list", length = length(groups) * 2)
names(pglmm_list) <- paste(rep(groups, each = 2), c("dN", "dS"), sep = ".")

for(group in groups[-9]){
  pglmm_list[[paste(group, "dS", sep = ".")]] <- x[[paste(group, "dS", sep = ".")]]
  pglmm_list[[paste(group, "dN", sep = ".")]] <- x[[paste(group, "dN", sep = ".")]]
}

# Prior specification for the model
prior <- list(R=list(V=diag(1),nu=0.002), 
              G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*1000)))

# Loop over each group
for (group in groups) {
  # Choose data based on group
  if(group %in% unique(subrate$ThermoMode)){
    subdata <- subrate %>% filter(ThermoMode == group)
  }
  if(group == "Non-long migrants"){
    subdata <- subrate %>% filter(Group == "Birds", Migration!="Long Migratory")
  }
  if(group == "Non migrants"){
    subdata <- subrate %>% filter(Group == "Birds", Migration=="Resident")
  }
  if(group %in% unique(subrate$Group)){
    subdata <- subrate %>% filter(Group == group)
  }
  
  # Scale and fit dS model
  fit_ds <- MCMCglmm(scale(log(dS)) ~ scale(abs(Lat)), prior=prior,
                     random = ~Species,
                     data = subdata, verbose = TRUE, 
                     ginverse = list(Species = inverseA(phylo)$Ainv), 
                     nitt = 60000, burnin = 10000, thin = 25, family = c("gaussian"))
  
  # Scale and fit dN model
  fit_dn <- MCMCglmm(scale(log(dN)) ~ scale(abs(Lat)), prior=prior,
                     random = ~Species,
                     data = subdata, verbose = TRUE, 
                     ginverse = list(Species = inverseA(phylo)$Ainv),
                     nitt = 60000, burnin = 10000, thin = 25, family = c("gaussian"))
  
  # Store models in list
  pglmm_list[[paste(group, "dS", sep = ".")]] <- fit_ds
  pglmm_list[[paste(group, "dN", sep = ".")]] <- fit_dn
  
  print(group)
}

#save(pglmm_list, file = "./Outputs/Data/PGLMM_latitudinal_pattern.rdata")
load(file = "./Outputs/Data/PGLMM_latitudinal_pattern.rdata")

summary_ds <-NULL
summary_dn <-NULL
for(group in groups){
  fit_ds <- pglmm_list[[paste(group, "dS", sep = ".")]]
  fit_dn <- pglmm_list[[paste(group, "dN", sep = ".")]]
  
  sum_ds <- summary(fit_ds)
  sum_dn <- summary(fit_dn)
  
  ds_fixed <- as.data.frame(sum_ds$solutions)
  ds_random <- as.data.frame(sum_ds$Gcovariances)%>%
    mutate(pMCMC=NA)
  
  ds_fixed <- data.frame(Var=row.names(ds_fixed), ds_fixed)
  ds_random <- data.frame(Var=row.names(ds_random), ds_random)
  summary_ds <- rbind(summary_ds, rbind(ds_fixed, ds_random)%>%mutate(Group=group))
  
  dn_fixed <- as.data.frame(sum_dn$solutions)
  dn_random <- as.data.frame(sum_dn$Gcovariances)%>%
    mutate(pMCMC=NA)
  
  dn_fixed <- data.frame(Var=row.names(dn_fixed), dn_fixed)
  dn_random <- data.frame(Var=row.names(dn_random), dn_random)
  summary_dn <- rbind(summary_dn, rbind(dn_fixed, dn_random)%>%mutate(Group=group))
}

write.csv(summary_ds,"summary_ds.csv", row.names = F)
write.csv(summary_dn,"summary_dn.csv", row.names = F)


# Extract parameter values from PGLMMs
pglmm_out <- NULL
for(i in 1:length(groups)){
  group <- groups[i]
  fit.out <- rbind(summary(pglmm_list[[paste(group, "dS", sep = ".")]])$solutions[2,], 
                   summary(pglmm_list[[paste(group, "dN", sep = ".")]])$solutions[2,])
  fit.out <- as.data.frame(fit.out) %>%
    mutate(Subrate = c("dS", "dN"),
           Group = group,
           ThermoMode = ifelse(Group %in% c("Birds", "Mammals", "Endotherms","Non-long migrants", "Non migrants"), "Endotherms", "Ectotherms"))
  pglmm_out <- rbind(pglmm_out, fit.out)
}




#Plot effect size
f2a1 <-pglmm_out%>%
  mutate(Group=ifelse(Group=="Non migrants", "Residents", Group))%>%
  mutate(sig=ifelse(pMCMC<0.05, "1", "0"),
         Group=factor(Group, levels=c("Endotherms","Ectotherms","Residents","Non-long migrants", "Birds","Mammals","Reptiles", "Amphibians","Fishes")))%>%
  filter(Subrate=="dS")%>%
  ggplot(aes(x=Group, y=post.mean, ymin=`l-95% CI`, ymax=`u-95% CI`, color=ThermoMode, shape=sig))+
  geom_pointrange(size=0.4)+
  scale_shape_manual(labels=c("0", "1"), values = c(1,19))+
  scale_color_manual(values = c("#2a6aaf", "#d3292f"))+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  labs(y="Effects of latitude", x="", title = "dS", tag="a")+
  theme_classic()+
  coord_flip()+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-0.3, 0.18),breaks = c(-0.3,-0.15,0,0.15), labels=c(-0.3,-0.15,0,0.15))+
  theme(axis.text = element_text(size=8),
        axis.title = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = 'none',
        plot.tag = element_text(size=10, face = "bold"),
        plot.tag.position = c(0.05, 1),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

f2a2 <-pglmm_out%>%
  mutate(Group=ifelse(Group=="Non migrants", "Residents", Group))%>%
  mutate(sig=ifelse(pMCMC<0.05, "1", "0"),
         Group=factor(Group, levels=c("Endotherms","Ectotherms","Residents","Non-long migrants", "Birds","Mammals","Reptiles", "Amphibians","Fishes")))%>%
  filter(Subrate=="dN")%>%
  ggplot(aes(x=Group, y=post.mean, ymin=`l-95% CI`, ymax=`u-95% CI`, color=ThermoMode, shape=sig))+
  geom_pointrange(size=0.4)+
  scale_shape_manual(labels=c("0", "1"), values = c(1,19))+
  scale_color_manual(values = c("#2a6aaf", "#d3292f"))+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  labs(y="Effects of latitude", x="", title = "dN")+
  theme_classic()+
  coord_flip()+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-0.3, 0.15),breaks = c(-0.3,-0.15,0,0.15), labels=c(-0.3,-0.15,0,0.15))+
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_blank(),
        axis.title = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = 'none',
        plot.tag = element_text(size=10, face = "bold"),
        plot.tag.position = c(0.05, 1),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))
cowplot::plot_grid(f2a1,f2a2, nrow = 1, align="hv")


#########################
#compare dS of migrants and residents
migrant <- subrate%>%
  filter(Group=="Birds")%>%
  mutate(Migration2=ifelse(Migration=="Resident", "Non-migrants", "Migrants"))

migrant%>%
  ggplot(aes(x=Migration, y=log10(dS)))+
  geom_boxplot()+
  ggpubr::stat_compare_means(method = "kruskal.test")+
  scale_y_continuous(limits = c(-8.6, -7.2),breaks = seq(-8.6, -7.2, 0.2))+
  theme_classic()+
  labs(x="", y=expression(log[10](dS)))+
  guides(y=guide_axis(cap='upper'))+
  theme(plot.title = element_text(hjust = 0.5, size=9),
        axis.title = element_text(size=9),
        axis.text = element_text(size=8),
        strip.text = element_text(size=9),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))


#dS corrected for body mass
pgls <- gls(log(dS)~log(BodyMass), migrant, correlation = corMartins(value=0.017, phy=drop.tip(phylo,setdiff(phylo$tip.label,migrant$Species)), form=~Species))
summary(pgls)

#data for plot
ds.corrected.mass <- data.frame(dS.Residual=residuals(pgls), Migration=migrant$Migration, Migration2=migrant$Migration2)
kruskal(ds.corrected.mass$dS.Residual, ds.corrected.mass$Migration, p.adj = "bonferroni")$groups

my_comparisons1 <- list( c("Long Migratory", "Short Migratory"), c("Short Migratory", "Resident"), c("Long Migratory", "Resident"))

f2b <- ds.corrected.mass%>%
  mutate(Migration=factor(Migration, levels = c("Long Migratory", "Short Migratory", "Resident")))%>%
  ggplot(aes(x=Migration, y=dS.Residual, colour=Migration))+
  geom_boxplot()+
  theme_classic()+
  labs(x="", y=expression(dS[Residuals]), tag="b")+
  guides(y=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-1.3,2.0), breaks = seq(-1.0,2.0,0.5))+
  scale_x_discrete(labels = c("Long-dist.\n migrants", "Short-dist.\n migrants", "Resident\n birds"))+
  ggpubr::stat_compare_means(comparisons = my_comparisons1, label = "p.signif", method="wilcox.test", size=2.5, label.x = c(1.5,2.5,2), label.y = c(1.5,1.5,1.7))+ 
  theme(plot.title = element_text(hjust = 0.5, size=8),
        axis.title = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        strip.text = element_text(size=8),
        legend.position = "none",
        plot.tag = element_text(size=10, face = "bold"),
        plot.tag.position = c(0.05, 1),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))


my_comparisons2 <- list( c("Migrants", "Non-migrants"))

f2 <- ds.corrected.mass%>%
  ggplot(aes(x=Migration2, y=dS.Residual, colour=Migration2))+
  geom_boxplot()+
  theme_classic()+
  labs(x="", y=expression(dS[Residuals]))+
  guides(y=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-1.3,2.0), breaks = seq(-1.5,2.0,0.5))+
  ggpubr::stat_compare_means(comparisons = my_comparisons2, method="wilcox.test", size=2.5)+ 
  theme(plot.title = element_text(hjust = 0.5, size=9),
        axis.title = element_text(size=9),
        axis.text = element_text(size=9),
        strip.text = element_text(size=9),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

#############################################################
k=8.62*10^-5
#Datasets of birds from Fristoe et al. (PNAS, 2015) (https://doi.org/10.1073/pnas.1521662112)
BMR.phy <- read.tree(paste0(workdir, "/DataFiles/BMR/BMR_phy.tre"))

BMR <- read.csv(paste0(workdir, "/DataFiles/BMR/BMR_birds.csv"))%>%
  filter(Species.Map %in% BMR.phy$tip.label)%>%
  mutate(B=BMR/Mass, ln.mass=log(Mass), ln.B=log(B), ln.BMR= log(BMR), Inv.kTb=1/((Tb+273.15)*k))

####################################################################
#Fitting logarithm of metabolic rate by body mass and body temperature 
#using Gillooly's model in endotherms
pgls.mass.b <-gls(ln.B~ln.mass+Migration+Inv.kTb, BMR, correlation=corPagel(1, BMR.phy, form = ~Species.Map))
summary(pgls.mass.b)



#Estimating residual of basal metabolic rate
endo.res <- data.frame(Ta=BMR$MeanTemp, BMR.res=residuals(pgls.mass.b), Migration=BMR$Migration, Class=BMR$Class)
lm1 <- lm(BMR.res~Ta, endo.res)
summary(lm1)

f2c <- endo.res%>%
  ggplot(aes(x=Ta, y=BMR.res,colour=Class))+
  geom_point(size=1)+
  scale_shape_manual(values=c(1,4))+
  geom_smooth(method="lm", size=0.8)+
  scale_color_manual(values=c("#d3292f"))+
  labs(x=expression("Ambient temperature (" * degree*"C)"), 
       y=expression("Mass specific BMR"[residuals]),tag="c")+
  #scale_y_continuous(limits = c(-1.2, 1.2),breaks = seq(-1, 1, 0.5))+
  theme_classic()+
  guides(y=guide_axis(cap='upper'))+
  theme(plot.title = element_text(hjust = 0.5, size=8),
        axis.text = element_text(size=8),
        axis.title = element_text(size=8),
        plot.tag = element_text(size=10, face = "bold"),
        plot.tag.position = c(0.05, 1),
        legend.position="none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

cowplot::plot_grid(f2a1,f2a2, f2b,f2c, nrow=1, align="hv")
ggsave(filename="./Outputs/MainFigures/Fig2abc.pdf", height=2.5, width=10.27)



################################################################################
#2. Analyzing latitudinal gradients in molecular rates across species using
# a Phylogenetic Generalized Linear Mixed Model (PGLMM) that accounts for 
# phylogenetic relatedness under the best traits evolution models as a random effect.

#function to fit the best model of traits evolution in the time tree
fitTraitEvolModel <- function(trait, data, phy){
  #match traits
  #trait=match.arg(trait, c("dS", "dN"))
  
  #prepare data
  dat <- data[,trait]
  names(dat)<-data$Species
  
  # define set of models to compare
  models=c("BM", "OU", "EB")
  
  ## ESTIMATING measurement error ##
  aic.se=numeric(length(models))
  lnl.se=numeric(length(models))
  opt.se=numeric(length(models))
  sigsq.se=numeric(length(models))
  for(m in 1:length(models)){
    tmp=fitContinuous(phy, dat=dat, SE=NA, model=models[m], ncores=2)
    aic.se[m]=tmp$opt$aicc
    lnl.se[m]=tmp$opt$lnL
    opt.se[m]=as.numeric(tmp$opt[1])
    sigsq.se[m]=as.numeric(tmp$opt[2])
  }
  
  ## ASSUMING no measurement error ##
  aic=numeric(length(models))
  lnl=numeric(length(models))
  opt=numeric(length(models))
  sigsq=numeric(length(models))
  
  for(m in 1:length(models)){
    tmp=fitContinuous(phy,dat,SE=0,model=models[m], ncores=2)
    aic[m]=tmp$opt$aicc
    lnl[m]=tmp$opt$lnL
    opt[m]=as.numeric(tmp$opt[1])
    sigsq[m]=as.numeric(tmp$opt[2])
  }
  
  ## COMPARE AIC ##
  names(aic.se)<-names(lnl.se)<-names(aic)<-names(lnl)<-names(opt.se)<-names(opt)<-models
  delta_aic<-function(x) x-x[which(x==min(x))]
  # no measurement error
  daic=delta_aic(aic)
  # measurement error
  daic.se=delta_aic(aic.se)
  
  res_aicc= rbind(aic, aic.se, daic, daic.se, opt, opt.se, sigsq, sigsq.se)
  rownames(res_aicc)=c("AICc","AICc_SE","dAICc", "dAICc_SE", "opt", "opt.se", "sigsq", "sigsq.se")
  
  
  row_col <- which(res_aicc[1:2,] == min(res_aicc[1:2,]), arr.ind = TRUE)
  best_model <- colnames(res_aicc)[row_col[2]]
  opt <- res_aicc[row_col[1]+4,best_model]
  sigsq <- res_aicc[row_col[1]+6,best_model]
  
  return(data.frame(best_model, opt, sigsq))
}



##########
#Fit the best models
groups <- c(unique(subrate$Group), unique(subrate$ThermoMode), "Non-long migrants", "Non-migrants")

# Creating a list to store PGLMM objects
trait_models <- vector("list", length = length(groups) * 2)
names(trait_models) <- paste(rep(groups, each = 2), c("dN", "dS"), sep = ".")


for (group in groups) {
  # Choose data based on group
  if(group %in% unique(subrate$ThermoMode)){
    subdata <- subrate %>% filter(ThermoMode == group)%>%mutate(dS=as.numeric(scale(log(dS))), dN=as.numeric(scale(log(dN))))
  }
  if(group == "Non-long migrants"){
    subdata <- subrate %>% filter(Group == "Birds", Migration!="Long Migratory")%>%mutate(dS=as.numeric(scale(log(dS))), dN=as.numeric(scale(log(dN))))
  }
  if(group == "Non-migrants"){
    subdata <- subrate %>% filter(Group == "Birds", Migration=="Resident")%>%mutate(dS=as.numeric(scale(log(dS))), dN=as.numeric(scale(log(dN))))
  }
  if(group %in% unique(subrate$Group)){
    subdata <- subrate %>% filter(Group == group)%>%mutate(dS=as.numeric(scale(log(dS))), dN=as.numeric(scale(log(dN))))
  }
  
  #phylogeny
  phy <- drop.tip(phylo, setdiff(phylo$tip.label, subdata$Species))
  
  #Model selection
  trait_models[[paste(group, "dS", sep = ".")]] <- fitTraitEvolModel(trait="dS", data=subdata, phy = phy)
  trait_models[[paste(group, "dN", sep = ".")]] <- fitTraitEvolModel(trait="dN", data=subdata, phy = phy)

}

#save(trait_models, file="./Outputs/Data/trait_models.rdata")
load(file="./Outputs/Data/trait_models.rdata")


trait_models%>%
  bind_rows()%>%
  mutate(Group=rep(groups, each=2), SubRate=rep(c("dN","dS"), length(groups)))

##############################
#Fit the PGLMM using the best model of traits evolution model
groups <- c(unique(subrate$Group), unique(subrate$ThermoMode), "Non-long migrants","Non-migrants")

pglmm_list <- vector("list", length = length(groups) * 2)
names(pglmm_list) <- paste(rep(groups, each = 2), c("dN", "dS"), sep = ".")

# Prior specification for the model
prior <- list(R=list(V=diag(1),nu=0.002), 
              G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*1000)))


# Loop over each group
for (group in groups) {
  # Choose data based on group
  if(group %in% unique(subrate$ThermoMode)){
    subdata <- subrate %>% filter(ThermoMode == group)
  }
  if(group == "Non-migrants"){
    subdata <- subrate %>% filter(Group == "Birds", Migration=="Resident")
  }
  if(group == "Non-long migrants"){
    subdata <- subrate %>% filter(Group == "Birds", Migration!="Long Migratory")
  }
  if(group %in% unique(subrate$Group)){
    subdata <- subrate %>% filter(Group == group)
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
  
  
  # Scale and fit dS model
  fit_ds <- MCMCglmm(scale(log(dS)) ~ scale(abs(Lat)), prior=prior,
                     random = ~Species,
                     data = subdata, verbose = TRUE, 
                     ginverse = list(Species = inverseA(phy_ds)$Ainv), 
                     nitt = 60000, burnin = 10000, thin = 25, family = c("gaussian"))
  
  # Scale and fit dN model
  fit_dn <- MCMCglmm(scale(log(dN)) ~ scale(abs(Lat)), prior=prior,
                     random = ~Species,
                     data = subdata, verbose = TRUE, 
                     ginverse = list(Species = inverseA(phy_dn)$Ainv),
                     nitt = 60000, burnin = 10000, thin = 25, family = c("gaussian"))
  
  # Store models in list
  pglmm_list[[paste(group, "dS", sep = ".")]] <- fit_ds
  pglmm_list[[paste(group, "dN", sep = ".")]] <- fit_dn
  
  print(group)
}

#save(pglmm_list, file = "./Outputs/Data/PGLMM_latitudinal_pattern_best_traits_model.rdata")
load(file = "./Outputs/Data/PGLMM_latitudinal_pattern_best_traits_model.rdata")

summary_ds <-NULL
summary_dn <-NULL
for(group in groups){
  fit_ds <- pglmm_list[[paste(group, "dS", sep = ".")]]
  fit_dn <- pglmm_list[[paste(group, "dN", sep = ".")]]
  
  sum_ds <- summary(fit_ds)
  sum_dn <- summary(fit_dn)
  
  ds_fixed <- as.data.frame(sum_ds$solutions)
  ds_random <- as.data.frame(sum_ds$Gcovariances)%>%
    mutate(pMCMC=NA)
  
  ds_fixed <- data.frame(Var=row.names(ds_fixed), ds_fixed)
  ds_random <- data.frame(Var=row.names(ds_random), ds_random)
  summary_ds <- rbind(summary_ds, rbind(ds_fixed, ds_random)%>%mutate(Group=group))
  
  dn_fixed <- as.data.frame(sum_dn$solutions)
  dn_random <- as.data.frame(sum_dn$Gcovariances)%>%
    mutate(pMCMC=NA)
  
  dn_fixed <- data.frame(Var=row.names(dn_fixed), dn_fixed)
  dn_random <- data.frame(Var=row.names(dn_random), dn_random)
  summary_dn <- rbind(summary_dn, rbind(dn_fixed, dn_random)%>%mutate(Group=group))
}

write.csv(summary_ds,"summary_ds.csv", row.names = F)
write.csv(summary_dn,"summary_dn.csv", row.names = F)



# Extract parameter values from PGLMMs
pglmm_out <- NULL
for(i in 1:length(groups)){
  group <- groups[i]
  fit.out <- rbind(summary(pglmm_list[[paste(group, "dS", sep = ".")]])$solutions[2,], 
                   summary(pglmm_list[[paste(group, "dN", sep = ".")]])$solutions[2,])
  fit.out <- as.data.frame(fit.out) %>%
    mutate(Subrate = c("dS", "dN"),
           Group = group,
           ThermoMode = ifelse(Group %in% c("Birds", "Mammals", "Endotherms","Non-long migrants", "Non-migrants"), "Endotherms", "Ectotherms"))
  pglmm_out <- rbind(pglmm_out, fit.out)
}


#Plot effect size
f1 <-pglmm_out%>%
  mutate(Group=ifelse(Group=="Non-migrants", "Residents", Group))%>%
  mutate(sig=ifelse(pMCMC<0.05, "1", "0"),
         Group=factor(Group, levels=c("Endotherms","Ectotherms","Residents", "Non-long migrants", "Birds","Mammals","Reptiles", "Amphibians","Fishes")))%>%
  filter(Subrate=="dS")%>%
  ggplot(aes(x=Group, y=post.mean, ymin=`l-95% CI`, ymax=`u-95% CI`, color=ThermoMode, shape=sig))+
  geom_pointrange(size=0.4)+
  scale_shape_manual(labels=c("0", "1"), values = c(1,19))+
  scale_color_manual(values = c("#2a6aaf", "#d3292f"))+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  labs(y="Effects of latitude", x="", title = "dS")+
  theme_classic()+
  coord_flip()+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-0.3, 0.18),breaks = c(-0.3,-0.15,0,0.15), labels=c(-0.3,-0.15,0,0.15))+
  theme(axis.text = element_text(size=7),
        axis.title = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        legend.position = 'none',
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f2 <-pglmm_out%>%
  mutate(Group=ifelse(Group=="Non-migrants", "Residents", Group))%>%
  mutate(sig=ifelse(pMCMC<0.05, "1", "0"),
         Group=factor(Group, levels=c("Endotherms","Ectotherms","Residents", "Non-long migrants", "Birds","Mammals","Reptiles", "Amphibians","Fishes")))%>%
  filter(Subrate=="dN")%>%
  ggplot(aes(x=Group, y=post.mean, ymin=`l-95% CI`, ymax=`u-95% CI`, color=ThermoMode, shape=sig))+
  geom_pointrange(size=0.4)+
  scale_shape_manual(labels=c("0", "1"), values = c(1,19))+
  scale_color_manual(values = c("#2a6aaf", "#d3292f"))+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  labs(y="Effects of latitude", x="", title = "dN")+
  theme_classic()+
  coord_flip()+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-0.3, 0.15),breaks = c(-0.3,-0.15,0,0.15), labels=c(-0.3,-0.15,0,0.15))+
  theme(axis.text = element_text(size=7),
        axis.title = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        legend.position = 'none',
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))
cowplot::plot_grid(f1,f2, nrow = 2, align="hv")


########################################################
#3. Analyzing latitudinal gradients in molecular rates across species using
# PGLMMs that accounts for phylogenetic relatedness, habitat and lineages 
# used for estimating molecular rates as random effects.
groups <- c(unique(subrate$Group), unique(subrate$ThermoMode), "Non-long migrants", "Non-migrants")

pglmm_list <- vector("list", length = length(groups) * 2)
names(pglmm_list) <- paste(rep(groups, each = 2), c("dN", "dS"), sep = ".")

# Prior specification for the model
prior <- list(R=list(V=diag(1),nu=0.002), 
              G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*1000),
                     G2=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*1000),
                     G3=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*1000)))

# Loop over each group
for (group in groups) {
  # Choose data based on group
  if(group %in% unique(subrate$ThermoMode)){
    subdata <- subrate %>% filter(ThermoMode == group)
  }
  if(group == "Non-migrants"){
    subdata <- subrate %>% filter(Group == "Birds", Migration=="Resident")
  }
  if(group == "Non-long migrants"){
    subdata <- subrate %>% filter(Group == "Birds", Migration!="Long Migratory")
  }
  if(group %in% unique(subrate$Group)){
    subdata <- subrate %>% filter(Group == group)
  }
  
  # Scale and fit dS model
  fit_ds <- MCMCglmm(scale(log(dS)) ~ scale(abs(Lat)),prior=prior,
                     random = ~Species+Habitat+Clade,
                     data = subdata, verbose = TRUE, 
                     ginverse = list(Species = inverseA(phylo)$Ainv), 
                     nitt = 60000, burnin = 10000, thin = 25, family = c("gaussian"))
  
  # Scale and fit dN model
  fit_dn <- MCMCglmm(scale(log(dN)) ~ scale(abs(Lat)), prior=prior,
                     random = ~Species+Habitat+Clade,
                     data = subdata, verbose = TRUE, 
                     ginverse = list(Species = inverseA(phylo)$Ainv),
                     nitt = 60000, burnin = 10000, thin = 25, family = c("gaussian"))
  
  # Store models in list
  pglmm_list[[paste(group, "dS", sep = ".")]] <- fit_ds
  pglmm_list[[paste(group, "dN", sep = ".")]] <- fit_dn
  
  print(group)
}


#save(pglmm_list, file = "./Outputs/Data/PGLMM_latitudinal_pattern_3randoms.rdata")
load(file = "./Outputs/Data/PGLMM_latitudinal_pattern_3randoms.rdata")

summary_ds <-NULL
summary_dn <-NULL
for(group in groups){
  fit_ds <- pglmm_list[[paste(group, "dS", sep = ".")]]
  fit_dn <- pglmm_list[[paste(group, "dN", sep = ".")]]
  
  sum_ds <- summary(fit_ds)
  sum_dn <- summary(fit_dn)
  
  ds_fixed <- as.data.frame(sum_ds$solutions)
  ds_random <- as.data.frame(sum_ds$Gcovariances)%>%
    mutate(pMCMC=NA)
  
  ds_fixed <- data.frame(Var=row.names(ds_fixed), ds_fixed)
  ds_random <- data.frame(Var=row.names(ds_random), ds_random)
  summary_ds <- rbind(summary_ds, rbind(ds_fixed, ds_random)%>%mutate(Group=group))
  
  dn_fixed <- as.data.frame(sum_dn$solutions)
  dn_random <- as.data.frame(sum_dn$Gcovariances)%>%
    mutate(pMCMC=NA)
  
  dn_fixed <- data.frame(Var=row.names(dn_fixed), dn_fixed)
  dn_random <- data.frame(Var=row.names(dn_random), dn_random)
  summary_dn <- rbind(summary_dn, rbind(dn_fixed, dn_random)%>%mutate(Group=group))
}

write.csv(summary_ds,"summary_ds.csv", row.names = F)
write.csv(summary_dn,"summary_dn.csv", row.names = F)

# Extract parameter values from PGLMMs
pglmm_out <- NULL
for(i in 1:length(groups)){
  group <- groups[i]
  fit.out <- rbind(summary(pglmm_list[[paste(group, "dS", sep = ".")]])$solutions[2,], 
                   summary(pglmm_list[[paste(group, "dN", sep = ".")]])$solutions[2,])
  fit.out <- as.data.frame(fit.out) %>%
    mutate(Subrate = c("dS", "dN"),
           Group = group,
           ThermoMode = ifelse(Group %in% c("Birds", "Mammals", "Endotherms","Non-long migrants", "Non-migrants"), "Endotherms", "Ectotherms"))
  pglmm_out <- rbind(pglmm_out, fit.out)
}


#Plot effect size
f3 <-pglmm_out%>%
  mutate(Group=ifelse(Group=="Non-migrants", "Residents", Group))%>%
  mutate(sig=ifelse(pMCMC<0.05, "1", "0"),
         Group=factor(Group, levels=c("Endotherms","Ectotherms","Residents", "Non-long migrants", "Birds","Mammals","Reptiles", "Amphibians","Fishes")))%>%
  filter(Subrate=="dS")%>%
  ggplot(aes(x=Group, y=post.mean, ymin=`l-95% CI`, ymax=`u-95% CI`, color=ThermoMode, shape=sig))+
  geom_pointrange(size=0.4)+
  scale_shape_manual(labels=c("0", "1"), values = c(1,19))+
  scale_color_manual(values = c("#2a6aaf", "#d3292f"))+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  labs(y="Effects of latitude", x="", title = "dS")+
  theme_classic()+
  coord_flip()+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-0.3, 0.18),breaks = c(-0.3,-0.15,0,0.15), labels=c(-0.3,-0.15,0,0.15))+
  theme(axis.text = element_text(size=7),
        axis.title = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        legend.position = 'none',
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f4 <-pglmm_out%>%
  mutate(Group=ifelse(Group=="Non-migrants", "Residents", Group))%>%
  mutate(sig=ifelse(pMCMC<0.05, "1", "0"),
         Group=factor(Group, levels=c("Endotherms","Ectotherms","Residents", "Non-long migrants", "Birds","Mammals","Reptiles", "Amphibians","Fishes")))%>%
  filter(Subrate=="dN")%>%
  ggplot(aes(x=Group, y=post.mean, ymin=`l-95% CI`, ymax=`u-95% CI`, color=ThermoMode, shape=sig))+
  geom_pointrange(size=0.4)+
  scale_shape_manual(labels=c("0", "1"), values = c(1,19))+
  scale_color_manual(values = c("#2a6aaf", "#d3292f"))+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  labs(y="Effects of latitude", x="", title = "dN")+
  theme_classic()+
  coord_flip()+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-0.3, 0.15),breaks = c(-0.3,-0.15,0,0.15), labels=c(-0.3,-0.15,0,0.15))+
  theme(axis.text = element_text(size=7),
        axis.title = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        legend.position = 'none',
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))
cowplot::plot_grid(f3,f4, nrow = 2, align="hv")

#############
#Fig.S2 Scatterplots show the relationships between molecular rates and 
#absolute midpoint latitude at the species level
f5 <- subrate%>%
  ggplot(aes(x=abs(Lat), y=log10(dS), colour=ThermoMode))+
  geom_point(size=0.5, alpha=0.5)+
  geom_smooth(method="lm", se=TRUE, size=0.5)+
  scale_colour_manual(values = c("#2a6aaf", "#d3292f"))+
  labs(x="", y=expression(log[10](dS)))+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-9.1, -7),breaks = seq(-9.0, -7, 0.5))+
  scale_x_continuous(limits = c(0, 90),breaks = seq(0,90,30))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=9),
        legend.position = "none",
        axis.title = element_text(size=9),
        axis.text = element_text(size=8),
        strip.text = element_text(size=9),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

f6 <- subrate%>%
  ggplot(aes(x=abs(Lat), y=log10(dN), colour=ThermoMode))+
  geom_point(size=0.5, alpha=0.5)+
  geom_smooth(method="lm", se=TRUE, size=0.5)+
  scale_colour_manual(values = c("#2a6aaf", "#d3292f"))+
  labs(x="", y=expression(log[10](dN)))+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-10.4, -8.8),breaks = seq(-10.4, -8.8, 0.4))+
  scale_x_continuous(limits = c(0, 80),breaks = seq(0,80,20))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=9),
        axis.title = element_text(size=9),
        axis.text = element_text(size=8),
        legend.position = c(0.8,0.9),
        legend.key.size = unit(0.3, "cm"),
        legend.title = element_blank(),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

cowplot::plot_grid(f5,f6, align="hv")
ggsave(filename = "./Outputs/Supplementary/Extended Fig2.pdf", width=8.27, height = 4.2)



f1 <- subrate%>%
  filter(Group=="Fishes")%>%
  ggplot(aes(x=abs(Lat), y=log10(dS)))+
  geom_point(size=0.5, alpha=0.5, colour="#2a6aaf")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#2a6aaf")+
  labs(x="", y=expression(log[10](dS)), title="Fishes")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-9.1, -7),breaks = seq(-9.0, -7, 0.5))+
  scale_x_continuous(limits = c(0, 90),breaks = seq(0,90,30))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=9),
        axis.title = element_text(size=9),
        axis.text = element_text(size=8),
        strip.text = element_text(size=9),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

f2 <- subrate%>%
  filter(Group=="Amphibians")%>%
  ggplot(aes(x=abs(Lat), y=log10(dS)))+
  geom_point(size=0.5, alpha=0.5, colour="#2a6aaf")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#2a6aaf")+
  labs(x="", y="", title="Amphibians")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-9, -7.75),breaks = seq(-9, -7.75, 0.25))+
  #scale_x_continuous(limits = c(0, 90),breaks = seq(0,90,30))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=9),
        axis.title = element_text(size=9),
        axis.text = element_text(size=8),
        strip.text = element_text(size=9),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))


f3 <- subrate%>%
  filter(Group=="Reptiles")%>%
  ggplot(aes(x=abs(Lat), y=log10(dS)))+
  geom_point(size=0.5, alpha=0.5, colour="#2a6aaf")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#2a6aaf")+
  labs(x="", y="", title="Reptiles")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-9.5, -7),breaks = seq(-9.5, -7, 0.5))+
  scale_x_continuous(limits = c(0, 60),breaks = seq(0,60,20))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=9),
        axis.title = element_text(size=9),
        axis.text = element_text(size=8),
        strip.text = element_text(size=9),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

f4 <- subrate%>%
  filter(Group=="Mammals")%>%
  ggplot(aes(x=abs(Lat), y=log10(dS)))+
  geom_point(size=0.5, alpha=0.5, colour="#d3292f")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#d3292f")+
  labs(x="", y="", title="Mammals")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-8.6, -7.2),breaks = seq(-8.6, -7.2, 0.4))+
  scale_x_continuous(limits = c(0, 80),breaks = seq(0,80,20))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=9),
        axis.title = element_text(size=9),
        axis.text = element_text(size=8),
        strip.text = element_text(size=9),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

f5 <- subrate%>%
  filter(Group=="Birds")%>%
  ggplot(aes(x=abs(Lat), y=log10(dS)))+
  geom_point(size=0.5, alpha=0.5, colour="#d3292f")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#d3292f")+
  labs(x="", y="", title="Birds")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  #scale_y_continuous(limits = c(-20, -16),breaks = seq(-20, -16, 1))+
  scale_x_continuous(limits = c(0, 80),breaks = seq(0,80,20))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=9),
        axis.title = element_text(size=9),
        axis.text = element_text(size=8),
        strip.text = element_text(size=9),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

f6 <- subrate%>%
  filter(Group=="Fishes")%>%
  ggplot(aes(x=abs(Lat), y=log10(dN)))+
  geom_point(size=0.5, alpha=0.5, colour="#2a6aaf")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#2a6aaf")+
  labs(x="Absolute Latitude", y=expression(log[10](dN)))+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-10.5, -9),breaks = seq(-10.5, -9, 0.5))+
  #scale_x_continuous(limits = c(0, 90),breaks = seq(0,90,30))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=9),
        axis.title = element_text(size=9),
        axis.text = element_text(size=8),
        strip.text = element_text(size=9),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

f7 <- subrate%>%
  filter(Group=="Amphibians")%>%
  ggplot(aes(x=abs(Lat), y=log10(dN)))+
  geom_point(size=0.5, alpha=0.5, colour="#2a6aaf")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#2a6aaf")+
  labs(x="Absolute Latitude", y="")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-10.7, -9),breaks = seq(-10.5, -9, 0.5))+
  #scale_x_continuous(limits = c(0, 90),breaks = seq(0,90,30))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=9),
        axis.title = element_text(size=9),
        axis.text = element_text(size=8),
        strip.text = element_text(size=9),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))


f8 <- subrate%>%
  filter(Group=="Reptiles")%>%
  ggplot(aes(x=abs(Lat), y=log10(dN)))+
  geom_point(size=0.5, alpha=0.5, colour="#2a6aaf")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#2a6aaf")+
  labs(x="Absolute Latitude", y="")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  #scale_y_continuous(limits = c(-23, -21),breaks = seq(-23, -21, 0.5))+
  scale_x_continuous(limits = c(0, 60),breaks = seq(0,60,20))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=9),
        axis.title = element_text(size=9),
        axis.text = element_text(size=8),
        strip.text = element_text(size=9),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

f9 <- subrate%>%
  filter(Group=="Mammals")%>%
  ggplot(aes(x=abs(Lat), y=log10(dN)))+
  geom_point(size=0.5, alpha=0.5, colour="#d3292f")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#d3292f")+
  labs(x="Absolute Latitude", y="")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-9.8, -8.8),breaks = seq(-9.8, -8.8, 0.2))+
  #scale_x_continuous(limits = c(0, 60),breaks = seq(0,60,20))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=9),
        axis.title = element_text(size=9),
        axis.text = element_text(size=8),
        strip.text = element_text(size=9),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

f10 <- subrate%>%
  filter(Group=="Birds")%>%
  ggplot(aes(x=abs(Lat), y=log10(dN)))+
  geom_point(size=0.5, alpha=0.5, colour="#d3292f")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#d3292f")+
  labs(x="Absolute Latitude", y="")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-10, -9),breaks = seq(-10, -9, 0.2))+
  #scale_x_continuous(limits = c(0, 60),breaks = seq(0,60,20))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=9),
        axis.title = element_text(size=9),
        axis.text = element_text(size=8),
        strip.text = element_text(size=9),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

cowplot::plot_grid(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10, nrow = 2, align = "hv")
ggsave(filename = "./Outputs/Supplementary/Extended Fig3.pdf", width=8.27, height = 3.6)




#############################################################
#2.Latitudinal gradients in molecular rates across space
rm(list=ls())
gc()
workdir <- "/Users/Tianlong/VertMolRate"
setwd(workdir)
#Loading packages
library(raster)
library(tidyverse)
library(RColorBrewer)
library(spdep)
library(sf)
library(spatialreg)
source(paste0(workdir, "/Scripts/source_functions.r"))

#input spatial join data and grids
load("./DataFiles/SpatialJoinFiles/grids.rdata")
load("./DataFiles/SpatialJoinFiles/all_spatial_join_grid.rdata")
load("./DataFiles/SpatialJoinFiles/all_spatial_join_ecoregion.rdata")




#input GIS poly and raster layers
marine.raster <- raster("./DataFiles/GISLayers/marine.tif")
terrestrial.raster <- raster("./DataFiles/GISLayers/terrestrial.tif")
outline.poly <- read_sf("./DataFiles/GISLayers/outline_moll.shp")
country.poly <- read_sf("./DataFiles/GISLayers/country.poly.merged.shp")
ecos.poly <- read.csv("./DataFiles/GISLayers/ecos.poly.csv")
terrestrial.ecoregions <- read.csv("./DataFiles/GISLayers/terrestrial.ecoregions.csv")
marine.ecoregions <- read.csv("./DataFiles/GISLayers/marine.ecoregions.csv")


#estimated mean substitution rate in grids
fresh.fishes <-mean_subrate_grid(spatial_join=all_spatial_join_grid%>%filter(Group=="Fishes"&(Habitat!="Marine")),
                               grids_poly=grids, Species=5, habitat="Terrestrial")
marine.fishes <- mean_subrate_grid(spatial_join=all_spatial_join_grid%>%filter(Group=="Fishes"&(Habitat!="Freshwater")),
                                 grids_poly=grids, Species=5, habitat="Marine")
terr.reptiles <- mean_subrate_grid(spatial_join=all_spatial_join_grid%>%filter(Group=="Reptiles"&Habitat!="Marine"),
                                   grids_poly=grids, Species=4, habitat="Terrestrial")
marine.reptiles <- mean_subrate_grid(spatial_join=all_spatial_join_grid%>%filter(Group=="Reptiles"&Habitat!="Terrestrial"),
                                     grids_poly=grids, Species=4, habitat="Marine")
amphibians <- mean_subrate_grid(spatial_join=all_spatial_join_grid%>%filter(Group=="Amphibians"),
                              grids_poly=grids, Species=4, habitat="Terrestrial")
terr.birds <- mean_subrate_grid(spatial_join=all_spatial_join_grid%>%filter(Group=="Birds"&Habitat!="Marine"),
                               grids_poly=grids, Species=5, habitat="Terrestrial")
marine.birds <- mean_subrate_grid(spatial_join=all_spatial_join_grid%>%filter(Group=="Birds"&Habitat!="Terrestrial"),
                                 grids_poly=grids, Species=5, habitat="Marine")
terr.mammals <- mean_subrate_grid(spatial_join=all_spatial_join_grid%>%filter(Group=="Mammals"&Habitat!="Marine"),
                                 grids_poly=grids, Species=5, habitat="Terrestrial")
marine.mammals <- mean_subrate_grid(spatial_join=all_spatial_join_grid%>%filter(Group=="Mammals"&Habitat!="Terrestrial"),
                                   grids_poly=grids, Species=5, habitat="Marine")
terr.ectotherm <- mean_subrate_grid(spatial_join=all_spatial_join_grid%>%filter((Group %in% c("Amphibians", "Reptiles","Fishes"))&Habitat!="Marine"),
                                    grids_poly=grids, Species=5, habitat="Terrestrial")
marine.ectotherm <- mean_subrate_grid(spatial_join=all_spatial_join_grid%>%filter((Group %in% c("Amphibians", "Reptiles","Fishes"))&(Habitat=="Marine"|Habitat=="Marine&Freshwater"|Habitat=="Terrestrial&Marine")),
                                      grids_poly=grids, Species=5, habitat="Marine")
terr.endotherm <- mean_subrate_grid(spatial_join=all_spatial_join_grid%>%filter((Group %in% c("Mammals", "Birds"))&Habitat!="Marine"),
                                    grids_poly=grids, Species=5, habitat="Terrestrial")
marine.endotherm <- mean_subrate_grid(spatial_join=all_spatial_join_grid%>%filter((Group %in% c("Mammals", "Birds"))&Habitat!="Terrestrial"),
                                      grids_poly=grids, Species=5, habitat="Marine")

fishes <- bind_rows(list(Marine=marine.fishes, Terrestrial=fresh.fishes), .id = "Habitat")
reptiles <- bind_rows(list(Marine=marine.reptiles, Terrestrial=terr.reptiles), .id = "Habitat")
amphibians <- amphibians%>%mutate(Habitat="Terrestrial")
birds <- bind_rows(list(Marine=marine.birds, Terrestrial=terr.birds), .id = "Habitat")
mammals <- bind_rows(list(Marine=marine.mammals, Terrestrial=terr.mammals), .id = "Habitat")

ectotherm <- bind_rows(list(Marine=marine.ectotherm, Terrestrial=terr.ectotherm), .id = "Habitat")
endotherm <- bind_rows(list(Marine=marine.endotherm, Terrestrial=terr.endotherm), .id = "Habitat")


all.groups <- bind_rows(list(Fishes=fishes,
                             Amphibians=amphibians,
                             Reptiles=reptiles,
                             Birds=birds,
                             Mammals=mammals), .id="Group")%>%
  mutate(Group=factor(Group, levels = c("Fishes","Amphibians", "Reptiles","Mammals","Birds")))


#estimated mean subrate in ecoregions
fresh.fishes.ecos <- mean_subrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter(Group=="Fishes"&System=="Terrestrial"),
                                          Species=5, habitat="Terrestrial")%>%filter(!is.na(SR))
marine.fishes.ecos <- mean_subrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter(Group=="Fishes"&System=="Marine"),
                                           Species=5, habitat="Marine")%>%filter(!is.na(SR))
terr.reptiles.ecos <- mean_subrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter(Group=="Reptiles"&System=="Terrestrial"),
                                             Species=3, habitat="Terrestrial")%>%filter(!is.na(SR))
marine.reptiles.ecos <- mean_subrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter(Group=="Reptiles"&System=="Marine"),
                                               Species=3, habitat="Marine")%>%filter(!is.na(SR))
amphibians.ecos <- mean_subrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter(Group=="Amphibians"&System=="Terrestrial"),
                                        Species=3, habitat="Terrestrial")%>%filter(!is.na(SR))

terr.birds.ecos <- mean_subrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter(Group=="Birds"&System=="Terrestrial"),
                                         Species=5, habitat="Terrestrial")%>%filter(!is.na(SR))
marine.birds.ecos <- mean_subrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter(Group=="Birds"&System=="Marine"),
                                           Species=5, habitat="Marine")%>%filter(!is.na(SR))
terr.mammals.ecos <- mean_subrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter(Group=="Mammal"&System=="Terrestrial"),
                                           Species=5, habitat="Terrestrial")%>%filter(!is.na(SR))

marine.mammals.ecos <- mean_subrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter(Group=="Mammal"&System=="Marine"),
                                             Species=5, habitat="Marine")%>%filter(!is.na(SR))

terr.ectotherm.ecos <- mean_subrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter((Group %in% c("Amphibians", "Reptiles","Fishes"))&System=="Terrestrial"),
                                              Species=5, habitat="Terrestrial")%>%filter(!is.na(SR))
marine.ectotherm.ecos <- mean_subrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter((Group %in% c("Amphibians", "Reptiles","Fishes"))&System=="Marine"),
                                                Species=5, habitat="Marine")%>%filter(!is.na(SR))

terr.endotherm.ecos <- mean_subrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter((Group %in% c("Birds", "Mammal"))&System=="Terrestrial"),
                                              Species=5, habitat="Terrestrial")%>%filter(!is.na(SR))

marine.endotherm.ecos <- mean_subrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter((Group %in% c("Birds", "Mammal"))&System=="Marine"),
                                                Species=5, habitat="Marine")%>%filter(!is.na(SR))

fishes.ecos <- bind_rows(list(Marine=marine.fishes.ecos, Terrestrial=fresh.fishes.ecos), .id = "Habitat")
reptiles.ecos <- bind_rows(list(Marine=marine.reptiles.ecos, Terrestrial=terr.reptiles.ecos), .id = "Habitat")
amphibians.ecos <- amphibians.ecos%>%mutate(Habitat="Terrestrial")
birds.ecos <- bind_rows(list(Marine=marine.birds.ecos, Terrestrial=terr.birds.ecos), .id = "Habitat")
mammals.ecos <- bind_rows(list(Marine=marine.mammals.ecos, Terrestrial=terr.mammals.ecos), .id = "Habitat")

ectotherm.ecos <- bind_rows(list(Marine=marine.ectotherm.ecos, Terrestrial=terr.ectotherm.ecos), .id = "Habitat")
endotherm.ecos <- bind_rows(list(Marine=marine.endotherm.ecos, Terrestrial=terr.endotherm.ecos), .id = "Habitat")




all.groups.ecos <- bind_rows(list(Fishes=fishes.ecos, Amphibians=amphibians.ecos,
                                  Reptiles=reptiles.ecos%>%filter(Habitat=="Terrestrial"),Birds=birds.ecos%>%filter(Habitat=="Terrestrial"), Mammals=mammals.ecos%>%filter(Habitat=="Terrestrial")), .id="Group")%>%
  mutate(Group=factor(Group, levels = c("Fishes","Amphibians", "Reptiles","Mammals","Birds")))


mean.subrate.ecos <- bind_rows(list(Fishes=fishes.ecos,
                                    Amphibians=amphibians.ecos,
                                    Reptiles=reptiles.ecos%>%filter(Habitat=="Terrestrial"),
                                    Birds=birds.ecos%>%filter(Habitat=="Terrestrial"), 
                                    Mammals=mammals.ecos%>%filter(Habitat=="Terrestrial"), 
                                    Endotherms=endotherm.ecos%>%filter(Habitat=="Terrestrial"), 
                                    Ectotherms=ectotherm.ecos), .id="Group")%>%
  mutate(Group=factor(Group, levels = c("Fishes","Amphibians", "Reptiles","Mammals","Birds", "Ectotherms", "Endotherms")))


#############################################################
#Spatial simultaneous autoregressive (SAR) models
############################################################
fit.vars <- c("gm.ds", "gm.dn",  "mid.ds", "mid.dn")
sar.summary.out <- NULL

for(id in unique(mean.subrate.ecos$Group)){
  out <- as.list(rep(NA, length(fit.vars)))
  names(out) <- fit.vars
  
  
  mean.subrate <- mean.subrate.ecos%>%filter(Group==id)%>%mutate(abs.lat=abs(Lat))
  
  nc.coords <- cbind(mean.subrate$Lon, mean.subrate$Lat)
  
  nc.5nn <- knearneigh(nc.coords, k=5, longlat = TRUE)
  nc.5nn.nb <- knn2nb(nc.5nn)
  
  #fit SAR model
  for (i in 1:length(fit.vars)) {
    var <- fit.vars[i]
    
    # fit ols
    fit.ols <- lm(scale(mean.subrate[,var,drop=T]) ~ scale(mean.subrate[,"abs.lat",drop=T]))
    
    # build SAR models with different combinations of distance and neighbourhood style
    neiStyle <- c('W','B','S', 'C', 'U', 'minmax')
    AICvec <- numeric(length = length(neiStyle))
    for (k in 1:length(neiStyle)) {
      nlw <- nb2listw(nc.5nn.nb, style=neiStyle[k], zero.policy=TRUE)
      sar_e <- lagsarlm(fit.ols, listw = nlw, na.action = na.omit, zero.policy = TRUE, tol.solve = 1.0e-25)
      AICvec[k] <- AIC(sar_e)			
    }
    bestStyle <- neiStyle[which.min(AICvec)]
    nlw <- nb2listw(nc.5nn.nb, style=bestStyle, zero.policy=TRUE)
    
    fit.sar <- lagsarlm(fit.ols, listw = nlw, na.action = na.omit, zero.policy = TRUE, tol.solve = 1.0e-25)
    mtest1 <- moran.test(resid(fit.ols), listw=nlw, na.action = na.omit, zero.policy = TRUE)
    mtest2 <- moran.test(resid(fit.sar), listw=nlw, na.action = na.omit, zero.policy = TRUE)
    aic.ols <- AIC(fit.ols)
    aic.sar <- AIC(fit.sar)
    res <- list(fit.ols = fit.ols, fit.sar = fit.sar, aic.ols = aic.ols, aic.sar = aic.sar, moran.ols = mtest1, moran.sar = mtest2)
    out[[i]] <- res
  }
  var.sp <- strsplit(gsub("","",fit.vars), split="\\.")
  names(var.sp) <- paste0("a", 1:4)
  
  var.sp<-var.sp%>%bind_rows()%>%t()%>%as.data.frame()%>%
    mutate(V1=ifelse(V1=="gm", "Geometric mean", ifelse(V1=="am", "Arithmetic mean", "Median")),
           V2=ifelse(V2=="ds", "dS", "dN"))
  
  # Summary
  x <- numeric(length(fit.vars))
  dff <- data.frame(Group=id, Avg.Value =var.sp[,1], SubRate = var.sp[,2], AIC.OLS=x, AIC.SAR=x, dAIC=x, 
                    OLS.Slope=x, OLS.Slope.SD=x, OLS.Pvalue=x, SAR.Slope=x, SAR.Slope.SD=x, SAR.Pvalue=x, 
                    Moran.OLS = x, Moran.OLS.Pvalue = x, Moran.SAR = x, Moran.SAR.Pvalue = x, stringsAsFactors=F)
  
  for (i in 1:length(out)) {
    
    fres <- out[[i]]
    
    dff$AIC.OLS[i] <- fres$aic.ols
    dff$AIC.SAR[i] <- fres$aic.sar
    dff$dAIC[i] <- fres$aic.ols - fres$aic.sar
    dff$OLS.Slope[i] <- summary(fres$fit.ols)$coef[2,1]
    dff$OLS.Slope.SD[i] <- summary(fres$fit.ols)$coef[2,2]
    dff$OLS.Pvalue[i] <- summary(fres$fit.ols)$coef[2,4]
    dff$SAR.Slope[i] <- summary(fres$fit.sar)$Coef[2,1]
    dff$SAR.Slope.SD[i] <- summary(fres$fit.sar)$Coef[2,2]
    dff$SAR.Pvalue[i] <- summary(fres$fit.sar)$Coef[2,4]
    dff$Moran.OLS[i] <- fres$moran.ols$estimate[1]
    dff$Moran.OLS.Pvalue[i] <- fres$moran.ols$p.value
    dff$Moran.SAR[i] <- fres$moran.sar$estimate[1]
    dff$Moran.SAR.Pvalue[i] <- fres$moran.sar$p.value
  }
  sar.summary.out <- rbind(sar.summary.out, dff)
  print(id)
}

#
save(sar.summary.out, file = "./Outputs/Data/SAR_latitudinal_pattern_assemblages.rdata")
load(file = "./Outputs/Data/SAR_latitudinal_pattern_assemblages.rdata")

write.csv(sar.summary.out, "sar.summary.out.csv")


########################################
source(paste0(workdir, "/SourceFunctions/source_functions.r"))

pdf("./Outputs/MainFigures/Fig2df.pdf", height = 2.8, width=3.8)
par(mar=c(0.6,0.6,0.6,0.6), mfrow=c(2,2))
plot_subrate_grids(mean.rate.grids=ectotherm%>%mutate(dS=gm.ds)%>%filter(!is.na(SR)), 
                   rate.type = "dS", geo.type="Global", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
mtext('dS', side=2, at=-3.0e+05, cex=0.6, line=-0.1, font=2)
title(main="Ectotherms", cex.main=0.8, line=-0.3)

plot_subrate_grids(mean.rate.grids=endotherm%>%mutate(dS=gm.ds)%>%filter(Habitat!="Marine"), 
                   rate.type = "dS", geo.type="Terrestrial", colramp=c("#ffffff", "#f4cdb9","#f47155","#d22428", "#660b11"))
#mtext('Endotherms', side=2, at=-3.0e+05, cex=0.6, line=-0.1, font=2)
title(main="Endotherms", cex.main=0.8, line=-0.3)

plot_subrate_grids(mean.rate.grids=ectotherm%>%mutate(dN=gm.dn), 
                   rate.type = "dN", geo.type="Global", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
title(main="Ectotherms", cex.main=0.8, line=-0.3)
mtext('dN', side=2, at=-3.0e+05, cex=0.6, line=-0.1, font=2)


plot_subrate_grids(mean.rate.grids=endotherm%>%mutate(dN=gm.dn)%>%filter(Habitat!="Marine"), 
                   rate.type = "dN", geo.type="Terrestrial", colramp=c("#ffffff", "#f4cdb9","#f47155","#d22428", "#660b11"))
title(main="Endotherms", cex.main=0.8, line=-0.3)
dev.off()


#Fig.2eg 
f1<-ectotherm.ecos%>%
  filter(!is.na(gm.ds))%>%
  ggplot(aes(x=abs(Lat), y=gm.ds*100))+
  geom_point(alpha=0.7, size=0.7, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y=expression("dS (" * 10^-8 * ")"), title="Ectotherms")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = 60, y = 1, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f2<-ectotherm.ecos%>%
  filter(!is.na(gm.dn))%>%
  ggplot(aes(x=abs(Lat), y=gm.dn*10000))+
  geom_point(alpha=0.7, size=0.7, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y=expression("dN (" * 10^-10 * ")"), title="Ectotherms")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = 60, y = 3, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))


f3<-endotherm.ecos%>%
  filter(!is.na(gm.ds), Habitat!="Marine")%>%
  ggplot(aes(x=abs(Lat), y=gm.ds*100))+
  geom_point(alpha=0.7, size=0.7, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="Latitude", y=expression("dS (" * 10^-8 * ")"), title="Endotherms")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = 60, y = 0.9, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f4<-endotherm.ecos%>%
  filter(!is.na(gm.dn), Habitat!="Marine")%>%
  ggplot(aes(x=abs(Lat), y=gm.dn*10000))+
  geom_point(alpha=0.7, size=0.7, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="Latitude", y=expression("dN (" * 10^-10 * ")"), title="Endotherms")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = 60, y = 2.85, label = expression(italic(P)[SAR] == 0.81), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

cowplot::plot_grid(f1,f3,f2,f4, align = "hv", nrow = 2)
#ggsave("./Outputs/MainFigures/Fig2c.pdf", height = 3.6, width=2.7)
ggsave("./Outputs/MainFigures/Fig2eg.pdf", height = 3.6, width=3.4)

###########################
#Extended Fig.5a Latitudinal gradients in molecular rates at assemblage level for each class.
##########################
pdf("./Outputs/Supplementary/Extended Fig5a.pdf", width=8.27, height = 2.3)
par(mar=c(0.6,0.6,0.6,0.6), mfrow=c(2,5))
plot_subrate_grids(mean.rate.grids=fishes%>%mutate(dS=gm.ds), rate.type = "dS", 
                   geo.type="Global", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
mtext('Geometric mean dS', side=2, at=-3e+05, cex=0.6, line=-0.3, font=2)
#mtext("a", side=3,font = 2, at=-1.5e+07, cex=0.7, line=-0.8)
title(main="Fishes", cex.main=0.8, line=-0.3)

plot_subrate_grids(mean.rate.grids=amphibians%>%mutate(dS=gm.ds), rate.type = "dS", 
                   geo.type="Terrestrial", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
mtext('', side=2, at=-3.0e+05, cex=0.6, line=-0.3, font=2)
#mtext("b", side=3,font = 2, at=-1.5e+07, cex=0.7, line=-0.8)
title(main="Amphibians", cex.main=0.8, line=-0.3)

plot_subrate_grids(mean.rate.grids=terr.reptiles%>%mutate(dS=gm.ds), rate.type = "dS", 
                   geo.type="Terrestrial", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
mtext('', side=2, at=-3.0e+05, cex=0.6, line=-0.3, font=2)
#mtext("c", side=3,font = 2, at=-1.5e+07, cex=0.7, line=-0.8)
title(main="Reptiles", cex.main=0.8, line=-0.3)

plot_subrate_grids(mean.rate.grids=terr.mammals%>%mutate(dS=gm.ds), rate.type = "dS", 
                   geo.type="Terrestrial", colramp=c("#ffffff", "#f4cdb9","#f47155","#d22428", "#660b11"))
mtext('', side=2, at=-3.0e+05, cex=0.6, line=-0.3, font=2)
#mtext("d", side=3,font = 2, at=-1.5e+07, cex=0.7, line=-0.8)
title(main="Mammals", cex.main=0.8, line=-0.3)

plot_subrate_grids(mean.rate.grids=terr.birds%>%mutate(dS=gm.ds), rate.type = "dS", 
                   geo.type="Terrestrial", colramp=c("#ffffff", "#f4cdb9","#f47155","#d22428", "#660b11"))
mtext('', side=2, at=-3.0e+05, cex=0.6, line=-0.3, font=2)
#mtext("e", side=3,font = 2, at=-1.5e+07, cex=0.7, line=-0.8)
title(main="Birds", cex.main=0.8, line=-0.3)

plot_subrate_grids(mean.rate.grids=fishes%>%mutate(dN=gm.dn), rate.type = "dN", 
                   geo.type="Global", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
mtext('Geometric mean dN', side=2, at=-3e+05, cex=0.6, line=-0.3, font=2)
title(main="Fishes", cex.main=0.8, line=-0.3)

plot_subrate_grids(mean.rate.grids=amphibians%>%mutate(dN=gm.dn), rate.type = "dN", 
                   geo.type="Terrestrial", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
title(main="Amphibians", cex.main=0.8, line=-0.3)

plot_subrate_grids(mean.rate.grids=terr.reptiles%>%mutate(dN=gm.dn), rate.type = "dN", 
                   geo.type="Terrestrial", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
title(main="Reptiles", cex.main=0.8, line=-0.3)

plot_subrate_grids(mean.rate.grids=terr.mammals%>%mutate(dN=gm.dn), rate.type = "dN", 
                   geo.type="Terrestrial", colramp=c("#ffffff", "#f4cdb9","#f47155","#d22428", "#660b11"))
title(main="Mammals", cex.main=0.8, line=-0.3)

plot_subrate_grids(mean.rate.grids=terr.birds%>%mutate(dN=gm.dn), rate.type = "dN", 
                   geo.type="Terrestrial", colramp=c("#ffffff", "#f4cdb9","#f47155","#d22428", "#660b11"))
title(main="Birds", cex.main=0.8, line=-0.3)
dev.off()


##############
#Extended Fig.5b Latitudinal gradients in molecular rates at assemblage level for each class.
f1 <- mean.subrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Fishes")%>%
  ggplot(aes(x=abs(Lat), y=gm.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y=expression("dS (" * 10^-8 * ")"), title="Fishes")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = 60, y = 1, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

f2 <- mean.subrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Amphibians", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=gm.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y="", title="Amphibians")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = 40, y = 1, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))


f3 <- mean.subrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Reptiles", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=gm.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y="", title="Reptiles")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = 40, y = 1.2, label = expression(italic(P)[SAR] == 0.0023), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

f4 <- mean.subrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Mammals", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=gm.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="Latitude", y="", title="Mammals")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.9,2.1), breaks = seq(0.9,2.1,0.2))+
  annotate("text", x = 60, y = 2.1, label = expression(italic(P)[SAR] == 0.96), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

f5 <- mean.subrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Birds", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=gm.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="Latitude", y="", title="Birds")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.8,1.5), breaks = seq(0.8,1.5,0.1))+
  annotate("text", x = 60, y = 1.5, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

f6 <- mean.subrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Fishes")%>%
  ggplot(aes(x=abs(Lat), y=gm.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y=expression("dN (" * 10^-10 * ")"), title="Fishes")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = 60, y = 2.8, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

f7 <- mean.subrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Amphibians", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=gm.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y="", title="Amphibians")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.6,3), breaks = seq(0.6,3,0.4))+
  annotate("text", x = 40, y = 3, label = expression(italic(P)[SAR] == 0.015), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

f8 <- mean.subrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Reptiles", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=gm.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y="", title="Reptiles")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(2.2,3.6), breaks = seq(2.2,3.6,0.2))+
  annotate("text", x = 40, y = 3.6, label = expression(italic(P)[SAR] == 0.031), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

f9 <- mean.subrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Mammals", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=gm.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="Latitude", y="", title="Mammals")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(3.9,5.5), breaks = seq(4,5.5,0.5))+
  annotate("text", x = 60, y = 5.5, label = expression(italic(P)[SAR] == 0.0042), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

f10 <- mean.subrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Birds", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=gm.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="Latitude", y="", title="Birds")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = 60, y = 3.3, label = expression(italic(P)[SAR] == 0.71), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

cowplot::plot_grid(f1,f2, f3, f4, f5, 
                   f6,f7, f8, f9, f10,
                   nrow = 2, align="hv")

ggsave("./Outputs/Supplementary/Extended Fig5b.pdf", width=8.27, height=3.2)

########################
#Extended Fig.6 Plot arithmetic mean values of molecular rates in cell grids for endotherms and ectotherms
#Extended Fig.6a
pdf("./Outputs/Supplementary/Extended Fig6a.pdf", height =1.4, width= 8.27)
par(mar=c(0.6,0.6,0.6,0.6), mfrow=c(1,4))
plot_subrate_grids(mean.rate.grids=ectotherm%>%mutate(dS=mid.ds)%>%filter(!is.na(SR)), 
                   rate.type = "dS", geo.type="Global", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
mtext('Ectotherms', side=2, at=-3.0e+05, cex=0.6, line=-0.1, font=2)
#mtext("a", side=3,font = 2, at=-1.5e+07, cex=0.7, line=-0.8)
title(main="Median dS", cex.main=0.8, line=-0.3)

plot_subrate_grids(mean.rate.grids=ectotherm%>%mutate(dN=mid.dn), 
                   rate.type = "dN", geo.type="Global", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
title(main="Median dN", cex.main=0.8, line=-0.3)
mtext('Ectotherms', side=2, at=-3.0e+05, cex=0.6, line=-0.1, font=2)
#mtext("b", side=3,font = 2, at=-1.5e+07, cex=0.7, line=-0.8)

plot_subrate_grids(mean.rate.grids=endotherm%>%mutate(dS=mid.ds)%>%filter(Habitat!="Marine"), 
                   rate.type = "dS", geo.type="Terrestrial", colramp=c("#ffffff", "#f4cdb9","#f47155","#d22428", "#660b11"))
mtext('Endotherms', side=2, at=-3.0e+05, cex=0.6, line=-0.1, font=2)
#mtext("c", side=3,font = 2, at=-1.5e+07, cex=0.7, line=-0.8)
title(main="Median dS", cex.main=0.8, line=-0.3)

plot_subrate_grids(mean.rate.grids=endotherm%>%mutate(dN=mid.dn)%>%filter(Habitat!="Marine"), 
                   rate.type = "dN", geo.type="Terrestrial", colramp=c("#ffffff", "#f4cdb9","#f47155","#d22428", "#660b11"))
title(main="Median dN", cex.main=0.8, line=-0.3)
mtext('Endotherms', side=2, at=-3.0e+05, cex=0.6, line=-0.1, font=2)
#mtext("d", side=3,font = 2, at=-1.5e+07, cex=0.7, line=-0.8)
dev.off()

#Extended Fig.6b 
f1<-ectotherm.ecos%>%
  filter(!is.na(mid.ds))%>%
  ggplot(aes(x=abs(Lat), y=mid.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y=expression("Median dS (" * 10^-8 * ")"), title="Ectotherms")+
  scale_y_continuous(limits = c(0.4, 1.2),breaks = seq(0.4, 1.2, 0.2))+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = 60, y = 1.2, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=9),
        axis.text = element_text(size=9),
        plot.title = element_text(hjust = 0.5, size=9),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f2<-ectotherm.ecos%>%
  filter(!is.na(mid.dn))%>%
  ggplot(aes(x=abs(Lat), y=mid.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y=expression("Median dN (" * 10^-10 * ")"), title="Ectotherms")+
  scale_y_continuous(limits = c(1.4, 3.5),breaks = seq(1.5, 3.5, 1))+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = 60, y = 3.5, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=9),
        axis.text = element_text(size=9),
        plot.title = element_text(hjust = 0.5, size=9),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))


f3<-endotherm.ecos%>%
  filter(!is.na(mid.ds), Habitat!="Marine")%>%
  ggplot(aes(x=abs(Lat), y=mid.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="Latitude", y=expression("Median dS (" * 10^-8 * ")"), title="Endotherms")+
  scale_y_continuous(limits = c(0.8, 1.6),breaks = seq(0.8, 1.6, 0.2))+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = 60, y = 1.6, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=9),
        axis.text = element_text(size=9),
        plot.title = element_text(hjust = 0.5, size=9),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f4<-endotherm.ecos%>%
  filter(!is.na(mid.dn), Habitat!="Marine")%>%
  ggplot(aes(x=abs(Lat), y=mid.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="Latitude", y=expression("Median dN (" * 10^-10 * ")"), title="Endotherms")+
  scale_y_continuous(limits = c(2.5, 4),breaks = seq(2.5, 4, 0.5))+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = 60, y = 4, label = expression(italic(P)[SAR] == 0.23), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=9),
        axis.text = element_text(size=9),
        plot.title = element_text(hjust = 0.5, size=9),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

cowplot::plot_grid(f1,f2,f3,f4, align = "hv", nrow = 1)
ggsave("./Outputs/Supplementary/Extended Fig6b.pdf", height =2, width= 8.27)


##################
###########################
#Extended Fig.7a Latitudinal gradients in molecular rates at assemblage level for each class.
##########################
pdf("./Outputs/Supplementary/Extended Fig7a.pdf", width=8.27, height = 2.3)
par(mar=c(0.6,0.6,0.6,0.6), mfrow=c(2,5))
plot_subrate_grids(mean.rate.grids=fishes%>%mutate(dS=mid.ds), rate.type = "dS", 
                   geo.type="Global", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
mtext('Median dS', side=2, at=-3e+05, cex=0.6, line=-0.3, font=2)
#mtext("a", side=3,font = 2, at=-1.5e+07, cex=0.7, line=-0.8)
title(main="Fishes", cex.main=0.8, line=-0.3)

plot_subrate_grids(mean.rate.grids=amphibians%>%mutate(dS=mid.ds), rate.type = "dS", 
                   geo.type="Terrestrial", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
mtext('', side=2, at=-3.0e+05, cex=0.6, line=-0.3, font=2)
#mtext("b", side=3,font = 2, at=-1.5e+07, cex=0.7, line=-0.8)
title(main="Amphibians", cex.main=0.8, line=-0.3)

plot_subrate_grids(mean.rate.grids=terr.reptiles%>%mutate(dS=mid.ds), rate.type = "dS", 
                   geo.type="Terrestrial", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
mtext('', side=2, at=-3.0e+05, cex=0.6, line=-0.3, font=2)
#mtext("c", side=3,font = 2, at=-1.5e+07, cex=0.7, line=-0.8)
title(main="Reptiles", cex.main=0.8, line=-0.3)

plot_subrate_grids(mean.rate.grids=terr.mammals%>%mutate(dS=mid.ds), rate.type = "dS", 
                   geo.type="Terrestrial", colramp=c("#ffffff", "#f4cdb9","#f47155","#d22428", "#660b11"))
mtext('', side=2, at=-3.0e+05, cex=0.6, line=-0.3, font=2)
#mtext("d", side=3,font = 2, at=-1.5e+07, cex=0.7, line=-0.8)
title(main="Mammals", cex.main=0.8, line=-0.3)

plot_subrate_grids(mean.rate.grids=terr.birds%>%mutate(dS=mid.ds), rate.type = "dS", 
                   geo.type="Terrestrial", colramp=c("#ffffff", "#f4cdb9","#f47155","#d22428", "#660b11"))
mtext('', side=2, at=-3.0e+05, cex=0.6, line=-0.3, font=2)
#mtext("e", side=3,font = 2, at=-1.5e+07, cex=0.7, line=-0.8)
title(main="Birds", cex.main=0.8, line=-0.3)

plot_subrate_grids(mean.rate.grids=fishes%>%mutate(dN=mid.dn), rate.type = "dN", 
                   geo.type="Global", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
mtext('Median dN', side=2, at=-3e+05, cex=0.6, line=-0.3, font=2)
title(main="Fishes", cex.main=0.8, line=-0.3)

plot_subrate_grids(mean.rate.grids=amphibians%>%mutate(dN=mid.dn), rate.type = "dN", 
                   geo.type="Terrestrial", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
title(main="Amphibians", cex.main=0.8, line=-0.3)

plot_subrate_grids(mean.rate.grids=terr.reptiles%>%mutate(dN=mid.dn), rate.type = "dN", 
                   geo.type="Terrestrial", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
title(main="Reptiles", cex.main=0.8, line=-0.3)

plot_subrate_grids(mean.rate.grids=terr.mammals%>%mutate(dN=mid.dn), rate.type = "dN", 
                   geo.type="Terrestrial", colramp=c("#ffffff", "#f4cdb9","#f47155","#d22428", "#660b11"))
title(main="Mammals", cex.main=0.8, line=-0.3)

plot_subrate_grids(mean.rate.grids=terr.birds%>%mutate(dN=mid.dn), rate.type = "dN", 
                   geo.type="Terrestrial", colramp=c("#ffffff", "#f4cdb9","#f47155","#d22428", "#660b11"))
title(main="Birds", cex.main=0.8, line=-0.3)
dev.off()


##############
#Extended Fig.7b Latitudinal gradients in molecular rates at assemblage level for each class.
f1 <- mean.subrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Fishes")%>%
  ggplot(aes(x=abs(Lat), y=mid.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y=expression("Median dS (" * 10^-8 * ")"), title="Fishes")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.5,1.1), breaks = seq(0.5,1.1,0.1))+
  annotate("text", x = 60, y = 1.1, label = expression(italic(P)[SAR] == 0.009), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f2 <- mean.subrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Amphibians", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=mid.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y="", title="Amphibians")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0,1.8), breaks = seq(0,1.8,0.3))+
  annotate("text", x = 40, y = 1.8, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f3 <- mean.subrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Reptiles", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=mid.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y="", title="Reptiles")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.3,1.5), breaks = seq(0.3,1.5,0.3))+
  annotate("text", x = 40, y = 1.5, label = expression(italic(P)[SAR] == 0.0013), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f4 <- mean.subrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Mammals", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=mid.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="Latitude", y="", title="Mammals")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.9,2.3), breaks = seq(0.9,2.3,0.2))+
  annotate("text", x = 60, y = 2.3, label = expression(italic(P)[SAR] == 0.48), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f5 <- mean.subrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Birds", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=mid.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="Latitude", y="", title="Birds")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.8,1.5), breaks = seq(0.8,1.5,0.1))+
  annotate("text", x = 60, y = 1.5, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f6 <- mean.subrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Fishes")%>%
  ggplot(aes(x=abs(Lat), y=mid.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y=expression("Median dN (" * 10^-10 * ")"), title="Fishes")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(1.4,2.9), breaks = seq(1.4,2.9,0.3))+
  annotate("text", x = 60, y = 2.9, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f7 <- mean.subrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Amphibians", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=mid.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y="", title="Amphibians")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.6,3.5), breaks = seq(0.6,3.4,0.4))+
  annotate("text", x = 40, y = 3.4, label = expression(italic(P)[SAR] == 0.013), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f8 <- mean.subrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Reptiles", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=mid.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y="", title="Reptiles")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(2,4), breaks = seq(2,4,0.5))+
  annotate("text", x = 40, y = 4, label = expression(italic(P)[SAR] == 0.022), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f9 <- mean.subrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Mammals", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=mid.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="Latitude", y="", title="Mammals")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(4,5.8), breaks = seq(4,5.8,0.3))+
  annotate("text", x = 60, y = 5.8, label = expression(italic(P)[SAR] == 0.002), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f10 <- mean.subrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Birds", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=mid.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="Latitude", y="", title="Birds")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(2.6,3.4), breaks = seq(2.6,3.4,0.2))+
  annotate("text", x = 60, y = 3.4, label = expression(italic(P)[SAR] == 0.29), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

cowplot::plot_grid(f1,f2, f3, f4, f5, 
                   f6,f7, f8, f9, f10,
                   nrow = 2, align="hv")

ggsave("./Outputs/Supplementary/Extended Fig7b.pdf", width=8.27, height=3.2)

