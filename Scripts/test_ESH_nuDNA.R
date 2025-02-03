################################
# R script to test evolutionary speed hypothesis (ESH) using molecular rates of nuclear genes
# Author: Tianlong Cai
# Email: caitianlong@westlake.edu.cn


##################################################################
#Part I: Analyzing the variation in molecular rates across different species and taxonomic groups.
###################################################################

#1. Analyzing how molecular rates varied at latitude and ambient temperature across species
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
# Input data for molecular rates
workdir <- "/Users/Tianlong/VertMolRate"
molrate <- read.csv(paste0(workdir, "/DataFiles/MolEvolRate/molrate_nudna.csv"))

# Input phylogeny
phylo <- read.tree(paste0(workdir, "/DataFiles/trees/phy_all_sampled_nudna.tre"))

# Check whether the phylogeny is ultrametric
is.ultrametric(phylo)

# If the tree is not ultrametric, force it to be ultrametric
phylo <- force.ultrametric(phylo, method = "nnls")  # Use nnls or extend

# 1.1 Phylogenetic Generalized Linear Mixed Models (PGLMMs) 
# PGLMMs account for phylogenetic relatedness as a random effect

# Define groups
groups <- c("Ectotherms", "Endotherms", "Fishes", "Amphibians", "Reptiles", "Mammals", "Birds")

# Initialize lists to store results for PGLMMs
pglmm_molrate_lat <- vector("list", length = length(groups) * 3)
names(pglmm_molrate_lat) <- paste(rep(groups, each = 3), c("dN", "dS", "dNdS"), sep = ".")

pglmm_molrate_temp <- vector("list", length = length(groups) * 3)
names(pglmm_molrate_temp) <- paste(rep(groups, each = 3), c("dN", "dS", "dNdS"), sep = ".")

# Prior specification for the model
prior <- list(R = list(V = diag(1), nu = 0.002), 
              G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1) * 1000)))

# Loop through each group to fit PGLMMs
for (group in groups) {
  
  # Select data for the current group
  if (group %in% unique(molrate$ThermoMode)) {
    subdata <- molrate %>% filter(ThermoMode == group)
  }
  if (group %in% unique(molrate$Group)) {
    subdata <- molrate %>% filter(Group == group)
  }
  
  # Fit models for dS, dN, and dNdS based on latitude and temperature
  fit_lat <- MCMCglmm(scale(log(dS)) ~ scale(abs(Lat)), prior = prior,
                      random = ~Species, data = subdata, verbose = TRUE, 
                      ginverse = list(Species = inverseA(phylo)$Ainv), 
                      nitt = 60000, burnin = 10000, thin = 25, family = c("gaussian"))
  
  fit_temp <- MCMCglmm(scale(log(dS)) ~ scale(AnnualTemp), prior = prior,
                       random = ~Species, data = subdata, verbose = TRUE, 
                       ginverse = list(Species = inverseA(phylo)$Ainv), 
                       nitt = 60000, burnin = 10000, thin = 25, family = c("gaussian"))
  
  # dN model
  fit_dn_lat <- MCMCglmm(scale(log(dN)) ~ scale(abs(Lat)), prior = prior,
                         random = ~Species, data = subdata, verbose = TRUE, 
                         ginverse = list(Species = inverseA(phylo)$Ainv),
                         nitt = 60000, burnin = 10000, thin = 25, family = c("gaussian"))
  
  fit_dn_temp <- MCMCglmm(scale(log(dN)) ~ scale(AnnualTemp), prior = prior,
                          random = ~Species, data = subdata, verbose = TRUE, 
                          ginverse = list(Species = inverseA(phylo)$Ainv), 
                          nitt = 60000, burnin = 10000, thin = 25, family = c("gaussian"))
  
  # dNdS model
  fit_dnds_lat <- MCMCglmm(scale(log(dNdS)) ~ scale(abs(Lat)), prior = prior,
                           random = ~Species, data = subdata, verbose = TRUE, 
                           ginverse = list(Species = inverseA(phylo)$Ainv),
                           nitt = 60000, burnin = 10000, thin = 25, family = c("gaussian"))
  
  fit_dnds_temp <- MCMCglmm(scale(log(dNdS)) ~ scale(AnnualTemp), prior = prior,
                            random = ~Species, data = subdata, verbose = TRUE, 
                            ginverse = list(Species = inverseA(phylo)$Ainv),
                            nitt = 60000, burnin = 10000, thin = 25, family = c("gaussian"))
  
  # Store model results in respective lists
  pglmm_molrate_lat[[paste(group, "dS", sep = ".")]] <- fit_lat
  pglmm_molrate_lat[[paste(group, "dN", sep = ".")]] <- fit_dn_lat
  pglmm_molrate_lat[[paste(group, "dNdS", sep = ".")]] <- fit_dnds_lat
  
  pglmm_molrate_temp[[paste(group, "dS", sep = ".")]] <- fit_temp
  pglmm_molrate_temp[[paste(group, "dN", sep = ".")]] <- fit_dn_temp
  pglmm_molrate_temp[[paste(group, "dNdS", sep = ".")]] <- fit_dnds_temp
  
  print(group)  # Print group name to track progress
}

# Optionally, save the results (uncomment to save)
# save(pglmm_molrate_lat, file = "./Outputs/Data/pglmm_molrate_latitudinal_pattern_nudna.rdata")
# save(pglmm_molrate_temp, file = "./Outputs/Data/pglmm_molrate_temperature_pattern_nudna.rdata")

# Load saved results (uncomment to load)
load(file = "./Outputs/Data/pglmm_molrate_latitudinal_pattern_nudna.rdata")
load(file = "./Outputs/Data/pglmm_molrate_temperature_pattern_nudna.rdata")


#########################################################################################
# Function to summarize PGLMMs of dN, dS, and dNdS for all groups
summary_pglmm_groups <- function(groups, molrate, pglmms) {
  fixed <- lapply(paste(groups, molrate, sep = "."), function(x) summary(pglmms[[x]])$solutions)
  random <- lapply(paste(groups, molrate, sep = "."), function(x) summary(pglmms[[x]])$Gcovariances)
  pglmm.out <- lapply(1:length(groups), 
                      function(x) cbind(data.frame(MolRate = molrate, Group = groups[x], 
                                                   Var = c(row.names(fixed[[x]]), row.names(random[[x]]))), 
                                        data.frame(rbind(fixed[[x]], cbind(random[[x]], pMCMC = NA)))))
  
  summary_out <- do.call(rbind, pglmm.out)
  return(summary_out)
}

# Extract parameter values from PGLMMs for both latitude and temperature
summary_dn_lat <- summary_pglmm_groups(groups = groups, molrate = "dN", pglmms = pglmm_molrate_lat)     
summary_ds_lat <- summary_pglmm_groups(groups = groups, molrate = "dS", pglmms = pglmm_molrate_lat)
summary_dnds_lat <- summary_pglmm_groups(groups = groups, molrate = "dNdS", pglmms = pglmm_molrate_lat)

summary_dn_temp <- summary_pglmm_groups(groups = groups, molrate = "dN", pglmms = pglmm_molrate_temp)  
summary_ds_temp <- summary_pglmm_groups(groups = groups, molrate = "dS", pglmms = pglmm_molrate_temp)
summary_dnds_temp <- summary_pglmm_groups(groups = groups, molrate = "dNdS", pglmms = pglmm_molrate_temp)

# All parameter values for PGLMM results (dN, dS, dNdS for both Latitude and Temperature)
pglmm_out <- rbind(
  bind_rows("Latitude" = summary_ds_lat, "AnnualTemp" = summary_ds_temp, .id = "Predictors"),
  bind_rows("Latitude" = summary_dn_lat, "AnnualTemp" = summary_dn_temp, .id = "Predictors"),
  bind_rows("Latitude" = summary_dnds_lat, "AnnualTemp" = summary_dnds_temp, .id = "Predictors")
)

###########################################################################################
# Function to extract MCMC slopes of PGLMMs
# This function extracts the slopes (posterior means) from the PGLMM results
extract_slope_pglmm <- function(groups, molrate, pglmms) {
  slopes <- lapply(1:length(groups), function(x) {
    data.frame(Group = groups[x], Slope = pglmms[[paste(groups[x], molrate, sep = ".")]]$Sol[, 2])
  })
  slopes <- do.call(rbind, slopes)
  colnames(slopes) <- c("Group", "Slope")
  return(slopes)
}

# Define groups for plotting
plot_groups <- c("Endotherms", "Birds", "Mammals", "Ectotherms", "Reptiles", "Amphibians", "Fishes")

# Extract slopes for different molecular rates (dS, dN, dNdS) for both Latitude and Temperature
slope_lat <- extract_slope_pglmm(groups = plot_groups, molrate = "dS", pglmms = pglmm_molrate_lat)
slope_dn_lat <- extract_slope_pglmm(groups = plot_groups, molrate = "dN", pglmms = pglmm_molrate_lat)
slope_dnds_lat <- extract_slope_pglmm(groups = plot_groups, molrate = "dNdS", pglmms = pglmm_molrate_lat)

slope_temp <- extract_slope_pglmm(groups = plot_groups, molrate = "dS", pglmms = pglmm_molrate_temp)
slope_dn_temp <- extract_slope_pglmm(groups = plot_groups, molrate = "dN", pglmms = pglmm_molrate_temp)
slope_dnds_temp <- extract_slope_pglmm(groups = plot_groups, molrate = "dNdS", pglmms = pglmm_molrate_temp)

###########################################################################################
# Perform Kruskal-Wallis test for each molecular rate and predictor (Latitude, Temperature)
# Adjust p-values using Bonferroni correction
comp_lat <- agricolae::kruskal(slope_lat$Slope, slope_lat$Group, p.adj = "bonferroni")
comp_dn_lat <- agricolae::kruskal(slope_dn_lat$Slope, slope_lat$Group, p.adj = "bonferroni")
comp_dnds_lat <- agricolae::kruskal(slope_dnds_lat$Slope, slope_lat$Group, p.adj = "bonferroni")
comp_temp <- agricolae::kruskal(slope_temp$Slope, slope_temp$Group, p.adj = "bonferroni")
comp_dn_temp <- agricolae::kruskal(slope_dn_temp$Slope, slope_temp$Group, p.adj = "bonferroni")
comp_dnds_temp <- agricolae::kruskal(slope_dnds_temp$Slope, slope_temp$Group, p.adj = "bonferroni")

# Create significance labels for the plot based on the Kruskal-Wallis results
sig.label <- data.frame(
  Group = plot_groups, x = 1:7, 
  ds_lat = comp_lat$groups[plot_groups, "groups"],
  dn_lat = comp_dn_lat$groups[plot_groups, "groups"],
  dnds_lat = comp_dnds_lat$groups[plot_groups, "groups"],
  ds_temp = comp_temp$groups[plot_groups, "groups"],
  dn_temp = comp_dn_temp$groups[plot_groups, "groups"],
  dnds_temp = comp_dnds_temp$groups[plot_groups, "groups"],
  ThermoMode = c(rep("Ectotherms", 4), rep("Endotherms", 3))
)

#Plot effect size
f3a <-pglmm_out%>%
  filter(Group %in% plot_groups)%>%
  mutate(ThermoMode = ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"))%>%
  filter(Var=="scale(AnnualTemp)")%>%
  mutate(sig=ifelse(pMCMC<0.05, "1", "0"),
         Group=factor(Group, levels=c("Endotherms","Birds","Mammals","Ectotherms","Reptiles", "Amphibians","Fishes")))%>%
  filter(MolRate=="dS")%>%
  ggplot(aes(x=Group, y=post.mean, ymin=l.95..CI, ymax=u.95..CI, color=ThermoMode, shape=sig))+
  geom_pointrange(size=0.3)+
  scale_shape_manual(labels=c("0", "1"), values = c(1,19))+
  scale_color_manual(values = c("#2a6aaf", "#d3292f"))+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  labs(y="", x="", title = "dS", tag="A")+
  theme_classic()+
  coord_flip()+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-0.2, 0.4),breaks = seq(-0.2,0.4,0.2), labels=c(-0.2,0,0.2,0.4))+
  theme(axis.text = element_text(size=7, color = "black"),
        axis.title = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        legend.position = 'none',
        plot.tag = element_text(size=8, face = "bold"),
        plot.tag.position = c(0.05, 1),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))+
  geom_text(data=sig.label, aes(x=x, y=0.4, label=ds_temp), size=2.3, inherit.aes = FALSE)

f3b <-pglmm_out%>%
  filter(Group %in% plot_groups)%>%
  mutate(ThermoMode = ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"))%>%
  filter(Var=="scale(AnnualTemp)")%>%
  mutate(sig=ifelse(pMCMC<0.05, "1", "0"),
         Group=factor(Group, levels=c("Endotherms","Birds","Mammals","Ectotherms","Reptiles", "Amphibians","Fishes")))%>%
  filter(MolRate=="dN")%>%
  ggplot(aes(x=Group, y=post.mean, ymin=l.95..CI, ymax=u.95..CI, color=ThermoMode, shape=sig))+
  geom_pointrange(size=0.3)+
  scale_shape_manual(labels=c("0", "1"), values = c(1,19))+
  scale_color_manual(values = c("#2a6aaf", "#d3292f"))+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  labs(y="", x="", title = "dN", tag="B")+
  theme_classic()+
  coord_flip()+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-0.4, 0.4),breaks = seq(-0.4,0.4,0.2), labels=seq(-0.4,0.4,0.2))+
  theme(axis.text = element_text(size=7,color = "black"),
        axis.title = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        legend.position = 'none',
        plot.tag = element_text(size=8, face = "bold"),
        plot.tag.position = c(0.05, 1),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))+
  geom_text(data=sig.label, aes(x=x, y=0.4, label=dn_temp), size=2.3, inherit.aes = FALSE)



f3c <-pglmm_out%>%
  filter(Group %in% plot_groups)%>%
  mutate(ThermoMode = ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"))%>%
  filter(Var=="scale(AnnualTemp)")%>%
  mutate(sig=ifelse(pMCMC<0.05, "1", "0"),
         Group=factor(Group, levels=c("Endotherms","Birds","Mammals","Ectotherms","Reptiles", "Amphibians","Fishes")))%>%
  filter(MolRate=="dNdS")%>%
  ggplot(aes(x=Group, y=post.mean, ymin=l.95..CI, ymax=u.95..CI, color=ThermoMode, shape=sig))+
  geom_pointrange(size=0.3)+
  scale_shape_manual(labels=c("0", "1"), values = c(1,19))+
  scale_color_manual(values = c("#2a6aaf", "#d3292f"))+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  labs(y="Effects of temperature", x="", title = "dN/dS", tag="C")+
  theme_classic()+
  coord_flip()+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-0.6, 0.3),breaks = seq(-0.6,0.3,0.3), labels=seq(-0.6,0.3,0.3))+
  theme(axis.text = element_text(size=7, color = "black"),
        axis.title = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        legend.position = 'none',
        plot.tag = element_text(size=8, face = "bold"),
        plot.tag.position = c(0.05, 1),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))+
  geom_text(data=sig.label, aes(x=x, y=0.3, label=dnds_temp), size=2.3, inherit.aes = FALSE)

# Combine the plots in a grid
cowplot::plot_grid(f3a,f3b, f3c, nrow = 3, align="hv")

# Save the figure as a PDF
ggsave(filename="./Outputs/MainFigures/Fig3ABC.pdf", height=4.3, width=1.8)


####################################################
#Extended Fig.S12
#Plot effect size
fs12a <-pglmm_out%>%
  filter(Group %in% plot_groups)%>%
  mutate(ThermoMode = ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"))%>%
  filter(Var=="scale(abs(Lat))")%>%
  mutate(sig=ifelse(pMCMC<0.05, "1", "0"),
         Group=factor(Group, levels=c("Endotherms","Birds","Mammals","Ectotherms","Reptiles", "Amphibians","Fishes")))%>%
  filter(MolRate=="dS")%>%
  ggplot(aes(x=Group, y=post.mean, ymin=l.95..CI, ymax=u.95..CI, color=ThermoMode, shape=sig))+
  geom_pointrange(size=0.3)+
  scale_shape_manual(labels=c("0", "1"), values = c(1,19))+
  scale_color_manual(values = c("#2a6aaf", "#d3292f"))+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  labs(y="Effects of latitude", x="", title = "dS", tag="A")+
  theme_classic()+
  coord_flip()+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  #scale_y_continuous(limits = c(-0.3, 0.3),breaks = seq(-0.3, 0.3,0.15), labels=c(-0.3, -0.15, 0,0.15,0.3))+
  theme(axis.text = element_text(size=7),
        axis.title = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = 'none',
        plot.tag = element_text(size=8, face = "bold"),
        plot.tag.position = c(0.05, 1),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))+
  geom_text(data=sig.label, aes(x=x, y=0.3, label=ds_lat), size=2.5, inherit.aes = FALSE)

fs12b <-pglmm_out%>%
  filter(Group %in% plot_groups)%>%
  mutate(ThermoMode = ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"))%>%
  filter(Var=="scale(abs(Lat))")%>%
  mutate(sig=ifelse(pMCMC<0.05, "1", "0"),
         Group=factor(Group, levels=c("Endotherms","Birds","Mammals","Ectotherms","Reptiles", "Amphibians","Fishes")))%>%
  filter(MolRate=="dN")%>%
  ggplot(aes(x=Group, y=post.mean, ymin=l.95..CI, ymax=u.95..CI, color=ThermoMode, shape=sig))+
  geom_pointrange(size=0.3)+
  scale_shape_manual(labels=c("0", "1"), values = c(1,19))+
  scale_color_manual(values = c("#2a6aaf", "#d3292f"))+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  labs(y="Effects of latitude", x="", title = "dN", tag="B")+
  theme_classic()+
  coord_flip()+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-0.3, 0.3),breaks = seq(-0.3, 0.3,0.15), labels=c(-0.3, -0.15, 0,0.15,0.3))+
  theme(axis.text = element_text(size=7),
        axis.title = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = 'none',
        plot.tag = element_text(size=8, face = "bold"),
        plot.tag.position = c(0.05, 1),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))+
  geom_text(data=sig.label, aes(x=x, y=0.3, label=dn_lat), size=2.5, inherit.aes = FALSE)



fs12c <-pglmm_out%>%
  filter(Group %in% plot_groups)%>%
  mutate(ThermoMode = ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"))%>%
  filter(Var=="scale(abs(Lat))")%>%
  mutate(sig=ifelse(pMCMC<0.05, "1", "0"),
         Group=factor(Group, levels=c("Endotherms","Birds","Mammals","Ectotherms","Reptiles", "Amphibians","Fishes")))%>%
  filter(MolRate=="dNdS")%>%
  ggplot(aes(x=Group, y=post.mean, ymin=l.95..CI, ymax=u.95..CI, color=ThermoMode, shape=sig))+
  geom_pointrange(size=0.3)+
  scale_shape_manual(labels=c("0", "1"), values = c(1,19))+
  scale_color_manual(values = c("#2a6aaf", "#d3292f"))+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  labs(y="Effects of latitude", x="", title = "dN/dS", tag="C")+
  theme_classic()+
  coord_flip()+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-0.2, 0.4),breaks = seq(-0.2, 0.4,0.2), labels=seq(-0.2, 0.4,0.2))+
  theme(axis.text = element_text(size=7),
        axis.title = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = 'none',
        plot.tag = element_text(size=8, face = "bold"),
        plot.tag.position = c(0.05, 1),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))+
  geom_text(data=sig.label, aes(x=x, y=0.4, label=dnds_lat), size=2.5, inherit.aes = FALSE)


cowplot::plot_grid(fs12a,fs12b, fs12c, nrow = 1, align="hv")
ggsave(filename="./Outputs/Supplementary/Extended Fig12.pdf", height=2.3, width=8.27)

################################################################################
# 1.2 PGLMMs account for phylogenetic relatedness under the best traits evolution models.
##########
#Fit the best models
groups <- c("Ectotherms", "Endotherms", "Fishes", "Amphibians", "Reptiles", "Mammals","Birds")

# Creating a list to store PGLMM objects
trait_models <- vector("list", length = length(groups) * 3)
names(trait_models) <- paste(rep(groups, each = 3), c("dN", "dS", "dNdS"), sep = ".")


for (group in groups) {
  # Choose data based on group
  if(group %in% unique(molrate$ThermoMode)){
    subdata <- molrate %>% filter(ThermoMode == group)%>%mutate(dS=as.numeric(scale(log(dS))), dN=as.numeric(scale(log(dN))))
  }
  if(group %in% unique(molrate$Group)){
    subdata <- molrate %>% filter(Group == group)%>%mutate(dS=as.numeric(scale(log(dS))), dN=as.numeric(scale(log(dN))))
  }
  
  #phylogeny
  phy <- drop.tip(phylo, setdiff(phylo$tip.label, subdata$Species))
  
  #Model selection
  trait_models[[paste(group, "dS", sep = ".")]] <- fitTraitEvolModel(trait="dS", data=subdata, phy = phy)
  trait_models[[paste(group, "dN", sep = ".")]] <- fitTraitEvolModel(trait="dN", data=subdata, phy = phy)
  trait_models[[paste(group, "dNdS", sep = ".")]] <- fitTraitEvolModel(trait="dNdS", data=subdata, phy = phy)
}

#save(trait_models, file="./Outputs/Data/best_trait_evol_models_nudna.rdata")
load(file="./Outputs/Data/best_trait_evol_models_nudna.rdata")

trait_models %>%bind_rows()%>%
  mutate(Group=rep(groups, each=3), molrate=rep(c("dN", "dS", "dNdS"), length(groups)))

##############################
#Fit the PGLMMs using the best trait evolution model
pglmm_molrate_lat_best_model <- vector("list", length = length(groups) * 3)
names(pglmm_molrate_lat_best_model) <- paste(rep(groups, each = 3), c("dN", "dS", "dNdS"), sep = ".")

pglmm_molrate_temp_best_model <- vector("list", length = length(groups) * 3)
names(pglmm_molrate_temp_best_model) <- paste(rep(groups, each = 3), c("dN", "dS", "dNdS"), sep = ".")

# Prior specification for the model
prior <- list(R=list(V=diag(1),nu=0.002), 
              G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*1000)))


# Loop over each group
for (group in groups) {
  # Choose data based on group
  if(group %in% unique(molrate$ThermoMode)){
    subdata <- molrate %>% filter(ThermoMode == group)
  }
  
  if(group %in% unique(molrate$Group)){
    subdata <- molrate %>% filter(Group == group)
  }
  
  phy <- drop.tip(phylo, setdiff(phylo$tip.label, subdata$Species))
  
  #The best traits evolution model of dS
  model.ds <- trait_models[[paste(group, "dS", sep = ".")]]
  
  #The best traits evolution model of dN
  model.dn <- trait_models[[paste(group, "dN", sep = ".")]]
  
  #The best traits evolution model of dNdS
  model.dnds <- trait_models[[paste(group, "dNdS", sep = ".")]]
  
  #rescale tree
  if(model.ds$best_model=="OU"){phy <- phytools::rescale(phy, model=model.ds$best_model, alpha=model.ds$opt, sigsq=model.ds$sigsq)}
  if(model.ds$best_model=="BM"){phy <- phytools::rescale(phy, model=model.ds$best_model, sigsq=model.ds$sigsq)}
  if(model.ds$best_model=="EB"){phy <- phytools::rescale(phy, model=model.ds$best_model, a=model.ds$opt, sigsq=model.ds$sigsq)}
  
  #rescale tree
  if(model.dn$best_model=="OU"){phy_dn <- phytools::rescale(phy, model=model.dn$best_model, alpha=model.dn$opt, sigsq=model.dn$sigsq)}
  if(model.dn$best_model=="BM"){phy_dn <- phytools::rescale(phy, model=model.dn$best_model, sigsq=model.dn$sigsq)}
  if(model.dn$best_model=="EB"){phy_dn <- phytools::rescale(phy, model=model.dn$best_model, a=model.dn$opt, sigsq=model.dn$sigsq)}
  
  #rescale tree
  if(model.dnds$best_model=="OU"){phy_dnds <- phytools::rescale(phy, model=model.dnds$best_model, alpha=model.dnds$opt, sigsq=model.dnds$sigsq)}
  if(model.dnds$best_model=="BM"){phy_dnds <- phytools::rescale(phy, model=model.dnds$best_model, sigsq=model.dnds$sigsq)}
  if(model.dnds$best_model=="EB"){phy_dnds <- phytools::rescale(phy, model=model.dnds$best_model, a=model.dnds$opt, sigsq=model.dnds$sigsq)}
  
  
  #modify edge.length in case of small values
  for (i in 1:length(phy$edge.length)) {
    if (phy$edge.length[i] < 1e-16) {
      phy$edge.length[i] = 1e-16
    }
  }
  for (i in 1:length(phy_dn$edge.length)) {
    if (phy_dn$edge.length[i] < 1e-16) {
      phy_dn$edge.length[i] = 1e-16
    }
  }
  for (i in 1:length(phy_dnds$edge.length)) {
    if (phy_dnds$edge.length[i] < 1e-16) {
      phy_dnds$edge.length[i] = 1e-16
    }
  }
  
  
  # Scale and fit dS model
  fit_lat <- MCMCglmm(scale(log(dS)) ~ scale(abs(Lat)), prior=prior,
                      random = ~Species,
                      data = subdata, verbose = TRUE, 
                      ginverse = list(Species = inverseA(phy)$Ainv), 
                      nitt = 60000, burnin = 10000, thin = 25, family = c("gaussian"))
  fit_temp <- MCMCglmm(scale(log(dS)) ~ scale(AnnualTemp), prior=prior,
                       random = ~Species,
                       data = subdata, verbose = TRUE, 
                       ginverse = list(Species = inverseA(phy)$Ainv), 
                       nitt = 60000, burnin = 10000, thin = 25, family = c("gaussian"))
  
  # Scale and fit dN model
  fit_dn_lat <- MCMCglmm(scale(log(dN)) ~ scale(abs(Lat)), prior=prior,
                         random = ~Species,
                         data = subdata, verbose = TRUE, 
                         ginverse = list(Species = inverseA(phy_dn)$Ainv),
                         nitt = 60000, burnin = 10000, thin = 25, family = c("gaussian"))
  
  fit_dn_temp <- MCMCglmm(scale(log(dN)) ~ scale(AnnualTemp), prior=prior,
                          random = ~Species,
                          data = subdata, verbose = TRUE, 
                          ginverse = list(Species = inverseA(phy_dn)$Ainv),
                          nitt = 60000, burnin = 10000, thin = 25, family = c("gaussian"))
  
  # Scale and fit dN/dS model
  fit_dnds_lat <- MCMCglmm(scale(log(dNdS)) ~ scale(abs(Lat)), prior=prior,
                           random = ~Species,
                           data = subdata, verbose = TRUE, 
                           ginverse = list(Species = inverseA(phy_dn)$Ainv),
                           nitt = 60000, burnin = 10000, thin = 25, family = c("gaussian"))
  fit_dnds_temp <- MCMCglmm(scale(log(dNdS)) ~ scale(AnnualTemp), prior=prior,
                            random = ~Species,
                            data = subdata, verbose = TRUE, 
                            ginverse = list(Species = inverseA(phy_dn)$Ainv),
                            nitt = 60000, burnin = 10000, thin = 25, family = c("gaussian"))
  
  # Store models in list
  pglmm_molrate_lat_best_model[[paste(group, "dS", sep = ".")]] <- fit_lat
  pglmm_molrate_lat_best_model[[paste(group, "dN", sep = ".")]] <- fit_dn_lat
  pglmm_molrate_lat_best_model[[paste(group, "dNdS", sep = ".")]] <- fit_dnds_lat
  
  pglmm_molrate_temp_best_model[[paste(group, "dS", sep = ".")]] <- fit_temp
  pglmm_molrate_temp_best_model[[paste(group, "dN", sep = ".")]] <- fit_dn_temp
  pglmm_molrate_temp_best_model[[paste(group, "dNdS", sep = ".")]] <- fit_dnds_temp
  
  print(group)
}

#save(pglmm_molrate_lat_best_model, file = "./Outputs/Data/pglmm_molrate_latitudinal_pattern_best_model_nudna.rdata")
#save(pglmm_molrate_temp_best_model, file = "./Outputs/Data/pglmm_molrate_temperature_pattern_best_model_nudna.rdata")

load(file = "./Outputs/Data/pglmm_molrate_latitudinal_pattern_best_model_nudna.rdata")
load(file = "./Outputs/Data/pglmm_molrate_temperature_pattern_best_model_nudna.rdata")

# Extract parameter values from PGLMMs
summary_dn_lat_best_model <- summary_pglmm_groups(groups=groups, molrate="dN", pglmms=pglmm_molrate_lat_best_model)
summary_lat_best_model <- summary_pglmm_groups(groups=groups, molrate="dS", pglmms=pglmm_molrate_lat_best_model)
summary_dnds_lat_best_model <- summary_pglmm_groups(groups=groups, molrate="dNdS", pglmms=pglmm_molrate_lat_best_model)

summary_dn_temp_best_model <- summary_pglmm_groups(groups=groups, molrate="dN", pglmms=pglmm_molrate_temp_best_model)
summary_temp_best_model <- summary_pglmm_groups(groups=groups, molrate="dS", pglmms=pglmm_molrate_temp_best_model)
summary_dnds_temp_best_model <- summary_pglmm_groups(groups=groups, molrate="dNdS", pglmms=pglmm_molrate_temp_best_model)

pglmm_out_best_model <- rbind(bind_rows("Latitude" = summary_lat_best_model, "AnnualTemp" = summary_temp_best_model, .id = "Predictors"),
                              bind_rows("Latitude" = summary_dn_lat_best_model, "AnnualTemp" = summary_dn_temp_best_model, .id = "Predictors"),
                              bind_rows("Latitude" = summary_dnds_lat_best_model, "AnnualTemp" = summary_dnds_temp_best_model, .id = "Predictors"))
#write.csv(pglmm_out_best_model, "pglmm_out.csv")
########################################################
#1.3. PGLMMs account for phylogenetic relatedness, habitat and lineages as random effects.
groups <- c(unique(molrate$Group), unique(molrate$ThermoMode))

pglmm_molrate_lat_3randoms <- vector("list", length = length(groups) * 3)
names(pglmm_molrate_lat_3randoms) <- paste(rep(groups, each = 3), c("dN", "dS", "dNdS"), sep = ".")

pglmm_molrate_temp_3randoms <- vector("list", length = length(groups) * 3)
names(pglmm_molrate_temp_3randoms) <- paste(rep(groups, each = 3), c("dN", "dS", "dNdS"), sep = ".")

# Prior specification for the model
prior <- list(R=list(V=diag(1),nu=0.002), 
              G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*1000),
                     G2=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*1000),
                     G3=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*1000)))

# Loop over each group
for (group in groups) {
  # Choose data based on group
  if(group %in% unique(molrate$ThermoMode)){
    subdata <- molrate %>% filter(ThermoMode == group)
  }
  
  if(group %in% unique(molrate$Group)){
    subdata <- molrate %>% filter(Group == group)
  }
  
  # Scale and fit dS model
  fit_lat <- MCMCglmm(scale(log(dS)) ~ scale(abs(Lat)),prior=prior,
                      random = ~Species+Habitat+Clade,
                      data = subdata, verbose = TRUE, 
                      ginverse = list(Species = inverseA(phylo)$Ainv), 
                      nitt = 60000, burnin = 10000, thin = 25, family = c("gaussian"))
  fit_temp <- MCMCglmm(scale(log(dS)) ~ scale(AnnualTemp),prior=prior,
                       random = ~Species+Habitat+Clade,
                       data = subdata, verbose = TRUE, 
                       ginverse = list(Species = inverseA(phylo)$Ainv), 
                       nitt = 60000, burnin = 10000, thin = 25, family = c("gaussian"))
  
  # Scale and fit dN model
  fit_dn_lat <- MCMCglmm(scale(log(dN)) ~ scale(abs(Lat)), prior=prior,
                         random = ~Species+Habitat+Clade,
                         data = subdata, verbose = TRUE, 
                         ginverse = list(Species = inverseA(phylo)$Ainv),
                         nitt = 60000, burnin = 10000, thin = 25, family = c("gaussian"))
  fit_dn_temp <- MCMCglmm(scale(log(dN)) ~ scale(AnnualTemp), prior=prior,
                          random = ~Species+Habitat+Clade,
                          data = subdata, verbose = TRUE, 
                          ginverse = list(Species = inverseA(phylo)$Ainv),
                          nitt = 60000, burnin = 10000, thin = 25, family = c("gaussian"))
  
  # Scale and fit dN/dS model
  fit_dnds_lat <- MCMCglmm(scale(log(dNdS)) ~ scale(abs(Lat)), prior=prior,
                           random = ~Species+Habitat+Clade,
                           data = subdata, verbose = TRUE, 
                           ginverse = list(Species = inverseA(phylo)$Ainv),
                           nitt = 60000, burnin = 10000, thin = 25, family = c("gaussian"))
  fit_dnds_temp <- MCMCglmm(scale(log(dNdS)) ~ scale(AnnualTemp), prior=prior,
                            random = ~Species+Habitat+Clade,
                            data = subdata, verbose = TRUE, 
                            ginverse = list(Species = inverseA(phylo)$Ainv),
                            nitt = 60000, burnin = 10000, thin = 25, family = c("gaussian"))
  
  # Store models in list
  pglmm_molrate_lat_3randoms[[paste(group, "dS", sep = ".")]] <- fit_lat
  pglmm_molrate_lat_3randoms[[paste(group, "dN", sep = ".")]] <- fit_dn_lat
  pglmm_molrate_lat_3randoms[[paste(group, "dNdS", sep = ".")]] <- fit_dnds_lat
  
  pglmm_molrate_temp_3randoms[[paste(group, "dS", sep = ".")]] <- fit_temp
  pglmm_molrate_temp_3randoms[[paste(group, "dN", sep = ".")]] <- fit_dn_temp
  pglmm_molrate_temp_3randoms[[paste(group, "dNdS", sep = ".")]] <- fit_dnds_temp
  
  
  print(group)
}


#save(pglmm_molrate_lat_3randoms, file = "./Outputs/Data/pglmm_molrate_latitudinal_pattern_3randoms_nudna.rdata")
#save(pglmm_molrate_temp_3randoms, file = "./Outputs/Data/pglmm_molrate_temperature_pattern_3randoms_nudna.rdata")

load(file = "./Outputs/Data/pglmm_molrate_latitudinal_pattern_3randoms_nudna.rdata")
load(file = "./Outputs/Data/pglmm_molrate_temperature_pattern_3randoms_nudna.rdata")

# Extract parameter values from PGLMMs
summary_dn_lat_3randoms <- summary_pglmm_groups(groups=groups, molrate="dN", pglmms=pglmm_molrate_lat_3randoms)
summary_lat_3randoms <- summary_pglmm_groups(groups=groups, molrate="dS", pglmms=pglmm_molrate_lat_3randoms)
summary_dnds_lat_3randoms <- summary_pglmm_groups(groups=groups, molrate="dNdS", pglmms=pglmm_molrate_lat_3randoms)

summary_dn_temp_3randoms <- summary_pglmm_groups(groups=groups, molrate="dN", pglmms=pglmm_molrate_temp_3randoms)
summary_temp_3randoms <- summary_pglmm_groups(groups=groups, molrate="dS", pglmms=pglmm_molrate_temp_3randoms)
summary_dnds_temp_3randoms <- summary_pglmm_groups(groups=groups, molrate="dNdS", pglmms=pglmm_molrate_temp_3randoms)

pglmm_out_3randoms <- rbind(bind_rows("Latitude" = summary_lat_3randoms, "AnnualTemp" = summary_temp_3randoms, .id = "Predictors"),
                            bind_rows("Latitude" = summary_dn_lat_3randoms, "AnnualTemp" = summary_dn_temp_3randoms, .id = "Predictors"),
                            bind_rows("Latitude" = summary_dnds_lat_3randoms, "AnnualTemp" = summary_dnds_temp_3randoms, .id = "Predictors"))

#write.csv(pglmm_out_3randoms, "pglmm_out_3randoms.csv")

#################################
#2. Analyzing how molecular rates vary with latitude and temperature at assemblage level
rm(list=ls())
gc()
workdir <- "/Users/Tianlong/VertMolRate"
setwd(workdir)
library(raster)
library(tidyverse)
library(RColorBrewer)
library(spdep)
library(sf)
library(spatialreg)

source(paste0(workdir, "/Scripts/source_functions.r"))

#input spatial join data and grids
load("./DataFiles/SpatialJoinFiles/grids.rdata")
load("./DataFiles/SpatialJoinFiles/all_spatial_join_grid_nudna.rdata")
load("./DataFiles/SpatialJoinFiles/all_spatial_join_ecoregion_nudna.rdata")


#input GIS poly and raster layers
marine.raster <- raster("./DataFiles/GISLayers/marine.tif")
terrestrial.raster <- raster("./DataFiles/GISLayers/terrestrial.tif")
outline.poly <- read_sf("./DataFiles/GISLayers/outline_moll.shp")
country.poly <- read_sf("./DataFiles/GISLayers/country.poly.merged.shp")
ecos.poly <- read.csv("./DataFiles/GISLayers/ecos.poly.csv")
terrestrial.ecoregions <- read.csv("./DataFiles/GISLayers/terrestrial.ecoregions.csv")
marine.ecoregions <- read.csv("./DataFiles/GISLayers/marine.ecoregions.csv")


#estimated mean substitution rate in grids
fresh.fishes <-mean_molrate_grid(spatial_join=all_spatial_join_grid%>%filter(Group=="Fishes"&(Habitat!="Marine")),
                                 grids_poly=grids, Species=5, habitat="Terrestrial")
marine.fishes <- mean_molrate_grid(spatial_join=all_spatial_join_grid%>%filter(Group=="Fishes"&(Habitat!="Freshwater")),
                                   grids_poly=grids, Species=5, habitat="Marine")
terr.reptiles <- mean_molrate_grid(spatial_join=all_spatial_join_grid%>%filter(Group=="Reptiles"&Habitat!="Marine"),
                                   grids_poly=grids, Species=4, habitat="Terrestrial")
marine.reptiles <- mean_molrate_grid(spatial_join=all_spatial_join_grid%>%filter(Group=="Reptiles"&Habitat!="Terrestrial"&Habitat!="Freshwater"),
                                     grids_poly=grids, Species=4, habitat="Marine")
amphibians <- mean_molrate_grid(spatial_join=all_spatial_join_grid%>%filter(Group=="Amphibians"),
                                grids_poly=grids, Species=4, habitat="Terrestrial")
terr.birds <- mean_molrate_grid(spatial_join=all_spatial_join_grid%>%filter(Group=="Birds"&Habitat!="Marine"),
                                grids_poly=grids, Species=5, habitat="Terrestrial")
marine.birds <- mean_molrate_grid(spatial_join=all_spatial_join_grid%>%filter(Group=="Birds"&Habitat!="Terrestrial"&Habitat!="Freshwater"),
                                  grids_poly=grids, Species=5, habitat="Marine")
terr.mammals <- mean_molrate_grid(spatial_join=all_spatial_join_grid%>%filter(Group=="Mammals"&Habitat!="Marine"),
                                  grids_poly=grids, Species=5, habitat="Terrestrial")
marine.mammals <- mean_molrate_grid(spatial_join=all_spatial_join_grid%>%filter(Group=="Mammals"&Habitat!="Terrestrial"&Habitat!="Freshwater"),
                                    grids_poly=grids, Species=5, habitat="Marine")
terr.ectotherm <- mean_molrate_grid(spatial_join=all_spatial_join_grid%>%filter((Group %in% c("Amphibians", "Reptiles","Fishes"))&Habitat!="Marine"),
                                    grids_poly=grids, Species=5, habitat="Terrestrial")
marine.ectotherm <- mean_molrate_grid(spatial_join=all_spatial_join_grid%>%filter((Group %in% c("Amphibians", "Reptiles","Fishes"))&Habitat!="Terrestrial"&Habitat!="Freshwater"),
                                      grids_poly=grids, Species=5, habitat="Marine")
terr.endotherm <- mean_molrate_grid(spatial_join=all_spatial_join_grid%>%filter((Group %in% c("Mammals", "Birds"))&Habitat!="Marine"),
                                    grids_poly=grids, Species=5, habitat="Terrestrial")
marine.endotherm <- mean_molrate_grid(spatial_join=all_spatial_join_grid%>%filter((Group %in% c("Mammals", "Birds"))&Habitat!="Terrestrial"&Habitat!="Freshwater"),
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



#estimated mean molrate in ecoregions
fresh.fishes.ecos <- mean_molrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter(Group=="Fishes"&Habitat!="Marine"),
                                            Species=5, habitat="Terrestrial")%>%filter(!is.na(SR))
marine.fishes.ecos <- mean_molrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter(Group=="Fishes"&Habitat!="Freshwater"),
                                             Species=5, habitat="Marine")%>%filter(!is.na(SR))
terr.reptiles.ecos <- mean_molrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter(Group=="Reptiles"&Habitat!="Marine"),
                                             Species=3, habitat="Terrestrial")%>%filter(!is.na(SR))
marine.reptiles.ecos <- mean_molrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter(Group=="Reptiles"&Habitat!="Terrestrial"&Habitat!="Freshwater"),
                                               Species=3, habitat="Marine")%>%filter(!is.na(SR))
amphibians.ecos <- mean_molrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter(Group=="Amphibians"),
                                          Species=3, habitat="Terrestrial")%>%filter(!is.na(SR))

terr.birds.ecos <- mean_molrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter(Group=="Birds"&Habitat!="Marine"),
                                          Species=5, habitat="Terrestrial")%>%filter(!is.na(SR))
marine.birds.ecos <- mean_molrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter(Group=="Birds"&Habitat!="Terrestrial"&Habitat!="Freshwater"),
                                            Species=5, habitat="Marine")%>%filter(!is.na(SR))
terr.mammals.ecos <- mean_molrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter(Group=="Mammals"&Habitat!="Marine"),
                                            Species=5, habitat="Terrestrial")%>%filter(!is.na(SR))

marine.mammals.ecos <- mean_molrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter(Group=="Mammals"&Habitat!="Terrestrial"&Habitat!="Freshwater"),
                                              Species=5, habitat="Marine")%>%filter(!is.na(SR))

terr.ectotherm.ecos <- mean_molrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter((Group %in% c("Amphibians", "Reptiles","Fishes"))&Habitat!="Marine"),
                                              Species=5, habitat="Terrestrial")%>%filter(!is.na(SR))
marine.ectotherm.ecos <- mean_molrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter((Group %in% c("Amphibians", "Reptiles","Fishes"))&Habitat!="Terrestrial"&Habitat!="Freshwater"),
                                                Species=5, habitat="Marine")%>%filter(!is.na(SR))

terr.endotherm.ecos <- mean_molrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter((Group %in% c("Birds", "Mammals"))&Habitat!="Marine"),
                                              Species=5, habitat="Terrestrial")%>%filter(!is.na(SR))

marine.endotherm.ecos <- mean_molrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter((Group %in% c("Birds", "Mammals"))&Habitat!="Terrestrial"&Habitat!="Freshwater"),
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


mean.molrate.ecos <- bind_rows(list(Fishes=fishes.ecos,
                                    Amphibians=amphibians.ecos,
                                    Reptiles=reptiles.ecos%>%filter(Habitat=="Terrestrial"),
                                    Birds=birds.ecos%>%filter(Habitat=="Terrestrial"), 
                                    Mammals=mammals.ecos%>%filter(Habitat=="Terrestrial"), 
                                    Endotherms=endotherm.ecos%>%filter(Habitat=="Terrestrial"), 
                                    Ectotherms=ectotherm.ecos), .id="Group")%>%
  mutate(Group=factor(Group, levels = c("Fishes","Amphibians", "Reptiles","Mammals","Birds", "Ectotherms", "Endotherms")))


#############################################################
#2.1 Spatial simultaneous autoregressive (SAR) models to examine relationship between 
#molecular rate and latitude/temperature at assemblages
############################################################
fit.vars <- c("gm.ds", "gm.dn", "gm.dnds", "mid.ds", "mid.dn", "mid.dnds")

#latitude
sar.summary.out.lat <- NULL

for(id in unique(mean.molrate.ecos$Group)){
  out <- as.list(rep(NA, length(fit.vars)))
  names(out) <- fit.vars
  
  
  mean.molrate <- mean.molrate.ecos%>%filter(Group==id)%>%mutate(abs.lat=abs(Lat))
  
  nc.coords <- cbind(mean.molrate$Lon, mean.molrate$Lat)
  
  nc.5nn <- knearneigh(nc.coords, k=5, longlat = TRUE)
  nc.5nn.nb <- knn2nb(nc.5nn)
  
  #fit SAR model
  for (i in 1:length(fit.vars)) {
    var <- fit.vars[i]
    
    # fit ols
    fit.ols <- lm(scale(mean.molrate[,var,drop=T]) ~ scale(mean.molrate[,"abs.lat",drop=T]))
    
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
  names(var.sp) <- paste0("a", 1:6)
  
  var.sp<-var.sp%>%bind_rows()%>%t()%>%as.data.frame()%>%
    mutate(V1=ifelse(V1=="gm", "Geometric mean", ifelse(V1=="am", "Arithmetic mean", "Median")),
           V2=ifelse(V2=="ds", "dS", ifelse(V2=="dn", "dN", "dNdS")))
  
  # Summary
  x <- numeric(length(fit.vars))
  dff <- data.frame(Group=id, Avg.Value =var.sp[,1], molrate = var.sp[,2], AIC.OLS=x, AIC.SAR=x, dAIC=x, 
                    OLS.Slope=x, OLS.Slope.SE=x, OLS.Pvalue=x, SAR.Slope=x, SAR.Slope.SE=x, SAR.Pvalue=x, 
                    Moran.OLS = x, Moran.OLS.Pvalue = x, Moran.SAR = x, Moran.SAR.Pvalue = x, stringsAsFactors=F)
  
  for (i in 1:length(out)) {
    
    fres <- out[[i]]
    
    dff$AIC.OLS[i] <- fres$aic.ols
    dff$AIC.SAR[i] <- fres$aic.sar
    dff$dAIC[i] <- fres$aic.ols - fres$aic.sar
    dff$OLS.Slope[i] <- summary(fres$fit.ols)$coef[2,1]
    dff$OLS.Slope.SE[i] <- summary(fres$fit.ols)$coef[2,2]
    dff$OLS.Pvalue[i] <- summary(fres$fit.ols)$coef[2,4]
    dff$SAR.Slope[i] <- summary(fres$fit.sar)$Coef[2,1]
    dff$SAR.Slope.SE[i] <- summary(fres$fit.sar)$Coef[2,2]
    dff$SAR.Pvalue[i] <- summary(fres$fit.sar)$Coef[2,4]
    dff$Moran.OLS[i] <- fres$moran.ols$estimate[1]
    dff$Moran.OLS.Pvalue[i] <- fres$moran.ols$p.value
    dff$Moran.SAR[i] <- fres$moran.sar$estimate[1]
    dff$Moran.SAR.Pvalue[i] <- fres$moran.sar$p.value
  }
  sar.summary.out.lat <- rbind(sar.summary.out.lat, dff)
  print(id)
}

#temperature
sar.summary.out.temp <- NULL

for(id in unique(mean.molrate.ecos$Group)){
  out <- as.list(rep(NA, length(fit.vars)))
  names(out) <- fit.vars
  
  
  mean.molrate <- mean.molrate.ecos%>%filter(Group==id)
  
  nc.coords <- cbind(mean.molrate$Lon, mean.molrate$Lat)
  
  nc.5nn <- knearneigh(nc.coords, k=5, longlat = TRUE)
  nc.5nn.nb <- knn2nb(nc.5nn)
  
  #fit SAR model
  for (i in 1:length(fit.vars)) {
    var <- fit.vars[i]
    
    # fit ols
    fit.ols <- lm(scale(mean.molrate[,var,drop=T]) ~ scale(mean.molrate[,"AnnualTemp",drop=T]))
    
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
  names(var.sp) <- paste0("a", 1:6)
  
  var.sp<-var.sp%>%bind_rows()%>%t()%>%as.data.frame()%>%
    mutate(V1=ifelse(V1=="gm", "Geometric mean", ifelse(V1=="am", "Arithmetic mean", "Median")),
           V2=ifelse(V2=="ds", "dS", ifelse(V2=="dn", "dN", "dNdS")))
  
  # Summary
  x <- numeric(length(fit.vars))
  dff <- data.frame(Group=id, Avg.Value =var.sp[,1], molrate = var.sp[,2], AIC.OLS=x, AIC.SAR=x, dAIC=x, 
                    OLS.Slope=x, OLS.Slope.SE=x, OLS.Pvalue=x, SAR.Slope=x, SAR.Slope.SE=x, SAR.Pvalue=x, 
                    Moran.OLS = x, Moran.OLS.Pvalue = x, Moran.SAR = x, Moran.SAR.Pvalue = x, stringsAsFactors=F)
  
  for (i in 1:length(out)) {
    
    fres <- out[[i]]
    
    dff$AIC.OLS[i] <- fres$aic.ols
    dff$AIC.SAR[i] <- fres$aic.sar
    dff$dAIC[i] <- fres$aic.ols - fres$aic.sar
    dff$OLS.Slope[i] <- summary(fres$fit.ols)$coef[2,1]
    dff$OLS.Slope.SE[i] <- summary(fres$fit.ols)$coef[2,2]
    dff$OLS.Pvalue[i] <- summary(fres$fit.ols)$coef[2,4]
    dff$SAR.Slope[i] <- summary(fres$fit.sar)$Coef[2,1]
    dff$SAR.Slope.SE[i] <- summary(fres$fit.sar)$Coef[2,2]
    dff$SAR.Pvalue[i] <- summary(fres$fit.sar)$Coef[2,4]
    dff$Moran.OLS[i] <- fres$moran.ols$estimate[1]
    dff$Moran.OLS.Pvalue[i] <- fres$moran.ols$p.value
    dff$Moran.SAR[i] <- fres$moran.sar$estimate[1]
    dff$Moran.SAR.Pvalue[i] <- fres$moran.sar$p.value
  }
  sar.summary.out.temp <- rbind(sar.summary.out.temp, dff)
  print(id)
}

sar.summary.out <-bind_rows(Temperature=sar.summary.out.temp, Latitude=sar.summary.out.lat, .id="Var")


#save results
#save(sar.summary.out, file = "./Outputs/Data/SAR_molrate_pattern_assemblages_nudna.rdata")
load(file = "./Outputs/Data/SAR_molrate_pattern_assemblages_nudna.rdata")

#write.csv(sar.summary.out,"sar.summary.out.csv")


########################################
#Fig.3DFH Molecular rate at grids
pdf("./Outputs/MainFigures/Fig3DFH.pdf", height = 4.2, width=3.6)
par(mar=c(0.6,0.6,0.6,0.6), mfrow=c(3,2))
plot_molrate_grids(mean.rate.grids=ectotherm%>%mutate(dS=gm.ds)%>%filter(!is.na(SR)), 
                   rate.type = "dS", geo.type="Global", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
mtext('dS', side=2, at=-3.0e+05, cex=0.6, line=-0.1, font=2)
mtext("D", side=3,font = 2, at=-1.5e+07, cex=0.6, line=-0.8)
title(main="Ectotherms", cex.main=0.8, line=-0.3)

plot_molrate_grids(mean.rate.grids=endotherm%>%mutate(dS=gm.ds)%>%filter(Habitat!="Marine"), 
                   rate.type = "dS", geo.type="Terrestrial", colramp=c("#ffffff", "#f4cdb9","#f47155","#d22428", "#660b11"))
#mtext('Endotherms', side=2, at=-3.0e+05, cex=0.6, line=-0.1, font=2)
title(main="Endotherms", cex.main=0.8, line=-0.3)

plot_molrate_grids(mean.rate.grids=ectotherm%>%mutate(dN=gm.dn), 
                   rate.type = "dN", geo.type="Global", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
#title(main="Ectotherms", cex.main=0.8, line=-0.3)
mtext("F", side=3,font = 2, at=-1.5e+07, cex=0.6, line=-0.8)
mtext('dN', side=2, at=-3.0e+05, cex=0.6, line=-0.1, font=2)


plot_molrate_grids(mean.rate.grids=endotherm%>%mutate(dN=gm.dn)%>%filter(Habitat!="Marine"), 
                   rate.type = "dN", geo.type="Terrestrial", colramp=c("#ffffff", "#f4cdb9","#f47155","#d22428", "#660b11"))
#title(main="Endotherms", cex.main=0.8, line=-0.3)

#dnds
plot_molrate_grids(mean.rate.grids=ectotherm%>%mutate(dNdS=gm.dnds), 
                   rate.type = "dNdS", geo.type="Global", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
#title(main="Ectotherms", cex.main=0.8, line=-0.3)
mtext("H", side=3,font = 2, at=-1.5e+07, cex=0.6, line=-0.8)
mtext('dN/dS', side=2, at=-3.0e+05, cex=0.6, line=-0.1, font=2)


plot_molrate_grids(mean.rate.grids=endotherm%>%mutate(dNdS=gm.dnds)%>%filter(Habitat!="Marine"), 
                   rate.type = "dNdS", geo.type="Terrestrial", colramp=c("#ffffff", "#f4cdb9","#f47155","#d22428", "#660b11"))
#title(main="Endotherms", cex.main=0.8, line=-0.3)
dev.off()

#P values
sar.summary.out%>%filter(Avg.Value=="Geometric mean")%>%filter(Group%in% c("Ectotherms","Endotherms"))%>%
  dplyr::select(1:4,"SAR.Pvalue")


#Fig.3EGI Plots show relationship between molecular rates and latitude at ecoregions 
f3e1<-ectotherm.ecos%>%
  filter(!is.na(gm.ds))%>%
  ggplot(aes(x=AnnualTemp, y=gm.ds*10^9))+
  geom_point(alpha=0.7, size=0.7, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="", y=expression("dS (" * 10^-9 * ")"~(sub/site/year)), title="Ectotherms", tag="E")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = -5, y = 2, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.tag = element_text(size=8, face = "bold"),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f3e2<-endotherm.ecos%>%
  filter(!is.na(gm.ds), Habitat!="Marine")%>%
  ggplot(aes(x=AnnualTemp, y=gm.ds*10^9))+
  geom_point(alpha=0.7, size=0.7, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="", y="", title="Endotherms")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = -5, y = 0.7, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f3g1<-ectotherm.ecos%>%
  filter(!is.na(gm.dn))%>%
  ggplot(aes(x=AnnualTemp, y=gm.dn*10^10))+
  geom_point(alpha=0.7, size=0.7, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="", y=expression("dN (" * 10^-10 * ")"~(sub/site/year)), tag="G")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = -5, y = 1.5, label = expression(italic(P)[SAR] == 0.25), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.tag = element_text(size=8, face = "bold"),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f3g2<-endotherm.ecos%>%
  filter(!is.na(gm.dn), Habitat!="Marine")%>%
  ggplot(aes(x=AnnualTemp, y=gm.dn*10^10))+
  geom_point(alpha=0.7, size=0.7, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="", y="")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = -5, y = 1, label = expression(italic(P)[SAR] == 0.002), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f3i1<-ectotherm.ecos%>%
  filter(!is.na(gm.dnds))%>%
  ggplot(aes(x=AnnualTemp, y=gm.dnds))+
  geom_point(alpha=0.7, size=0.7, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x=expression(Temperature~(degree*C)), y="dN/dS", tag="I")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = -5, y = 0.16, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.tag = element_text(size=8, face = "bold"),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))



f3i2<-endotherm.ecos%>%
  filter(!is.na(gm.dnds), Habitat!="Marine")%>%
  ggplot(aes(x=AnnualTemp, y=gm.dnds))+
  geom_point(alpha=0.7, size=0.7, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x=expression(Temperature~(degree*C)), y="")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits=c(0.12, 0.16), breaks = seq(0.12, 0.16, 0.01))+
  annotate("text", x = -5, y = 0.16, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

cowplot::plot_grid(f3e1,f3e2, f3g1, f3g2, f3i1, f3i2, align = "hv", nrow = 3)
ggsave("./Outputs/MainFigures/Fig3EGI.pdf", height = 4.83, width=3.4)

#cowplot::plot_grid(f3a,f3b,f3c, f3e1, f3g1, f3i1, f3e2, f3g2, f3i2, align = "hv", byrow=F, nrow = 3)
#ggsave("./Outputs/MainFigures/Fig3.pdf", height = 4.83, width=6)

#Fig.3EGI Plots show relationship between molecular rates and latitude at ecoregions 
f3e1<-ectotherm.ecos%>%
  filter(!is.na(gm.ds))%>%
  ggplot(aes(x=abs(Lat), y=gm.ds*100))+
  geom_point(alpha=0.7, size=0.7, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="", y=expression("dS (" * 10^-8 * ")"~(sub/site/year)), title="Ectotherms", tag="E")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = 60, y = 1, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.tag = element_text(size=8, face = "bold"),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f3e2<-endotherm.ecos%>%
  filter(!is.na(gm.ds), Habitat!="Marine")%>%
  ggplot(aes(x=abs(Lat), y=gm.ds*100))+
  geom_point(alpha=0.7, size=0.7, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="", y="", title="Endotherms")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = 60, y = 0.9, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f3g1<-ectotherm.ecos%>%
  filter(!is.na(gm.dn))%>%
  ggplot(aes(x=abs(Lat), y=gm.dn*10000))+
  geom_point(alpha=0.7, size=0.7, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="", y=expression("dN (" * 10^-10 * ")"~(sub/site/year)), tag="G")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = 60, y = 3, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.tag = element_text(size=8, face = "bold"),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f3g2<-endotherm.ecos%>%
  filter(!is.na(gm.dn), Habitat!="Marine")%>%
  ggplot(aes(x=abs(Lat), y=gm.dn*10000))+
  geom_point(alpha=0.7, size=0.7, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="", y="")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = 60, y = 4, label = expression(italic(P)[SAR] == 0.88), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f3i1<-ectotherm.ecos%>%
  filter(!is.na(gm.dnds))%>%
  ggplot(aes(x=abs(Lat), y=gm.dnds))+
  geom_point(alpha=0.7, size=0.7, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y="dN/dS", tag="I")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = 60, y = 0.05, label = expression(italic(P)[SAR] == 0.42), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.tag = element_text(size=8, face = "bold"),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))



f3i2<-endotherm.ecos%>%
  filter(!is.na(gm.dnds), Habitat!="Marine")%>%
  ggplot(aes(x=abs(Lat), y=gm.dnds))+
  geom_point(alpha=0.7, size=0.7, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="Latitude", y="")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits=c(0.022, 0.035), breaks = seq(0.02, 0.035, 0.005))+
  annotate("text", x = 60, y = 0.035, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

cowplot::plot_grid(f3e1,f3e2, f3g1, f3g2, f3i1, f3i2, align = "hv", nrow = 3)
ggsave("./Outputs/MainFigures/Fig3DFH.pdf", height = 5, width=3.4)


##################################################################
#Part II: R script to predict molecular rates using multiple PGLMMs
###################################################################
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
molrate <- read.csv(paste0(workdir, "/DataFiles/MolEvolRate/MolRate_nudna.csv"))%>%
  mutate(AnnualTemp=ifelse(ThermoMode=="Ectotherms"&AnnualTemp<0.01, 0.01, AnnualTemp))


row.names(molrate) <- molrate$Species

#input phylogeny
phylo <- read.tree(paste0(workdir, "/DataFiles/trees/phy_all_sampled_nudna.tre"))

# If the tree is not ultrametric, force it to be ultrametric
phylo <- force.ultrametric(phylo, method = "nnls")  # Use nnls or extend


#get most recent common ancestor (node) of species pair
Ainv <- inverseA(phylo)$Ainv


############################################################
#1 PGLMMs for each class, endotherm and ectotherms using all data
#################################################################
#1.1  Multiple pglmms account for random effects of phylogenetic signal.
pred1 <- c("AnnualTemp", "BodyMass", "Longevity","MaturityAge","Fecundity")
pred2 <- c("AnnualTemp", "BodyMass", "Longevity","MaturityAge","Fecundity","dS")

prior <- list(R=list(V=diag(1),nu=0.002), 
              G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*1000)))

groups <- c("Ectotherms", "Endotherms", "Fishes", "Amphibians", "Reptiles", "Mammals","Birds")

multi_pglmm_list <- vector("list", length = length(groups) * 3)
names(multi_pglmm_list) <- paste(rep(groups, each = 3), c("dN", "dS", "dNdS"), sep = ".")

# Prior specification for the model
for (group in groups) {
  # Choose data based on group
  if(group %in% unique(molrate$ThermoMode)){
    subdata <- molrate%>% 
      filter(ThermoMode == group)%>%
      mutate(dS=scale(log(dS)), dN=scale(log(dN)), dNdS=scale(log(dNdS)), AnnualTemp=scale(AnnualTemp),
             BodyMass=scale(log(BodyMass)), Fecundity=scale(log(Fecundity)), 
             MaturityAge=scale(log(MaturityAge)), Longevity=scale(log(Longevity)))
    
  }else{
    subdata <- molrate %>% 
      filter(Group == group)%>%
      mutate(dS=scale(log(dS)), dN=scale(log(dN)), dNdS=scale(log(dNdS)), AnnualTemp=scale(AnnualTemp),
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
  fit_dnds <- MCMCglmm(formula(paste("dNdS~", paste(pred1, collapse="+"))), 
                       random=~Species, prior=prior,
                       ginverse = list(Species=inverseA(phylo)$Ainv),
                       data=subdata,verbose=TRUE, 
                       nitt=60000, burnin=10000, thin=25, family = c("gaussian"))
  # Store models in list
  multi_pglmm_list[[paste(group, "dS", sep = ".")]] <- fit_ds
  multi_pglmm_list[[paste(group, "dN", sep = ".")]] <- fit_dn
  multi_pglmm_list[[paste(group, "dNdS", sep = ".")]] <- fit_dnds
  
  print(group)
}


#save(multi_pglmm_list, file="./Outputs/Data/multiple_pglmm_nudna.rdata")
load(file="./Outputs/Data/multiple_pglmm_nudna.rdata")


#########
#function to
summary_multi_pglmm_groups <- function(groups, molrate, pglmms){
  fixed <- lapply(paste0(groups, ".", molrate), function(x) summary(pglmms[[x]])$solutions)
  random <- lapply(paste0(groups, ".", molrate), function(x) summary(pglmms[[x]])$Gcovariances)
  r2 <- lapply(paste0(groups, ".", molrate), function(x) r.squaredMCMCglmm(pglmms[[x]]))
  names(r2)<- paste0(groups, ".", molrate)
  
  r2.ci <- lapply(paste0(groups, ".", molrate), function(x) 
    matrix(c(r2[[x]]$R2m, r2[[x]]$R2m.CIl, r2[[x]]$R2m.CIu, nrow(pglmms[[x]]$Sol),
             r2[[x]]$R2c, r2[[x]]$R2c.CIl, r2[[x]]$R2c.CIu, nrow(pglmms[[x]]$Sol),
             r2[[x]]$R2r, r2[[x]]$R2r.CIl, r2[[x]]$R2r.CIu, nrow(pglmms[[x]]$Sol)), 
           nrow=3, ncol=4, byrow = T))
  r2.ci <- lapply(r2.ci, function(mat) {
    row.names(mat) <- c("Rm2", "Rc2", "Rr2")
    colnames(mat) <- colnames(random[[1]]) 
    return(mat)
  })
  
  
  pglmm_summary_out <- lapply(1:length(groups), function(x)
    cbind(data.frame(MolRate=molrate, Group=groups[x], Var=c(row.names(fixed[[x]]), row.names(random[[x]]), row.names(r2.ci[[x]]))), 
          data.frame(rbind(fixed[[x]], cbind(random[[x]], pMCMC=NA), cbind(r2.ci[[x]], pMCMC=NA)))))
  
  return(do.call(rbind, pglmm_summary_out))
}

pglmm_summary_out <- rbind(summary_multi_pglmm_groups(groups=groups, molrate="dS", multi_pglmm_list),
                           summary_multi_pglmm_groups(groups=groups, molrate="dN", multi_pglmm_list),
                           summary_multi_pglmm_groups(groups=groups, molrate="dNdS", multi_pglmm_list))

#write.csv(pglmm_summary_out,"pglmm_summary_out.csv")

#Fig.3JKL
f3j1 <- pglmm_summary_out%>%
  filter(Var != "(Intercept)", !is.na(pMCMC))%>%
  mutate(ThermoMode=ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"),
         sig=ifelse(pMCMC<0.05, "1", "0"))%>%
  filter(MolRate=="dS")%>%
  filter(Group =="Ectotherms")%>%
  ggplot(aes(x=Var, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=l.95..CI, ymax=u.95..CI), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="dS", title = "Ectotherms", tag="J")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  #scale_y_continuous(limits = c(-0.45,0.21), breaks = seq(-0.4,0.2,0.2))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=7, color = "black"),
        axis.text = element_text(size=7, color = "black"),
        plot.tag = element_text(size=8, face = "bold"),
        strip.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()


f3j2 <- pglmm_summary_out%>%
  filter(Var != "(Intercept)", !is.na(pMCMC))%>%
  mutate(ThermoMode=ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"),
         sig=ifelse(pMCMC<0.05, "1", "0"))%>%
  filter(MolRate=="dS")%>%
  filter(Group =="Endotherms")%>%
  ggplot(aes(x=Var, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=l.95..CI, ymax=u.95..CI), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="dS", title = "Endotherms")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  #scale_y_continuous(limits = c(-0.46,0.21), breaks = seq(-0.4,0.2,0.2))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=7, color = "black"),
        axis.text = element_text(size=7, color = "black"),
        plot.tag = element_text(size=8, face = "bold"),
        strip.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()



f3k1 <- pglmm_summary_out%>%
  filter(Var != "(Intercept)", !is.na(pMCMC))%>%
  mutate(ThermoMode=ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"),
         sig=ifelse(pMCMC<0.05, "1", "0"))%>%
  filter(MolRate=="dN")%>%
  filter(Group =="Ectotherms")%>%
  ggplot(aes(x=Var, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=l.95..CI, ymax=u.95..CI), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="dN", title = "Ectotherms", tag="K")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.2,1.2), breaks = seq(0,1.2,0.4))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=7, color = "black"),
        axis.text = element_text(size=7, color = "black"),
        plot.tag = element_text(size=8, face = "bold"),
        strip.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()


f3k2 <- pglmm_summary_out%>%
  filter(Var != "(Intercept)", !is.na(pMCMC))%>%
  mutate(ThermoMode=ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"),
         sig=ifelse(pMCMC<0.05, "1", "0"))%>%
  filter(MolRate=="dN")%>%
  filter(Group =="Endotherms")%>%
  ggplot(aes(x=Var, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=l.95..CI, ymax=u.95..CI), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="dN", title = "Endotherms")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.2,1.2), breaks = seq(0,1.2,0.4))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=7, color = "black"),
        axis.text = element_text(size=7, color = "black"),
        plot.tag = element_text(size=8, face = "bold"),
        strip.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

f3l1 <- pglmm_summary_out%>%
  filter(Var != "(Intercept)", !is.na(pMCMC))%>%
  mutate(ThermoMode=ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"),
         sig=ifelse(pMCMC<0.05, "1", "0"))%>%
  filter(MolRate=="dNdS")%>%
  filter(Group =="Ectotherms")%>%
  ggplot(aes(x=Var, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=l.95..CI, ymax=u.95..CI), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="dN/dS", title = "Ectotherms", tag="L")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  #scale_y_continuous(limits = c(-0.22,0.4), breaks = seq(-0.2,0.4,0.2))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=7, color = "black"),
        axis.text = element_text(size=7, color = "black"),
        plot.tag = element_text(size=8, face = "bold"),
        strip.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()


f3l2 <- pglmm_summary_out%>%
  filter(Var != "(Intercept)", !is.na(pMCMC))%>%
  mutate(ThermoMode=ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"),
         sig=ifelse(pMCMC<0.05, "1", "0"))%>%
  filter(MolRate=="dNdS")%>%
  filter(Group =="Endotherms")%>%
  ggplot(aes(x=Var, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=l.95..CI, ymax=u.95..CI), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="dN/dS", title = "Endotherms")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  #scale_y_continuous(limits = c(-0.21,0.42), breaks = seq(-0.2,0.4,0.2))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=7, color = "black"),
        axis.text = element_text(size=7, color = "black"),
        plot.tag = element_text(size=8, face = "bold"),
        strip.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()




f3j1+f3j2+f3k1+f3k2+f3l1+f3l2+
  plot_layout(ncol = 6, nrow = 1, byrow = F)

#ggsave(filename="./Outputs/MainFigures/Fig3JKL.pdf", height=1.8, width=11)



#Extended Fig.S13
fs13.1 <- pglmm_summary_out%>%
  filter(Var != "(Intercept)", !is.na(pMCMC))%>%
  mutate(ThermoMode=ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"),
         sig=ifelse(pMCMC<0.05, "1", "0"))%>%
  filter(MolRate=="dS")%>%
  filter(Group =="Fishes")%>%
  ggplot(aes(x=Var, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=l.95..CI, ymax=u.95..CI), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="dS", title = "Fishes")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.5,0.5), breaks = seq(-0.5,0.5,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8, color="black"),
        axis.text = element_text(size=8, color="black"),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

fs13.2 <- pglmm_summary_out%>%
  filter(Var != "(Intercept)", !is.na(pMCMC))%>%
  mutate(ThermoMode=ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"),
         sig=ifelse(pMCMC<0.05, "1", "0"))%>%
  filter(MolRate=="dS")%>%
  filter(Group =="Amphibians")%>%
  ggplot(aes(x=Var, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=l.95..CI, ymax=u.95..CI), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="", title = "Amphibians")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.6,0.6), breaks = c(-0.6,-0.4,-0.2, 0, 0.2, 0.4,0.6))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8, color="black"),
        axis.text.x = element_text(size=8, color="black"),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

fs13.3 <- pglmm_summary_out%>%
  filter(Var != "(Intercept)", !is.na(pMCMC))%>%
  mutate(ThermoMode=ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"),
         sig=ifelse(pMCMC<0.05, "1", "0"))%>%
  filter(MolRate=="dS")%>%
  filter(Group =="Reptiles")%>%
  ggplot(aes(x=Var, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=l.95..CI, ymax=u.95..CI), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="", title = "Reptiles")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.5,0.5), breaks = seq(-0.5,0.5,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8, color="black"),
        axis.text.x = element_text(size=8, color="black"),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

fs13.4 <- pglmm_summary_out%>%
  filter(Var != "(Intercept)", !is.na(pMCMC))%>%
  mutate(ThermoMode=ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"),
         sig=ifelse(pMCMC<0.05, "1", "0"))%>%
  filter(MolRate=="dS")%>%
  filter(Group =="Mammals")%>%
  ggplot(aes(x=Var, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=l.95..CI, ymax=u.95..CI), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="", title = "Mammals")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.63,0.5), breaks = seq(-0.5,0.5,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8, color="black"),
        axis.text.x = element_text(size=8, color="black"),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

fs13.5 <- pglmm_summary_out%>%
  filter(Var != "(Intercept)", !is.na(pMCMC))%>%
  mutate(ThermoMode=ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"),
         sig=ifelse(pMCMC<0.05, "1", "0"))%>%
  filter(MolRate=="dS")%>%
  filter(Group =="Birds")%>%
  ggplot(aes(x=Var, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=l.95..CI, ymax=u.95..CI), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="", title = "Birds")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.5,0.5), breaks = seq(-0.5,0.5,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8, color="black"),
        axis.text.x = element_text(size=8, color="black"),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

fs13.6 <- pglmm_summary_out%>%
  filter(Var != "(Intercept)", !is.na(pMCMC))%>%
  mutate(ThermoMode=ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"),
         sig=ifelse(pMCMC<0.05, "1", "0"))%>%
  filter(MolRate=="dN")%>%
  filter(Group =="Fishes")%>%
  ggplot(aes(x=Var, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=l.95..CI, ymax=u.95..CI), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="dN", title = "Fishes")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.5,1), breaks = seq(-0.5,1,0.5))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8, color="black"),
        axis.text = element_text(size=8, color="black"),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

fs13.7 <- pglmm_summary_out%>%
  filter(Var != "(Intercept)", !is.na(pMCMC))%>%
  mutate(ThermoMode=ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"),
         sig=ifelse(pMCMC<0.05, "1", "0"))%>%
  filter(MolRate=="dN")%>%
  filter(Group =="Amphibians")%>%
  ggplot(aes(x=Var, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=l.95..CI, ymax=u.95..CI), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="", title = "Amphibians")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  #scale_y_continuous(limits = c(-0.5,0.75), breaks = seq(-0.5,0.75,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8, color="black"),
        axis.text.x = element_text(size=8, color="black"),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

fs13.8 <- pglmm_summary_out%>%
  filter(Var != "(Intercept)", !is.na(pMCMC))%>%
  mutate(ThermoMode=ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"),
         sig=ifelse(pMCMC<0.05, "1", "0"))%>%
  filter(MolRate=="dN")%>%
  filter(Group =="Reptiles")%>%
  ggplot(aes(x=Var, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=l.95..CI, ymax=u.95..CI), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="", title = "Reptiles")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.3,1.2), breaks = seq(-0.3,1.2,0.3))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8, color="black"),
        axis.text.x = element_text(size=8, color="black"),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

fs13.9 <- pglmm_summary_out%>%
  filter(Var != "(Intercept)", !is.na(pMCMC))%>%
  mutate(ThermoMode=ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"),
         sig=ifelse(pMCMC<0.05, "1", "0"))%>%
  filter(MolRate=="dN")%>%
  filter(Group =="Mammals")%>%
  ggplot(aes(x=Var, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=l.95..CI, ymax=u.95..CI), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="", title = "Mammals")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  #scale_y_continuous(limits = c(-0.5,0.75), breaks = seq(-0.5,0.75,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8, color="black"),
        axis.text.x = element_text(size=8, color="black"),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

fs13.10 <- pglmm_summary_out%>%
  filter(Var != "(Intercept)", !is.na(pMCMC))%>%
  mutate(ThermoMode=ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"),
         sig=ifelse(pMCMC<0.05, "1", "0"))%>%
  filter(MolRate=="dN")%>%
  filter(Group =="Birds")%>%
  ggplot(aes(x=Var, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=l.95..CI, ymax=u.95..CI), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="", title = "Birds")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  #scale_y_continuous(limits = c(-0.5,0.75), breaks = seq(-0.5,0.75,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8, color="black"),
        axis.text.x = element_text(size=8, color="black"),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

fs13.11 <- pglmm_summary_out%>%
  filter(Var != "(Intercept)", !is.na(pMCMC))%>%
  mutate(ThermoMode=ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"),
         sig=ifelse(pMCMC<0.05, "1", "0"))%>%
  filter(MolRate=="dNdS")%>%
  filter(Group =="Fishes")%>%
  ggplot(aes(x=Var, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=l.95..CI, ymax=u.95..CI), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="dNdS", title = "Fishes")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  #scale_y_continuous(limits = c(-0.5,0.75), breaks = seq(-0.5,0.75,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8, color="black"),
        axis.text = element_text(size=8, color="black"),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

fs13.12 <- pglmm_summary_out%>%
  filter(Var != "(Intercept)", !is.na(pMCMC))%>%
  mutate(ThermoMode=ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"),
         sig=ifelse(pMCMC<0.05, "1", "0"))%>%
  filter(MolRate=="dNdS")%>%
  filter(Group =="Amphibians")%>%
  ggplot(aes(x=Var, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=l.95..CI, ymax=u.95..CI), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="", title = "Amphibians")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.6,0.6), breaks = seq(-0.6,0.6,0.3))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8, color="black"),
        axis.text.x = element_text(size=8, color="black"),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

fs13.13 <- pglmm_summary_out%>%
  filter(Var != "(Intercept)", !is.na(pMCMC))%>%
  mutate(ThermoMode=ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"),
         sig=ifelse(pMCMC<0.05, "1", "0"))%>%
  filter(MolRate=="dNdS")%>%
  filter(Group =="Reptiles")%>%
  ggplot(aes(x=Var, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=l.95..CI, ymax=u.95..CI), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="", title = "Reptiles")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  #scale_y_continuous(limits = c(-0.5,0.75), breaks = seq(-0.5,0.75,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8, color="black"),
        axis.text.x = element_text(size=8, color="black"),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

fs13.14 <- pglmm_summary_out%>%
  filter(Var != "(Intercept)", !is.na(pMCMC))%>%
  mutate(ThermoMode=ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"),
         sig=ifelse(pMCMC<0.05, "1", "0"))%>%
  filter(MolRate=="dNdS")%>%
  filter(Group =="Mammals")%>%
  ggplot(aes(x=Var, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=l.95..CI, ymax=u.95..CI), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="", title = "Mammals")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  #scale_y_continuous(limits = c(-0.5,0.75), breaks = seq(-0.5,0.75,0.25))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8, color="black"),
        axis.text.x = element_text(size=8, color="black"),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

fs13.15 <- pglmm_summary_out%>%
  filter(Var != "(Intercept)", !is.na(pMCMC))%>%
  mutate(ThermoMode=ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"),
         sig=ifelse(pMCMC<0.05, "1", "0"))%>%
  filter(MolRate=="dNdS")%>%
  filter(Group =="Birds")%>%
  ggplot(aes(x=Var, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=l.95..CI, ymax=u.95..CI), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="", title = "Birds")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.4,0.8), breaks = seq(-0.4,0.8,0.4))+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8, color="black"),
        axis.text.x = element_text(size=8, color="black"),
        axis.text.y = element_blank(),
        strip.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))+
  coord_flip()

fs13.1+fs13.2+fs13.3+fs13.4+fs13.5+
  fs13.6+fs13.7+fs13.8+fs13.9+fs13.10+
  fs13.11+fs13.12+fs13.13+fs13.14+fs13.15+
  plot_layout(ncol = 5, nrow = 3)

ggsave(filename="./Outputs/Supplementary/Extended Fig13.pdf", height=5.5, width=8.27)



#################################################
#1.2 Multiple PGLMMs account for random effects of phylogenetic signal under the best-fit evolutionary model
load(file="./Outputs/Data/best_trait_evol_models_nudna.rdata")

pred1 <- c("AnnualTemp", "BodyMass", "Longevity","MaturityAge","Fecundity")
pred2 <- c("AnnualTemp", "BodyMass", "Longevity","MaturityAge","Fecundity","dS")

prior <- list(R=list(V=diag(1),nu=0.002), 
              G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*1000)))

groups <- c("Ectotherms", "Endotherms", "Fishes", "Amphibians", "Reptiles", "Mammals","Birds")

multi_pglmm_list_best_model <- vector("list", length = length(groups) * 3)
names(multi_pglmm_list_best_model) <- paste(rep(groups, each = 3), c("dN", "dS", "dNdS"), sep = ".")

# Prior specification for the model
for (group in groups) {
  # Choose data based on group
  if(group %in% unique(molrate$ThermoMode)){
    subdata <- molrate%>% 
      filter(ThermoMode == group)%>%
      mutate(dS=scale(log(dS)), dN=scale(log(dN)),dNdS=scale(log(dNdS)), AnnualTemp=scale(AnnualTemp),
             BodyMass=scale(log(BodyMass)), Fecundity=scale(log(Fecundity)), 
             MaturityAge=scale(log(MaturityAge)), Longevity=scale(log(Longevity)))
    
  }else{
    subdata <- molrate %>% 
      filter(Group == group)%>%
      mutate(dS=scale(log(dS)), dN=scale(log(dN)),dNdS=scale(log(dNdS)), AnnualTemp=scale(AnnualTemp),
             BodyMass=scale(log(BodyMass)), Fecundity=scale(log(Fecundity)), 
             MaturityAge=scale(log(MaturityAge)), Longevity=scale(log(Longevity)))
  }
  
  phy <- drop.tip(phylo, setdiff(phylo$tip.label, subdata$Species))
  
  #The best traits evolution model of dS
  model.ds <- best_trait_evol_models[[paste(group, "dS", sep = ".")]]
  
  #The best traits evolution model of dN
  model.dn <- best_trait_evol_models[[paste(group, "dN", sep = ".")]]
  
  #The best traits evolution model of dN
  model.dnds <- best_trait_evol_models[[paste(group, "dNdS", sep = ".")]]
  
  #rescale tree
  if(model.ds$best_model=="OU"){phy_ds <- phytools::rescale(phy, model=model.ds$best_model, alpha=model.ds$opt, sigsq=model.ds$sigsq)}
  if(model.ds$best_model=="BM"){phy_ds <- phytools::rescale(phy, model=model.ds$best_model, sigsq=model.ds$sigsq)}
  if(model.ds$best_model=="EB"){phy_ds <- phytools::rescale(phy, model=model.ds$best_model, a=model.ds$opt, sigsq=model.ds$sigsq)}
  
  #rescale tree
  if(model.dn$best_model=="OU"){phy_dn <- phytools::rescale(phy, model=model.dn$best_model, alpha=model.dn$opt, sigsq=model.dn$sigsq)}
  if(model.dn$best_model=="BM"){phy_dn <- phytools::rescale(phy, model=model.dn$best_model, sigsq=model.dn$sigsq)}
  if(model.dn$best_model=="EB"){phy_dn <- phytools::rescale(phy, model=model.dn$best_model, a=model.dn$opt, sigsq=model.dn$sigsq)}
  
  #rescale tree
  if(model.dnds$best_model=="OU"){phy_dnds <- phytools::rescale(phy, model=model.dnds$best_model, alpha=model.dnds$opt, sigsq=model.dnds$sigsq)}
  if(model.dnds$best_model=="BM"){phy_dnds <- phytools::rescale(phy, model=model.dnds$best_model, sigsq=model.dnds$sigsq)}
  if(model.dnds$best_model=="EB"){phy_dnds <- phytools::rescale(phy, model=model.dnds$best_model, a=model.dnds$opt, sigsq=model.dnds$sigsq)}
  
  
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
  for (i in 1:length(phy_dnds$edge.length)) {
    if (phy_dnds$edge.length[i] < 1e-16) {
      phy_dnds$edge.length[i] = 1e-16
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
  fit_dnds <- MCMCglmm(formula(paste("dNdS~", paste(pred1, collapse="+"))), 
                       random=~Species, prior=prior,
                       ginverse = list(Species=inverseA(phy_dnds)$Ainv),
                       data=subdata,verbose=TRUE, 
                       nitt=60000, burnin=10000, thin=25, family = c("gaussian"))
  
  # Store models in list
  multi_pglmm_list_best_model[[paste(group, "dS", sep = ".")]] <- fit_ds
  multi_pglmm_list_best_model[[paste(group, "dN", sep = ".")]] <- fit_dn
  multi_pglmm_list_best_model[[paste(group, "dNdS", sep = ".")]] <- fit_dnds
  print(group)
}

#save(multi_pglmm_list_best_model, file="./Outputs/Data/multiple_pglmm_best_traits_model_nudna.rdata")
load(file="./Outputs/Data/multiple_pglmm_best_traits_model_nudna.rdata")


pglmm_summary_out_best_model <- rbind(summary_multi_pglmm_groups(groups=groups, molrate="dS", multi_pglmm_list_best_model),
                                      summary_multi_pglmm_groups(groups=groups, molrate="dN", multi_pglmm_list_best_model),
                                      summary_multi_pglmm_groups(groups=groups, molrate="dNdS", multi_pglmm_list_best_model))

#write.csv(pglmm_summary_out_best_model,"pglmm_summary_out_best_model.csv")

################################################
#1.3 Multiple PGLMMa account for randoms effects of phylogenetic signal, habitat and lineage for estimating molecular rates.
groups <- c("Ectotherms", "Endotherms", "Fishes", "Amphibians", "Reptiles", "Mammals","Birds")

multi_pglmm_list_3randoms <- vector("list", length = length(groups) * 3)
names(multi_pglmm_list_3randoms) <- paste(rep(groups, each = 3), c("dN", "dS", "dNdS"), sep = ".")

prior <- list(R=list(V=diag(1),nu=0.002), 
              G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*1000),
                     G2=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*1000),
                     G3=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*1000)))

for (group in groups) {
  # Choose data based on group
  if(group %in% unique(molrate$ThermoMode)){
    subdata <- molrate%>% 
      filter(ThermoMode == group)%>%
      mutate(dS=scale(log(dS)), dN=scale(log(dN)),dNdS=scale(log(dNdS)), AnnualTemp=scale(AnnualTemp),
             BodyMass=scale(log(BodyMass)), Fecundity=scale(log(Fecundity)), 
             MaturityAge=scale(log(MaturityAge)), Longevity=scale(log(Longevity)))
    
  }else{
    subdata <- molrate %>% 
      filter(Group == group)%>%
      mutate(dS=scale(log(dS)), dN=scale(log(dN)),dNdS=scale(log(dNdS)), AnnualTemp=scale(AnnualTemp),
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
  fit_dnds <- MCMCglmm(formula(paste("dNdS~", paste(pred1, collapse="+"))), 
                       random=~Species+Habitat+Clade, prior=prior,
                       ginverse = list(Species=inverseA(phylo)$Ainv),
                       data=subdata,verbose=TRUE, 
                       nitt=60000, burnin=10000, thin=25, family = c("gaussian"))
  
  # Store models in list
  multi_pglmm_list_3randoms[[paste(group, "dS", sep = ".")]] <- fit_ds
  multi_pglmm_list_3randoms[[paste(group, "dN", sep = ".")]] <- fit_dn
  multi_pglmm_list_3randoms[[paste(group, "dNdS", sep = ".")]] <- fit_dnds
  print(group)
}
#save(multi_pglmm_list_3randoms, file="./Outputs/Data/multiple_pglmm_three_random_effects_nudna.rdata")
load(file="./Outputs/Data/multiple_pglmm_three_random_effects_nudna.rdata")

pglmm_summary_out_3randoms <- rbind(summary_multi_pglmm_groups(groups=groups, molrate="dS", multi_pglmm_list_3randoms),
                                    summary_multi_pglmm_groups(groups=groups, molrate="dN", multi_pglmm_list_3randoms),
                                    summary_multi_pglmm_groups(groups=groups, molrate="dNdS", multi_pglmm_list_3randoms))

#write.csv(pglmm_summary_out_3randoms, "pglmm_summary_out_3randoms.csv")


#######################
#Part III. R script to examine relationships between diversification rates and molecular rates
#########################
library(tidyverse)
library(MCMCglmm)
library(MuMIn)
library(patchwork)
rm(list=ls())
workdir <- "/Users/Tianlong/VertMolRate"
setwd(workdir)
sisters <- read.csv("./DataFiles/sisters_family/sisters_nudna.csv", row.names = 1)

groups <- c(unique(sisters$Group))

ols_list <- vector("list", length = length(groups) * 3)
names(ols_list) <- paste(rep(groups, each = 3), c("dN", "dS", "dNdS"), sep = ".")


for(group in groups){
  if(group %in% unique(sisters$ThermoMode)){
    subdata <- sisters%>% 
      filter(ThermoMode == group)
  }else{
    subdata <- sisters %>% 
      filter(Group == group)
  }
  
  fit_dn <- lm(diff.dn~diff.spp+0, subdata)
  fit_ds <- lm(diff.ds~diff.spp+0, subdata)
  fit_dnds <- lm(diff.dnds~diff.spp+0, subdata)
  
  ols_list[[paste(group, "dS", sep = ".")]] <- fit_ds
  ols_list[[paste(group, "dN", sep = ".")]] <- fit_dn
  ols_list[[paste(group, "dNdS", sep = ".")]] <- fit_dnds
}

ols_out <- NULL
for(i in 1:length(groups)){
  group <- groups[i]
  fit.ds.out<- summary(ols_list[[paste(group, "dS", sep = ".")]])$coefficients
  fit.dn.out<- summary(ols_list[[paste(group, "dN", sep = ".")]])$coefficients
  fit.dnds.out<- summary(ols_list[[paste(group, "dNdS", sep = ".")]])$coefficients
  
  fit.out <- data.frame(MolRate=c("dS", "dN", "dNdS"), rbind(fit.ds.out, fit.dn.out, fit.dnds.out))%>%
    rename(Slope=Estimate, Slope.SD=Std..Error, t=t.value, p=Pr...t..)%>%
    mutate(Group=group, ThermoMode=ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"))%>%
    select(ThermoMode, Group, MolRate, Slope, Slope.SD, t, p)
  ols_out <- rbind(ols_out, fit.out)
}

ols_out

#Extended Fig.S14
f1 <- sisters%>%
  filter(Group=="Fishes")%>%
  ggplot(aes(x=diff.spp, y=diff.ds, colour=Group))+
  geom_point(size=1)+
  geom_smooth(method="lm", formula = y~x+0)+
  geom_errorbar(aes(ymin=diff.ds-diff.ds.sd, ymax=diff.ds+diff.ds.sd), width=.2, linewidth=0.2)+
  scale_color_manual(values = c("#2a6aaf"))+
  scale_x_continuous(limits = c(0,6), breaks = seq(0,6,2))+
  scale_y_continuous(limits = c(-1,2), breaks = seq(-1,2,0.5))+
  labs(x="Diff. in ln(Clade Size)", y="Diff. in ln(dS)", title="Fishes")+
  annotate("text", x=4,y=2, label="p = 0.88", size=2.5)+
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
  ggplot(aes(x=diff.spp, y=diff.ds, colour=Group))+
  geom_point(size=1)+
  geom_smooth(method="lm", formula = y~x+0)+
  geom_errorbar(aes(ymin=diff.ds-diff.ds.sd, ymax=diff.ds+diff.ds.sd), width=.2, linewidth=0.2)+
  scale_color_manual(values = c("#2a6aaf"))+
  scale_x_continuous(limits = c(0,6), breaks = seq(0,6,2))+
  scale_y_continuous(limits = c(-0.75,0.5), breaks = seq(-0.75,0.5,0.25))+
  labs(x="Diff. in ln(Clade Size)", y="", title="Amphibians")+
  annotate("text", x=4,y=0.5, label="p < 0.001", size=2.5)+
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
  ggplot(aes(x=diff.spp, y=diff.ds, colour=Group))+
  geom_point(size=1)+
  geom_smooth(method="lm", formula = y~x+0)+
  geom_errorbar(aes(ymin=diff.ds-diff.ds.sd, ymax=diff.ds+diff.ds.sd), width=.2, linewidth=0.2)+
  scale_color_manual(values = c("#2a6aaf"))+
  scale_x_continuous(limits = c(0,6), breaks = seq(0,6,2))+
  scale_y_continuous(limits = c(-0.5,0.5), breaks = seq(-0.5,0.5,0.25))+
  labs(x="Diff. in ln(Clade Size)", y="", title="Reptiles")+
  annotate("text", x=3,y=0.5, label="p = 0.02", size=2.5)+
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
  ggplot(aes(x=diff.spp, y=diff.ds, colour=Group))+
  geom_point(size=1)+
  geom_smooth(method="lm", formula = y~x+0)+
  geom_errorbar(aes(ymin=diff.ds-diff.ds.sd, ymax=diff.ds+diff.ds.sd), width=.2, linewidth=0.2)+
  scale_color_manual(values = c("#d3292f"))+
  scale_x_continuous(limits = c(0,4), breaks = seq(0,4,1))+
  scale_y_continuous(limits = c(-0.5,1), breaks = seq(-0.5,1,0.5))+
  labs(x="Diff. in ln(Clade Size)", y="", title="Mammals")+
  annotate("text", x=3,y=1, label="p = 0.33", size=2.5)+
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
  ggplot(aes(x=diff.spp, y=diff.ds, colour=Group))+
  geom_point(size=1)+
  geom_smooth(method="lm", formula = y~x+0)+
  geom_errorbar(aes(ymin=diff.ds-diff.ds.sd, ymax=diff.ds+diff.ds.sd), width=.2, linewidth=0.2)+
  scale_color_manual(values = c("#d3292f"))+
  scale_x_continuous(limits = c(0,6), breaks = seq(0,6,2))+
  scale_y_continuous(limits = c(-0.8,1.2), breaks = seq(-0.8,1.2,0.4))+
  labs(x="Diff. in ln(Clade Size)", y="", title="Birds")+
  annotate("text", x=4,y=1.2, label="p = 0.01", size=2.5)+
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
  ggplot(aes(x=diff.spp, y=diff.dn, colour=Group))+
  geom_point(size=1)+
  geom_smooth(method="lm", formula = y~x+0)+
  geom_errorbar(aes(ymin=diff.dn-diff.dn.sd, ymax=diff.dn+diff.dn.sd), width=.2, linewidth=0.2)+
  scale_color_manual(values = c("#2a6aaf"))+
  scale_x_continuous(limits = c(0,6), breaks = seq(0,6,2))+
  scale_y_continuous(limits = c(-1.1,2), breaks = seq(-1,2,0.5))+
  labs(x="Diff. in ln(Clade Size)", y="Diff. in ln(dN)", title="Fishes")+
  annotate("text", x=4,y=2, label="p = 0.66", size=2.5)+
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
  ggplot(aes(x=diff.spp, y=diff.dn, colour=Group))+
  geom_point(size=1)+
  geom_smooth(method="lm", formula = y~x+0)+
  geom_errorbar(aes(ymin=diff.dn-diff.dn.sd, ymax=diff.dn+diff.dn.sd), width=.2, linewidth=0.2)+
  scale_color_manual(values = c("#2a6aaf"))+
  scale_x_continuous(limits = c(0,6), breaks = seq(0,6,2))+
  scale_y_continuous(limits = c(-0.4,0.6), breaks = seq(-0.4,0.6,0.2))+
  labs(x="Diff. in ln(Clade Size)", y="", title="Amphibians")+
  annotate("text", x=4,y=0.6, label="p < 0.001", size=2.5)+
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
  ggplot(aes(x=diff.spp, y=diff.dn, colour=Group))+
  geom_point(size=1)+
  geom_smooth(method="lm", formula = y~x+0)+
  geom_errorbar(aes(ymin=diff.dn-diff.dn.sd, ymax=diff.dn+diff.dn.sd), width=.2, linewidth=0.2)+
  scale_color_manual(values = c("#2a6aaf"))+
  scale_x_continuous(limits = c(0,6), breaks = seq(0,6,2))+
  scale_y_continuous(limits = c(-0.6,0.6), breaks = round(seq(-0.6,0.6,0.2),2))+
  labs(x="Diff. in ln(Clade Size)", y="", title="Reptiles")+
  annotate("text", x=3,y=0.6, label="p = 0.045", size=2.5)+
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
  ggplot(aes(x=diff.spp, y=diff.dn, colour=Group))+
  geom_point(size=1)+
  geom_smooth(method="lm", formula = y~x+0)+
  geom_errorbar(aes(ymin=diff.dn-diff.dn.sd, ymax=diff.dn+diff.dn.sd), width=.2, linewidth=0.2)+
  scale_color_manual(values = c("#d3292f"))+
  scale_x_continuous(limits = c(0,4), breaks = seq(0,4,1))+
  scale_y_continuous(limits = c(-1,1), breaks = seq(-1,1,0.5))+
  labs(x="Diff. in ln(Clade Size)", y="", title="Mammals")+
  annotate("text", x=3,y=1, label="p = 0.77", size=2.5)+
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
  ggplot(aes(x=diff.spp, y=diff.dn, colour=Group))+
  geom_point(size=1)+
  geom_smooth(method="lm", formula = y~x+0)+
  geom_errorbar(aes(ymin=diff.dn-diff.dn.sd, ymax=diff.dn+diff.dn.sd), width=.2, linewidth=0.2)+
  scale_color_manual(values = c("#d3292f"))+
  scale_x_continuous(limits = c(0,6), breaks = seq(0,6,2))+
  scale_y_continuous(limits = c(-0.6,0.8), breaks = round(seq(-0.6,0.8,0.2),2))+
  labs(x="Diff. in ln(Clade Size)", y="", title="Birds")+
  annotate("text", x=4,y=0.8, label="p = 0.03", size=2.5)+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))


f11 <- sisters%>%
  filter(Group=="Fishes")%>%
  ggplot(aes(x=diff.spp, y=diff.dnds, colour=Group))+
  geom_point(size=1)+
  geom_smooth(method="lm", formula = y~x+0)+
  geom_errorbar(aes(ymin=diff.dnds-diff.dnds.sd, ymax=diff.dnds+diff.dnds.sd), width=.2, linewidth=0.2)+
  scale_color_manual(values = c("#2a6aaf"))+
  scale_x_continuous(limits = c(0,6), breaks = seq(0,6,2))+
  scale_y_continuous(limits = c(-0.6,0.6), breaks = seq(-0.6,0.6,0.3))+
  labs(x="Diff. in ln(Clade Size)", y="Diff. in ln(dN/dS)", title="Fishes")+
  annotate("text", x=4,y=0.6, label="p = 0.56", size=2.5)+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f12 <- sisters%>%
  filter(Group=="Amphibians")%>%
  ggplot(aes(x=diff.spp, y=diff.dnds, colour=Group))+
  geom_point(size=1)+
  geom_smooth(method="lm", formula = y~x+0)+
  geom_errorbar(aes(ymin=diff.dnds-diff.dnds.sd, ymax=diff.dnds+diff.dnds.sd), width=.2, linewidth=0.2)+
  scale_color_manual(values = c("#2a6aaf"))+
  scale_x_continuous(limits = c(0,6), breaks = seq(0,6,2))+
  scale_y_continuous(limits = c(-0.4,0.4), breaks = seq(-0.4,0.4,0.2))+
  labs(x="Diff. in ln(Clade Size)", y="", title="Amphibians")+
  annotate("text", x=4,y=0.4, label="p = 0.11", size=2.5)+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f13 <- sisters%>%
  filter(Group=="Reptiles")%>%
  ggplot(aes(x=diff.spp, y=diff.dnds, colour=Group))+
  geom_point(size=1)+
  geom_smooth(method="lm", formula = y~x+0)+
  geom_errorbar(aes(ymin=diff.dnds-diff.dnds.sd, ymax=diff.dnds+diff.dnds.sd), width=.2, linewidth=0.2)+
  scale_color_manual(values = c("#2a6aaf"))+
  scale_x_continuous(limits = c(0,6), breaks = seq(0,6,2))+
  scale_y_continuous(limits = c(-0.4,0.22), breaks = round(seq(-0.4,0.2,0.2),2))+
  labs(x="Diff. in ln(Clade Size)", y="", title="Reptiles")+
  annotate("text", x=3,y=0.22, label="p = 0.46", size=2.5)+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f14 <- sisters%>%
  filter(Group=="Mammals")%>%
  ggplot(aes(x=diff.spp, y=diff.dnds, colour=Group))+
  geom_point(size=1)+
  geom_smooth(method="lm", formula = y~x+0)+
  geom_errorbar(aes(ymin=diff.dnds-diff.dnds.sd, ymax=diff.dnds+diff.dnds.sd), width=.2, linewidth=0.2)+
  scale_color_manual(values = c("#d3292f"))+
  scale_x_continuous(limits = c(0,6), breaks = seq(0,6,2))+
  scale_y_continuous(limits = c(-0.4,0.2), breaks = seq(-0.4,0.2,0.2))+
  labs(x="Diff. in ln(Clade Size)", y="", title="Mammals")+
  annotate("text", x=3,y=0.2, label="p = 0.009", size=2.5)+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))
f15 <- sisters%>%
  filter(Group=="Birds")%>%
  ggplot(aes(x=diff.spp, y=diff.dnds, colour=Group))+
  geom_point(size=1)+
  geom_smooth(method="lm", formula = y~x+0)+
  geom_errorbar(aes(ymin=diff.dnds-diff.dnds.sd, ymax=diff.dnds+diff.dnds.sd), width=.2, linewidth=0.2)+
  scale_color_manual(values = c("#d3292f"))+
  scale_x_continuous(limits = c(0,6), breaks = seq(0,6,2))+
  scale_y_continuous(limits = c(-0.3,0.3), breaks = round(seq(-0.3,0.3,0.1),2))+
  labs(x="Diff. in ln(Clade Size)", y="", title="Birds")+
  annotate("text", x=4,y=0.3, label="p = 0.06", size=2.5)+
  theme_classic()+
  guides(y=guide_axis(cap="upper"), x=guide_axis(cap="upper"))+
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f11+f12+f13+f14+f15+
  plot_layout(ncol = 5, nrow = 3)

ggsave(filename="./Outputs/Supplementary/Extended Fig14.pdf", height=4.2, width=8.27)
