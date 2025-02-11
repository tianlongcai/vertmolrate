################################
# R script to test evolutionary speed hypothesis (ESH) using molecular rates of mitochondrial genes
# Author: Tianlong Cai
# Email: caitianlong@westlake.edu.cn

##################################################################
#Part I: Analyzing the variation in molecular rates across different species and taxonomic groups.
###################################################################

############################################
rm(list=ls())
gc()
workdir <- "/Users/tianlong/VertMolRate"
setwd(workdir)
library(tidyverse)
library(PupillometryR)
library(cowplot)
library(RColorBrewer)
library(pastecs)
library(phytools)
library(MCMCglmm)
library(geiger)
library(plotrix)

#####################################################################
#Fig1. Plot substitution rates in tips and speciation rate along phylogeny
#functions
source(paste0(workdir, "/Scripts/source_functions.r"))
plot2.phylo <- diversitree:::plot2.phylo

#Input phylogeny
all_vert_phy <- read.tree(paste0(workdir, "/DataFiles/trees/phy_all_sampled_mtdna_outgroup.tre"))

#Input substitution rates
molrate <- read.csv(paste0(workdir, "/DataFiles/MolEvolRate/molrate_mtDNA.csv"))

#Check matched tip names in substitution rates
row.names(molrate) <- molrate$Species
unmatched_species <- name.check(all_vert_phy, molrate)$tree_not_data

#Dropping unmatched tips
sampled_vert_phy<-drop.tip(all_vert_phy, unmatched_species[-which(unmatched_species=="Outgroup")])


molrate<-molrate[sampled_vert_phy$tip.label[-which(sampled_vert_phy$tip.label=="Outgroup")],]%>%
  mutate(MajorGroup=toupper(MajorGroup))

#Basic statistics for raw data
group1 <- molrate%>%
  group_by(ThermoMode)%>%
  summarize(
    n = length(Species),
    dN_mean = mean(dN),
    dN_median = median(dN),
    dN_min = min(dN),
    dN_max = max(dN),
    dN_sd = sd(dN),
    dN_se = dN_sd / sqrt(n),
    dN_Q0.025 = quantile(dN, 0.025),
    dN_Q0.975 = quantile(dN, 0.975),
    dN_95CI_Lower = unlist(t.test(dN)[4])[1],
    dN_95CI_Upper = unlist(t.test(dN)[4])[2],
    dS_mean = mean(dS),
    dS_median = median(dS),
    dS_min = min(dS),
    dS_max = max(dS),
    dS_sd = sd(dS),
    dS_se = dS_sd / sqrt(n),
    dS_Q0.025 = quantile(dS, 0.025),
    dS_Q0.975 = quantile(dS, 0.975),
    dS_95CI_Lower = unlist(t.test(dS)[4])[1],
    dS_95CI_Upper = unlist(t.test(dS)[4])[2],
    dNdS_mean = mean(dNdS),
    dNdS_median = median(dNdS),
    dNdS_min = min(dNdS),
    dNdS_max = max(dNdS),
    dNdS_sd = sd(dNdS),
    dNdS_se = dNdS_sd / sqrt(n),
    dNdS_Q0.025 = quantile(dNdS, 0.025),
    dNdS_Q0.975 = quantile(dNdS, 0.975),
    dNdS_95CI_Lower = unlist(t.test(dNdS)[4])[1],
    dNdS_95CI_Upper = unlist(t.test(dNdS)[4])[2])%>%
  rename(Group=ThermoMode)

group2 <- molrate%>%
  group_by(Group)%>%
  summarize(
    n = length(Species),
    dN_mean = mean(dN),
    dN_median = median(dN),
    dN_min = min(dN),
    dN_max = max(dN),
    dN_sd = sd(dN),
    dN_se = dN_sd / sqrt(n),
    dN_Q0.025 = quantile(dN, 0.025),
    dN_Q0.975 = quantile(dN, 0.975),
    dN_95CI_Lower = unlist(t.test(dN)[4])[1],
    dN_95CI_Upper = unlist(t.test(dN)[4])[2],
    dS_mean = mean(dS),
    dS_median = median(dS),
    dS_min = min(dS),
    dS_max = max(dS),
    dS_sd = sd(dS),
    dS_se = dS_sd / sqrt(n),
    dS_Q0.025 = quantile(dS, 0.025),
    dS_Q0.975 = quantile(dS, 0.975),
    dS_95CI_Lower = unlist(t.test(dS)[4])[1],
    dS_95CI_Upper = unlist(t.test(dS)[4])[2],
    dNdS_mean = mean(dNdS),
    dNdS_median = median(dNdS),
    dNdS_min = min(dNdS),
    dNdS_max = max(dNdS),
    dNdS_sd = sd(dNdS),
    dNdS_se = dNdS_sd / sqrt(n),
    dNdS_Q0.025 = quantile(dNdS, 0.025),
    dNdS_Q0.975 = quantile(dNdS, 0.975),
    dNdS_95CI_Lower = unlist(t.test(dNdS)[4])[1],
    dNdS_95CI_Upper = unlist(t.test(dNdS)[4])[2])

bas_stat <- rbind(group2, group1)

#################################################
#Extract node for calibrations
#The calibration points displayed on the phylogenetic tree are fewer than the actual fossil calibration points used for time estimation 
calibrations<- read.csv("./DataFiles/calibrations.csv")

for(i in 1:nrow(calibrations)) {
  sp1 <- calibrations$Sp1[i]
  sp2 <- calibrations$Sp2[i]
  if(!sp1 %in% sampled_vert_phy$tip.label){
    calibrations[calibrations$Sp1==sp1,"In_data"] <- FALSE
  }
  
  if(!sp2 %in% sampled_vert_phy$tip.label){
    calibrations[calibrations$Sp2==sp2,"In_data"] <- FALSE
  }
}

calibrations <- calibrations%>%
  mutate(Clade=paste0("Calibration", 1:173))%>%
  filter(In_data==TRUE)


fossils <- calibrations$Clade
fossils_node <- matrix(NA, length(fossils), 2)
colnames(fossils_node) <- c("Clade", "Fossil.Node")

for(i in 1:length(fossils)){
  clade <- fossils[i]
  species <- c(calibrations[i,]$Sp1, calibrations[i,]$Sp2)
  node <- getMRCA(sampled_vert_phy, species)
  fossils_node[i,1] <- clade
  fossils_node[i,2] <- node
}
fossils_node <- as.data.frame(fossils_node)
fossils_node$Fossil.Node <- as.numeric(fossils_node$Fossil.Node)

#Extract node for each clade
clades <- unique(molrate$Clade.Label)
clade_nodes <- matrix(NA, length(clades), 2)
colnames(clade_nodes) <- c("Clade", "Clade.Node")

for(i in 1:length(clades)){
  clade <- clades[i]
  species <- molrate%>%filter(Clade.Label==clade)%>%pull(Species)
  node <- getMRCA(sampled_vert_phy, species)
  clade_nodes[i,1] <- clade
  clade_nodes[i,2] <- node
}
clade_nodes <- as.data.frame(clade_nodes)
clade_nodes$Clade.Node <- as.numeric(clade_nodes$Clade.Node)

#Extract node for each major group
groups <- unique(molrate$MajorGroup)
group_nodes <- matrix(NA, length(groups), 2)
colnames(group_nodes) <- c("MajorGroup", "Group.Node")

for(i in 1:length(groups)){
  group <- groups[i]
  species <- molrate%>%filter(MajorGroup==group)%>%pull(Species)
  node <- getMRCA(sampled_vert_phy, species)
  group_nodes[i,1] <- group
  group_nodes[i,2] <- node
}
group_nodes <- as.data.frame(group_nodes)
group_nodes$Group.Node <- as.numeric(group_nodes$Group.Node)

#Merged clade nodes and groups nodes into molrate
molrate <- molrate%>%
  #inner_join(clade_nodes, by=c("Clade.Label"="Clade"))%>%
  inner_join(group_nodes, by="MajorGroup")

#scale dN and dS to plot bars
molrate$dS <- (molrate$dS/max(molrate$dS))*60
molrate$dN <- (molrate$dN/max(molrate$dN))*60
molrate$dNdS <- (molrate$dNdS/max(molrate$dNdS))*60

#define colors of molrate bars
molrate_bar_colour <- c("#93c47d", "#fdb96b", "#2985cc")


#Extracted molrate to plotting bars
molrate_plot <- as.matrix(molrate[,c("dS","dN", "dNdS")])
row.names(molrate_plot) <- molrate$Species


pdf("./Outputs/MainFigures/Fig1.pdf", width = 8.27, height = 8.27)
#plotting phylogeny using different colors to show diversification rate
par(xpd = TRUE)
plot2.phylo(sampled_vert_phy, show.tip.label=FALSE,  cex=0.05, 
            label.offset=0, type="f", edge.width=0.2, no.margin=FALSE, 
            root.edge=TRUE, x.lim=c(-860, 860))
for (j in 1:3) ring(molrate_plot[, j], sampled_vert_phy, offset = j*60 - 60 +2, col = molrate_bar_colour[j])
legend(110, 190, legend=c("dS", "dN", "dN/dS"),bty='n',
       col=molrate_bar_colour, pch=15, cex=0.4)

###########################################################
nodes<-clade_nodes[,2]
clade.labels<-as.character(1:length(clades))

#clade.col <- clade_colors[-1]
#names(clade.col) <- clade.labels
clade.col <- rep(brewer.pal(n = 7, name = 'Greys')[c(3,7)], 43)
names(clade.col) <- clade.labels
group.col <-brewer.pal(n = 9, name = 'Greys')[c(5,8,6,9,9,9,7)]
#group.col <- c("#2a6aaf", "#2a6aaf", "#d3292f", "#2a6aaf", "#2a6aaf", "#2a6aaf", "#d3292f")
#group.col <-brewer.pal(n = 7, name = 'Blues')[c(3,6,4,7,7,7,5)]
#display.brewer.all()
for(i in 1:length(fossils)) {
  arc.cladelabels(text=fossils[i], node=fossils_node$Fossil.Node[i], lab.offset=1.44, ln.offset=1.44, col = "red", lwd=0, point.cex=0.3, point.lwd=0.2, cex=0.1)
} 

for(i in 1:7) {
  arc.cladelabels(text=group_nodes$MajorGroup[i],node=group_nodes$Group.Node[i], lab.offset=1.44, ln.offset=1.44, col = group.col[i], lwd=8, point.cex=0, cex=0.4, text.color = "white")
}

for(i in 1:length(clades)) {
  arc.cladelabels(text=clade.labels[i],node=nodes[i], lab.offset=1.47, ln.offset=1.47, col = clade.col[i], lwd=3, point.cex=0, point.lwd=0.2, cex=0.4)
} 



#Time scale
lines<-rev(seq(0,450,50))
for (i in 1:length(lines)){
  draw.circle(0,0, radius=lines[i],col=NA, lty=2, lwd=0.5, border=grey(0.8,alpha=0.6))
}

text(y=-9, x=seq(0,450,50), labels=rev(seq(0,450,50)), cex=0.3)

dev.off()


#################################################################
#Extended Fig.1: Raincloud plot showing the variation in molecular rates across the five classes of vertebrates.

#input data
molrate <- read.csv(paste0(workdir, "/DataFiles/MolEvolRate/molrateDNA.csv"))%>%
  mutate(Group=factor(Group, levels = c("Fishes","Amphibians", "Reptiles","Mammals","Birds")))%>%
  mutate(dS=dS*10^8, dN=dN*10^10)

range(molrate$dS)
range(molrate$dN)
range(molrate$dNdS)

#Kruskal Wallis test and multiple comparison of groups
comp.ds <- agricolae::kruskal(molrate$dS, molrate$Group, p.adj = "bonferroni")
comp.dn <- agricolae::kruskal(molrate$dN, molrate$Group, p.adj = "bonferroni")
comp.dnds <- agricolae::kruskal(molrate$dNdS, molrate$Group, p.adj = "bonferroni")

#Data for plots
#Extended Fig.S2
groups <- c("Fishes","Amphibians", "Reptiles","Mammals","Birds")
y1 <- molrate%>%group_by(Group)%>%summarise(y=max(dS))%>%pull(y)+0.3
y2 <- molrate%>%group_by(Group)%>%summarise(y=max(dN))%>%pull(y)+0.4
y3 <- molrate%>%group_by(Group)%>%summarise(y=max(dNdS))%>%pull(y)+0.02

label <- data.frame(Group,groups, x=(1:5)+0.2, y1=y1, y2=y2, y3=y3,
                    label.ds=comp.ds$groups[groups,"groups"],
                    label.dn=comp.dn$groups[groups,"groups"],
                    label.dnds=comp.dnds$groups[groups,"groups"],
                    ThermoMode = c("Ectotherms", "Ectotherms", "Ectotherms", "Endotherms","Endotherms"))

f1 <- ggplot(data = molrate, aes(y = dS, x = Group, colour=ThermoMode, fill=ThermoMode)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), adjust =1) +
  geom_point(aes(y = dS), position = position_jitter(width = 0.12), alpha = 0.1, size = 0.5, shape=1) +
  geom_boxplot(width = 0.2, outlier.shape = NA, lwd=0.5, fill = NA) +
  scale_fill_manual(values = c("#2a6aaf", "#d3292f"))+
  scale_color_manual(values = c("#2a6aaf", "#d3292f"))+
  scale_x_discrete(limits=c("Fishes","Amphibians", "Reptiles","Mammals","Birds"))+
  scale_y_continuous(limits = c(0,7), breaks = seq(0,7,1))+
  theme_classic()+
  labs(x="", y = expression("dS (" * 10^-8 * " sub/site/year)"))+
  guides(x=guide_axis(cap='upper'), y=guide_axis(cap='upper'))+
  theme(axis.title = element_text(size=9),
        axis.text = element_text(size = 9),
        legend.position = 'none',
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))+
  geom_text(data=label, aes(x=x, y=y1, label=label.ds), size=2.5)+
  coord_flip()


f2 <-ggplot(data = molrate, aes(y = dN, x = Group, colour=ThermoMode, fill=ThermoMode)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), adjust =1) +
  geom_point(aes(y = dN), position = position_jitter(width = 0.12), alpha = 0.1, size = 0.5, shape=1) +
  geom_boxplot(width = .2, outlier.shape = NA, lwd=0.5, fill=NA) +
  scale_fill_manual(values = c("#2a6aaf", "#d3292f"))+
  scale_color_manual(values = c("#2a6aaf", "#d3292f"))+
  scale_x_discrete(limits=c("Fishes","Amphibians", "Reptiles","Mammals","Birds"))+
  theme_classic()+
  labs(x = "", y = expression("dN (" * 10^-10 * " sub/site/year)"))+
  scale_y_continuous(limits = c(0,16), breaks = seq(0,16,4))+
  guides(x=guide_axis(cap='upper'), y=guide_axis(cap='upper'))+
  theme(axis.title = element_text(size=9),
        axis.text = element_text(size = 9),
        legend.position = 'none',
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))+
  geom_text(data=label, aes(x=x, y=y2, label=label.ds), size=2.5)+
  coord_flip()

f3 <-ggplot(data = molrate, aes(y = dNdS, x = Group, colour=ThermoMode, fill=ThermoMode)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), adjust =1) +
  geom_point(aes(y = dNdS), position = position_jitter(width = 0.12), alpha = 0.1, size = 0.5, shape=1) +
  geom_boxplot(width = .2, outlier.shape = NA, lwd=0.5, fill=NA) +
  scale_fill_manual(values = c("#2a6aaf", "#d3292f"))+
  scale_color_manual(values = c("#2a6aaf", "#d3292f"))+
  scale_x_discrete(limits=c("Fishes","Amphibians", "Reptiles","Mammals","Birds"))+
  theme_classic()+
  labs(x = "", y = expression("dN/dS"))+
  #scale_y_continuous(limits = c(0,12), breaks = seq(0,12,3))+
  guides(x=guide_axis(cap='upper'), y=guide_axis(cap='upper'))+
  theme(axis.title = element_text(size=9),
        axis.text = element_text(size = 9),
        legend.position = 'none',
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))+
  geom_text(data=label, aes(x=x, y=y3, label=label.dnds), size=2.5)+
  coord_flip()
cowplot::plot_grid(f1,f2, f3, nrow=1, align = "hv")
ggsave("./Outputs/Supplementary/Fig.S2.pdf", width = 8.3, height = 2.8)


######################################################
#Extended Fig.S5 Phylogenetic signal test of dS and dN for each class of vertebrates
#input phylogeny
phy <- read.tree(paste0(workdir, "/DataFiles/trees/phy_all_sampled_mtdna.tre"))

#force ultrametric phylogeny
phy <- force.ultrametric(phy, method="nnls")#nnls or extend

#Groups
groups <- unique(molrate$Group)

phylo.signal <- matrix(NA, 5, 6)
colnames(phylo.signal) <- c("ds.lambda", "ds.p", "dn.lambda", "dn.p", "dnds.lambda", "dnds.p")
row.names(phylo.signal) <- groups

for(i in 1:length(groups)){
  fgroup <- groups[i]
  
  #extract data
  fdat <- molrate%>%filter(Group==fgroup)
  
  row.names(fdat) <- fdat$Species
  
  #trait data with corresponding tree tips
  ds <- fdat$dS; names(ds)<-fdat$Species
  dn <- fdat$dN; names(dn)<-fdat$Species
  dnds <- fdat$dNdS; names(dnds)<-fdat$Species
  
  #phylogeny corresponding to trait data
  tree <- drop.tip(phy, name.check(phy, fdat)$tree_not_data)
  
  #compute phylogenetic signal using Pagel's lambda
  lambda.ds <- phylosig(tree, ds, method="lambda", test=TRUE)
  lambda.dn <- phylosig(tree, dn, method="lambda", test=TRUE)
  lambda.dnds <- phylosig(tree, dnds, method="lambda", test=TRUE)
  
  #matrix of lambda and p values
  phylo.signal[i,] <- c(lambda.ds$lambda, lambda.ds$P, lambda.dn$lambda, lambda.dn$P, lambda.dnds$lambda, lambda.dnds$P)
  print(fgroup)
}

#save data
#save(phylo.signal, file = "./Outputs/Data/phylogenetic_signal_mt.rdata")
load(file = "./Outputs/Data/phylogenetic_signal_mt.rdata")

phylo.signal <- as.data.frame(phylo.signal)
phylo.signal$Group <- row.names(phylo.signal)


phylo.signal <- phylo.signal%>%
  mutate(sig.dn=ifelse(dn.p<0.001,"***", ifelse(dn.p<0.01, "**", ifelse(dn.p<0.05, "*", "ns"))),
         sig.ds=ifelse(ds.p<0.001,"***", ifelse(ds.p<0.01, "**", ifelse(ds.p<0.05, "*", "ns"))),
         sig.dnds=ifelse(dnds.p<0.001,"***", ifelse(dnds.p<0.01, "**", ifelse(dnds.p<0.05, "*", "ns"))),
         Group=factor(Group, levels = c("Fishes","Amphibians", "Reptiles","Mammals","Birds")),
         ThermoMode=ifelse(Group %in% c("Birds","Mammals"), "Endotherms", "Ectotherms"))%>%
  rename(dS=ds.lambda, dN=dn.lambda, dNdS=dnds.lambda)%>%
  select(ThermoMode, Group, dS, dN, dNdS, sig.ds, sig.dn, sig.dnds)

#Extended Fig.S4 
fs5a<-phylo.signal%>%
  ggplot(aes(x = Group, y = dS, fill=ThermoMode)) +
  geom_bar(stat = "identity", position = "dodge")+
  scale_fill_manual(values=c("#2a6aaf", "#d3292f"))+
  geom_text(
    aes(label = sig.ds),
    colour = "black", size = 2.5,
    position = position_dodge(.9))+
  theme_classic()+
  scale_y_continuous(limits = c(0, 1),breaks = seq(0, 1, 0.25))+
  labs(x="", y=expression(paste("dS (" * lambda["Pagel"] *")")))+
  guides(y=guide_axis(cap='upper'))+
  theme(legend.title = element_blank(),
        legend.key.width = unit(0.3, "cm"),
        legend.key.height = unit(0.3, "cm"),
        legend.position = "none",
        legend.text = element_text(size=7),
        legend.background = element_blank(),
        axis.title = element_text(size=9),
        axis.text.y = element_text(size=9),
        axis.text.x = element_text(angle = 45, hjust = 1, size=9),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))


fs5b<-phylo.signal%>%
  ggplot(aes(x = Group, y = dN, fill=ThermoMode)) +
  geom_bar(stat = "identity", position = "dodge")+
  geom_text(
    aes(label = sig.dn),
    colour = "black", size = 2.5,
    position = position_dodge(.9))+
  scale_fill_manual(values=c("#2a6aaf", "#d3292f"))+
  scale_y_continuous(limits = c(0, 1),breaks = seq(0, 1, 0.25))+
  theme_classic()+
  labs(x="", y=expression(paste("dN (" * lambda["Pagel"] *")")))+
  guides(y=guide_axis(cap='upper'))+
  theme(legend.title = element_blank(),
        legend.key.width = unit(0.3, "cm"),
        legend.key.height = unit(0.3, "cm"),
        legend.position = "none",
        legend.text = element_text(size=7),
        legend.background = element_blank(),
        axis.title = element_text(size=9),
        axis.text.y = element_text(size=9),
        axis.text.x = element_text(angle = 45, hjust = 1, size=9),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs5c<-phylo.signal%>%
  ggplot(aes(x = Group, y = dNdS, fill=ThermoMode)) +
  geom_bar(stat = "identity", position = "dodge")+
  geom_text(
    aes(label = sig.dnds),
    colour = "black", size = 2.5,
    position = position_dodge(.9))+
  scale_fill_manual(values=c("#2a6aaf", "#d3292f"))+
  scale_y_continuous(limits = c(0, 1),breaks = seq(0, 1, 0.25))+
  theme_classic()+
  labs(x="", y=expression(paste("dN/dS (" * lambda["Pagel"] *")")))+
  guides(y=guide_axis(cap='upper'))+
  theme(legend.title = element_blank(),
        legend.key.width = unit(0.3, "cm"),
        legend.key.height = unit(0.3, "cm"),
        legend.position = "none",
        legend.text = element_text(size=7),
        legend.background = element_blank(),
        axis.title = element_text(size=9),
        axis.text.y = element_text(size=9),
        axis.text.x = element_text(angle = 45, hjust = 1, size=9),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))
plot_grid(fs5a, fs5b, fs5c, nrow = 1, align = "hv")
ggsave("./Outputs/Supplementary/Fig.S5.pdf", width = 8.3, height = 3.6)


##################################################################
#Part II: Analyzing how molecular rates vary with latitude and temperature at species and assemblage level
###################################################################
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
workdir <- "/Users/Tianlong/VertMolRate"
molrate <- read.csv(paste0(workdir, "/DataFiles/MolEvolRate/molrate_mtDNA.csv"))

# Load input phylogeny
phylo <- read.tree(paste0(workdir, "/DataFiles/trees/phy_all_sampled_mtdna.tre"))

# Check whether the phylogeny is ultrametric
is.ultrametric(phylo)

# Force the tree to be ultrametric
phylo <- force.ultrametric(phylo, method = "nnls")  # Use nnls or extend

#########################################################################################
#1. Analyzing how molecular rates varied at latitude and ambient temperature across species
#############################################################################################
# 1.1 PGLMMs account for phylogenetic relatedness as a random effect.

#groups
groups <- c("Ectotherms", "Endotherms", "Fishes", "Amphibians", "Reptiles", "Mammals","Birds", "Non-long migrants", "Non-migrants")

#Define two lists to store results of pglmms
pglmm_molrate_lat <- vector("list", length = length(groups) * 3)
names(pglmm_molrate_lat) <- paste(rep(groups, each = 3), c("dN", "dS", "dNdS"), sep = ".")

pglmm_molrate_temp <- vector("list", length = length(groups) * 3)
names(pglmm_molrate_temp) <- paste(rep(groups, each = 3), c("dN", "dS", "dNdS"), sep = ".")


# Prior specification for the model
prior <- list(R=list(V=diag(1),nu=0.002), 
              G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*1000)))

# Loop over each group
for (group in groups) {
  # Choose data based on group
  if(group %in% unique(molrate$ThermoMode)){
    subdata <- molrate %>% filter(ThermoMode == group)
  }
  if(group == "Non-long migrants"){
    subdata <- molrate %>% filter(Group == "Birds", Migration!="Long Migratory")
  }
  if(group == "Non-migrants"){
    subdata <- molrate %>% filter(Group == "Birds", Migration=="Resident")
  }
  if(group %in% unique(molrate$Group)){
    subdata <- molrate %>% filter(Group == group)
  }
  
  # Scale and fit dS model
  fit_lat <- MCMCglmm(scale(log(dS)) ~ scale(abs(Lat)), prior=prior,
                      random = ~Species, data = subdata, verbose = TRUE, 
                      ginverse = list(Species = inverseA(phylo)$Ainv), 
                      nitt = 60000, burnin = 10000, thin = 25, family = c("gaussian"))
  
  fit_temp <- MCMCglmm(scale(log(dS)) ~ scale(AnnualTemp), prior=prior,
                       random = ~Species, data = subdata, verbose = TRUE, 
                       ginverse = list(Species = inverseA(phylo)$Ainv), 
                       nitt = 60000, burnin = 10000, thin = 25, family = c("gaussian"))
  
  # Scale and fit dN model
  fit_dn_lat <- MCMCglmm(scale(log(dN)) ~ scale(abs(Lat)), prior=prior,
                         random = ~Species, data = subdata, verbose = TRUE, 
                         ginverse = list(Species = inverseA(phylo)$Ainv),
                         nitt = 60000, burnin = 10000, thin = 25,  family = c("gaussian"))
  
  fit_dn_temp <- MCMCglmm(scale(log(dN)) ~ scale(AnnualTemp), prior=prior,
                          random = ~Species, data = subdata, verbose = TRUE, 
                          ginverse = list(Species = inverseA(phylo)$Ainv), 
                          nitt = 60000, burnin = 10000, thin = 25,family = c("gaussian"))
  
  # Scale and fit dN/dS model
  fit_dnds_lat <- MCMCglmm(scale(log(dNdS)) ~ scale(abs(Lat)), prior=prior,
                           random = ~Species, data = subdata, verbose = TRUE, 
                           ginverse = list(Species = inverseA(phylo)$Ainv),
                           nitt = 60000, burnin = 10000, thin = 25, family = c("gaussian"))
  
  fit_dnds_temp <- MCMCglmm(scale(log(dNdS)) ~ scale(AnnualTemp), prior=prior,
                            random = ~Species, data = subdata, verbose = TRUE, 
                            ginverse = list(Species = inverseA(phylo)$Ainv),
                            nitt = 60000, burnin = 10000, thin = 25, family = c("gaussian"))
  
  # Store models in list
  pglmm_molrate_lat[[paste(group, "dS", sep = ".")]] <- fit_lat
  pglmm_molrate_lat[[paste(group, "dN", sep = ".")]] <- fit_dn_lat
  pglmm_molrate_lat[[paste(group, "dNdS", sep = ".")]] <- fit_dnds_lat
  
  pglmm_molrate_temp[[paste(group, "dS", sep = ".")]] <- fit_temp
  pglmm_molrate_temp[[paste(group, "dN", sep = ".")]] <- fit_dn_temp
  pglmm_molrate_temp[[paste(group, "dNdS", sep = ".")]] <- fit_dnds_temp
  
  print(group)
}

#save(pglmm_molrate_lat, file = "./Outputs/Data/pglmm_molrate_latitudinal_pattern_mtdna.rdata")
#save(pglmm_molrate_temp, file = "./Outputs/Data/pglmm_molrate_temperature_pattern_mtdna.rdata")
load(file = "./Outputs/Data/pglmm_molrate_latitudinal_pattern_mtdna.rdata")
load(file = "./Outputs/Data/pglmm_molrate_temperature_pattern_mtdna.rdata")


#a function to summary pglmms of dN, dS and dN/dS for all groups
summary_pglmm_groups <- function(groups, molrate, pglmms){
  fixed <- lapply(paste(groups, molrate, sep = "."), function (x) summary(pglmms[[x]])$solutions)
  random <- lapply(paste(groups, molrate, sep = "."), function (x) summary(pglmms[[x]])$Gcovariances)
  pglmm.out <- lapply(1:length(groups), 
                      function(x) cbind(data.frame(MolRate=molrate, Group=groups[x], Var=c(row.names(fixed[[x]]), row.names(random[[x]]))), 
                                        data.frame(rbind(fixed[[x]], cbind(random[[x]], pMCMC=NA)))))
  
  summary_out <- do.call(rbind, pglmm.out)
  return(summary_out)
}



# Extract parameter values from PGLMMs
summary_dn_lat <- summary_pglmm_groups(groups=groups, molrate="dN", pglmms=pglmm_molrate_lat)     
summary_ds_lat <- summary_pglmm_groups(groups=groups, molrate="dS", pglmms=pglmm_molrate_lat)
summary_dnds_lat <- summary_pglmm_groups(groups=groups, molrate="dNdS", pglmms=pglmm_molrate_lat)

summary_dn_temp <- summary_pglmm_groups(groups=groups, molrate="dN", pglmms=pglmm_molrate_temp)  
summary_ds_temp <- summary_pglmm_groups(groups=groups, molrate="dS", pglmms=pglmm_molrate_temp)
summary_dnds_temp <- summary_pglmm_groups(groups=groups, molrate="dNdS", pglmms=pglmm_molrate_temp)

# All parameter values
pglmm_out <- rbind(bind_rows("Latitude" = summary_ds_lat, "AnnualTemp" = summary_ds_temp, .id = "Predictors"),
                   bind_rows("Latitude" = summary_dn_lat, "AnnualTemp" = summary_dn_temp, .id = "Predictors"),
                   bind_rows("Latitude" = summary_dnds_lat, "AnnualTemp" = summary_dnds_temp, .id = "Predictors"))


#A function to extact MCMC slopes of PGLMMS
extarct_slope_pglmm <- function(groups, molrate, pglmms) {
  slopes <- lapply(1:length(groups), function (x) data.frame(Group=groups[x], Slope=pglmms[[paste(groups[x], molrate, sep = ".")]]$Sol[,2]))
  slopes <- do.call(rbind, slopes)
  colnames(slopes) <- c("Group", "Slope")
  return(slopes)
}

#Groups for plot
plot_groups <- c("Endotherms","Birds","Mammals","Ectotherms","Reptiles", "Amphibians","Fishes")

slope_lat <- extarct_slope_pglmm(groups=plot_groups, molrate="dS", pglmms=pglmm_molrate_lat)
slope_dn_lat <- extarct_slope_pglmm(groups=plot_groups, molrate="dN", pglmms=pglmm_molrate_lat)
slope_dnds_lat <- extarct_slope_pglmm(groups=plot_groups, molrate="dNdS", pglmms=pglmm_molrate_lat)

slope_temp <- extarct_slope_pglmm(groups=plot_groups, molrate="dS", pglmms=pglmm_molrate_temp)
slope_dn_temp <- extarct_slope_pglmm(groups=plot_groups, molrate="dN", pglmms=pglmm_molrate_temp)
slope_dnds_temp <- extarct_slope_pglmm(groups=plot_groups, molrate="dNdS", pglmms=pglmm_molrate_temp)



#Kruskal Wallis test and multiple comparison of groups
comp_lat <- agricolae::kruskal(slope_lat$Slope, slope_lat$Group, p.adj = "bonferroni")
comp_dn_lat <- agricolae::kruskal(slope_dn_lat$Slope, slope_lat$Group, p.adj = "bonferroni")
comp_dnds_lat <- agricolae::kruskal(slope_dnds_lat$Slope, slope_lat$Group, p.adj = "bonferroni")
comp_temp <- agricolae::kruskal(slope_temp$Slope, slope_temp$Group, p.adj = "bonferroni")
comp_dn_temp <- agricolae::kruskal(slope_dn_temp$Slope, slope_temp$Group, p.adj = "bonferroni")
comp_dnds_temp <- agricolae::kruskal(slope_dnds_temp$Slope, slope_temp$Group, p.adj = "bonferroni")



sig.label <- data.frame(Group=plot_groups, x=1:7, 
                        ds_lat=comp_lat$groups[plot_groups,"groups"],
                        dn_lat=comp_dn_lat$groups[plot_groups,"groups"],
                        dnds_lat=comp_dnds_lat$groups[plot_groups,"groups"],
                        ds_temp=comp_temp$groups[plot_groups,"groups"],
                        dn_temp=comp_dn_temp$groups[plot_groups,"groups"],
                        dnds_temp=comp_dnds_temp$groups[plot_groups,"groups"],
                        ThermoMode = c(rep("Ectotherms", 4), rep("Endotherms",3)))

#Plot effect size
f2a <-pglmm_out%>%
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
  labs(y="Effects of temperature", x="", title = "dS", tag="A")+
  theme_classic()+
  coord_flip()+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-0.2, 0.4),breaks = seq(-0.2,0.4,0.2), labels=c(-0.2,0,0.2,0.4))+
  theme(axis.text = element_text(size=7, color = "black"),
        axis.title = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = 'none',
        plot.tag = element_text(size=8, face = "bold"),
        plot.tag.position = c(0.05, 1),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))+
  geom_text(data=sig.label, aes(x=x, y=0.4, label=ds_temp), size=2.3, inherit.aes = FALSE)

f2b <-pglmm_out%>%
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
  labs(y="Effects of temperature", x="", title = "dN", tag="B")+
  theme_classic()+
  coord_flip()+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-0.2, 0.4),breaks = seq(-0.2,0.4,0.2), labels=seq(-0.2,0.4,0.2))+
  theme(axis.text = element_text(size=7,color = "black"),
        axis.title = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = 'none',
        plot.tag = element_text(size=8, face = "bold"),
        plot.tag.position = c(0.05, 1),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))+
  geom_text(data=sig.label, aes(x=x, y=0.4, label=dn_temp), size=2.3, inherit.aes = FALSE)



f2c <-pglmm_out%>%
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
  scale_y_continuous(limits = c(-0.3, 0.3),breaks = seq(-0.3,0.3,0.15), labels=seq(-0.3,0.3,0.15))+
  theme(axis.text = element_text(size=7, color = "black"),
        axis.title = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = 'none',
        plot.tag = element_text(size=8, face = "bold"),
        plot.tag.position = c(0.05, 1),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))+
  geom_text(data=sig.label, aes(x=x, y=0.3, label=dnds_temp), size=2.3, inherit.aes = FALSE)


cowplot::plot_grid(f2a,f2b, f2c, nrow = 3, align="hv")

ggsave(filename="./Outputs/MainFigures/Fig2abc.pdf", height=4.3, width=1.8)

#Plot effect size to show how molecular rates varied with latitude
#Extended Fig.S6
fs6a <-pglmm_out%>%
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
  scale_y_continuous(limits = c(-0.3, 0.15),breaks = seq(-0.3,0.15,0.15), labels=c(-0.3,-0.15,0,0.15))+
  theme(axis.text = element_text(size=7),
        axis.title = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = 'none',
        plot.tag = element_text(size=8, face = "bold"),
        plot.tag.position = c(0.05, 1),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))+
  geom_text(data=sig.label, aes(x=x, y=0.15, label=ds_lat), size=2.5, inherit.aes = FALSE)

fs6b <-pglmm_out%>%
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
  scale_y_continuous(limits = c(-0.3, 0.15),breaks = seq(-0.3,0.15,0.15), labels=c(-0.3,-0.15,0,0.15))+
  theme(axis.text = element_text(size=7),
        axis.title = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = 'none',
        plot.tag = element_text(size=8, face = "bold"),
        plot.tag.position = c(0.05, 1),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))+
  geom_text(data=sig.label, aes(x=x, y=0.15, label=dn_lat), size=2.5, inherit.aes = FALSE)



fs6c <-pglmm_out%>%
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
  scale_y_continuous(limits = c(-0.3, 0.15),breaks = seq(-0.3,0.15,0.15), labels=c(-0.3,-0.15,0,0.15))+
  theme(axis.text = element_text(size=7),
        axis.title = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = 'none',
        plot.tag = element_text(size=8, face = "bold"),
        plot.tag.position = c(0.05, 1),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))+
  geom_text(data=sig.label, aes(x=x, y=0.15, label=dnds_lat), size=2.5, inherit.aes = FALSE)


cowplot::plot_grid(fs6a,fs6b, fs6c, nrow = 1, align="hv")
ggsave(filename="./Outputs/Supplementary/Fig.S6.pdf", height=2.3, width=8.27)

####################################################################################
# 1.2 PGLMMs using data remove relative molecular rates less than 0.4 sub/site to reduce the impacts of substitution saturation
#groups
groups <- c("Ectotherms", "Endotherms", "Fishes", "Amphibians", "Reptiles", "Mammals","Birds", "Non-long migrants", "Non-migrants")

#Define two lists to store results of pglmms
pglmm_molrate_lat_0.4 <- vector("list", length = length(groups) * 3)
names(pglmm_molrate_lat_0.4) <- paste(rep(groups, each = 3), c("dN", "dS", "dNdS"), sep = ".")

pglmm_molrate_temp_0.4 <- vector("list", length = length(groups) * 3)
names(pglmm_molrate_temp_0.4) <- paste(rep(groups, each = 3), c("dN", "dS", "dNdS"), sep = ".")


# Prior specification for the model
prior <- list(R=list(V=diag(1),nu=0.002), 
              G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*1000)))

# Loop over each group
for (group in groups) {
  # Choose data based on group
  if(group %in% unique(molrate$ThermoMode)){
    subdata <- molrate %>% filter(dS*Tip.Age*1000000<=0.4) %>% filter(ThermoMode == group)
  }
  if(group == "Non-long migrants"){
    subdata <- molrate %>% filter(dS*Tip.Age*1000000<=0.4) %>% filter(Group == "Birds", Migration!="Long Migratory")
  }
  if(group == "Non-migrants"){
    subdata <- molrate %>% filter(dS*Tip.Age*1000000<=0.4) %>% filter(Group == "Birds", Migration=="Resident")
  }
  if(group %in% unique(molrate$Group)){
    subdata <- molrate %>% filter(dS*Tip.Age*1000000<=0.4) %>% filter(Group == group)
  }
  
  # Scale and fit dS model
  fit_lat <- MCMCglmm(scale(log(dS)) ~ scale(abs(Lat)), prior=prior,
                      random = ~Species, data = subdata, verbose = TRUE, 
                      ginverse = list(Species = inverseA(phylo)$Ainv), 
                      nitt = 60000, burnin = 10000, thin = 25, family = c("gaussian"))
  
  fit_temp <- MCMCglmm(scale(log(dS)) ~ scale(AnnualTemp), prior=prior,
                       random = ~Species, data = subdata, verbose = TRUE, 
                       ginverse = list(Species = inverseA(phylo)$Ainv), 
                       nitt = 60000, burnin = 10000, thin = 25, family = c("gaussian"))
  
  # Scale and fit dN model
  fit_dn_lat <- MCMCglmm(scale(log(dN)) ~ scale(abs(Lat)), prior=prior,
                         random = ~Species, data = subdata, verbose = TRUE, 
                         ginverse = list(Species = inverseA(phylo)$Ainv),
                         nitt = 60000, burnin = 10000, thin = 25,  family = c("gaussian"))
  
  fit_dn_temp <- MCMCglmm(scale(log(dN)) ~ scale(AnnualTemp), prior=prior,
                          random = ~Species, data = subdata, verbose = TRUE, 
                          ginverse = list(Species = inverseA(phylo)$Ainv), 
                          nitt = 60000, burnin = 10000, thin = 25,family = c("gaussian"))
  
  # Scale and fit dN/dS model
  fit_dnds_lat <- MCMCglmm(scale(log(dNdS)) ~ scale(abs(Lat)), prior=prior,
                           random = ~Species, data = subdata, verbose = TRUE, 
                           ginverse = list(Species = inverseA(phylo)$Ainv),
                           nitt = 60000, burnin = 10000, thin = 25, family = c("gaussian"))
  
  fit_dnds_temp <- MCMCglmm(scale(log(dNdS)) ~ scale(AnnualTemp), prior=prior,
                            random = ~Species, data = subdata, verbose = TRUE, 
                            ginverse = list(Species = inverseA(phylo)$Ainv),
                            nitt = 60000, burnin = 10000, thin = 25, family = c("gaussian"))
  
  # Store models in list
  pglmm_molrate_lat_0.4[[paste(group, "dS", sep = ".")]] <- fit_lat
  pglmm_molrate_lat_0.4[[paste(group, "dN", sep = ".")]] <- fit_dn_lat
  pglmm_molrate_lat_0.4[[paste(group, "dNdS", sep = ".")]] <- fit_dnds_lat
  
  pglmm_molrate_temp_0.4[[paste(group, "dS", sep = ".")]] <- fit_temp
  pglmm_molrate_temp_0.4[[paste(group, "dN", sep = ".")]] <- fit_dn_temp
  pglmm_molrate_temp_0.4[[paste(group, "dNdS", sep = ".")]] <- fit_dnds_temp
  
  print(group)
}

#save(pglmm_molrate_lat_0.4, file = "./Outputs/Data/pglmm_molrate_latitudinal_pattern_mtdna_0.4.rdata")
#save(pglmm_molrate_temp_0.4, file = "./Outputs/Data/pglmm_molrate_temperature_pattern_mtdna_0.4.rdata")
load(file = "./Outputs/Data/pglmm_molrate_latitudinal_pattern_mtdna_0.4.rdata")
load(file = "./Outputs/Data/pglmm_molrate_temperature_pattern_mtdna_0.4.rdata")


# Extract parameter values from PGLMMs
summary_dn_lat_0.4 <- summary_pglmm_groups(groups=groups, molrate="dN", pglmms=pglmm_molrate_lat_0.4)     
summary_lat_0.4 <- summary_pglmm_groups(groups=groups, molrate="dS", pglmms=pglmm_molrate_lat_0.4)
summary_dnds_lat_0.4 <- summary_pglmm_groups(groups=groups, molrate="dNdS", pglmms=pglmm_molrate_lat_0.4)

summary_dn_temp_0.4 <- summary_pglmm_groups(groups=groups, molrate="dN", pglmms=pglmm_molrate_temp_0.4)  
summary_temp_0.4 <- summary_pglmm_groups(groups=groups, molrate="dS", pglmms=pglmm_molrate_temp_0.4)
summary_dnds_temp_0.4 <- summary_pglmm_groups(groups=groups, molrate="dNdS", pglmms=pglmm_molrate_temp_0.4)

# All parameter values
pglmm_out_0.4 <- rbind(bind_rows("Latitude" = summary_lat_0.4, "AnnualTemp" = summary_temp_0.4, .id = "Predictors"),
                       bind_rows("Latitude" = summary_dn_lat_0.4, "AnnualTemp" = summary_dn_temp_0.4, .id = "Predictors"),
                       bind_rows("Latitude" = summary_dnds_lat_0.4, "AnnualTemp" = summary_dnds_temp_0.4, .id = "Predictors"))

#write.csv(pglmm_out_0.4, "pglmm_out.csv")


################################################################################
# 1.3 PGLMMs account for phylogenetic relatedness under the best traits evolution models.
##########
#Fit the best models
groups <- c("Ectotherms", "Endotherms", "Fishes", "Amphibians", "Reptiles", "Mammals","Birds", "Non-long migrants", "Non-migrants")

# Creating a list to store PGLMM objects
best_trait_evol_models <- vector("list", length = length(groups) * 3)
names(best_trait_evol_models) <- paste(rep(groups, each = 3), c("dN", "dS", "dNdS"), sep = ".")


for (group in groups) {
  # Choose data based on group
  if(group %in% unique(molrate$ThermoMode)){
    subdata <- molrate %>% filter(ThermoMode == group)%>%mutate(dS=as.numeric(scale(log(dS))), dN=as.numeric(scale(log(dN))))
  }
  if(group == "Non-long migrants"){
    subdata <- molrate %>% filter(Group == "Birds", Migration!="Long Migratory")%>%mutate(dS=as.numeric(scale(log(dS))), dN=as.numeric(scale(log(dN))))
  }
  if(group == "Non-migrants"){
    subdata <- molrate %>% filter(Group == "Birds", Migration=="Resident")%>%mutate(dS=as.numeric(scale(log(dS))), dN=as.numeric(scale(log(dN))))
  }
  if(group %in% unique(molrate$Group)){
    subdata <- molrate %>% filter(Group == group)%>%mutate(dS=as.numeric(scale(log(dS))), dN=as.numeric(scale(log(dN))))
  }
  
  #phylogeny
  phy <- drop.tip(phylo, setdiff(phylo$tip.label, subdata$Species))
  
  #Model selection
  best_trait_evol_models[[paste(group, "dS", sep = ".")]] <- fitTraitEvolModel(trait="dS", data=subdata, phy = phy)
  best_trait_evol_models[[paste(group, "dN", sep = ".")]] <- fitTraitEvolModel(trait="dN", data=subdata, phy = phy)
  best_trait_evol_models[[paste(group, "dNdS", sep = ".")]] <- fitTraitEvolModel(trait="dNdS", data=subdata, phy = phy)
}

#save(best_trait_evol_models, file="./Outputs/Data/best_trait_evol_models_mtdna.rdata")
load(file="./Outputs/Data/best_trait_evol_models_mtdna.rdata")


best_trait_evol_models%>%bind_rows()%>%
  mutate(Group=rep(groups, each=3), molrate=rep(c("dN","dS", "dNdS"), length(groups)))

##############################
#Fit the PGLMM using the best model of traits evolution model
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
  if(group == "Non-migrants"){
    subdata <- molrate %>% filter(Group == "Birds", Migration=="Resident")
  }
  if(group == "Non-long migrants"){
    subdata <- molrate %>% filter(Group == "Birds", Migration!="Long Migratory")
  }
  if(group %in% unique(molrate$Group)){
    subdata <- molrate %>% filter(Group == group)
  }
  
  phy <- drop.tip(phylo, setdiff(phylo$tip.label, subdata$Species))
  
  #The best traits evolution model of dS
  model.ds <- best_trait_evol_models[[paste(group, "dS", sep = ".")]]
  
  #The best traits evolution model of dN
  model.dn <- best_trait_evol_models[[paste(group, "dN", sep = ".")]]
  
  #The best traits evolution model of dNdS
  model.dnds <- best_trait_evol_models[[paste(group, "dNdS", sep = ".")]]
  
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

#save(pglmm_molrate_lat_best_model, file = "./Outputs/Data/pglmm_molrate_latitudinal_pattern_best_model_mtdna.rdata")
#save(pglmm_molrate_temp_best_model, file = "./Outputs/Data/pglmm_molrate_temperature_pattern_best_model_mtdna.rdata")

load(file = "./Outputs/Data/pglmm_molrate_latitudinal_pattern_best_model_mtdna.rdata")
load(file = "./Outputs/Data/pglmm_molrate_temperature_pattern_best_model_mtdna.rdata")

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
#1.4. PGLMMs account for phylogenetic relatedness, habitat and lineages as random effects.
groups <- c(unique(molrate$Group), unique(molrate$ThermoMode), "Non-long migrants", "Non-migrants")

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
  if(group == "Non-migrants"){
    subdata <- molrate %>% filter(Group == "Birds", Migration=="Resident")
  }
  if(group == "Non-long migrants"){
    subdata <- molrate %>% filter(Group == "Birds", Migration!="Long Migratory")
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


#save(pglmm_molrate_lat_3randoms, file = "./Outputs/Data/pglmm_molrate_latitudinal_pattern_3randoms_mtdna.rdata")
#save(pglmm_molrate_temp_3randoms, file = "./Outputs/Data/pglmm_molrate_temperature_pattern_3randoms_mtdna.rdata")

load(file = "./Outputs/Data/pglmm_molrate_latitudinal_pattern_3randoms_mtdna.rdata")
load(file = "./Outputs/Data/pglmm_molrate_temperature_pattern_3randoms_mtdna.rdata")

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


#############
#2. Plots to show relationships between molecular rates and absolute midpoint latitude/annual temperature at the species level
#Extended Fig.S3 Plots for endotherms and ectotherms
fs3a <- molrate%>%
  ggplot(aes(x=abs(Lat), y=log10(dS), colour=ThermoMode))+
  geom_point(size=0.5, alpha=0.5)+
  geom_smooth(method="lm", se=TRUE, size=0.5)+
  scale_colour_manual(values = c("#2a6aaf", "#d3292f"))+
  labs(x="Absolute Latitude", y=expression(log[10](dS)~(substitution/site/year)))+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-9.1, -7),breaks = seq(-9.0, -7, 0.5))+
  scale_x_continuous(limits = c(0, 90),breaks = seq(0,90,30))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=9),
        legend.position = "none",
        axis.title = element_text(size=9),
        axis.text = element_text(size=8, color = "black"),
        strip.text = element_text(size=9),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs3b <- molrate%>%
  ggplot(aes(x=abs(Lat), y=log10(dN), colour=ThermoMode))+
  geom_point(size=0.5, alpha=0.5)+
  geom_smooth(method="lm", se=TRUE, size=0.5)+
  scale_colour_manual(values = c("#2a6aaf", "#d3292f"))+
  labs(x="Absolute Latitude", y=expression(log[10](dN)~(substitution/site/year)))+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-10.4, -8.8),breaks = seq(-10.4, -8.8, 0.4))+
  scale_x_continuous(limits = c(0, 80),breaks = seq(0,80,20))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=9),
        axis.title = element_text(size=9),
        axis.text = element_text(size=8, color = "black"),
        legend.position = "none",
        legend.key.size = unit(0.3, "cm"),
        legend.title = element_blank(),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs3c <- molrate%>%
  ggplot(aes(x=abs(Lat), y=log10(dNdS), colour=ThermoMode))+
  geom_point(size=0.5, alpha=0.5)+
  geom_smooth(method="lm", se=TRUE, size=0.5)+
  scale_colour_manual(values = c("#2a6aaf", "#d3292f"))+
  labs(x="Absolute Latitude", y=expression(log[10](dN/dS)))+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-2.5, -0.5),breaks = seq(-2.5, -0.5, 0.5))+
  scale_x_continuous(limits = c(0, 80),breaks = seq(0,80,20))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=9),
        axis.title = element_text(size=9),
        axis.text = element_text(size=8, color = "black"),
        legend.position = c(0.8,0.9),
        legend.key.size = unit(0.3, "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.background = element_blank(),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs3d <- molrate%>%
  ggplot(aes(x=AnnualTemp, y=log10(dS), colour=ThermoMode))+
  geom_point(size=0.5, alpha=0.5)+
  geom_smooth(method="lm", se=TRUE, size=0.5)+
  scale_colour_manual(values = c("#2a6aaf", "#d3292f"))+
  labs(x=expression("Mean Annual Temperature"~(degree * C)), y=expression(log[10](dS)~(substitution/site/year)))+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-9.1, -7),breaks = seq(-9.0, -7, 0.5))+
  scale_x_continuous(limits = c(-15, 30),breaks = seq(-15, 30,5))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=9),
        legend.position = "none",
        axis.title = element_text(size=9),
        axis.text = element_text(size=8, color = "black"),
        strip.text = element_text(size=9),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs3e <- molrate%>%
  ggplot(aes(x=AnnualTemp, y=log10(dN), colour=ThermoMode))+
  geom_point(size=0.5, alpha=0.5)+
  geom_smooth(method="lm", se=TRUE, size=0.5)+
  scale_colour_manual(values = c("#2a6aaf", "#d3292f"))+
  labs(x=expression("Mean Annual Temperature"~(degree * C)), y=expression(log[10](dN)~(substitution/site/year)))+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-10.4, -8.8),breaks = seq(-10.4, -8.8, 0.4))+
  scale_x_continuous(limits = c(-15, 30),breaks = seq(-15, 30,5))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=9),
        axis.title = element_text(size=9),
        axis.text = element_text(size=8, color = "black"),
        legend.position = "none",
        legend.key.size = unit(0.3, "cm"),
        legend.title = element_blank(),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs3f <- molrate%>%
  ggplot(aes(x=AnnualTemp, y=log10(dNdS), colour=ThermoMode))+
  geom_point(size=0.5, alpha=0.5)+
  geom_smooth(method="lm", se=TRUE, size=0.5)+
  scale_colour_manual(values = c("#2a6aaf", "#d3292f"))+
  labs(x=expression("Mean Annual Temperature"~(degree * C)), y=expression(log[10](dN/dS)))+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-2.5, -0.5),breaks = seq(-2.5, -0.5, 0.5))+
  scale_x_continuous(limits = c(-15, 30),breaks = seq(-15, 30,5))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=9),
        axis.title = element_text(size=9),
        axis.text = element_text(size=8, color = "black"),
        legend.position = c(0.8,0.9),
        legend.key.size = unit(0.3, "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.background = element_blank(),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))
cowplot::plot_grid(fs3a,fs3b, fs3c, fs3d,fs3e, fs3f, nrow=2, align="hv")
ggsave(filename = "./Outputs/Supplementary/Fig.S3.pdf", width=8.27, height = 5)


#Extended Fig.S4 Plots for each class of vertebrates
fs41 <- molrate%>%
  filter(Group=="Fishes")%>%
  ggplot(aes(x=abs(Lat), y=log10(dS)))+
  geom_point(size=0.5, alpha=0.5, colour="#2a6aaf")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#2a6aaf")+
  labs(x="", y=expression(log[10](dS)~(substitution/site/year)), title="Fishes")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-9.1, -7),breaks = seq(-9.0, -7, 0.5))+
  scale_x_continuous(limits = c(0, 80),breaks = seq(0,80,20))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=8),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        strip.text = element_text(size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs42 <- molrate%>%
  filter(Group=="Amphibians")%>%
  ggplot(aes(x=abs(Lat), y=log10(dS)))+
  geom_point(size=0.5, alpha=0.5, colour="#2a6aaf")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#2a6aaf")+
  labs(x="", y="", title="Amphibians")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-9, -7.75),breaks = seq(-9, -7.75, 0.25))+
  #scale_x_continuous(limits = c(0, 90),breaks = seq(0,90,30))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=8),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        strip.text = element_text(size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))


fs43 <- molrate%>%
  filter(Group=="Reptiles")%>%
  ggplot(aes(x=abs(Lat), y=log10(dS)))+
  geom_point(size=0.5, alpha=0.5, colour="#2a6aaf")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#2a6aaf")+
  labs(x="", y="", title="Reptiles")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-9, -7.5),breaks = seq(-9, -7.5, 0.25))+
  scale_x_continuous(limits = c(0, 60),breaks = seq(0,60,20))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=8),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        strip.text = element_text(size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs44 <- molrate%>%
  filter(Group=="Mammals")%>%
  ggplot(aes(x=abs(Lat), y=log10(dS)))+
  geom_point(size=0.5, alpha=0.5, colour="#d3292f")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#d3292f")+
  labs(x="", y="", title="Mammals")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-8.6, -7.2),breaks = seq(-8.6, -7.2, 0.2))+
  scale_x_continuous(limits = c(0, 80),breaks = seq(0,80,20))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=8),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        strip.text = element_text(size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs45 <- molrate%>%
  filter(Group=="Birds")%>%
  ggplot(aes(x=abs(Lat), y=log10(dS)))+
  geom_point(size=0.5, alpha=0.5, colour="#d3292f")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#d3292f")+
  labs(x="", y="", title="Birds")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-8.5, -7.5),breaks = seq(-8.5, -7.5, 0.25))+
  scale_x_continuous(limits = c(0, 80),breaks = seq(0,80,20))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=8),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        strip.text = element_text(size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs46 <- molrate%>%
  filter(Group=="Fishes")%>%
  ggplot(aes(x=abs(Lat), y=log10(dN)))+
  geom_point(size=0.5, alpha=0.5, colour="#2a6aaf")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#2a6aaf")+
  labs(x="", y=expression(log[10](dN)~(substitution/site/year)))+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-10.5, -9),breaks = seq(-10.5, -9, 0.5))+
  #scale_x_continuous(limits = c(0, 90),breaks = seq(0,90,30))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=8),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        strip.text = element_text(size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs47 <- molrate%>%
  filter(Group=="Amphibians")%>%
  ggplot(aes(x=abs(Lat), y=log10(dN)))+
  geom_point(size=0.5, alpha=0.5, colour="#2a6aaf")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#2a6aaf")+
  labs(x="", y="")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-10.7, -9),breaks = seq(-10.5, -9, 0.5))+
  #scale_x_continuous(limits = c(0, 90),breaks = seq(0,90,30))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=8),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        strip.text = element_text(size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))


fs48 <- molrate%>%
  filter(Group=="Reptiles")%>%
  ggplot(aes(x=abs(Lat), y=log10(dN)))+
  geom_point(size=0.5, alpha=0.5, colour="#2a6aaf")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#2a6aaf")+
  labs(x="", y="")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-10, -9.2),breaks = seq(-10, -9.2, 0.2))+
  scale_x_continuous(limits = c(0, 60),breaks = seq(0,60,20))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=8),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        strip.text = element_text(size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs49 <- molrate%>%
  filter(Group=="Mammals")%>%
  ggplot(aes(x=abs(Lat), y=log10(dN)))+
  geom_point(size=0.5, alpha=0.5, colour="#d3292f")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#d3292f")+
  labs(x="", y="")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-9.8, -8.8),breaks = seq(-9.8, -8.8, 0.2))+
  #scale_x_continuous(limits = c(0, 60),breaks = seq(0,60,20))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=8),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        strip.text = element_text(size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs410 <- molrate%>%
  filter(Group=="Birds")%>%
  ggplot(aes(x=abs(Lat), y=log10(dN)))+
  geom_point(size=0.5, alpha=0.5, colour="#d3292f")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#d3292f")+
  labs(x="", y="")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-10, -9),breaks = seq(-10, -9, 0.2))+
  #scale_x_continuous(limits = c(0, 60),breaks = seq(0,60,20))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=8),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        strip.text = element_text(size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs411 <- molrate%>%
  filter(Group=="Fishes")%>%
  ggplot(aes(x=abs(Lat), y=log10(dNdS)))+
  geom_point(size=0.5, alpha=0.5, colour="#2a6aaf")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#2a6aaf")+
  labs(x="Absolute Latitude", y=expression(log[10](dN/dS)))+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-2.5, -0.5),breaks = seq(-2.5, -0.5, 0.5))+
  #scale_x_continuous(limits = c(0, 90),breaks = seq(0,90,30))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=8),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        strip.text = element_text(size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs412 <- molrate%>%
  filter(Group=="Amphibians")%>%
  ggplot(aes(x=abs(Lat), y=log10(dNdS)))+
  geom_point(size=0.5, alpha=0.5, colour="#2a6aaf")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#2a6aaf")+
  labs(x="Absolute Latitude", y="")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-2.1, -0.9),breaks = seq(-2.1, -0.9, 0.3))+
  #scale_x_continuous(limits = c(0, 90),breaks = seq(0,90,30))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=8),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        strip.text = element_text(size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))


fs413 <- molrate%>%
  filter(Group=="Reptiles")%>%
  ggplot(aes(x=abs(Lat), y=log10(dNdS)))+
  geom_point(size=0.5, alpha=0.5, colour="#2a6aaf")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#2a6aaf")+
  labs(x="Absolute Latitude", y="")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-2, -0.5),breaks = seq(-2, -0.5, 0.5))+
  scale_x_continuous(limits = c(0, 60),breaks = seq(0,60,20))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=8),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        strip.text = element_text(size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs414 <- molrate%>%
  filter(Group=="Mammals")%>%
  ggplot(aes(x=abs(Lat), y=log10(dNdS)))+
  geom_point(size=0.5, alpha=0.5, colour="#d3292f")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#d3292f")+
  labs(x="Absolute Latitude", y="")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  #scale_y_continuous(limits = c(-9.8, -8.8),breaks = seq(-9.8, -8.8, 0.2))+
  #scale_x_continuous(limits = c(0, 60),breaks = seq(0,60,20))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=8),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        strip.text = element_text(size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs415 <- molrate%>%
  filter(Group=="Birds")%>%
  ggplot(aes(x=abs(Lat), y=log10(dNdS)))+
  geom_point(size=0.5, alpha=0.5, colour="#d3292f")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#d3292f")+
  labs(x="Absolute Latitude", y="")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-2.2, -0.8),breaks = seq(-2.2, -0.8, 0.2))+
  #scale_x_continuous(limits = c(0, 60),breaks = seq(0,60,20))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=8),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        strip.text = element_text(size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs416 <- molrate%>%
  filter(Group=="Fishes")%>%
  ggplot(aes(x=AnnualTemp, y=log10(dS)))+
  geom_point(size=0.5, alpha=0.5, colour="#2a6aaf")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#2a6aaf")+
  labs(x="", y=expression(log[10](dS)~(substitution/site/year)), title="Fishes")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-9.1, -7),breaks = seq(-9.0, -7, 0.5))+
  scale_x_continuous(limits = c(-10, 30),breaks = seq(-10, 30,10))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=8),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        strip.text = element_text(size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs417 <- molrate%>%
  filter(Group=="Amphibians")%>%
  ggplot(aes(x=AnnualTemp, y=log10(dS)))+
  geom_point(size=0.5, alpha=0.5, colour="#2a6aaf")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#2a6aaf")+
  labs(x="", y="", title="Amphibians")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-9, -7.75),breaks = seq(-9, -7.75, 0.25))+
  scale_x_continuous(limits = c(-10, 30),breaks = seq(-10, 30,10))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=8),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        strip.text = element_text(size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))


fs418 <- molrate%>%
  filter(Group=="Reptiles")%>%
  ggplot(aes(x=AnnualTemp, y=log10(dS)))+
  geom_point(size=0.5, alpha=0.5, colour="#2a6aaf")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#2a6aaf")+
  labs(x="", y="", title="Reptiles")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-9, -7.5),breaks = seq(-9, -7.5, 0.25))+
  scale_x_continuous(limits = c(-10, 30),breaks = seq(-10, 30,10))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=8),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        strip.text = element_text(size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs419 <- molrate%>%
  filter(Group=="Mammals")%>%
  ggplot(aes(x=AnnualTemp, y=log10(dS)))+
  geom_point(size=0.5, alpha=0.5, colour="#d3292f")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#d3292f")+
  labs(x="", y="", title="Mammals")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-8.6, -7.2),breaks = seq(-8.6, -7.2, 0.2))+
  scale_x_continuous(limits = c(-10, 30),breaks = seq(-10, 30,10))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=8),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        strip.text = element_text(size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs420 <- molrate%>%
  filter(Group=="Birds")%>%
  ggplot(aes(x=AnnualTemp, y=log10(dS)))+
  geom_point(size=0.5, alpha=0.5, colour="#d3292f")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#d3292f")+
  labs(x="", y="", title="Birds")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-8.5, -7.5),breaks = seq(-8.5, -7.5, 0.25))+
  scale_x_continuous(limits = c(-10, 30),breaks = seq(-10, 30,10))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=8),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        strip.text = element_text(size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs421 <- molrate%>%
  filter(Group=="Fishes")%>%
  ggplot(aes(x=AnnualTemp, y=log10(dN)))+
  geom_point(size=0.5, alpha=0.5, colour="#2a6aaf")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#2a6aaf")+
  labs(x="", y=expression(log[10](dN)~(substitution/site/year)))+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-10.5, -9),breaks = seq(-10.5, -9, 0.5))+
  scale_x_continuous(limits = c(-10, 30),breaks = seq(-10, 30,10))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=8),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        strip.text = element_text(size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs422 <- molrate%>%
  filter(Group=="Amphibians")%>%
  ggplot(aes(x=AnnualTemp, y=log10(dN)))+
  geom_point(size=0.5, alpha=0.5, colour="#2a6aaf")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#2a6aaf")+
  labs(x="", y="")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-10.7, -9),breaks = seq(-10.5, -9, 0.5))+
  scale_x_continuous(limits = c(-10, 30),breaks = seq(-10, 30,10))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=8),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        strip.text = element_text(size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))


fs423 <- molrate%>%
  filter(Group=="Reptiles")%>%
  ggplot(aes(x=AnnualTemp, y=log10(dN)))+
  geom_point(size=0.5, alpha=0.5, colour="#2a6aaf")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#2a6aaf")+
  labs(x="", y="")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-10, -9.2),breaks = seq(-10, -9.2, 0.2))+
  scale_x_continuous(limits = c(-10, 30),breaks = seq(-10, 30,10))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=8),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        strip.text = element_text(size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs424 <- molrate%>%
  filter(Group=="Mammals")%>%
  ggplot(aes(x=AnnualTemp, y=log10(dN)))+
  geom_point(size=0.5, alpha=0.5, colour="#d3292f")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#d3292f")+
  labs(x="", y="")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-9.8, -8.8),breaks = seq(-9.8, -8.8, 0.2))+
  scale_x_continuous(limits = c(-10, 30),breaks = seq(-10, 30,10))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=8),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        strip.text = element_text(size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs425 <- molrate%>%
  filter(Group=="Birds")%>%
  ggplot(aes(x=AnnualTemp, y=log10(dN)))+
  geom_point(size=0.5, alpha=0.5, colour="#d3292f")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#d3292f")+
  labs(x="", y="")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-10, -9),breaks = seq(-10, -9, 0.2))+
  scale_x_continuous(limits = c(-10, 30),breaks = seq(-10, 30,10))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=8),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        strip.text = element_text(size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs426 <- molrate%>%
  filter(Group=="Fishes")%>%
  ggplot(aes(x=AnnualTemp, y=log10(dNdS)))+
  geom_point(size=0.5, alpha=0.5, colour="#2a6aaf")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#2a6aaf")+
  labs(x=expression("Mean Annual temperature"~(degree * C)), y=expression(log[10](dN/dS)))+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-2.5, -0.5),breaks = seq(-2.5, -0.5, 0.5))+
  scale_x_continuous(limits = c(-10, 30),breaks = seq(-10, 30,10))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=8),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        strip.text = element_text(size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs427 <- molrate%>%
  filter(Group=="Amphibians")%>%
  ggplot(aes(x=AnnualTemp, y=log10(dNdS)))+
  geom_point(size=0.5, alpha=0.5, colour="#2a6aaf")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#2a6aaf")+
  labs(x=expression("Mean Annual temperature"~(degree * C)), y="")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-2.1, -0.9),breaks = seq(-2.1, -0.9, 0.3))+
  scale_x_continuous(limits = c(-10, 30),breaks = seq(-10, 30,10))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=8),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        strip.text = element_text(size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))


fs428 <- molrate%>%
  filter(Group=="Reptiles")%>%
  ggplot(aes(x=AnnualTemp, y=log10(dNdS)))+
  geom_point(size=0.5, alpha=0.5, colour="#2a6aaf")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#2a6aaf")+
  labs(x=expression("Mean Annual temperature"~(degree * C)), y="")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-2, -0.5),breaks = seq(-2, -0.5, 0.5))+
  scale_x_continuous(limits = c(-10, 30),breaks = seq(-10, 30,10))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=8),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        strip.text = element_text(size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs429 <- molrate%>%
  filter(Group=="Mammals")%>%
  ggplot(aes(x=AnnualTemp, y=log10(dNdS)))+
  geom_point(size=0.5, alpha=0.5, colour="#d3292f")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#d3292f")+
  labs(x=expression("Mean Annual temperature"~(degree * C)), y="")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  #scale_y_continuous(limits = c(-9.8, -8.8),breaks = seq(-9.8, -8.8, 0.2))+
  scale_x_continuous(limits = c(-10, 30),breaks = seq(-10, 30,10))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=8),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        strip.text = element_text(size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs430 <- molrate%>%
  filter(Group=="Birds")%>%
  ggplot(aes(x=AnnualTemp, y=log10(dNdS)))+
  geom_point(size=0.5, alpha=0.5, colour="#d3292f")+
  geom_smooth(method="lm", se=TRUE, size=0.5, colour="#d3292f")+
  labs(x=expression("Mean Annual temperature"~(degree * C)), y="")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-2.2, -0.8),breaks = seq(-2.2, -0.8, 0.2))+
  scale_x_continuous(limits = c(-10, 30),breaks = seq(-10, 30,10))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=8),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        strip.text = element_text(size=8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

cowplot::plot_grid(fs41,fs42,fs43,fs44,fs45,
                   fs46,fs47,fs48,fs49,fs410,
                   fs411,fs412,fs413,fs414,fs415,
                   fs416,fs417,fs418,fs419,fs420,
                   fs421,fs422,fs423,fs424,fs425,
                   fs426,fs427,fs428,fs429,fs430,nrow = 6, align = "hv")
ggsave(filename = "./Outputs/Supplementary/Fig.S4.pdf", width=8.27, height = 10.0)




#############################################################
#2.Analyzing how molecular rates vary with latitude and temperature at assemblage level
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
load("./DataFiles/SpatialJoinFiles/all_spatial_join_grid_mtdna.rdata")
load("./DataFiles/SpatialJoinFiles/all_spatial_join_ecoregion_mtdna.rdata")

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
marine.reptiles <- mean_molrate_grid(spatial_join=all_spatial_join_grid%>%filter(Group=="Reptiles"&Habitat!="Terrestrial"),
                                     grids_poly=grids, Species=4, habitat="Marine")
amphibians <- mean_molrate_grid(spatial_join=all_spatial_join_grid%>%filter(Group=="Amphibians"),
                                grids_poly=grids, Species=4, habitat="Terrestrial")
terr.birds <- mean_molrate_grid(spatial_join=all_spatial_join_grid%>%filter(Group=="Birds"&Habitat!="Marine"),
                                grids_poly=grids, Species=5, habitat="Terrestrial")
marine.birds <- mean_molrate_grid(spatial_join=all_spatial_join_grid%>%filter(Group=="Birds"&Habitat!="Terrestrial"),
                                  grids_poly=grids, Species=5, habitat="Marine")
terr.mammals <- mean_molrate_grid(spatial_join=all_spatial_join_grid%>%filter(Group=="Mammals"&Habitat!="Marine"),
                                  grids_poly=grids, Species=5, habitat="Terrestrial")
marine.mammals <- mean_molrate_grid(spatial_join=all_spatial_join_grid%>%filter(Group=="Mammals"&Habitat!="Terrestrial"),
                                    grids_poly=grids, Species=5, habitat="Marine")
terr.ectotherm <- mean_molrate_grid(spatial_join=all_spatial_join_grid%>%filter((Group %in% c("Amphibians", "Reptiles","Fishes"))&Habitat!="Marine"),
                                    grids_poly=grids, Species=5, habitat="Terrestrial")
marine.ectotherm <- mean_molrate_grid(spatial_join=all_spatial_join_grid%>%filter((Group %in% c("Amphibians", "Reptiles","Fishes"))&(Habitat=="Marine"|Habitat=="Marine&Freshwater"|Habitat=="Terrestrial&Marine")),
                                      grids_poly=grids, Species=5, habitat="Marine")
terr.endotherm <- mean_molrate_grid(spatial_join=all_spatial_join_grid%>%filter((Group %in% c("Mammals", "Birds"))&Habitat!="Marine"),
                                    grids_poly=grids, Species=5, habitat="Terrestrial")
marine.endotherm <- mean_molrate_grid(spatial_join=all_spatial_join_grid%>%filter((Group %in% c("Mammals", "Birds"))&Habitat!="Terrestrial"),
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
fresh.fishes.ecos <- mean_molrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter(Group=="Fishes"&System=="Terrestrial"),
                                            Species=5, habitat="Terrestrial")%>%filter(!is.na(SR))
marine.fishes.ecos <- mean_molrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter(Group=="Fishes"&System=="Marine"),
                                             Species=5, habitat="Marine")%>%filter(!is.na(SR))
terr.reptiles.ecos <- mean_molrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter(Group=="Reptiles"&System=="Terrestrial"),
                                             Species=3, habitat="Terrestrial")%>%filter(!is.na(SR))
marine.reptiles.ecos <- mean_molrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter(Group=="Reptiles"&System=="Marine"),
                                               Species=3, habitat="Marine")%>%filter(!is.na(SR))
amphibians.ecos <- mean_molrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter(Group=="Amphibians"&System=="Terrestrial"),
                                          Species=3, habitat="Terrestrial")%>%filter(!is.na(SR))

terr.birds.ecos <- mean_molrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter(Group=="Birds"&System=="Terrestrial"),
                                          Species=5, habitat="Terrestrial")%>%filter(!is.na(SR))
marine.birds.ecos <- mean_molrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter(Group=="Birds"&System=="Marine"),
                                            Species=5, habitat="Marine")%>%filter(!is.na(SR))
terr.mammals.ecos <- mean_molrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter(Group=="Mammal"&System=="Terrestrial"),
                                            Species=5, habitat="Terrestrial")%>%filter(!is.na(SR))

marine.mammals.ecos <- mean_molrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter(Group=="Mammal"&System=="Marine"),
                                              Species=5, habitat="Marine")%>%filter(!is.na(SR))

terr.ectotherm.ecos <- mean_molrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter((Group %in% c("Amphibians", "Reptiles","Fishes"))&System=="Terrestrial"),
                                              Species=5, habitat="Terrestrial")%>%filter(!is.na(SR))
marine.ectotherm.ecos <- mean_molrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter((Group %in% c("Amphibians", "Reptiles","Fishes"))&System=="Marine"),
                                                Species=5, habitat="Marine")%>%filter(!is.na(SR))

terr.endotherm.ecos <- mean_molrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter((Group %in% c("Birds", "Mammal"))&System=="Terrestrial"),
                                              Species=5, habitat="Terrestrial")%>%filter(!is.na(SR))

marine.endotherm.ecos <- mean_molrate_ecoregion(spatial_join=all_spatial_join_ecoregion%>%filter((Group %in% c("Birds", "Mammal"))&System=="Marine"),
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
#2.1. Spatial simultaneous autoregressive (SAR) models to examine relationship between 
#molecular rate and latitude at assemblages
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
#save(sar.summary.out, file = "./Outputs/Data/SAR_molrate_pattern_assemblages_mtdna.rdata")
load(file = "./Outputs/Data/SAR_molrate_pattern_assemblages_mtdna.rdata")

sar.summary.out

#write.csv(sar.summary.out,"sar.summary.out.csv")

########################################
#Fig.2DFH Molecular rate at grids
pdf("./Outputs/MainFigures/Fig2DFH.pdf", height = 4.2, width=3.6)
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


#Fig.2EGI Plots show relationship between molecular rates and latitude at ecoregions 
f2e1<-ectotherm.ecos%>%
  filter(!is.na(gm.ds))%>%
  ggplot(aes(x=AnnualTemp, y=gm.ds*100))+
  geom_point(alpha=0.7, size=0.7, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="", y=expression("dS (" * 10^-8 * ")"~(sub/site/year)), title="Ectotherms", tag="E")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = -5, y = 1, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.tag = element_text(size=8, face = "bold"),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f2e2<-endotherm.ecos%>%
  filter(!is.na(gm.ds), Habitat!="Marine")%>%
  ggplot(aes(x=AnnualTemp, y=gm.ds*100))+
  geom_point(alpha=0.7, size=0.7, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="", y="", title="Endotherms")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = -5, y = 1, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f2g1<-ectotherm.ecos%>%
  filter(!is.na(gm.dn))%>%
  ggplot(aes(x=AnnualTemp, y=gm.dn*10000))+
  geom_point(alpha=0.7, size=0.7, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="", y=expression("dN (" * 10^-10 * ")"~(sub/site/year)), tag="G")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = -5, y = 3, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.tag = element_text(size=8, face = "bold"),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f2g2<-endotherm.ecos%>%
  filter(!is.na(gm.dn), Habitat!="Marine")%>%
  ggplot(aes(x=AnnualTemp, y=gm.dn*10000))+
  geom_point(alpha=0.7, size=0.7, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="", y="")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = -5, y = 4, label = expression(italic(P)[SAR] == 0.06), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

f2i1<-ectotherm.ecos%>%
  filter(!is.na(gm.dnds))%>%
  ggplot(aes(x=AnnualTemp, y=gm.dnds))+
  geom_point(alpha=0.7, size=0.7, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x=expression(Temperature~(degree*C)), y="dN/dS", tag="I")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = -5, y = 0.05, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.tag = element_text(size=8, face = "bold"),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))



f2i2<-endotherm.ecos%>%
  filter(!is.na(gm.dnds), Habitat!="Marine")%>%
  ggplot(aes(x=AnnualTemp, y=gm.dnds))+
  geom_point(alpha=0.7, size=0.7, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x=expression(Temperature~(degree*C)), y="")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits=c(0.02, 0.04), breaks = seq(0.02, 0.04, 0.01))+
  annotate("text", x = -5, y = 0.04, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

cowplot::plot_grid(f2e1,f2e2, f2g1, f2g2, f2i1, f2i2, align = "hv", nrow = 3)
ggsave("./Outputs/MainFigures/Fig2EGI.pdf", height = 4.83, width=3.4)

#cowplot::plot_grid(f2a,f2b,f2c, f2e1, f2g1, f2i1, f2e2, f2g2, f2i2, align = "hv", byrow=F, nrow = 3)
#ggsave("./Outputs/MainFigures/Fig2dfh.pdf", height = 4.83, width=6)

#Fig.2EGI Plots show relationship between molecular rates and latitude at ecoregions 
f2e1<-ectotherm.ecos%>%
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

f2e2<-endotherm.ecos%>%
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

f2g1<-ectotherm.ecos%>%
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

f2g2<-endotherm.ecos%>%
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

f2i1<-ectotherm.ecos%>%
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



f2i2<-endotherm.ecos%>%
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

cowplot::plot_grid(f2e1,f2e2, f2g1, f2g2, f2i1, f2i2, align = "hv", nrow = 3)
ggsave("./Outputs/MainFigures/Fig2DFH.pdf", height = 5, width=3.4)

###########################
#Extended Fig.S7a Latitudinal gradients in molecular rates at assemblage level for each class.
##########################
pdf("./Outputs/Supplementary/Fig.S7ABC.pdf", width=8.27, height = 4)
par(mar=c(0.6,0.6,0.6,0.6), mfrow=c(3,5))
plot_molrate_grids(mean.rate.grids=fishes%>%mutate(dS=gm.ds), rate.type = "dS", 
                   geo.type="Global", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
mtext('Geometric mean dS', side=2, at=-3e+05, cex=0.5, line=-0.3, font=2)
mtext("A", side=3,font = 2, at=-1.5e+07, cex=0.6, line=-0.8)
title(main="Fishes", cex.main=0.8, line=-0.3)

plot_molrate_grids(mean.rate.grids=amphibians%>%mutate(dS=gm.ds), rate.type = "dS", 
                   geo.type="Terrestrial", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
mtext('', side=2, at=-3.0e+05, cex=0.6, line=-0.3, font=2)
#mtext("b", side=3,font = 2, at=-1.5e+07, cex=0.6, line=-0.8)
title(main="Amphibians", cex.main=0.8, line=-0.3)

plot_molrate_grids(mean.rate.grids=terr.reptiles%>%mutate(dS=gm.ds), rate.type = "dS", 
                   geo.type="Terrestrial", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
mtext('', side=2, at=-3.0e+05, cex=0.6, line=-0.3, font=2)
#mtext("c", side=3,font = 2, at=-1.5e+07, cex=0.6, line=-0.8)
title(main="Reptiles", cex.main=0.8, line=-0.3)

plot_molrate_grids(mean.rate.grids=terr.mammals%>%mutate(dS=gm.ds), rate.type = "dS", 
                   geo.type="Terrestrial", colramp=c("#ffffff", "#f4cdb9","#f47155","#d22428", "#660b11"))
mtext('', side=2, at=-3.0e+05, cex=0.6, line=-0.3, font=2)
#mtext("d", side=3,font = 2, at=-1.5e+07, cex=0.6, line=-0.8)
title(main="Mammals", cex.main=0.8, line=-0.3)

plot_molrate_grids(mean.rate.grids=terr.birds%>%mutate(dS=gm.ds), rate.type = "dS", 
                   geo.type="Terrestrial", colramp=c("#ffffff", "#f4cdb9","#f47155","#d22428", "#660b11"))
mtext('', side=2, at=-3.0e+05, cex=0.6, line=-0.3, font=2)
#mtext("e", side=3,font = 2, at=-1.5e+07, cex=0.6, line=-0.8)
title(main="Birds", cex.main=0.8, line=-0.3)

plot_molrate_grids(mean.rate.grids=fishes%>%mutate(dN=gm.dn), rate.type = "dN", 
                   geo.type="Global", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
mtext('Geometric mean dN', side=2, at=-3e+05, cex=0.5, line=-0.3, font=2)
mtext("B", side=3,font = 2, at=-1.5e+07, cex=0.6, line=-0.8)
title(main="Fishes", cex.main=0.8, line=-0.3)

plot_molrate_grids(mean.rate.grids=amphibians%>%mutate(dN=gm.dn), rate.type = "dN", 
                   geo.type="Terrestrial", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
title(main="Amphibians", cex.main=0.8, line=-0.3)

plot_molrate_grids(mean.rate.grids=terr.reptiles%>%mutate(dN=gm.dn), rate.type = "dN", 
                   geo.type="Terrestrial", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
title(main="Reptiles", cex.main=0.8, line=-0.3)

plot_molrate_grids(mean.rate.grids=terr.mammals%>%mutate(dN=gm.dn), rate.type = "dN", 
                   geo.type="Terrestrial", colramp=c("#ffffff", "#f4cdb9","#f47155","#d22428", "#660b11"))
title(main="Mammals", cex.main=0.8, line=-0.3)

plot_molrate_grids(mean.rate.grids=terr.birds%>%mutate(dN=gm.dn), rate.type = "dN", 
                   geo.type="Terrestrial", colramp=c("#ffffff", "#f4cdb9","#f47155","#d22428", "#660b11"))
title(main="Birds", cex.main=0.8, line=-0.3)

plot_molrate_grids(mean.rate.grids=fishes%>%mutate(dNdS=gm.dnds), rate.type = "dNdS", 
                   geo.type="Global", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
mtext('Geometric mean dN/dS', side=2, at=-3e+05, cex=0.5, line=-0.3, font=2)
mtext("C", side=3,font = 2, at=-1.5e+07, cex=0.6, line=-0.8)
title(main="Fishes", cex.main=0.8, line=-0.3)

plot_molrate_grids(mean.rate.grids=amphibians%>%mutate(dNdS=gm.dnds), rate.type = "dNdS", 
                   geo.type="Terrestrial", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
title(main="Amphibians", cex.main=0.8, line=-0.3)

plot_molrate_grids(mean.rate.grids=terr.reptiles%>%mutate(dNdS=gm.dnds), rate.type = "dNdS", 
                   geo.type="Terrestrial", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
title(main="Reptiles", cex.main=0.8, line=-0.3)

plot_molrate_grids(mean.rate.grids=terr.mammals%>%mutate(dNdS=gm.dnds), rate.type = "dNdS", 
                   geo.type="Terrestrial", colramp=c("#ffffff", "#f4cdb9","#f47155","#d22428", "#660b11"))
title(main="Mammals", cex.main=0.8, line=-0.3)

plot_molrate_grids(mean.rate.grids=terr.birds%>%mutate(dNdS=gm.dnds), rate.type = "dNdS", 
                   geo.type="Terrestrial", colramp=c("#ffffff", "#f4cdb9","#f47155","#d22428", "#660b11"))
title(main="Birds", cex.main=0.8, line=-0.3)
dev.off()


##############
#Extended Fig.S7d-f Latitudinal gradients in molecular rates at assemblage level for each class.
fs7d1 <- mean.molrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Fishes")%>%
  ggplot(aes(x=abs(Lat), y=gm.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y=expression("dS (" * 10^-8 * ")"~(sub/site/year)), title="Fishes", tag="D")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = 60, y = 1, label = expression(italic(P)[SAR] < 0.001), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        plot.tag = element_text(size=8, face = "bold"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs7d2 <- mean.molrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Amphibians", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=gm.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y="", title="Amphibians")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = 40, y = 1, label = expression(italic(P)[SAR] < 0.001), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))


fs7d3 <- mean.molrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Reptiles", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=gm.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y="", title="Reptiles")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = 40, y = 1.2, label = expression(italic(P)[SAR] == 0.0035), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs7d4 <- mean.molrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Mammals", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=gm.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="Latitude", y="", title="Mammals")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.9,2.1), breaks = seq(0.9,2.1,0.2))+
  annotate("text", x = 60, y = 2.1, label = expression(italic(P)[SAR] == 0.91), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs7d5 <- mean.molrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Birds", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=gm.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="Latitude", y="", title="Birds")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.8,1.5), breaks = seq(0.8,1.5,0.1))+
  annotate("text", x = 60, y = 1.5, label = expression(italic(P)[SAR] < 0.001), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs7e1 <- mean.molrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Fishes")%>%
  ggplot(aes(x=abs(Lat), y=gm.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y=expression("dN (" * 10^-10 * ")"~(sub/site/year)), title="Fishes", tag="E")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = 60, y = 2.8, label = expression(italic(P)[SAR] < 0.001), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.tag = element_text(size=8, face = "bold"),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs7e2 <- mean.molrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Amphibians", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=gm.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y="", title="Amphibians")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.6,3), breaks = seq(0.6,3,0.4))+
  annotate("text", x = 40, y = 3, label = expression(italic(P)[SAR] == 0.015), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs7e3 <- mean.molrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Reptiles", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=gm.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y="", title="Reptiles")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(2.2,3.6), breaks = seq(2.2,3.6,0.2))+
  annotate("text", x = 40, y = 3.6, label = expression(italic(P)[SAR] == 0.011), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs7e4 <- mean.molrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Mammals", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=gm.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="Latitude", y="", title="Mammals")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(3.9,5.5), breaks = seq(4,5.5,0.5))+
  annotate("text", x = 60, y = 5.5, label = expression(italic(P)[SAR] == 0.0012), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs7e5 <- mean.molrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Birds", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=gm.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="Latitude", y="", title="Birds")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = 60, y = 3.3, label = expression(italic(P)[SAR] == 0.55), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs7f1 <- mean.molrate.ecos%>%
  filter(!is.na(gm.dnds), Group=="Fishes")%>%
  ggplot(aes(x=abs(Lat), y=gm.dnds))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y="dN/dS", title="Fishes", tag="F")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits=c(0.02, 0.05), breaks = seq(0.02, 0.05, 0.01))+
  #scale_y_continuous(limits=c(0.02, 0.035), breaks = seq(0.02, 0.035, 0.005))+
  annotate("text", x = 60, y = 0.05, label = expression(italic(P)[SAR] == 0.014), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.tag = element_text(size=8, face = "bold"),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))


fs7f2 <- mean.molrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Amphibians", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=gm.dnds))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y="", title="Amphibians")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.01,0.05), breaks = seq(0.01,0.05,0.01))+
  annotate("text", x = 40, y = 0.05, label = expression(italic(P)[SAR] == 0.03), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs7f3 <- mean.molrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Reptiles", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=gm.dnds))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y="", title="Reptiles")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.02,0.07), breaks = seq(0.02,0.07,0.01))+
  annotate("text", x = 40, y = 0.07, label = expression(italic(P)[SAR] == 0.02), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs7f4 <- mean.molrate.ecos%>%
  filter(!is.na(gm.dnds), Group=="Mammals", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=gm.dnds))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="Latitude", y="", title="Mammals")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.02,0.05), breaks = seq(0.02,0.05,0.01))+
  annotate("text", x = 60, y = 0.05, label = expression(italic(P)[SAR] == 0.15), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs7f5 <- mean.molrate.ecos%>%
  filter(!is.na(gm.dnds), Group=="Birds", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=gm.dnds))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="Latitude", y="", title="Birds")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.02,0.035), breaks = seq(0.02,0.035,0.005))+
  annotate("text", x = 60, y = 0.035, label = expression(italic(P)[SAR] < 0.001), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))


#Extended Fig.S7g-i Latitudinal gradients in molecular rates at assemblage level for each class.
sar.summary.out%>%filter(Avg.Value=="Geometric mean")%>%filter(Var=="Temperature")%>%
  dplyr::select(1:4,"SAR.Pvalue")


fs7g1 <- mean.molrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Fishes")%>%
  ggplot(aes(x=AnnualTemp, y=gm.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x=expression(Temperature~(degree*C)), y=expression("dS (" * 10^-8 * ")"~(sub/site/year)), title="Fishes", tag="G")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = -5, y = 1, label = expression(italic(P)[SAR] < 0.001), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        plot.tag = element_text(size=8, face = "bold"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs7g2 <- mean.molrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Amphibians", Habitat=="Terrestrial")%>%
  ggplot(aes(x=AnnualTemp, y=gm.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x=expression(Temperature~(degree*C)), y="", title="Amphibians")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = -5, y = 1, label = expression(italic(P)[SAR] == 0.005), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))


fs7g3 <- mean.molrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Reptiles", Habitat=="Terrestrial")%>%
  ggplot(aes(x=AnnualTemp, y=gm.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x=expression(Temperature~(degree*C)), y="", title="Reptiles")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = -5, y = 1.2, label = expression(italic(P)[SAR] < 0.001), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs7g4 <- mean.molrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Mammals", Habitat=="Terrestrial")%>%
  ggplot(aes(x=AnnualTemp, y=gm.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x=expression(Temperature~(degree*C)), y="", title="Mammals")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.9,2.1), breaks = seq(0.9,2.1,0.2))+
  annotate("text", x = -5, y = 2.1, label = expression(italic(P)[SAR] == 0.27), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs7g5 <- mean.molrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Birds", Habitat=="Terrestrial")%>%
  ggplot(aes(x=AnnualTemp, y=gm.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x=expression(Temperature~(degree*C)), y="", title="Birds")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.8,1.5), breaks = seq(0.8,1.5,0.1))+
  annotate("text", x = -5, y = 1.5, label = expression(italic(P)[SAR] < 0.001), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs7h1 <- mean.molrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Fishes")%>%
  ggplot(aes(x=AnnualTemp, y=gm.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x=expression(Temperature~(degree*C)), y=expression("dN (" * 10^-10 * ")"~(sub/site/year)), title="Fishes", tag="H")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = -5, y = 2.8, label = expression(italic(P)[SAR] == 0.003), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.tag = element_text(size=8, face = "bold"),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs7h2 <- mean.molrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Amphibians", Habitat=="Terrestrial")%>%
  ggplot(aes(x=AnnualTemp, y=gm.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x=expression(Temperature~(degree*C)), y="", title="Amphibians")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.6,3), breaks = seq(0.6,3,0.4))+
  annotate("text", x = -5, y = 3, label = expression(italic(P)[SAR] == 0.024), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs7h3 <- mean.molrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Reptiles", Habitat=="Terrestrial")%>%
  ggplot(aes(x=AnnualTemp, y=gm.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x=expression(Temperature~(degree*C)), y="", title="Reptiles")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(2.2,3.6), breaks = seq(2.2,3.6,0.2))+
  annotate("text", x = -5, y = 3.6, label = expression(italic(P)[SAR] == 0.001), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs7h4 <- mean.molrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Mammals", Habitat=="Terrestrial")%>%
  ggplot(aes(x=AnnualTemp, y=gm.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x=expression(Temperature~(degree*C)), y="", title="Mammals")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(3.9,5.5), breaks = seq(4,5.5,0.5))+
  annotate("text", x = -5, y = 5.5, label = expression(italic(P)[SAR] == 0.06), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs7h5 <- mean.molrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Birds", Habitat=="Terrestrial")%>%
  ggplot(aes(x=AnnualTemp, y=gm.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x=expression(Temperature~(degree*C)), y="", title="Birds")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = -5, y = 3.3, label = expression(italic(P)[SAR] == 0.02), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs7i1 <- mean.molrate.ecos%>%
  filter(!is.na(gm.dnds), Group=="Fishes")%>%
  ggplot(aes(x=AnnualTemp, y=gm.dnds))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x=expression(Temperature~(degree*C)), y="dN/dS", title="Fishes", tag="I")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits=c(0.02, 0.05), breaks = seq(0.02, 0.05, 0.01))+
  #scale_y_continuous(limits=c(0.02, 0.035), breaks = seq(0.02, 0.035, 0.005))+
  annotate("text", x = -5, y = 0.05, label = expression(italic(P)[SAR] < 0.001), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.tag = element_text(size=8, face = "bold"),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))


fs7i2 <- mean.molrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Amphibians", Habitat=="Terrestrial")%>%
  ggplot(aes(x=AnnualTemp, y=gm.dnds))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x=expression(Temperature~(degree*C)), y="", title="Amphibians")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.01,0.05), breaks = seq(0.01,0.05,0.01))+
  annotate("text", x = -5, y = 0.05, label = expression(italic(P)[SAR] == 0.62), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs7i3 <- mean.molrate.ecos%>%
  filter(!is.na(gm.dn), Group=="Reptiles", Habitat=="Terrestrial")%>%
  ggplot(aes(x=AnnualTemp, y=gm.dnds))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x=expression(Temperature~(degree*C)), y="", title="Reptiles")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.02,0.07), breaks = seq(0.02,0.07,0.01))+
  annotate("text", x = -5, y = 0.07, label = expression(italic(P)[SAR] < 0.001), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs7i4 <- mean.molrate.ecos%>%
  filter(!is.na(gm.dnds), Group=="Mammals", Habitat=="Terrestrial")%>%
  ggplot(aes(x=AnnualTemp, y=gm.dnds))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x=expression(Temperature~(degree*C)), y="", title="Mammals")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.02,0.05), breaks = seq(0.02,0.05,0.01))+
  annotate("text", x = -5, y = 0.05, label = expression(italic(P)[SAR] == 0.99), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs7i5 <- mean.molrate.ecos%>%
  filter(!is.na(gm.dnds), Group=="Birds", Habitat=="Terrestrial")%>%
  ggplot(aes(x=AnnualTemp, y=gm.dnds))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x=expression(Temperature~(degree*C)), y="", title="Birds")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.02,0.035), breaks = seq(0.02,0.035,0.005))+
  annotate("text", x = -5, y = 0.035, label = expression(italic(P)[SAR] < 0.001), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))
cowplot::plot_grid(fs7d1,fs7d2, fs7d3, fs7d4, fs7d5, 
                   fs7e1,fs7e2, fs7e3, fs7e4, fs7e5,
                   fs7f1,fs7f2, fs7f3, fs7f4, fs7f5,
                   fs7g1,fs7g2, fs7g3, fs7g4, fs7g5, 
                   fs7h1,fs7h2, fs7h3, fs7h4, fs7h5,
                   fs7i1,fs7i2, fs7i3, fs7i4, fs7i5,
                   nrow = 6, align="hv")

ggsave("./Outputs/Supplementary/Fig.S7D-I.pdf", width=8.27, height=10)

########################
#Extended Fig.S8 Plot arithmetic mean values of molecular rates in cell grids for endotherms and ectotherms
#Extended Fig.S8ab
pdf("./Outputs/Supplementary/Fig.S8AB.pdf", height =4, width= 8.27)
par(mar=c(0.6,0.6,0.6,0.6), mfrow=c(2,3))
plot_molrate_grids(mean.rate.grids=ectotherm%>%mutate(dS=mid.ds)%>%filter(!is.na(SR)), 
                   rate.type = "dS", geo.type="Global", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
mtext('Ectotherms', side=2, at=-3.0e+05, cex=0.6, line=-0.1, font=2)
mtext("A", side=3,font = 2, at=-1.5e+07, cex=0.6, line=-0.8)
title(main="Median dS", cex.main=0.8, line=-0.3)

plot_molrate_grids(mean.rate.grids=ectotherm%>%mutate(dN=mid.dn), 
                   rate.type = "dN", geo.type="Global", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
title(main="Median dN", cex.main=0.8, line=-0.3)
mtext('Ectotherms', side=2, at=-3.0e+05, cex=0.6, line=-0.1, font=2)
#mtext("b", side=3,font = 2, at=-1.5e+07, cex=0.6, line=-0.8)

plot_molrate_grids(mean.rate.grids=ectotherm%>%mutate(dNdS=mid.dnds), 
                   rate.type = "dNdS", geo.type="Global", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
title(main="Median dN/dS", cex.main=0.8, line=-0.3)
mtext('Ectotherms', side=2, at=-3.0e+05, cex=0.6, line=-0.1, font=2)
#mtext("b", side=3,font = 2, at=-1.5e+07, cex=0.6, line=-0.8)

plot_molrate_grids(mean.rate.grids=endotherm%>%mutate(dS=mid.ds)%>%filter(Habitat!="Marine"), 
                   rate.type = "dS", geo.type="Terrestrial", colramp=c("#ffffff", "#f4cdb9","#f47155","#d22428", "#660b11"))
mtext('Endotherms', side=2, at=-3.0e+05, cex=0.6, line=-0.1, font=2)
mtext("B", side=3,font = 2, at=-1.5e+07, cex=0.6, line=-0.8)
title(main="Median dS", cex.main=0.8, line=-0.3)

plot_molrate_grids(mean.rate.grids=endotherm%>%mutate(dN=mid.dn)%>%filter(Habitat!="Marine"), 
                   rate.type = "dN", geo.type="Terrestrial", colramp=c("#ffffff", "#f4cdb9","#f47155","#d22428", "#660b11"))
mtext('Endotherms', side=2, at=-3.0e+05, cex=0.6, line=-0.1, font=2)
title(main="Median dN", cex.main=0.8, line=-0.3)
#mtext("d", side=3,font = 2, at=-1.5e+07, cex=0.6, line=-0.8)

plot_molrate_grids(mean.rate.grids=endotherm%>%mutate(dNdS=mid.dnds)%>%filter(Habitat!="Marine"), 
                   rate.type = "dNdS", geo.type="Terrestrial", colramp=c("#ffffff", "#f4cdb9","#f47155","#d22428", "#660b11"))
mtext('Endotherms', side=2, at=-3.0e+05, cex=0.6, line=-0.1, font=2)
title(main="Median dN/dS", cex.main=0.8, line=-0.3)
#mtext("d", side=3,font = 2, at=-1.5e+07, cex=0.6, line=-0.8)

dev.off()

#Extended Fig.S8CD 
fs8c1<-ectotherm.ecos%>%
  filter(!is.na(mid.ds))%>%
  ggplot(aes(x=abs(Lat), y=mid.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y=expression("Median dS (" * 10^-8 * ")"~(sub/site/year)), title="Ectotherms", tag="C")+
  scale_y_continuous(limits = c(0.4, 1.2),breaks = seq(0.4, 1.2, 0.2))+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = 60, y = 1.2, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.tag = element_text(size=8, face = "bold"),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

fs8c2<-ectotherm.ecos%>%
  filter(!is.na(mid.dn))%>%
  ggplot(aes(x=abs(Lat), y=mid.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y=expression("Median dN (" * 10^-10 * ")"~(sub/site/year)), title="Ectotherms")+
  scale_y_continuous(limits = c(1.4, 3.5),breaks = seq(1.5, 3.5, 1))+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = 60, y = 3.5, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

fs8c3<-ectotherm.ecos%>%
  filter(!is.na(mid.dnds))%>%
  ggplot(aes(x=abs(Lat), y=mid.dnds))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y="Median dN/dS", title="Ectotherms")+
  scale_y_continuous(limits = c(0.02, 0.05),breaks = seq(0, 0.05, 0.01))+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = 60, y = 0.05, label = expression(italic(P)[SAR] == 0.18), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))


fs8d1<-endotherm.ecos%>%
  filter(!is.na(mid.ds), Habitat!="Marine")%>%
  ggplot(aes(x=abs(Lat), y=mid.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="Latitude", y=expression("Median dS (" * 10^-8 * ")"~(sub/site/year)), title="Endotherms", tag="D")+
  scale_y_continuous(limits = c(0.8, 1.6),breaks = seq(0.8, 1.6, 0.2))+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = 60, y = 1.6, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.tag = element_text(size=8, face = "bold"),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

fs8d2<-endotherm.ecos%>%
  filter(!is.na(mid.dn), Habitat!="Marine")%>%
  ggplot(aes(x=abs(Lat), y=mid.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="Latitude", y=expression("Median dN (" * 10^-10 * ")"~(sub/site/year)), title="Endotherms")+
  scale_y_continuous(limits = c(2.5, 4),breaks = seq(2.5, 4, 0.5))+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = 60, y = 4, label = expression(italic(P)[SAR] == 0.15), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

fs8d3<-endotherm.ecos%>%
  filter(!is.na(mid.dn), Habitat!="Marine")%>%
  ggplot(aes(x=abs(Lat), y=mid.dnds))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="Latitude", y="Median dN/dS", title="Endotherms")+
  scale_y_continuous(limits = c(0.02, 0.035),breaks = seq(0.02, 0.035, 0.005))+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = 60, y = 0.035, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

sar.summary.out%>%filter(Avg.Value=="Median")%>%filter(Var=="Temperature")%>%
  dplyr::select(1:4,"SAR.Pvalue")


#Extended Fig.S8EF 
fs8e1<-ectotherm.ecos%>%
  filter(!is.na(mid.ds))%>%
  ggplot(aes(x=AnnualTemp, y=mid.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x=expression(Temperature~(degree*C)), y=expression("Median dS (" * 10^-8 * ")"~(sub/site/year)), title="Ectotherms", tag="E")+
  scale_y_continuous(limits = c(0.4, 1.2),breaks = seq(0.4, 1.2, 0.2))+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = -5, y = 1.2, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.tag = element_text(size=8, face = "bold"),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

fs8e2<-ectotherm.ecos%>%
  filter(!is.na(mid.dn))%>%
  ggplot(aes(x=AnnualTemp, y=mid.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x=expression(Temperature~(degree*C)), y=expression("Median dN (" * 10^-10 * ")"~(sub/site/year)), title="Ectotherms")+
  scale_y_continuous(limits = c(1.4, 3.5),breaks = seq(1.5, 3.5, 1))+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = -5, y = 3.5, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

fs8e3<-ectotherm.ecos%>%
  filter(!is.na(mid.dnds))%>%
  ggplot(aes(x=AnnualTemp, y=mid.dnds))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x=expression(Temperature~(degree*C)), y="Median dN/dS", title="Ectotherms")+
  scale_y_continuous(limits = c(0.02, 0.05),breaks = seq(0, 0.05, 0.01))+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = -5, y = 0.05, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))


fs8f1<-endotherm.ecos%>%
  filter(!is.na(mid.ds), Habitat!="Marine")%>%
  ggplot(aes(x=AnnualTemp, y=mid.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x=expression(Temperature~(degree*C)), y=expression("Median dS (" * 10^-8 * ")"~(sub/site/year)), title="Endotherms", tag="F")+
  scale_y_continuous(limits = c(0.8, 1.6),breaks = seq(0.8, 1.6, 0.2))+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = -5, y = 1.6, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.tag = element_text(size=8, face = "bold"),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

fs8f2<-endotherm.ecos%>%
  filter(!is.na(mid.dn), Habitat!="Marine")%>%
  ggplot(aes(x=AnnualTemp, y=mid.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x=expression(Temperature~(degree*C)), y=expression("Median dN (" * 10^-10 * ")"~(sub/site/year)), title="Endotherms")+
  scale_y_continuous(limits = c(2.5, 4),breaks = seq(2.5, 4, 0.5))+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = -5, y = 4, label = expression(italic(P)[SAR] == 0.003), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

fs8f3<-endotherm.ecos%>%
  filter(!is.na(mid.dn), Habitat!="Marine")%>%
  ggplot(aes(x=AnnualTemp, y=mid.dnds))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x=expression(Temperature~(degree*C)), y="Median dN/dS", title="Endotherms")+
  scale_y_continuous(limits = c(0.02, 0.035),breaks = seq(0.02, 0.035, 0.005))+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  annotate("text", x = -5, y = 0.035, label = expression(italic(P)[SAR] < 0.001), size=2.5)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))


cowplot::plot_grid(fs8c1,fs8c2,fs8c3,fs8d1,fs8d2,fs8d3,
                   fs8e1,fs8e2,fs8e3,fs8f1,fs8f2,fs8f3,
                   align = "hv", nrow = 4)
ggsave("./Outputs/Supplementary/Fig.S8C-F.pdf", height =10, width= 8.27)


###########################
#Extended Fig.S9abc Latitudinal gradients in molecular rates at assemblage level for each class.
##########################
pdf("./Outputs/Supplementary/Fig.S9ABC.pdf", width=8.27, height = 4)
par(mar=c(0.6,0.6,0.6,0.6), mfrow=c(3,5))
plot_molrate_grids(mean.rate.grids=fishes%>%mutate(dS=mid.ds), rate.type = "dS", 
                   geo.type="Global", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
mtext('Median dS', side=2, at=-3e+05, cex=0.5, line=-0.3, font=2)
mtext("A", side=3,font = 2, at=-1.5e+07, cex=0.6, line=-0.8)
title(main="Fishes", cex.main=0.8, line=-0.3)

plot_molrate_grids(mean.rate.grids=amphibians%>%mutate(dS=mid.ds), rate.type = "dS", 
                   geo.type="Terrestrial", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
mtext('', side=2, at=-3.0e+05, cex=0.6, line=-0.3, font=2)
#mtext("b", side=3,font = 2, at=-1.5e+07, cex=0.6, line=-0.8)
title(main="Amphibians", cex.main=0.8, line=-0.3)

plot_molrate_grids(mean.rate.grids=terr.reptiles%>%mutate(dS=mid.ds), rate.type = "dS", 
                   geo.type="Terrestrial", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
mtext('', side=2, at=-3.0e+05, cex=0.6, line=-0.3, font=2)
#mtext("c", side=3,font = 2, at=-1.5e+07, cex=0.6, line=-0.8)
title(main="Reptiles", cex.main=0.8, line=-0.3)

plot_molrate_grids(mean.rate.grids=terr.mammals%>%mutate(dS=mid.ds), rate.type = "dS", 
                   geo.type="Terrestrial", colramp=c("#ffffff", "#f4cdb9","#f47155","#d22428", "#660b11"))
mtext('', side=2, at=-3.0e+05, cex=0.6, line=-0.3, font=2)
#mtext("d", side=3,font = 2, at=-1.5e+07, cex=0.6, line=-0.8)
title(main="Mammals", cex.main=0.8, line=-0.3)

plot_molrate_grids(mean.rate.grids=terr.birds%>%mutate(dS=mid.ds), rate.type = "dS", 
                   geo.type="Terrestrial", colramp=c("#ffffff", "#f4cdb9","#f47155","#d22428", "#660b11"))
mtext('', side=2, at=-3.0e+05, cex=0.6, line=-0.3, font=2)
#mtext("e", side=3,font = 2, at=-1.5e+07, cex=0.6, line=-0.8)
title(main="Birds", cex.main=0.8, line=-0.3)

plot_molrate_grids(mean.rate.grids=fishes%>%mutate(dN=mid.dn), rate.type = "dN", 
                   geo.type="Global", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
mtext('Median dN', side=2, at=-3e+05, cex=0.5, line=-0.3, font=2)
mtext("B", side=3,font = 2, at=-1.5e+07, cex=0.6, line=-0.8)
title(main="Fishes", cex.main=0.8, line=-0.3)

plot_molrate_grids(mean.rate.grids=amphibians%>%mutate(dN=mid.dn), rate.type = "dN", 
                   geo.type="Terrestrial", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
title(main="Amphibians", cex.main=0.8, line=-0.3)

plot_molrate_grids(mean.rate.grids=terr.reptiles%>%mutate(dN=mid.dn), rate.type = "dN", 
                   geo.type="Terrestrial", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
title(main="Reptiles", cex.main=0.8, line=-0.3)

plot_molrate_grids(mean.rate.grids=terr.mammals%>%mutate(dN=mid.dn), rate.type = "dN", 
                   geo.type="Terrestrial", colramp=c("#ffffff", "#f4cdb9","#f47155","#d22428", "#660b11"))
title(main="Mammals", cex.main=0.8, line=-0.3)

plot_molrate_grids(mean.rate.grids=terr.birds%>%mutate(dN=mid.dn), rate.type = "dN", 
                   geo.type="Terrestrial", colramp=c("#ffffff", "#f4cdb9","#f47155","#d22428", "#660b11"))
title(main="Birds", cex.main=0.8, line=-0.3)

plot_molrate_grids(mean.rate.grids=fishes%>%mutate(dNdS=mid.dnds), rate.type = "dNdS", 
                   geo.type="Global", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
mtext('Median dN/dS', side=2, at=-3e+05, cex=0.5, line=-0.3, font=2)
mtext("C", side=3,font = 2, at=-1.5e+07, cex=0.6, line=-0.8)
title(main="Fishes", cex.main=0.8, line=-0.3)

plot_molrate_grids(mean.rate.grids=amphibians%>%mutate(dNdS=mid.dnds), rate.type = "dNdS", 
                   geo.type="Terrestrial", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
title(main="Amphibians", cex.main=0.8, line=-0.3)

plot_molrate_grids(mean.rate.grids=terr.reptiles%>%mutate(dNdS=mid.dnds), rate.type = "dNdS", 
                   geo.type="Terrestrial", colramp=c("#ffffff","#cfe0ea", "#95b9cd", "#4076ab", "#0c10a0"))
title(main="Reptiles", cex.main=0.8, line=-0.3)

plot_molrate_grids(mean.rate.grids=terr.mammals%>%mutate(dNdS=mid.dnds), rate.type = "dNdS", 
                   geo.type="Terrestrial", colramp=c("#ffffff", "#f4cdb9","#f47155","#d22428", "#660b11"))
title(main="Mammals", cex.main=0.8, line=-0.3)

plot_molrate_grids(mean.rate.grids=terr.birds%>%mutate(dNdS=mid.dnds), rate.type = "dNdS", 
                   geo.type="Terrestrial", colramp=c("#ffffff", "#f4cdb9","#f47155","#d22428", "#660b11"))
title(main="Birds", cex.main=0.8, line=-0.3)
dev.off()


##############
#Extended Fig.S9DEF Latitudinal gradients in molecular rates at assemblage level for each class.
fs9d1 <- mean.molrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Fishes")%>%
  ggplot(aes(x=abs(Lat), y=mid.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y=expression("Median dS (" * 10^-8 * ")"~(sub/site/year)), title="Fishes", tag="D")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.5,1.1), breaks = seq(0.5,1.1,0.1))+
  annotate("text", x = 60, y = 1.1, label = expression(italic(P)[SAR] == 0.004), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.tag = element_text(size=8, face = "bold"),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

fs9d2 <- mean.molrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Amphibians", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=mid.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y="", title="Amphibians")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0,1.8), breaks = seq(0,1.8,0.3))+
  annotate("text", x = 40, y = 1.8, label = expression(italic(P)[SAR] < 0.001), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

fs9d3 <- mean.molrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Reptiles", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=mid.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y="", title="Reptiles")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.3,1.5), breaks = seq(0.3,1.5,0.3))+
  annotate("text", x = 40, y = 1.5, label = expression(italic(P)[SAR] == 0.001), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

fs9d4 <- mean.molrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Mammals", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=mid.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="Latitude", y="", title="Mammals")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.9,2.3), breaks = seq(0.9,2.3,0.2))+
  annotate("text", x = 60, y = 2.3, label = expression(italic(P)[SAR] == 0.84), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

fs9d5 <- mean.molrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Birds", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=mid.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="Latitude", y="", title="Birds")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.8,1.5), breaks = seq(0.8,1.5,0.1))+
  annotate("text", x = 60, y = 1.5, label = expression(italic(P)[SAR] < 0.001), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

fs9e1 <- mean.molrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Fishes")%>%
  ggplot(aes(x=abs(Lat), y=mid.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y=expression("Median dN (" * 10^-10 * ")"~(sub/site/year)), title="Fishes", tag="E")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(1.4,2.9), breaks = seq(1.4,2.9,0.3))+
  annotate("text", x = 60, y = 2.9, label = expression(italic(P)[SAR] < 0.001), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.tag = element_text(size=8, face = "bold"),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

fs9e2 <- mean.molrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Amphibians", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=mid.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y="", title="Amphibians")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.6,3.5), breaks = seq(0.6,3.4,0.4))+
  annotate("text", x = 40, y = 3.4, label = expression(italic(P)[SAR] == 0.013), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

fs9e3 <- mean.molrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Reptiles", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=mid.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y="", title="Reptiles")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(2,4), breaks = seq(2,4,0.5))+
  annotate("text", x = 40, y = 4, label = expression(italic(P)[SAR] == 0.003), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

fs9e4 <- mean.molrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Mammals", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=mid.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="Latitude", y="", title="Mammals")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(4,5.8), breaks = seq(4,5.8,0.3))+
  annotate("text", x = 60, y = 5.8, label = expression(italic(P)[SAR] < 0.001), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

fs9e5 <- mean.molrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Birds", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=mid.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="Latitude", y="", title="Birds")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(2.6,3.4), breaks = seq(2.6,3.4,0.2))+
  annotate("text", x = 60, y = 3.4, label = expression(italic(P)[SAR] == 0.21), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

fs9f1 <- mean.molrate.ecos%>%
  filter(!is.na(mid.dnds), Group=="Fishes")%>%
  ggplot(aes(x=abs(Lat), y=mid.dnds))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y="Median dN/dS", title="Fishes", tag="F")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits=c(0.015, 0.05), breaks = seq(0.015, 0.05, 0.005))+
  #scale_y_continuous(limits=c(0.02, 0.035), breaks = seq(0.02, 0.035, 0.005))+
  annotate("text", x = 60, y = 0.05, label = expression(italic(P)[SAR] == 0.12), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.tag = element_text(size=8, face = "bold"),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))


fs9f2 <- mean.molrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Amphibians", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=mid.dnds))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y="", title="Amphibians")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.01,0.06), breaks = seq(0.01,0.06,0.01))+
  annotate("text", x = 40, y = 0.06, label = expression(italic(P)[SAR] == 0.72), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs9f3 <- mean.molrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Reptiles", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=mid.dnds))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x="Latitude", y="", title="Reptiles")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.02,0.07), breaks = seq(0.02,0.07,0.01))+
  annotate("text", x = 40, y = 0.07, label = expression(italic(P)[SAR] < 0.001), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs9f4 <- mean.molrate.ecos%>%
  filter(!is.na(mid.dnds), Group=="Mammals", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=mid.dnds))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="Latitude", y="", title="Mammals")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.02,0.04), breaks = seq(0.02,0.04,0.005))+
  annotate("text", x = 60, y = 0.04, label = expression(italic(P)[SAR] == 0.14), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs9f5 <- mean.molrate.ecos%>%
  filter(!is.na(mid.dnds), Group=="Birds", Habitat=="Terrestrial")%>%
  ggplot(aes(x=abs(Lat), y=mid.dnds))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x="Latitude", y="", title="Birds")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.02,0.03), breaks = seq(0.02,0.03,0.002))+
  annotate("text", x = 60, y = 0.03, label = expression(italic(P)[SAR] < 0.001), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

########################################
sar.summary.out%>%filter(Avg.Value=="Median", Var=="Temperature")%>%dplyr::select(1:4, SAR.Pvalue)
fs9g1 <- mean.molrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Fishes")%>%
  ggplot(aes(x=AnnualTemp, y=mid.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x=expression(Temperature~(degree*C)), y=expression("Median dS (" * 10^-8 * ")"~(sub/site/year)), title="Fishes", tag="G")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.5,1.1), breaks = seq(0.5,1.1,0.1))+
  annotate("text", x = -5, y = 1.1, label = expression(italic(P)[SAR] < 0.001), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.tag = element_text(size=8, face = "bold"),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

fs9g2 <- mean.molrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Amphibians", Habitat=="Terrestrial")%>%
  ggplot(aes(x=AnnualTemp, y=mid.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x=expression(Temperature~(degree*C)), y="", title="Amphibians")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0,1.8), breaks = seq(0,1.8,0.3))+
  annotate("text", x = -5, y = 1.8, label = expression(italic(P)[SAR] == 0.004), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

fs9g3 <- mean.molrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Reptiles", Habitat=="Terrestrial")%>%
  ggplot(aes(x=AnnualTemp, y=mid.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x=expression(Temperature~(degree*C)), y="", title="Reptiles")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.3,1.5), breaks = seq(0.3,1.5,0.3))+
  annotate("text", x = -5, y = 1.5, label = expression(italic(P)[SAR] < 0.001), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

fs9g4 <- mean.molrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Mammals", Habitat=="Terrestrial")%>%
  ggplot(aes(x=AnnualTemp, y=mid.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x=expression(Temperature~(degree*C)), y="", title="Mammals")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.9,2.3), breaks = seq(0.9,2.3,0.2))+
  annotate("text", x = -5, y = 2.3, label = expression(italic(P)[SAR] == 0.61), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

fs9g5 <- mean.molrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Birds", Habitat=="Terrestrial")%>%
  ggplot(aes(x=AnnualTemp, y=mid.ds*100))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x=expression(Temperature~(degree*C)), y="", title="Birds")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.8,1.5), breaks = seq(0.8,1.5,0.1))+
  annotate("text", x = -5, y = 1.5, label = expression(italic(P)[SAR] < 0.001), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

fs9h1 <- mean.molrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Fishes")%>%
  ggplot(aes(x=AnnualTemp, y=mid.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x=expression(Temperature~(degree*C)), y=expression("Median dN (" * 10^-10 * ")"~(sub/site/year)), title="Fishes", tag="H")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(1.4,2.9), breaks = seq(1.4,2.9,0.3))+
  annotate("text", x = -5, y = 2.9, label = expression(italic(P)[SAR] < 0.001), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.tag = element_text(size=8, face = "bold"),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

fs9h2 <- mean.molrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Amphibians", Habitat=="Terrestrial")%>%
  ggplot(aes(x=AnnualTemp, y=mid.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x=expression(Temperature~(degree*C)), y="", title="Amphibians")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.6,3.5), breaks = seq(0.6,3.4,0.4))+
  annotate("text", x = -5, y = 3.4, label = expression(italic(P)[SAR] == 0.018), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

fs9h3 <- mean.molrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Reptiles", Habitat=="Terrestrial")%>%
  ggplot(aes(x=AnnualTemp, y=mid.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x=expression(Temperature~(degree*C)), y="", title="Reptiles")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(2,4), breaks = seq(2,4,0.5))+
  annotate("text", x = -5, y = 4, label = expression(italic(P)[SAR] < 0.001), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

fs9h4 <- mean.molrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Mammals", Habitat=="Terrestrial")%>%
  ggplot(aes(x=AnnualTemp, y=mid.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x=expression(Temperature~(degree*C)), y="", title="Mammals")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(4,5.8), breaks = seq(4,5.8,0.3))+
  annotate("text", x = -5, y = 5.8, label = expression(italic(P)[SAR] == 0.008), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

fs9h5 <- mean.molrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Birds", Habitat=="Terrestrial")%>%
  ggplot(aes(x=AnnualTemp, y=mid.dn*10000))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x=expression(Temperature~(degree*C)), y="", title="Birds")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(2.6,3.4), breaks = seq(2.6,3.4,0.2))+
  annotate("text", x = -5, y = 3.4, label = expression(italic(P)[SAR] == 0.18), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))

fs9i1 <- mean.molrate.ecos%>%
  filter(!is.na(mid.dnds), Group=="Fishes")%>%
  ggplot(aes(x=AnnualTemp, y=mid.dnds))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x=expression(Temperature~(degree*C)), y="Median dN/dS", title="Fishes", tag="I")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits=c(0.015, 0.05), breaks = seq(0.015, 0.05, 0.005))+
  #scale_y_continuous(limits=c(0.02, 0.035), breaks = seq(0.02, 0.035, 0.005))+
  annotate("text", x = -5, y = 0.05, label = expression(italic(P)[SAR] < 0.001), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.tag = element_text(size=8, face = "bold"),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))


fs9i2 <- mean.molrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Amphibians", Habitat=="Terrestrial")%>%
  ggplot(aes(x=AnnualTemp, y=mid.dnds))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x=expression(Temperature~(degree*C)), y="", title="Amphibians")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.01,0.06), breaks = seq(0.01,0.06,0.01))+
  annotate("text", x = -5, y = 0.06, label = expression(italic(P)[SAR] == 0.86), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs9i3 <- mean.molrate.ecos%>%
  filter(!is.na(mid.dn), Group=="Reptiles", Habitat=="Terrestrial")%>%
  ggplot(aes(x=AnnualTemp, y=mid.dnds))+
  geom_point(alpha=0.7, size=0.8, colour="#2a6aaf")+
  geom_smooth(method="lm",  size=0.4, colour="#2a6aaf")+
  labs(x=expression(Temperature~(degree*C)), y="", title="Reptiles")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.02,0.07), breaks = seq(0.02,0.07,0.01))+
  annotate("text", x = -5, y = 0.07, label = expression(italic(P)[SAR] < 0.001), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs9i4 <- mean.molrate.ecos%>%
  filter(!is.na(mid.dnds), Group=="Mammals", Habitat=="Terrestrial")%>%
  ggplot(aes(x=AnnualTemp, y=mid.dnds))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x=expression(Temperature~(degree*C)), y="", title="Mammals")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.02,0.04), breaks = seq(0.02,0.04,0.005))+
  annotate("text", x = -5, y = 0.04, label = expression(italic(P)[SAR] == 0.81), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

fs9i5 <- mean.molrate.ecos%>%
  filter(!is.na(mid.dnds), Group=="Birds", Habitat=="Terrestrial")%>%
  ggplot(aes(x=AnnualTemp, y=mid.dnds))+
  geom_point(alpha=0.7, size=0.8, colour="#d3292f")+
  geom_smooth(method="lm",  size=0.4, colour="#d3292f")+
  labs(x=expression(Temperature~(degree*C)), y="", title="Birds")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0.02,0.03), breaks = seq(0.02,0.03,0.002))+
  annotate("text", x = -5, y = 0.03, label = expression(italic(P)[SAR] < 0.001), size=2.0)+
  theme_classic()+
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))

cowplot::plot_grid(fs9d1,fs9d2, fs9d3, fs9d4, fs9d5, 
                   fs9e1,fs9e2, fs9e3, fs9e4, fs9e5,
                   fs9f1,fs9f2, fs9f3, fs9f4, fs9f5,
                   fs9g1,fs9g2, fs9g3, fs9g4, fs9g5, 
                   fs9h1,fs9h2, fs9h3, fs9h4, fs9h5,
                   fs9i1,fs9i2, fs9i3, fs9i4, fs9i5,
                   nrow = 6, align="hv")

ggsave("./Outputs/Supplementary/Extended Fig8D-I.pdf", width=8.27, height=10)



##################################################################
#Part III: R script to predict molecular rates using multiple PGLMMs
###################################################################

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
molrate <- read.csv(paste0(workdir, "/DataFiles/MolEvolRate/MolRate_mtDNA.csv"))%>%
  mutate(AnnualTemp=ifelse(ThermoMode=="Ectotherms"&AnnualTemp<0.01, 0.01, AnnualTemp))


row.names(molrate) <- molrate$Species

#input phylogeny
phy <- read.tree(paste0(workdir, "/DataFiles/trees/phy_all_sampled_mtdna.tre"))
phylo <- force.ultrametric(phy, method="nnls")#nnls or extend


#get most recent common ancestor (node) of species pair
Ainv <- inverseA(phylo)$Ainv


####################################
#1. Correlations between molecular rates and life-history traits
#Extended Fig.S10
groups <- c("Fishes","Amphibians","Reptiles","Mammals","Birds")

pdf("./Outputs/Supplementary/Extended Fig.S10.pdf", width=9.27, height = 4.2)
layout(matrix(1:15, nrow = 3, ncol = 5, byrow = F))

for(group in groups){
  mydata<-molrate%>%
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
  
  mydata<-molrate%>%
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
  mydata<-molrate%>%
    filter(Group==group)%>%
    mutate(dNdS=log(dNdS),BodyMass=log(BodyMass), Fecundity=log(Fecundity), 
           MaturityAge=log(MaturityAge), Longevity=log(Longevity))%>%
    select(dNdS, BodyMass, Fecundity, MaturityAge, Longevity)
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
scale_molrate <- molrate %>%
  mutate(dS=scale(log(dS)), dN=scale(log(dN)),AnnualTemp=scale(AnnualTemp),
         BodyMass=scale(log(BodyMass)), Fecundity=scale(log(Fecundity)), 
         MaturityAge=scale(log(MaturityAge)), Longevity=scale(log(Longevity)))


#################################
prior <- list(R=list(V=diag(1),nu=0.002), 
              G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*1000)))


full.group.interact.ds<- MCMCglmm(dS~AnnualTemp*Group+BodyMass*Group+Fecundity*Group+MaturityAge*Group+Longevity*Group, 
                                  random=~Species, prior = prior, ginverse=list(Species=Ainv), 
                                  data=scale_molrate, verbose=TRUE,
                                  nitt=60000, burnin=10000, thin=25, family = c("gaussian"))

full.thermo.interact.ds<- MCMCglmm(dS~AnnualTemp*ThermoMode+BodyMass*ThermoMode+Fecundity*ThermoMode+MaturityAge*ThermoMode+Longevity*ThermoMode, 
                                   random=~Species, prior = prior, ginverse=list(Species=Ainv), 
                                   data=scale_molrate, verbose=TRUE,
                                   nitt=60000, burnin=10000, thin=25, family = c("gaussian"))

full.no.interact.ds <- MCMCglmm(dS~AnnualTemp+BodyMass+Fecundity+MaturityAge+Longevity, 
                                random=~Species, prior = prior, ginverse=list(Species=Ainv),
                                data=scale_molrate, verbose=TRUE,
                                nitt=60000, burnin=10000, thin=25, family = c("gaussian"))

fixed.only.group.interact.ds<- MCMCglmm(dS~AnnualTemp*Group+BodyMass*Group+Fecundity*Group+MaturityAge*Group+Longevity*Group,
                                        data=scale_molrate, verbose=TRUE,
                                        nitt=60000, burnin=10000, thin=25, family = c("gaussian"))

fixed.only.thermo.interact.ds<- MCMCglmm(dS~AnnualTemp*ThermoMode+BodyMass*ThermoMode+Fecundity*ThermoMode+MaturityAge*ThermoMode+Longevity*ThermoMode,
                                         data=scale_molrate, verbose=TRUE,
                                         nitt=60000, burnin=10000, thin=25, family = c("gaussian"))


fixed.only.no.interact.ds <- MCMCglmm(dS~AnnualTemp+BodyMass+Fecundity+MaturityAge+Longevity, 
                                      data=scale_molrate, verbose=TRUE, 
                                      nitt=60000, burnin=10000, thin=25, family = c("gaussian"))

phylo.only.ds<- MCMCglmm(dS~1,
                         random=~Species, prior = prior, ginverse=list(Species=Ainv),
                         data=scale_molrate, verbose=TRUE, 
                         nitt=60000, burnin=10000, thin=25, family = c("gaussian"))


full.group.interact.dn<- MCMCglmm(dN~dS*Group+AnnualTemp*Group+BodyMass*Group+Fecundity*Group+MaturityAge*Group+Longevity*Group, 
                                  random=~Species, prior = prior, ginverse=list(Species=Ainv), 
                                  data=scale_molrate, verbose=TRUE,
                                  nitt=60000, burnin=10000, thin=25, family = c("gaussian"))

full.thermo.interact.dn<- MCMCglmm(dN~dS*ThermoMode+AnnualTemp*ThermoMode+BodyMass*ThermoMode+Fecundity*ThermoMode+MaturityAge*ThermoMode+Longevity*ThermoMode, 
                                   random=~Species, prior = prior, ginverse=list(Species=Ainv), 
                                   data=scale_molrate, verbose=TRUE,
                                   nitt=60000, burnin=10000, thin=25, family = c("gaussian"))

full.no.interact.dn <- MCMCglmm(dN~dS+AnnualTemp+BodyMass+Fecundity+MaturityAge+Longevity, 
                                random=~Species, prior = prior, ginverse=list(Species=Ainv),
                                data=scale_molrate, verbose=TRUE,
                                nitt=60000, burnin=10000, thin=25, family = c("gaussian"))

fixed.only.group.interact.dn<- MCMCglmm(dN~dS*Group+AnnualTemp*Group+BodyMass*Group+Fecundity*Group+MaturityAge*Group+Longevity*Group,
                                        data=scale_molrate, verbose=TRUE,
                                        nitt=60000, burnin=10000, thin=25, family = c("gaussian"))

fixed.only.thermo.interact.dn<- MCMCglmm(dN~dS*ThermoMode+AnnualTemp*ThermoMode+BodyMass*ThermoMode+Fecundity*ThermoMode+MaturityAge*ThermoMode+Longevity*ThermoMode,
                                         data=scale_molrate, verbose=TRUE,
                                         nitt=60000, burnin=10000, thin=25, family = c("gaussian"))


fixed.only.no.interact.dn <- MCMCglmm(dN~dS+AnnualTemp+BodyMass+Fecundity+MaturityAge+Longevity, 
                                      data=scale_molrate, verbose=TRUE, 
                                      nitt=60000, burnin=10000, thin=25, family = c("gaussian"))

phylo.only.dn<- MCMCglmm(dN~1,
                         random=~Species, prior = prior, ginverse=list(Species=Ainv),
                         data=scale_molrate, verbose=TRUE, 
                         nitt=60000, burnin=10000, thin=25, family = c("gaussian"))

full.group.interact.dnds<- MCMCglmm(dNdS~AnnualTemp*Group+BodyMass*Group+Fecundity*Group+MaturityAge*Group+Longevity*Group, 
                                    random=~Species, prior = prior, ginverse=list(Species=Ainv), 
                                    data=scale_molrate, verbose=TRUE,
                                    nitt=60000, burnin=10000, thin=25, family = c("gaussian"))

full.thermo.interact.dnds<- MCMCglmm(dNdS~AnnualTemp*ThermoMode+BodyMass*ThermoMode+Fecundity*ThermoMode+MaturityAge*ThermoMode+Longevity*ThermoMode, 
                                     random=~Species, prior = prior, ginverse=list(Species=Ainv), 
                                     data=scale_molrate, verbose=TRUE,
                                     nitt=60000, burnin=10000, thin=25, family = c("gaussian"))

full.no.interact.dnds <- MCMCglmm(dNdS~AnnualTemp+BodyMass+Fecundity+MaturityAge+Longevity, 
                                  random=~Species, prior = prior, ginverse=list(Species=Ainv),
                                  data=scale_molrate, verbose=TRUE,
                                  nitt=60000, burnin=10000, thin=25, family = c("gaussian"))

fixed.only.group.interact.dnds<- MCMCglmm(dNdS~AnnualTemp*Group+BodyMass*Group+Fecundity*Group+MaturityAge*Group+Longevity*Group,
                                          data=scale_molrate, verbose=TRUE,
                                          nitt=60000, burnin=10000, thin=25, family = c("gaussian"))

fixed.only.thermo.interact.dnds<- MCMCglmm(dNdS~AnnualTemp*ThermoMode+BodyMass*ThermoMode+Fecundity*ThermoMode+MaturityAge*ThermoMode+Longevity*ThermoMode,
                                           data=scale_molrate, verbose=TRUE,
                                           nitt=60000, burnin=10000, thin=25, family = c("gaussian"))


fixed.only.no.interact.dnds <- MCMCglmm(dNdS~AnnualTemp+BodyMass+Fecundity+MaturityAge+Longevity, 
                                        data=scale_molrate, verbose=TRUE, 
                                        nitt=60000, burnin=10000, thin=25, family = c("gaussian"))

phylo.only.dnds<- MCMCglmm(dNdS~1,
                           random=~Species, prior = prior, ginverse=list(Species=Ainv),
                           data=scale_molrate, verbose=TRUE, 
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
                       phylo.only.dn=phylo.only.dn,
                       full.group.interact.dnds=full.group.interact.dnds,
                       full.thermo.interact.dnds=full.thermo.interact.dnds,
                       full.no.interact.dnds=full.no.interact.dnds,
                       fixed.only.group.interact.dnds=fixed.only.group.interact.dnds, 
                       fixed.only.thermo.interact.dnds=fixed.only.thermo.interact.dnds,
                       fixed.only.no.interact.dnds=fixed.only.no.interact.dnds,
                       phylo.only.dnds=phylo.only.dnds)


#save(pglmm_model_sel, file="./Outputs/Data/multiple_pglmm_model_selection_mtdna.rdata")
load(file="./Outputs/Data/multiple_pglmm_model_selection_mtdna.rdata")

full.group.interact.ds<-pglmm_model_sel[["full.group.interact.ds"]]
full.thermo.interact.ds<-pglmm_model_sel[["full.thermo.interact.ds"]]
full.no.interact.ds<-pglmm_model_sel[["full.no.interact.ds"]]
fixed.only.group.interact.ds<-pglmm_model_sel[["fixed.only.group.interact.ds"]]
fixed.only.thermo.interact.ds<-pglmm_model_sel[["fixed.only.thermo.interact.ds"]]
fixed.only.no.interact.ds<-pglmm_model_sel[["fixed.only.no.interact.ds"]]
phylo.only.ds<-pglmm_model_sel[["phylo.only.ds"]]

full.group.interact.dn<-pglmm_model_sel[["full.group.interact.dn"]]
full.thermo.interact.dn<-pglmm_model_sel[["full.thermo.interact.dn"]]
full.no.interact.dn<-pglmm_model_sel[["full.no.interact.dn"]]
fixed.only.group.interact.dn<-pglmm_model_sel[["fixed.only.group.interact.dn"]]
fixed.only.thermo.interact.dn<-pglmm_model_sel[["fixed.only.thermo.interact.dn"]]
fixed.only.no.interact.dn<-pglmm_model_sel[["fixed.only.no.interact.dn"]]
phylo.only.dn<-pglmm_model_sel[["phylo.only.dn"]]

full.group.interact.dnds<-pglmm_model_sel[["full.group.interact.dnds"]]
full.thermo.interact.dnds<-pglmm_model_sel[["full.thermo.interact.dnds"]]
full.no.interact.dnds<-pglmm_model_sel[["full.no.interact.dnds"]]
fixed.only.group.interact.dnds<-pglmm_model_sel[["fixed.only.group.interact.dnds"]]
fixed.only.thermo.interact.dnds<-pglmm_model_sel[["fixed.only.thermo.interact.dnds"]]
fixed.only.no.interact.dnds<-pglmm_model_sel[["fixed.only.no.interact.dnds"]]
phylo.only.dnds<-pglmm_model_sel[["phylo.only.dnds"]]


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
dnds <- model.sel(full.group.interact.dnds, full.thermo.interact.dnds, full.no.interact.dnds, fixed.only.group.interact.dnds, 
                  fixed.only.thermo.interact.dnds, fixed.only.no.interact.dnds, phylo.only.dnds,  rank=DIC,
                  extra = c(R2m=function(x) r.squaredMCMCglmm(x)$R2m,
                            R2m.sd=function(x) r.squaredMCMCglmm(x)$R2m.sd,
                            R2c=function(x) r.squaredMCMCglmm(x)$R2c,
                            R2c.sd=function(x) r.squaredMCMCglmm(x)$R2c.sd,
                            R2r=function(x) r.squaredMCMCglmm(x)$R2r,
                            R2r.sd=function(x) r.squaredMCMCglmm(x)$R2r.sd))

model_sel <- list(ds,dn,dnds)
#save(model_sel, file="./Outputs/Data/model_sel_mtdna.rdata")
load(file="./Outputs/Data/model_sel_mtdna.rdata")

ds <- model_sel[[1]]%>%as.data.frame()%>%dplyr::select(R2m:weight)
dn <- model_sel[[2]]%>%as.data.frame()%>%dplyr::select(R2m:weight)
dnds <- model_sel[[3]]%>%as.data.frame()%>%dplyr::select(R2m:weight)
write.csv(rbind(ds,dn,dnds), "model.sel.csv")

df1 <- data.frame(
  model = c("Fixed-only", "Fixed-Interact", "Best Model", "Phylo-only"),
  Mol.Rate = rep(c("dS", "dN", "dNdS"),each=4),
  R2c = c(ds["fixed.only.no.interact.ds","R2c"], 
          ds["fixed.only.group.interact.ds","R2c"], 
          ds["full.group.interact.ds","R2c"], 
          ds["phylo.only.ds", "R2c"],
          dn["fixed.only.no.interact.dn","R2c"], 
          dn["fixed.only.group.interact.dn","R2c"], 
          dn["full.group.interact.dn","R2c"],
          dn["phylo.only.dn", "R2c"],
          dnds["fixed.only.no.interact.dnds","R2c"], 
          dnds["fixed.only.group.interact.dnds","R2c"], 
          dnds["full.group.interact.dnds","R2c"],
          dnds["phylo.only.dnds", "R2c"]),
  
  R2c.sd = c(ds["fixed.only.no.interact.ds","R2c.sd"], 
             ds["fixed.only.group.interact.ds","R2c.sd"], 
             ds["full.group.interact.ds","R2c.sd"], 
             ds["phylo.only.ds", "R2c.sd"],
             dn["fixed.only.no.interact.dn","R2c.sd"], 
             dn["fixed.only.group.interact.dn","R2c.sd"], 
             dn["full.group.interact.dn","R2c.sd"],
             dn["phylo.only.dn", "R2c.sd"],
             dnds["fixed.only.no.interact.dnds","R2c.sd"], 
             dnds["fixed.only.group.interact.dnds","R2c.sd"], 
             dnds["full.group.interact.dnds","R2c.sd"],
             dnds["phylo.only.dnds", "R2c.sd"]))%>%
  mutate(model=factor(model, levels = c("Fixed-only", "Fixed-Interact", "Best Model", "Phylo-only")))%>%
  mutate(Mol.Rate=factor(Mol.Rate, levels=c("dS","dN", "dNdS")))

df2 <- data.frame(R2=rep(c("R2.Phy", "R2m"), 3),
                  Mean=c(ds[1,"R2r"], ds[1,"R2m"], dn[1,"R2r"], dn[1,"R2m"],dnds[1,"R2r"], dnds[1,"R2m"]),
                  SD=c(ds[1,"R2r.sd"], ds[1,"R2m.sd"], dn[1,"R2r.sd"], dn[1,"R2m.sd"],dnds[1,"R2r.sd"], dnds[1,"R2m.sd"]),
                  group=rep(c("dS", "dN", "dNdS"), each=2))%>%
  mutate(group=factor(group, levels=c("dS","dN", "dNdS")))


#Fig.2j-k
f2j <- ggplot(df1, aes(model, R2c, group=Mol.Rate, colour = Mol.Rate))+
  geom_line(size=0.4)+
  geom_point(size=1.2)+
  geom_errorbar(aes(ymin = R2c-R2c.sd, ymax = R2c+R2c.sd), width=.3, size=0.4)+
  scale_color_manual(values=c("#93c47d", "#fdb96b", "#2985cc"))+
  labs(x="", y=expression(R[c] ^  "2" ), title = "Model selection", tag="J")+
  guides(y=guide_axis(cap='upper'), x=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0,1))+
  theme_classic()+
  theme(axis.title = element_text(size=7, color = "black"),
        axis.text.y = element_text(size=7, color = "black"),
        plot.title = element_text(hjust = 0.5, size=8),
        plot.tag = element_text(size=8, face = "bold"),
        legend.text = element_text(size=7),
        legend.position = c(0.30, 0.8),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key.height = unit(0.2, 'cm'),
        legend.key.width = unit(0.4, 'cm'),
        axis.text.x = element_text(size=7, angle = 20, vjust = 1, hjust=0.8, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))


f2k<-df2%>%
  ggplot(aes(x=group, y=Mean, fill=R2)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar(aes(ymin = Mean-SD, ymax = Mean+SD), width=.3, position=position_dodge(0.9), size=0.2)+
  labs(x="", y=expression("R" ^ 2), title = "Best model", tag="K")+
  scale_fill_manual(labels=expression("R"["phy"]^2,"R"["fixed"]^2), values=c("#7570b3", "#e7298a"))+
  guides(y=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(0,1))+
  theme_classic()+
  theme(axis.title = element_text(size=7, color = "black"),
        axis.text = element_text(size=7, color = "black"),
        legend.text = element_text(size=7),
        plot.tag = element_text(size=8, face = "bold"),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.text.align = 0,
        legend.position = c(0.80, 0.9),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key.height = unit(0.1, 'cm'),
        legend.key.width = unit(0.1, 'cm'),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.line = element_line(color = "black", size = 0.2))


######################
summary_multi_pglmm <- function(molrate, pglmm){
  fixed <- summary(pglmm)$solutions
  random <- summary(pglmm)$Gcovariances
  r2 <- r.squaredMCMCglmm(pglmm)
  r2.ci <- matrix(c(r2$R2m, r2$R2m.CIl, r2$R2m.CIu, nrow(pglmm$Sol),
                    r2$R2c, r2$R2c.CIl, r2$R2c.CIu, nrow(pglmm$Sol),
                    r2$R2r, r2$R2r.CIl, r2$R2r.CIu, nrow(pglmm$Sol)), 
                  nrow=3, ncol=4, byrow = T)
  row.names(r2.ci) <- c("Rm2","Rc2","Rr2")
  colnames(r2.ci) <- colnames(random)
  
  pglmm_summary_out <- cbind(data.frame(MolRate=molrate, Var=c(row.names(fixed), row.names(random), row.names(r2.ci))), 
                             data.frame(rbind(fixed, cbind(random, pMCMC=NA), cbind(r2.ci, pMCMC=NA))))
  
  return(pglmm_summary_out)
}

pglmm_summary_out <- rbind(summary_multi_pglmm(molrate="dS", full.thermo.interact.ds),
                           summary_multi_pglmm(molrate="dN", full.thermo.interact.dn),
                           summary_multi_pglmm(molrate="dNdS", full.thermo.interact.dnds))


############################################################
#3.PGLMMs for each class, endotherm and ectotherms using all data
#################################################################
#3.1  Multiple pglmms account for random effects of phylogenetic signal.
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

#save(multi_pglmm_list, file="./Outputs/Data/multiple_pglmm_mtdna.rdata")
load(file="./Outputs/Data/multiple_pglmm_mtdna.rdata")


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


#Fig.2l-n
f2l1 <- pglmm_summary_out%>%
  filter(Var != "(Intercept)", !is.na(pMCMC))%>%
  mutate(ThermoMode=ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"),
         sig=ifelse(pMCMC<0.05, "1", "0"))%>%
  filter(MolRate=="dS")%>%
  filter(Group =="Ectotherms")%>%
  ggplot(aes(x=Var, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=l.95..CI, ymax=u.95..CI), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="dS", title = "Ectotherms", tag="L")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.45,0.21), breaks = seq(-0.4,0.2,0.2))+
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


f2l2 <- pglmm_summary_out%>%
  filter(Var != "(Intercept)", !is.na(pMCMC))%>%
  mutate(ThermoMode=ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"),
         sig=ifelse(pMCMC<0.05, "1", "0"))%>%
  filter(MolRate=="dS")%>%
  filter(Group =="Endotherms")%>%
  ggplot(aes(x=Var, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=l.95..CI, ymax=u.95..CI), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="dS", title = "Endotherms")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.46,0.21), breaks = seq(-0.4,0.2,0.2))+
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



f2m1 <- pglmm_summary_out%>%
  filter(Var != "(Intercept)", !is.na(pMCMC))%>%
  mutate(ThermoMode=ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"),
         sig=ifelse(pMCMC<0.05, "1", "0"))%>%
  filter(MolRate=="dN")%>%
  filter(Group =="Ectotherms")%>%
  ggplot(aes(x=Var, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=l.95..CI, ymax=u.95..CI), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="dN", title = "Ectotherms",tag="M")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.22,0.4), breaks = seq(-0.2,0.4,0.2))+
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


f2m2 <- pglmm_summary_out%>%
  filter(Var != "(Intercept)", !is.na(pMCMC))%>%
  mutate(ThermoMode=ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"),
         sig=ifelse(pMCMC<0.05, "1", "0"))%>%
  filter(MolRate=="dN")%>%
  filter(Group =="Endotherms")%>%
  ggplot(aes(x=Var, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=l.95..CI, ymax=u.95..CI), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="dN", title = "Endotherms")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp", "dS"))+
  scale_y_continuous(limits = c(-0.21,0.4), breaks = seq(-0.2,0.4,0.2))+
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

f2n1 <- pglmm_summary_out%>%
  filter(Var != "(Intercept)", !is.na(pMCMC))%>%
  mutate(ThermoMode=ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"),
         sig=ifelse(pMCMC<0.05, "1", "0"))%>%
  filter(MolRate=="dNdS")%>%
  filter(Group =="Ectotherms")%>%
  ggplot(aes(x=Var, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=l.95..CI, ymax=u.95..CI), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#2a6aaf"))+
  labs(y="Effect Size", x="dN/dS", title = "Ectotherms", tag="N")+
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


f2n2 <- pglmm_summary_out%>%
  filter(Var != "(Intercept)", !is.na(pMCMC))%>%
  mutate(ThermoMode=ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms"),
         sig=ifelse(pMCMC<0.05, "1", "0"))%>%
  filter(MolRate=="dNdS")%>%
  filter(Group =="Endotherms")%>%
  ggplot(aes(x=Var, y=post.mean, colour=ThermoMode))+
  geom_point(aes(shape=sig), size=2.5, show.legend = F)+
  scale_shape_manual(values = c(1,16))+
  geom_errorbar(aes(ymin=l.95..CI, ymax=u.95..CI), width=0, size=0.5)+
  geom_hline(yintercept = 0, linetype=2, colour="grey")+
  scale_color_manual(values = c("#d3292f"))+
  labs(y="Effect Size", x="dN/dS", title = "Endotherms")+
  scale_x_discrete(limits=c("Fecundity","MaturityAge","Longevity", "BodyMass","AnnualTemp"))+
  scale_y_continuous(limits = c(-0.21,0.42), breaks = seq(-0.2,0.4,0.2))+
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



cowplot::plot_grid(f2j,f2k,f2l1,f2l2,f2m1,f2m2,f2n1,f2n2, nrow=2, byrow = F, align = "hv")
#ggsave(filename="./Outputs/MainFigures/Fig2JKLMN.pdf", height=3.8, width=8.27)



#Extended Fig.S11
fs11.1 <- pglmm_summary_out%>%
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

fs11.2 <- pglmm_summary_out%>%
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

fs11.3 <- pglmm_summary_out%>%
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

fs11.4 <- pglmm_summary_out%>%
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

fs11.5 <- pglmm_summary_out%>%
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

fs11.6 <- pglmm_summary_out%>%
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

fs11.7 <- pglmm_summary_out%>%
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

fs11.8 <- pglmm_summary_out%>%
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

fs11.9 <- pglmm_summary_out%>%
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

fs11.10 <- pglmm_summary_out%>%
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

fs11.11 <- pglmm_summary_out%>%
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

fs11.12 <- pglmm_summary_out%>%
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

fs11.13 <- pglmm_summary_out%>%
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

fs11.14 <- pglmm_summary_out%>%
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

fs11.15 <- pglmm_summary_out%>%
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

fs11.1+fs11.2+fs11.3+fs11.4+fs11.5+
  fs11.6+fs11.7+fs11.8+fs11.9+fs11.10+
  fs11.11+fs11.12+fs11.13+fs11.14+fs11.15+
  plot_layout(ncol = 5, nrow = 3)

ggsave(filename="./Outputs/Supplementary/Fig.S11.pdf", height=5.5, width=8.27)



#################################################
#3.2 Multiple PGLMMs account for random effects of phylogenetic signal under the best-fit evolutionary model
load(file="./Outputs/Data/best_trait_evol_models_mtdna.rdata")

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

#save(multi_pglmm_list_best_model, file="./Outputs/Data/multiple_pglmm_best_traits_model_mtdna.rdata")
load(file="./Outputs/Data/multiple_pglmm_best_traits_model_mtdna.rdata")


pglmm_summary_out_best_model <- rbind(summary_multi_pglmm_groups(groups=groups, molrate="dS", multi_pglmm_list_best_model),
                                      summary_multi_pglmm_groups(groups=groups, molrate="dN", multi_pglmm_list_best_model),
                                      summary_multi_pglmm_groups(groups=groups, molrate="dNdS", multi_pglmm_list_best_model))

write.csv(pglmm_summary_out_best_model, "pglmm_summary_out_best_model.csv")


################################################
#Multiple PGLMMa account for randoms effects of phylogenetic signal, habitat and lineage for estimating molecular rates.
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

#save(multi_pglmm_list_3randoms, file="./Outputs/Data/multiple_pglmm_three_random_effects_mtdna.rdata")
load(file="./Outputs/Data/multiple_pglmm_three_random_effects_mtdna.rdata")

pglmm_summary_out_3randoms <- rbind(summary_multi_pglmm_groups(groups=groups, molrate="dS", multi_pglmm_list_3randoms),
                                    summary_multi_pglmm_groups(groups=groups, molrate="dN", multi_pglmm_list_3randoms),
                                    summary_multi_pglmm_groups(groups=groups, molrate="dNdS", multi_pglmm_list_3randoms))

write.csv(pglmm_summary_out_3randoms, "pglmm_summary_out_3randoms.csv")

##############
#Multiple PGLMM using data without missing traits
groups <- c("Ectotherms", "Endotherms", "Fishes", "Amphibians", "Reptiles", "Mammals","Birds")

multi_pglmm_list_no_missing <- vector("list", length = length(groups) * 3)
names(multi_pglmm_list_no_missing) <- paste(rep(groups, each = 3), c("dN", "dS", "dNdS"), sep = ".")


prior <- list(R=list(V=diag(1),nu=0.002), 
              G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*1000)))

for (group in groups) {
  # Choose data based on group
  if(group %in% unique(molrate$ThermoMode)){
    subdata <- molrate%>% 
      filter(ThermoMode == group)%>%
      filter(!is.na(Ref.Fecundity), !is.na(Ref.Longevity),!is.na(Ref.MaturityAge))%>%
      mutate(dS=scale(log(dS)), dN=scale(log(dN)), dNdS=scale(log(dNdS)), AnnualTemp=scale(AnnualTemp),
             BodyMass=scale(log(BodyMass)), Fecundity=scale(log(Fecundity)), 
             MaturityAge=scale(log(MaturityAge)), Longevity=scale(log(Longevity)))
    
  }else{
    subdata <- molrate %>% 
      filter(Group == group)%>%
      filter(!is.na(Ref.Fecundity), !is.na(Ref.Longevity),!is.na(Ref.MaturityAge))%>%
      mutate(dS=scale(log(dS)), dN=scale(log(dN)),dNdS=scale(log(dNdS)), AnnualTemp=scale(AnnualTemp),
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
  fit_dnds <- MCMCglmm(formula(paste("dN~", paste(pred1, collapse="+"))), 
                       random=~Species, prior=prior,
                       ginverse = list(Species=inverseA(phylo)$Ainv),
                       data=subdata,verbose=TRUE, 
                       nitt=60000, burnin=10000, thin=25, family = c("gaussian"))
  
  # Store models in list
  multi_pglmm_list_no_missing[[paste(group, "dS", sep = ".")]] <- fit_ds
  multi_pglmm_list_no_missing[[paste(group, "dN", sep = ".")]] <- fit_dn
  multi_pglmm_list_no_missing[[paste(group, "dNdS", sep = ".")]] <- fit_dnds
  
  print(group)
}

#save(multi_pglmm_list_no_missing, file="./Outputs/Data/mltiple_gplmm_without_missing_traits_mtdna.rdata")
load(file="./Outputs/Data/multiple_pglmm_without_missing_traits_mtdna.rdata")

pglmm_summary_out_no_missing <- rbind(summary_multi_pglmm_groups(groups=groups, molrate="dS", multi_pglmm_list_no_missing),
                                      summary_multi_pglmm_groups(groups=groups, molrate="dN", multi_pglmm_list_no_missing),
                                      summary_multi_pglmm_groups(groups=groups, molrate="dNdS", multi_pglmm_list_no_missing))
#write.csv(pglmm_summary_out_no_missing,"pglmm_summary_out_no_missing.csv")
################


#################################################################
#3.4  Multiple pglmms using data remove relative molecular rates less than 0.4 sub/site to reduce the impacts of substitution saturation
pred1 <- c("AnnualTemp", "BodyMass", "Longevity","MaturityAge","Fecundity")
pred2 <- c("AnnualTemp", "BodyMass", "Longevity","MaturityAge","Fecundity","dS")

prior <- list(R=list(V=diag(1),nu=0.002), 
              G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*1000)))

groups <- c("Ectotherms", "Endotherms", "Fishes", "Amphibians", "Reptiles", "Mammals","Birds")

multi_pglmm_list_0.4 <- vector("list", length = length(groups) * 3)
names(multi_pglmm_list_0.4) <- paste(rep(groups, each = 3), c("dN", "dS", "dNdS"), sep = ".")

# Prior specification for the model
for (group in groups) {
  # Choose data based on group
  if(group %in% unique(molrate$ThermoMode)){
    subdata <- molrate%>% 
      filter(dS*Tip.Age*1000000<=0.4)%>%
      filter(ThermoMode == group)%>%
      mutate(dS=scale(log(dS)), dN=scale(log(dN)),dNdS=scale(log(dNdS)), AnnualTemp=scale(AnnualTemp),
             BodyMass=scale(log(BodyMass)), Fecundity=scale(log(Fecundity)), 
             MaturityAge=scale(log(MaturityAge)), Longevity=scale(log(Longevity)))
    
  }else{
    subdata <- molrate %>% 
      filter(dS*Tip.Age*1000000<=0.4)%>%
      filter(Group == group)%>%
      mutate(dS=scale(log(dS)), dN=scale(log(dN)),dNdS=scale(log(dNdS)), AnnualTemp=scale(AnnualTemp),
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
  multi_pglmm_list_0.4[[paste(group, "dS", sep = ".")]] <- fit_ds
  multi_pglmm_list_0.4[[paste(group, "dN", sep = ".")]] <- fit_dn
  multi_pglmm_list_0.4[[paste(group, "dNdS", sep = ".")]] <- fit_dnds
  print(group)
}

#save(multi_pglmm_list_0.4, file="./Outputs/Data/multiple_pglmm_mtdna_0.4.rdata")
load(file="./Outputs/Data/multiple_pglmm_mtdna_0.4.rdata")



pglmm_summary_out_0.4 <- rbind(summary_multi_pglmm_groups(groups=groups, molrate="dS", multi_pglmm_list_0.4),
                               summary_multi_pglmm_groups(groups=groups, molrate="dN", multi_pglmm_list_0.4),
                               summary_multi_pglmm_groups(groups=groups, molrate="dNdS", multi_pglmm_list_0.4))


#write.csv(pglmm_summary_out_0.4,"pglmm_summary_out_0.4.csv")

##################################################################
#Part IV: Examining relationships between diversification rates and molecular rates
###################################################################

rm(list=ls())
gc()
#define work direction
workdir <- "/Users/tianlong/VertMolRate"
setwd(workdir)

# Load necessary libraries
library(tidyverse)   # For data manipulation and visualization
library(MCMCglmm)    # For mixed models (if required)
library(MuMIn)       # For model selection and multi-model inference
library(patchwork)   # For combining ggplot objects

# Load the data from CSV file
sisters <- read.csv("./DataFiles/sisters_family/sisters_mtdna.csv", row.names = 1)

# Identify unique groups (for example: Fishes, Amphibians, Reptiles, etc.)
groups <- c(unique(sisters$Group))

# Create a list to store model results (one for each group and molecular rate)
ols_list <- vector("list", length = length(groups) * 3)
names(ols_list) <- paste(rep(groups, each = 3), c("dN", "dS", "dNdS"), sep = ".")

# Loop through each group and fit linear models for dN, dS, and dNdS
for(group in groups){
  
  # Filter data by thermo-mode (Endotherms/Ectotherms) or group
  if(group %in% unique(sisters$ThermoMode)){
    subdata <- sisters %>% filter(ThermoMode == group)
  } else {
    subdata <- sisters %>% filter(Group == group)
  }
  
  # Fit linear models for each molecular rate (dN, dS, and dNdS) using clade size differences as predictors
  fit_dn <- lm(diff.dn ~ diff.spp + 0, subdata)  # dN model
  fit_ds <- lm(diff.ds ~ diff.spp + 0, subdata)  # dS model
  fit_dnds <- lm(diff.dnds ~ diff.spp + 0, subdata)  # dNdS model
  
  # Store model results for each molecular rate (dN, dS, dNdS)
  ols_list[[paste(group, "dS", sep = ".")]] <- fit_ds
  ols_list[[paste(group, "dN", sep = ".")]] <- fit_dn
  ols_list[[paste(group, "dNdS", sep = ".")]] <- fit_dnds
}

# Collect and summarize results from the linear models
ols_out <- NULL
for(i in 1:length(groups)){
  group <- groups[i]
  
  # Extract summary statistics for each model (dS, dN, dNdS)
  fit.ds.out <- summary(ols_list[[paste(group, "dS", sep = ".")]])$coefficients
  fit.dn.out <- summary(ols_list[[paste(group, "dN", sep = ".")]])$coefficients
  fit.dnds.out <- summary(ols_list[[paste(group, "dNdS", sep = ".")]])$coefficients
  
  # Combine the results into a data frame for easier interpretation
  fit.out <- data.frame(MolRate = c("dS", "dN", "dNdS"), 
                        rbind(fit.ds.out, fit.dn.out, fit.dnds.out)) %>%
    rename(Slope = Estimate, Slope.SD = Std..Error, t = t.value, p = Pr...t..) %>%
    mutate(Group = group, 
           ThermoMode = ifelse(Group %in% c("Birds", "Mammals", "Endotherms"), "Endotherms", "Ectotherms")) %>%
    select(ThermoMode, Group, MolRate, Slope, Slope.SD, t, p)
  
  # Append the results to the final output
  ols_out <- rbind(ols_out, fit.out)
}


#Extended Fig.S12
f1 <- sisters%>%
  filter(Group=="Fishes")%>%
  ggplot(aes(x=diff.spp, y=diff.ds, colour=Group))+
  geom_point(size=1)+
  geom_smooth(method="lm", formula = y~x+0)+
  geom_errorbar(aes(ymin=diff.ds-diff.ds.sd, ymax=diff.ds+diff.ds.sd), width=.2, linewidth=0.2)+
  scale_color_manual(values = c("#2a6aaf"))+
  scale_x_continuous(limits = c(0,6), breaks = seq(0,6,2))+
  scale_y_continuous(limits = c(-2,2), breaks = seq(-2,2,1))+
  labs(x="Diff. in ln(Clade Size)", y="Diff. in ln(dS)", title="Fishes")+
  annotate("text", x=4,y=2, label="p = 0.966", size=2.5)+
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
  scale_y_continuous(limits = c(-1,0.5), breaks = seq(-1,0.5,0.5))+
  labs(x="Diff. in ln(Clade Size)", y="", title="Amphibians")+
  annotate("text", x=4,y=0.5, label="p = 0.167", size=2.5)+
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
  scale_x_continuous(limits = c(0,4), breaks = seq(0,4,1))+
  scale_y_continuous(limits = c(-1.5,1), breaks = seq(-1.5,1,0.5))+
  labs(x="Diff. in ln(Clade Size)", y="", title="Reptiles")+
  annotate("text", x=3,y=1, label="p = 0.052", size=2.5)+
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
  scale_y_continuous(limits = c(-1,1), breaks = seq(-1,1,0.5))+
  labs(x="Diff. in ln(Clade Size)", y="", title="Mammals")+
  annotate("text", x=3,y=1, label="p = 0.927", size=2.5)+
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
  annotate("text", x=4,y=1.2, label="p = 0.007", size=2.5)+
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
  scale_y_continuous(limits = c(-2,2), breaks = seq(-2,2,1))+
  labs(x="Diff. in ln(Clade Size)", y="Diff. in ln(dN)", title="Fishes")+
  annotate("text", x=4,y=2, label="p = 0.699", size=2.5)+
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
  scale_y_continuous(limits = c(-1,0.5), breaks = seq(-1,0.5,0.5))+
  labs(x="Diff. in ln(Clade Size)", y="", title="Amphibians")+
  annotate("text", x=4,y=0.5, label="p = 0.389", size=2.5)+
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
  scale_x_continuous(limits = c(0,4), breaks = seq(0,4,1))+
  scale_y_continuous(limits = c(-0.5,1), breaks = seq(-0.5,1,0.5))+
  labs(x="Diff. in ln(Clade Size)", y="", title="Reptiles")+
  annotate("text", x=3,y=1, label="p = 0.602", size=2.5)+
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
  scale_y_continuous(limits = c(-1.5,1.5), breaks = seq(-1.5,1.5,0.5))+
  labs(x="Diff. in ln(Clade Size)", y="", title="Mammals")+
  annotate("text", x=3,y=1.5, label="p = 0.451", size=2.5)+
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
  scale_y_continuous(limits = c(-0.8,1.2), breaks = seq(-0.8,1.2,0.4))+
  labs(x="Diff. in ln(Clade Size)", y="", title="Birds")+
  annotate("text", x=4,y=1.2, label="p = 0.392", size=2.5)+
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
  #scale_y_continuous(limits = c(-0.6,0.6), breaks = seq(-0.6,0.6,0.3))+
  labs(x="Diff. in ln(Clade Size)", y="Diff. in ln(dN/dS)", title="Fishes")+
  annotate("text", x=4,y=1, label="p = 0.637", size=2.5)+
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
  #scale_y_continuous(limits = c(-0.4,0.4), breaks = seq(-0.4,0.4,0.2))+
  labs(x="Diff. in ln(Clade Size)", y="", title="Amphibians")+
  annotate("text", x=4,y=0.5, label="p = 0.835", size=2.5)+
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
  #scale_y_continuous(limits = c(-0.4,0.22), breaks = round(seq(-0.4,0.2,0.2),2))+
  labs(x="Diff. in ln(Clade Size)", y="", title="Reptiles")+
  annotate("text", x=3,y=1, label="p = 0.179", size=2.5)+
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
  #scale_y_continuous(limits = c(-0.4,0.2), breaks = seq(-0.4,0.2,0.2))+
  labs(x="Diff. in ln(Clade Size)", y="", title="Mammals")+
  annotate("text", x=3,y=1, label="p = 0.536", size=2.5)+
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
  #scale_y_continuous(limits = c(-0.3,0.3), breaks = round(seq(-0.3,0.3,0.1),2))+
  labs(x="Diff. in ln(Clade Size)", y="", title="Birds")+
  annotate("text", x=4,y=1, label="p = 0.217", size=2.5)+
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

ggsave(filename="./Outputs/Supplementary/Fig.S12.pdf", height=5, width=8.27)



###############################################
#Part V: 
#############################################################
library(tidyverse)
library(nlme)
library(ape)
library(phytools)
library(agricolae)
######
workdir <- "/Users/Tianlong/VertMolRate"
molrate <- read.csv(paste0(workdir, "/DataFiles/MolEvolRate/molrate_mtDNA.csv"))

# Load input phylogeny
phylo <- read.tree(paste0(workdir, "/DataFiles/trees/phy_all_sampled_mtdna.tre"))

# Check whether the phylogeny is ultrametric
is.ultrametric(phylo)

# Force the tree to be ultrametric
phylo <- force.ultrametric(phylo, method = "nnls")  # Use nnls or extend


##################################################
#5.1 Compare dS of migrants and residents
migrant <- molrate%>%
  filter(Group=="Birds")%>%
  mutate(Migration2=ifelse(Migration=="Resident", "Non-migrants", "Migrants"))

#dS corrected for body mass
pgls <- gls(log(dS)~log(BodyMass), migrant, correlation = corMartins(value=0.017, phy=drop.tip(phylo,setdiff(phylo$tip.label,migrant$Species)), form=~Species))
summary(pgls)

#data for plot
ds.corrected.mass <- data.frame(dS.Residual=residuals(pgls), Migration=migrant$Migration, Migration2=migrant$Migration2)
kruskal(ds.corrected.mass$dS.Residual, ds.corrected.mass$Migration, p.adj = "bonferroni")$groups

my_comparisons1 <- list( c("Long Migratory", "Short Migratory"), c("Short Migratory", "Resident"), c("Long Migratory", "Resident"))

#Extended Fig.S16
fs16a <- ds.corrected.mass%>%
  mutate(Migration=factor(Migration, levels = c("Long Migratory", "Short Migratory", "Resident")))%>%
  ggplot(aes(x=Migration, y=dS.Residual, colour=Migration))+
  geom_boxplot()+
  theme_classic()+
  labs(x="", y=expression(dS[Residuals]), tag="a")+
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

fs16b <- ds.corrected.mass%>%
  ggplot(aes(x=Migration2, y=dS.Residual, colour=Migration2))+
  geom_boxplot()+
  theme_classic()+
  labs(x="", y=expression(dS[Residuals]), tag="b")+
  guides(y=guide_axis(cap='upper'))+
  scale_y_continuous(limits = c(-1.3,2.0), breaks = seq(-1.5,2.0,0.5))+
  ggpubr::stat_compare_means(comparisons = my_comparisons2, method="wilcox.test", size=2.5)+ 
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

##############################################################
#5.2. Accounting for mass-specific BMR by body mass, body temperature, migratory status and phylogenetic relatedness
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

fs16c <- endo.res%>%
  ggplot(aes(x=Ta, y=BMR.res,colour=Class))+
  geom_point(size=1)+
  scale_shape_manual(values=c(1,4))+
  geom_smooth(method="lm", size=0.8)+
  scale_color_manual(values=c("#d3292f"))+
  labs(x=expression("Mean Annual Temperature" ~ ( degree * C)), 
       y=expression("Mass specific BMR"[residuals]),tag="c")+
  annotate("text", x = 0, y = -0.5, label = expression(italic(P)[OLS] < 0.001), size=2.5)+
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

cowplot::plot_grid(fs16a,fs16b, fs16c, nrow=1, align="hv")
ggsave(filename="./Outputs/Supplementary/Extended Fig.S16.pdf", height=2.8, width=8.27)




