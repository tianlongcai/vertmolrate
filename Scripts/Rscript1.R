################################
# R script to analyze variation of molecular rates across species and groups
# Author: Tianlong Cai
# Email: caitianlong@westlake.edu.cn


############################################
rm(list=ls())
gc()
workdir <- "/Users/tianlong/Desktop/VertMolRate"
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
####################################################################
#functions
source(paste0(workdir, "/SourceFunctions/source_functions.r"))
plot2.phylo <- diversitree:::plot2.phylo

#Input phylogeny
all_vert_phy <- read.tree(paste0(workdir, "/DataFiles/trees/phy_all_sampled_vertebrates_outgroup.tre"))

#Input substitution rates
subrate <- read.csv(paste0(workdir, "/DataFiles/MolEvolRate/VertMolRate.csv"))

#Check matched tip names in substitution rates
row.names(subrate) <- subrate$Species
unmatched_species <- name.check(all_vert_phy, subrate)$tree_not_data

#Dropping unmatched tips
sampled_vert_phy<-drop.tip(all_vert_phy, unmatched_species[-which(unmatched_species=="Outgroup")])


subrate<-subrate[sampled_vert_phy$tip.label[-which(sampled_vert_phy$tip.label=="Outgroup")],]%>%
  mutate(MajorGroup=toupper(MajorGroup))

#################################################
#Extract node for each clade
clades <- unique(subrate$Clade.Label)
clade_nodes <- matrix(NA, length(clades), 2)
colnames(clade_nodes) <- c("Clade", "Clade.Node")

for(i in 1:length(clades)){
  clade <- clades[i]
  species <- subrate%>%filter(Clade.Label==clade)%>%pull(Species)
  node <- getMRCA(sampled_vert_phy, species)
  clade_nodes[i,1] <- clade
  clade_nodes[i,2] <- node
}
clade_nodes <- as.data.frame(clade_nodes)
clade_nodes$Clade.Node <- as.numeric(clade_nodes$Clade.Node)


#Extract node for each major group
groups <- unique(subrate$MajorGroup)
group_nodes <- matrix(NA, length(groups), 2)
colnames(group_nodes) <- c("MajorGroup", "Group.Node")

for(i in 1:length(groups)){
  group <- groups[i]
  species <- subrate%>%filter(MajorGroup==group)%>%pull(Species)
  node <- getMRCA(sampled_vert_phy, species)
  group_nodes[i,1] <- group
  group_nodes[i,2] <- node
}
group_nodes <- as.data.frame(group_nodes)
group_nodes$Group.Node <- as.numeric(group_nodes$Group.Node)

#Merged clade nodes and groups nodes into subrate
subrate <- subrate%>%
  inner_join(clade_nodes, by=c("Clade.Label"="Clade"))%>%
  inner_join(group_nodes, by="MajorGroup")

#scale dN and dS to plot bars
subrate$dS <- (subrate$dS/max(subrate$dS))*60
subrate$dN <- (subrate$dN/max(subrate$dN))*60


#define colors of subrate bars
subrate_bar_colour <- c("#93c47d", "#fdb96b")


#Extracted subrate to plotting bars
subrate_plot <- as.matrix(subrate[,c("dS","dN")])
row.names(subrate_plot) <- subrate$Species


pdf("./Outputs/MainFigures/Fig1.pdf", width = 8.27, height = 8.27)
#plotting phylogeny using different colors to show diversification rate
par(xpd = TRUE)
plot2.phylo(sampled_vert_phy, show.tip.label=FALSE,  cex=0.05, 
            label.offset=0, type="f", edge.width=0.2, no.margin=FALSE, 
            root.edge=TRUE, x.lim=c(-860, 860))
for (j in 1:2) ring(subrate_plot[, j], sampled_vert_phy, offset = j*60 - 60 +2, col = subrate_bar_colour[j])
legend(110, 190, legend=c("dS", "dN"),bty='n',
       col=subrate_bar_colour, pch=15, cex=0.4)

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

for(i in 1:7) {
  arc.cladelabels(text=group_nodes$MajorGroup[i],node=group_nodes$Group.Node[i], lab.offset=1.33, ln.offset=1.33, col = group.col[i], lwd=8, point.cex=0, cex=0.4, text.color = "white")
}

for(i in 1:length(clades)) {
  arc.cladelabels(text=clade.labels[i],node=nodes[i], lab.offset=1.36, ln.offset=1.36, col = clade.col[i], lwd=3, point.cex=0.3, point.lwd=0.2, cex=0.4)
} 

#Time scale
lines<-rev(seq(0,450,50))
for (i in 1:length(lines)){
  draw.circle(0,0, radius=lines[i],col=NA, lty=2, lwd=0.5, border=grey(0.8,alpha=0.6))
}

text(y=-9, x=seq(0,450,50), labels=rev(seq(0,450,50)), cex=0.3)

dev.off()


#################################################################
#Extended Fig1. Raincloud plot to show molecular rates of five classes of vertebrates
#Kruskal Wallis test and multiple comparison of groups.
#input data
subrate <- read.csv(paste0(workdir, "/DataFiles/MolEvolRate/VertMolRate.csv"))%>%
  mutate(Group=factor(Group, levels = c("Fishes","Amphibians", "Reptiles","Mammals","Birds")))%>%
  mutate(dS=dS*10^8, dN=dN*10^10)

range(subrate$dS)
range(subrate$dN)

comp.ds <- agricolae::kruskal(subrate$dS, subrate$Group, p.adj = "bonferroni")
comp.dn <- agricolae::kruskal(subrate$dN, subrate$Group, p.adj = "bonferroni")

group <- c("Fishes","Amphibians", "Reptiles","Mammals","Birds")
y1 <- subrate%>%group_by(Group)%>%summarise(y=max(dS))%>%pull(y)+0.1
y2 <- subrate%>%group_by(Group)%>%summarise(y=max(dN))%>%pull(y)+0.2

label <- data.frame(x=(1:5)+0.2, 
                 y1=y1,
                 y2=y2,
                 label.ds=comp.ds$groups[match(group,row.names(comp.ds$groups)),"groups"],
                 label.dn=comp.dn$groups[match(group,row.names(comp.dn$groups)),"groups"],
                 ThermoMode = c("Ectotherms", "Ectotherms", "Ectotherms", "Endotherms","Endotherms"))
row.names(label) <- group


f1 <- ggplot(data = subrate, aes(y = dS, x = Group, colour=ThermoMode, fill=ThermoMode)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), adjust =1) +
  geom_point(aes(y = dS), position = position_jitter(width = 0.12), alpha = 0.1, size = 0.5, shape=1) +
  geom_boxplot(width = 0.2, outlier.shape = NA, lwd=0.5, fill = NA) +
  scale_fill_manual(values = c("#2a6aaf", "#d3292f"))+
  scale_color_manual(values = c("#2a6aaf", "#d3292f"))+
  scale_x_discrete(limits=c("Fishes","Amphibians", "Reptiles","Mammals","Birds"))+
  scale_y_continuous(limits = c(0,5), breaks = seq(0,5,1))+
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



f2 <-ggplot(data = subrate, aes(y = dN, x = Group, colour=ThermoMode, fill=ThermoMode)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), adjust =1) +
  geom_point(aes(y = dN), position = position_jitter(width = 0.12), alpha = 0.1, size = 0.5, shape=1) +
  geom_boxplot(width = .2, outlier.shape = NA, lwd=0.5, fill=NA) +
  scale_fill_manual(values = c("#2a6aaf", "#d3292f"))+
  scale_color_manual(values = c("#2a6aaf", "#d3292f"))+
  scale_x_discrete(limits=c("Fishes","Amphibians", "Reptiles","Mammals","Birds"))+
  theme_classic()+
  labs(x = "", y = expression("dN (" * 10^-10 * " sub/site/year)"))+
  scale_y_continuous(limits = c(0,12), breaks = seq(0,12,3))+
  guides(x=guide_axis(cap='upper'), y=guide_axis(cap='upper'))+
  theme(axis.title = element_text(size=9),
        axis.text = element_text(size = 9),
        legend.position = 'none',
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))+
  geom_text(data=label, aes(x=x, y=y2, label=label.ds), size=2.5)+
  coord_flip()

cowplot::plot_grid(f1,f2, nrow=1, align = "hv")
ggsave("./Outputs/Supplementary/Extended Fig1.pdf", width = 8.3, height = 4)


######################################################
#Extended Fig.4 Phylogenetic signal test of dS and dN for each class of vertebrates
#input phylogeny
phy <- read.tree(paste0(workdir, "/DataFiles/trees/phy_all_sampled_vertebrates.tre"))
#force ultrametric phylogeny
phy <- force.ultrametric(phy, method="nnls")#nnls or extend

groups <- unique(subrate$Group)

phylo.signal <- matrix(NA, 5, 4)
colnames(phylo.signal) <- c("ds.lambda", "ds.p", "dn.lambda", "dn.p")
row.names(phylo.signal) <- groups

for(i in 1:length(groups)){
  fgroup <- groups[i]
  
  #extract data
  fdat <- subrate%>%filter(Group==fgroup)
  
  row.names(fdat) <- fdat$Species
  
  #trait data with corresponding tree tips
  ds <- fdat$dS; names(ds)<-fdat$Species
  dn <- fdat$dN; names(dn)<-fdat$Species
  
  #phylogeny corresponding to trait data
  tree <- drop.tip(phy, name.check(phy, fdat)$tree_not_data)
  
  #compute phylogenetic signal using Pagel's lambda
  lambda.ds <- phylosig(tree, ds, method="lambda", test=TRUE)
  lambda.dn <- phylosig(tree, dn, method="lambda", test=TRUE)
  
  #matrix of lambda and p values
  phylo.signal[i,] <- c(lambda.ds$lambda, lambda.ds$P, lambda.dn$lambda, lambda.dn$P)
  print(fgroup)
}

#save data
save(phylo.signal, file = "./Outputs/Data/Phylogenetic_signal.rdata")
load(file = "./Outputs/Data/Phylogenetic_signal.rdata")

phylo.signal <- as.data.frame(phylo.signal)
phylo.signal$Group <- row.names(phylo.signal)


phylo.signal <- phylo.signal%>%
  mutate(sig.dn=ifelse(dn.p<0.001,"***", ifelse(dn.p<0.01, "**", ifelse(dn.p<0.05, "*", "ns"))),
         sig.ds=ifelse(ds.p<0.001,"***", ifelse(ds.p<0.01, "**", ifelse(ds.p<0.05, "*", "ns"))),
         Group=factor(Group, levels = c("Fishes","Amphibians", "Reptiles","Mammals","Birds")),
         ThermoMode=ifelse(Group %in% c("Birds","Mammals"), "Endotherms", "Ectotherms"))%>%
  rename(dS=ds.lambda, dN=dn.lambda)%>%
  select(ThermoMode, Group, dS, dN, sig.dn, sig.ds)

#Extended Fig.2 Phylogenetic signals of dN and dS
f3<-phylo.signal%>%
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
        axis.text = element_text(size=9),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))


f4<-phylo.signal%>%
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
        axis.text = element_text(size=9),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.2),
        axis.line = element_line(color = "black", linewidth = 0.2))
plot_grid(f3, f4, ncol = 2, align = "hv")
ggsave("./Outputs/Supplementary/Extended Fig4.pdf", width = 8.3, height = 4)


group1 <- subrate%>%
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
    dS_95CI_Upper = unlist(t.test(dS)[4])[2])%>%
  rename(Group=ThermoMode)

group2 <- subrate%>%
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
    dS_95CI_Upper = unlist(t.test(dS)[4])[2])

bas_stat <- rbind(group2, group1)
write.csv(bas_stat, "bas_stat.csv")



