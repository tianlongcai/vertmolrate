#modified "arc.cladelabels()" function in the phytools package to plot clade labels in fan tree
arc.cladelabels<-function(tree=NULL,text,node,ln.offset=1.02, type="fan", lwd=5, col="grey", point.cex=1, text.color="black",
                           point.pch=21, point.lwd=0.5, lab.offset=1.06,cex=1,orientation="curved",...){
  obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  if(obj$type!="fan") stop("method works only for type=\"fan\"")
  h<-max(sqrt(obj$xx^2+obj$yy^2))
  if(hasArg(mark.node)) mark.node<-list(...)$mark.node
  else mark.node<-TRUE
  if(mark.node) points(obj$xx[node],obj$yy[node],pch=point.pch,cex=point.cex,lwd=point.lwd,
                       bg="red")
  if(is.null(tree)){
    tree<-list(edge=obj$edge,tip.label=1:obj$Ntip,
               Nnode=obj$Nnode)
    class(tree)<-"phylo"
  }
  d<-getDescendants(tree,node)
  d<-sort(d[d<=Ntip(tree)])
  deg<-atan(obj$yy[d]/obj$xx[d])*180/pi
  ii<-intersect(which(obj$yy[d]>=0),which(obj$xx[d]<0))
  deg[ii]<-180+deg[ii]
  ii<-intersect(which(obj$yy[d]<0),which(obj$xx[d]<0))
  deg[ii]<-180+deg[ii]
  ii<-intersect(which(obj$yy[d]<0),which(obj$xx[d]>=0))
  deg[ii]<-360+deg[ii]
  draw.arc(x=0,y=0,radius=ln.offset*h,deg1=min(deg),
           deg2=max(deg), lwd = lwd, col=col,lend=1)
  if(orientation=="curved")
    arctext(text,radius=lab.offset*h,
            middle=mean(range(deg*pi/180)),cex=cex, col=text.color)
  else if(orientation=="horizontal"){
    x0<-lab.offset*cos(median(deg)*pi/180)*h
    y0<-lab.offset*sin(median(deg)*pi/180)*h
    text(x=x0,y=y0,label=text,
         adj=c(if(x0>=0) 0 else 1,if(y0>=0) 0 else 1),
         offset=0)
  }
}

#histgram of speciation rate with gradient colours
color.legend <- function (xl, yb, xr, yt, legend, rect.col, cex = 1, align = "lt", 
gradient = "x", border="grey", lwd=0.1,...) {

	oldcex <- par("cex")
	par(xpd = TRUE, cex = cex)
	gradient.rect(xl, yb, xr, yt, col = rect.col, nslices = length(rect.col), gradient = gradient, border=border, lwd=lwd)
	if (gradient == "x") {
		xsqueeze <- (xr - xl)/(2 * length(rect.col))
		textx <- seq(xl + xsqueeze, xr - xsqueeze, length.out = length(legend))
			if (match(align, "rb", 0)) {
				texty <- yb - 0.2 * strheight("O")
				textadj <- c(0.5, 1)
				}
			else {
				texty <- yt + 0.2 * strheight("O")
				textadj <- c(0.5, 0)
			}
		}
	else {
	ysqueeze <- (yt - yb)/(2 * length(rect.col))
	texty <- seq(yb + ysqueeze, yt - ysqueeze, length.out = length(legend))
	if (match(align, "rb", 0)) {
	textx <- xr + 0.2 * strwidth("O")
	textadj <- c(0, 0.5)
	}
	else {
	textx <- xl - 0.2 * strwidth("O")
	textadj <- c(1, 0.5)
	}
	}
	text(textx, texty, labels = legend, adj = textadj, ...)
	par(xpd = FALSE, cex = oldcex)
}

gradient.rect <- function (xleft, ybottom, xright, ytop, reds, greens, blues, col = NULL, nslices = 50, gradient = "x", border = par("fg"), lwd) {
	if (is.null(col)) 
	col <- color.gradient(reds, greens, blues, nslices)
	else nslices <- length(col)
	nrect <- max(unlist(lapply(list(xleft, ybottom, xright, ytop), 
	length)))
	if (nrect > 1) {
		if (length(xleft) < nrect) 
		xleft <- rep(xleft, length.out = nrect)
		if (length(ybottom) < nrect) 
		ybottom <- rep(ybottom, length.out = nrect)
		if (length(xright) < nrect) 
		xright <- rep(xright, length.out = nrect)
		if (length(ytop) < nrect) 
		ytop <- rep(ytop, length.out = nrect)
		for (i in 1:nrect) gradient.rect(xleft[i], ybottom[i], 
		xright[i], ytop[i], reds, greens, blues, col, nslices, 
		gradient, border = border)
		}
		else {
			if (gradient == "x") {
				xinc <- (xright - xleft)/nslices
				xlefts <- seq(xleft, xright - xinc, length = nslices)
				xrights <- xlefts + xinc
				rect(xlefts, ybottom, xrights, ytop, col = col, lty = 0, lwd=lwd)
				rect(xlefts[1], ybottom, xrights[nslices], ytop, 
				border = border, lwd=lwd)
				}
				else {
					yinc <- (ytop - ybottom)/nslices
					ybottoms <- seq(ybottom, ytop - yinc, length = nslices)
					ytops <- ybottoms + yinc
					rect(xleft, ybottoms, xright, ytops, col = col, lty = 0, lwd=lwd)
					rect(xleft, ybottoms[1], xright, ytops[nslices], border = border, lwd=lwd)
					}
				}
			invisible(col)
}


#Build Function to Return x and y axises label in central position in ggplot2 plots
rotatedAxisElementText = function(angle,position='x',size){
  angle = angle[1]; 
  position = position[1]
  positions = list(x=0,y=90,top=180,right=270)
  if(!position %in% names(positions))
    stop(sprintf("'position' must be one of [%s]",paste(names(positions),collapse=", ")),call.=FALSE)
  if(!is.numeric(angle))
    stop("'angle' must be numeric",call.=FALSE)
  rads  = (angle - positions[[ position ]])*pi/180
  hjust = 0.5*(1 - sin(rads))
  vjust = 0.5*(1 + cos(rads))
  element_text(angle=angle,vjust=vjust,hjust=hjust,size=size)
}


add.outgroup <- function(phy, root.node=NULL){ #from Sergio Sanchez https://stat.ethz.ch/pipermail/r-sig-phylo/2015-February/003886.html
  require(ape)
  if (is.null(phy$edge.length)){
    tt <- rtree(2)
    tt$edge.length <- NULL
    tt$tip.label <- c("outgroup","drop")
    phy$root.edge <- 1
    if (!is.null(root.node)){
      rphy <- root(phy, node=root.node)
      ot <- bind.tree(rphy, tt, position=1)
    } else {
      ot <- bind.tree(phy, tt, position=1)
    }
    ot <- bind.tree(phy, tt, position=1)
    res <- drop.tip(ot, "drop")
    return(res)
  } else if (is.ultrametric(phy)){
    th <- as.numeric(sort(branching.times(phy), decreasing=TRUE))[1]
    re <- th/5
    phy$root.edge <- re
    tt <- rtree(2)
    tt$edge.length <- c(0,0)
    tt$tip.label <- c("outgroup","drop")
    tt$root.edge <- th + re
    ot <- bind.tree(phy, tt, position=re)
    res <- drop.tip(ot, "drop")
    return(res)
  } else {
    tl <- max(phy$edge.length)
    re <- tl/3
    tt <- rtree(2)
    tt$edge.length <- c(0,0)
    tt$tip.label <- c("outgroup","drop")
    tt$root.edge <- tl
    if (!is.null(root.node)){
      rphy <- root(phy, node=root.node)
      rphy$root.edge <- re
      ot <- bind.tree(rphy, tt, position=re)
    } else {
      phy$root.edge <- re
      ot <- bind.tree(phy, tt, position=re)
    }
    res <- drop.tip(ot, "drop")
    return(res)
  }
}




r.squaredMCMCglmm <- function(fit){
  #number of VCVs
  NoVCV <- length(colnames(fit$VCV))
  
  #alternative with credible intervals
  #Number of samples
  MCMCSamples <- length(fit$VCV[,1])
  
  #Variances of fixed effects for all intervals
  vmVarFixed<-numeric(MCMCSamples)
  for(i in 1:MCMCSamples){
    Var<-var(as.vector(fit$Sol[i,] %*% t(fit$X)))
    vmVarFixed[i]<-Var
  }
  
  #Variances of random effects for all intervals
  if((NoVCV-1)==0){
    vmVarRandom <- 0
  } else{
    if((NoVCV-1)==1){
      vmVarRandom <-fit$VCV[,1]
    } else{
      vmVarRandom <-fit$VCV[,1]
      for(k in 2:(NoVCV-1)){
        vmVarRandom <- vmVarRandom + fit$VCV[,k]
      }
    }
  }
  
  #total variances for all intervals
  tVar <- vmVarFixed
  for(j in 1:NoVCV){
    tVar <- tVar+ fit$VCV[,j]
  }
  
  #######################
  #Marginal r squares
  R2ms <- vmVarFixed/tVar
  R2m <- round(mean(R2ms), 4)
  R2m.CI <- c(round(HPDinterval(R2ms), 4))
  R2m.sd <- sd(R2ms)
  
  #conditional r squares
  R2cs<-(vmVarFixed+vmVarRandom)/tVar
  R2c <- round(mean(R2cs), 4)
  R2c.CI <- c(round(HPDinterval(R2cs),4))
  R2c.sd <- sd(R2cs)
  
  #r square for random effects
  R2rs <- R2cs-R2ms
  R2r <- round(mean(R2rs), 4)
  R2r.CI <- c(round(HPDinterval(R2rs),4))
  R2r.sd <- sd(R2rs)
  
  
  r2 <- matrix(dimnames = list(NULL, c("R2m","R2m.CIl","R2m.CIu","R2m.sd","R2c","R2c.CIl","R2c.CIu","R2c.sd","R2r","R2r.CIl","R2r.CIu","R2r.sd")),
               c(R2m,R2m.CI,R2m.sd,R2c,R2c.CI, R2c.sd, R2r, R2r.CI, R2r.sd), ncol = 12, nrow = 1,byrow = T)
  r2 <- as.data.frame(r2)
  return(r2)
}

r.squaredMCMCglmm.phy <- function(fit, phylo.var=NULL){
  #number of VCVs
  NoVCV <- length(colnames(fit$VCV))
  
  #alternative with credible intervals
  #Number of samples
  MCMCSamples <- length(fit$VCV[,1])
  
  #Variances of fixed effects for all intervals
  vmVarFixed<-numeric(MCMCSamples)
  for(i in 1:MCMCSamples){
    Var<-var(as.vector(fit$Sol[i,] %*% t(fit$X)))
    vmVarFixed[i]<-Var
  }
  
  #Variances of random effects for all intervals
  if((NoVCV-1)==0){
    vmVarRandom <- 0
  } else{
    if((NoVCV-1)==1){
      vmVarRandom <-fit$VCV[,1]
    } else{
      vmVarRandom <-fit$VCV[,1]
      for(k in 2:(NoVCV-1)){
        vmVarRandom <- vmVarRandom + fit$VCV[,k]
      }
    }
  }
  
  #Variances of phylogeny for all intervals
  if(phylo.var %in% names(fit$VCV[1,])){
    if(length(phylo.var)==1){
      phyloVar <- fit$VCV[,phylo.var]
    }else{
      phyloVar <- rowSums(fit$VCV[,phylo.var])
    }
    
  }else{
    phyloVar <-0
  }
  
  
  #total variances for all intervals
  tVar <- vmVarFixed
  for(j in 1:NoVCV){
    tVar <- tVar+ fit$VCV[,j]
  }
  
  #######################
  #Marginal r squares
  R2ms <- vmVarFixed/tVar
  R2m <- round(median(R2ms), 4)
  R2m.CI <- c(round(HPDinterval(R2ms), 4))
  
  #conditional r squares
  R2cs<-(vmVarFixed+vmVarRandom)/tVar
  R2c <- round(median(R2cs), 4)
  R2c.CI <- c(round(HPDinterval(R2cs),4))
  
  #r square for random effects
  R2rs <- R2cs-R2ms
  R2r <- round(median(R2rs), 4)
  R2r.CI <- c(round(HPDinterval(R2rs),4))
  
  #r for phylogeny
  R2ps <- phyloVar/tVar
  R2p <- round(median(R2ps), 4)
  R2p.CI <- c(round(HPDinterval(R2ps),4))
  
  r2 <- matrix(dimnames = list(NULL, c("R2m","R2m.CIl","R2m.CIu","R2c","R2c.CIl","R2c.CIu","R2r","R2r.CIl","R2r.CIu", "R2p","R2p.CIl","R2p.CIu")),
               c(R2m,R2m.CI,R2c,R2c.CI, R2r, R2r.CI, R2p, R2p.CI), ncol = 12, nrow = 1,byrow = T)
  r2 <- as.data.frame(r2)
  return(r2)
}

#a function to put all summary of PGLMM to table
summary_table_pglmm <- function(fit){
  x <- summary(fit)
  r2 <- r.squaredMCMCglmm(fit)
  sum_table <- x$solutions%>%
    as.data.frame()%>%
    dplyr::select(-eff.samp)
  
  row.names(sum_table) <- c("Intercept", paste0("Fixed effect: ", row.names(sum_table)[-1]))
  
  G <- x$Gcovariances%>%
    as.data.frame()%>%
    mutate(pMCMC=NA)%>%
    dplyr::select(-eff.samp)
  row.names(G) <- paste0("Random effect: ",row.names(G))
  R <- x$Rcovariances%>%
    as.data.frame()%>%
    dplyr::select(-eff.samp)%>%
    mutate(pMCMC=NA)
  row.names(R) <- "Residual"
  
  MR2 <-rbind(c(as.numeric(r2[1:3]),NA), c(as.numeric(r2[4:6]),NA))%>%
    as.data.frame()%>%
    rename(post.mean=V1, `l-95% CI`=V2, `u-95% CI`=V3, pMCMC=V4)
  
  row.names(MR2) <- c("R2.m", "R2.c")
  
  DIC <- data.frame(post.mean=x$DIC)%>%
    mutate(`l-95% CI`=NA, `u-95% CI`=NA, pMCMC=NA)
  row.names(DIC) <- "DIC"
  
  sum_table <- rbind(sum_table, G, R, MR2, DIC)
  return(sum_table)
}


vif.MCMCglmm <- function (fit, intercept.columns = c(1)) {
  nF <- fit$Fixed$nfl
  v <- cov(as.matrix(fit$X[,1:nF]))
  nam <- colnames(fit$Sol[,1:nF])
  v <- v[-intercept.columns, -intercept.columns, drop = FALSE]
  nam <- nam[-intercept.columns]
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}




require("tidyverse")
#a function to transform extreme values of substitution rate
trans.extr.val <- function(x) {
  x.max <- quantile(x, 0.99, na.rm=T)
  x.min <- quantile(x, 0.01, na.rm=T)
  x <- ifelse(x>x.max, x.max, ifelse(x<x.min, x.min, x))
  return(x)
}

trans.extr <- function(x) {
  x.max <- quantile(x, 0.99, na.rm=T)
  x.min <- quantile(x, 0.01, na.rm=T)
  return(c(x.min,x.max))
}

#define a function to estimate mean evolution rate in grids
mean_subrate_grid <- function(spatial_join, grids_poly, Species=5, habitat="Terrestrial"){
  if(!habitat %in% c("Terrestrial", "Marine")){
    print("Warning: habitat must be Terrestrial or Marine!")
  }
  #Unweighted geometric mean substitution rate in grids
  geometric_mean_subrate <- spatial_join%>%
    group_by(GridID)%>%
    summarise(SR=length(Species), gm.dn=exp(mean(log(dN.Mid))), gm.ds=exp(mean(log(dS.Mid))))%>%
    filter(SR>=Species)%>%
    #filter(gm.dn>=trans.extr(gm.dn)[1], gm.dn<=trans.extr(gm.dn)[2], gm.ds >=trans.extr(gm.ds)[1], gm.ds <=trans.extr(gm.ds)[2])%>%
    mutate(gm.dn=trans.extr.val(gm.dn), gm.ds=trans.extr.val(gm.ds))%>%
    dplyr::select(GridID, SR, gm.ds, gm.dn)
  
  #Unweighted geometric mean substitution in grids
  arithmetic_mean_subrate <- spatial_join%>%
    group_by(GridID)%>%
    summarise(SR=length(Species), am.dn=mean(dN.Mid), am.ds=mean(dS.Mid))%>%
    filter(SR>=Species)%>%
    #filter(am.dn>=trans.extr(am.dn)[1], am.dn<=trans.extr(am.dn)[2], am.ds >=trans.extr(am.ds)[1], am.ds <=trans.extr(am.ds)[2])%>%
    mutate(am.dn=trans.extr.val(am.dn), am.ds=trans.extr.val(am.ds))%>%
    dplyr::select(GridID, am.ds, am.dn)
  
  #Middle value of substitution in grids
  mid_subrate <- spatial_join%>%
    group_by(GridID)%>%
    summarise(SR=length(Species), mid.dn=median(dN.Mid), mid.ds=median(dS.Mid))%>%
    filter(SR>=Species)%>%
    #filter(mid.dn>=trans.extr(mid.dn)[1], mid.dn<=trans.extr(mid.dn)[2], mid.ds >=trans.extr(mid.ds)[1], mid.ds <=trans.extr(mid.ds)[2])%>%
    mutate(mid.dn=trans.extr.val(mid.dn), mid.ds=trans.extr.val(mid.ds))%>%
    dplyr::select(GridID, mid.ds, mid.dn)
  
  if(habitat=="Terrestrial"){
    mean_subrate <- grids_poly%>%
      dplyr::select(GridID, Lat, Lon, Terrest)%>%
      left_join(geometric_mean_subrate, by="GridID")%>%
      left_join(arithmetic_mean_subrate, by="GridID")%>%
      left_join(mid_subrate, by="GridID")%>%
      filter(Terrest==1)%>%
      dplyr::select(-Terrest)
  }
  if(habitat=="Marine"){
    mean_subrate <- grids_poly%>%
      dplyr::select(GridID, Lat, Lon, Marine)%>%
      left_join(geometric_mean_subrate, by="GridID")%>%
      left_join(arithmetic_mean_subrate, by="GridID")%>%
      left_join(mid_subrate, by="GridID")%>%
      filter(Marine==1)%>%
      dplyr::select(-Marine)
  }
  
  return(mean_subrate)
}

#define a function to estimate mean evolution rate in ecoregions
mean_subrate_ecoregion <- function(spatial_join, Species=5, habitat="Terrestrial"){
  
  if(!habitat %in% c("Terrestrial", "Marine")){
    print("Warning: habitat must be Terrestrial or Marine!")
  }
  
  geometric_mean_subrate <-spatial_join%>%
    group_by(ECO_CODE)%>%
    summarise(SR=length(Species), gm.dn=exp(mean(log(dN.Mid))), gm.ds=exp(mean(log(dS.Mid))))%>%
    filter(SR>=Species)%>%
    mutate(gm.dn=trans.extr.val(gm.dn), gm.ds=trans.extr.val(gm.ds))%>%
    dplyr::select(ECO_CODE, SR, gm.dn, gm.ds)
  
  arithmetic_mean_subrate <- spatial_join%>%
    group_by(ECO_CODE)%>%
    summarise(SR=length(Species), am.dn=mean(dN.Mid), am.ds=mean(dS.Mid))%>%
    filter(SR>=Species)%>%
    mutate(am.dn=trans.extr.val(am.dn), am.ds=trans.extr.val(am.ds))%>%
    dplyr::select(ECO_CODE, am.dn, am.ds)
  
  mid_subrate <- spatial_join%>%
    group_by(ECO_CODE)%>%
    summarise(SR=length(Species), mid.dn=median(dN.Mid), mid.ds=median(dS.Mid))%>%
    filter(SR>=Species)%>%
    mutate(mid.dn=trans.extr.val(mid.dn), mid.ds=trans.extr.val(mid.ds))%>%
    dplyr::select(ECO_CODE, mid.ds, mid.dn)
  
  if(habitat=="Marine"){
    mean_subrate <- marine.ecoregions%>%
      as.data.frame()%>%
      dplyr::select(ECO_CODE, Lat, Lon)%>%
      left_join(geometric_mean_subrate, by="ECO_CODE")%>%
      left_join(arithmetic_mean_subrate, by="ECO_CODE")%>%
      left_join(mid_subrate, by="ECO_CODE")
  }
  if(habitat=="Terrestrial"){
    mean_subrate <- terrestrial.ecoregions%>%
      as.data.frame()%>%
      dplyr::select(ECO_CODE, Lat, Lon)%>%
      left_join(geometric_mean_subrate, by="ECO_CODE")%>%
      left_join(arithmetic_mean_subrate, by="ECO_CODE")%>%
      left_join(mid_subrate, by="ECO_CODE")
  }
  return(mean_subrate)
}




#a function to plotting average substitution rates in geographic grids
plot_subrate_grids <- function(mean.rate.grids, rate.type, geo.type, colramp, asp=1.2){
  if(!"Lon" %in% names(mean.rate.grids)& "Lat" %in% names(mean.rate.grids)){
    print("Warning: Lon and Lat need to be included in mean.rate.grids data frame!")
  }
  if(!rate.type %in% c("dN", "dS")){
    print("Warning: rate.type must be one of dN and dS!")
  }
  if(!rate.type %in%names(mean.rate.grids)){
    print("Warning: 'rate.type' needs to be included in mean.rate.grids data frame!")
  }
  
  #data frame to raster
  ras <- rasterFromXYZ(mean.rate.grids[,c("Lon", "Lat", rate.type)], res=1, crs=CRS('+proj=longlat +datum=WGS84'))
  #ras <- disaggregate(ras, 5)
  #transform coordinates to Mollweide
  ras_moll <- projectRaster(ras, crs = CRS('+proj=moll'))
  
  #colors for mapping
  #colramp <- c("#4C7AB5","#FBFC9F","#D93529")#bule to red
  #colramp <- c("#ffffff","#dcdce3","#a3a1c7","#69539b","#3d1476")#bule
  #colramp <- c("#ffffff", "#f4cdb9","#f47155","#d22428", "#660b11")#red
  
  colours <- colorRampPalette(colramp)
  
  if(rate.type=="dS"){
    x <- mean.rate.grids%>%filter(!is.na(dS))%>%pull(dS)
  }
  if(rate.type=="dN"){
    x <- mean.rate.grids%>%filter(!is.na(dN))%>%pull(dN)
  }
  
  cuts <- BAMMtools::getJenksBreaks(mean.rate.grids%>%filter(!is.na(rate.type))%>%pull(rate.type), 21, subset = NULL)
  #cuts <- quantile(x, seq(0, 1, 0.05))
  names(cuts) <- NULL
  
  #plotting raster
  if(geo.type=="Terrestrial"){
    image(terrestrial.raster, col='#e1e1e1', axes=FALSE, xlab='', ylab='', asp=asp)#for NULL data
    image(ras_moll, breaks=cuts, col=colours(20), add=TRUE)
    image(marine.raster, col='#e1e1e1', add=TRUE)
  }
  
  if(geo.type=="Global"){
    image(marine.raster, col='#e1e1e1', axes=FALSE, xlab='', ylab='', asp=asp)#for NULL data
    image(terrestrial.raster, col='#e1e1e1', add=TRUE)#for NULL data
    image(ras_moll, breaks=cuts, col=colours(20), add=TRUE)
  }
  
  
  plot(country.poly, col=NA, lwd=0.3, add=TRUE)
  plot(outline.poly, col=NA, lwd=0.5, add=TRUE)
  
  #added legend
  epm::addLegend(ras_moll, ramp=colramp, ncolors = 10, 
                 location=c(-6000000,6000000,-9900000,-9400000), 
                 direction="horizontal", nTicks=0,side = 1,
                 labelDist=-0.3, cex.axis=0.6)
}



