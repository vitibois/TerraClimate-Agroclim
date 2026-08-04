library(terra)
library(ggplot2)
library(data.table)

setwd("E:/Benjamin/_Recherches_Dijon/_AAA_TerraClimate/2023_WWV_ANALYSIS/")

  # Read data sets
dobs <- cbind(period="2001-2020",as.data.frame(fread("AgroIndices_obs_VGDB_Pts_v1.2.3.csv")))
dfut <- cbind(period="2042",as.data.frame(fread("AgroIndices_plus2C_VGDB_Pts_v1.2.3.csv")))
  # Merge obs and future data sets
SameNames <- names(dobs)[names(dobs) %in% names(dfut)]
d <- rbind(dobs[,SameNames],dfut[,SameNames])

length(unique(dobs$CN_REG))

# Change regions names
  # NAPA 
names(d)
grep("Napa",unique(d$CN_REG),value=T)
  # The régional AVA Nappa Valley is identified as "AVA Napa Valley"
d$CN_REG <- sub("US_Napa Valley", "US_AVA Napa Valley", d$CN_REG)
  # Other "finer" AVA within NAPA are renamed as "Napa Valley"
Napa <- c("US_Calistoga",
          "US_St Helena",
          "US_Rutherford",
          "US_Oakville",
          "US_Yountville",
          "US_Oak Knool District")
d$CN_REG[d$CN_REG %in% Napa] <- "US_Napa Valley"
  # Marlborough
d$CN_REG[d$CN_REG == "NZ_Wairau Valley"] <- "NZ_Marlborough"

  # Selected regions
myRegions <- read.csv("ReferencesRegions_v02_Mendoza.csv")
myRegions <- myRegions[myRegions$Include==T,]
mymatch <- match(d$CN_REG,myRegions$CN_REG)
d$Type <- myRegions$Type[mymatch]
d$Style <- myRegions$Style[mymatch]
d$Type[is.na(d$Type)]="Undertermined"
d$Style[is.na(d$Style)]="Undertermined"
summary(factor(d$Type)) ; summary(factor(d$Style))
  #Transformation en facteur
myRegions$CN <- unlist(lapply(strsplit(myRegions$CN_REG, split="_"), function(x) x[1]))
myRegions$REG <- unlist(lapply(strsplit(myRegions$CN_REG, split="_"), function(x) x[2]))
myRegions$CN <- factor(myRegions$CN, levels=unique(myRegions$CN))
myRegions$REG <- factor(myRegions$REG, levels=unique(myRegions$REG))

  # Identify Regions and Country
d$CN <- substr(d$CN_REG,1,2)
d$REG <- substr(d$CN_REG,4,nchar(d$CN_REG))
head(d)

  # We remove koppen A (Tropical) and T (Polar) climates
d <- d[d$KOPPEN_Avg %in% 4:27,]

  # Thresholds selected
th <- read.csv("Criterions.csv")
th <- th[!is.na(th$Rule),]
th$var <- paste0(th$Index,"_",th$Statistic.over.time)
myrules <- paste(paste(th$Index, th$Statistic.over.time, sep="_"),th$Sign,th$Rule)


# Classes en pluie correspondantes aux classes DI
di.classes <- c(-1000, -100, 50, 150, 200)
rr49.classes <- c(0,200,350,475)
di.lims <- c(-100,50,150)
di.lab.y <- c(-150,-25,100, 200)
hi.lims <- c(1500,1800,2100,2400,3000,4000)
hi.lab.y <- c(1200,1650,1950, 2250, 2700, 3400)
# hi.names <- c("Very\ncool","Cool","Temp.","Temp.\nwarm", "Warm", "Hot")
# di.names <- c("Very dry","Mod. dry","Sub-humid","Humid")
# myxlab <- "Huglin index [°C.days]"
# myylab <- "Dryness index [mm]"

hi.names <- c("Très\nfrais","Frais","Temp.","Temp.\nchaud", "Chaud", "Très chaud")
di.names <- c("Très Sec","Sec","Sub-humide","Humide")
myxlab <- "Indice de Huglin [°C.jours]"
myylab <- "Indice de sècheresse [mm]"
myzlab <- "Indice de stress thermique [°C]"



################# PLOTS ########################
library(ggplot2)

  # A function to add the rectangle to map the criterion space retained
makerect <- function(Type="Still"){
  rxmax <- th[th$Type==Type & th$var == xvar & th$Sign == "<","Rule"]
  if(length(rxmax)==0) rxmax = xlims[2]
  rxmin <- th[th$Type==Type & th$var == xvar & th$Sign == ">","Rule"]
  if(length(rxmin)==0) rxmin = xlims[1]
  rymax <- th[th$Type==Type & th$var == yvar & th$Sign == "<","Rule"]
  if(length(rymax)==0)  rymax = ylims[2]
  rymin <- th[th$Type==Type & th$var == yvar & th$Sign == ">","Rule"]
  if(length(rymin)==0)  rymin  = ylims[1]
  rect <- c(rxmin,rxmax,rymin,rymax)
  return(rect)
}

# 1 - Selected regions - Obs and plus 2 facets + global criterions area ----
  # A subset of selected regions
ds <- d[d$CN_REG %in% myRegions$CN_REG,]
ds$CN <- factor(as.character(ds$CN), levels=myRegions$CN)
  # A 10% sample of the all wine producing regions
ds$REG <- factor()'
set.seed(55) ; ech <- sample(1:nrow(d),round(nrow(d)*0.1)) ; dsamp <- d[ech,]
dsamp <- dsamp[,-grep("Type",names(dsamp))] # We don't want to discriminate the type of wine for the whole data set

#mywidth <- 10
mywidth <- 6
pdf("./out_plots/Mendoza_Climagram_GeneralClimateFeatures_Obs_and_Plus2C_RegionsSubSet.pdf", width = mywidth, height = 6)
  mytitle <- "Diseases risks and Fruit developement conditions"
  xvar = "DSP_Avg" ; xlab = "Disease period precipitation [mm]" ; 
  yvar = "FSP_Avg" ; ylab= "Fruit development period precipitation [mm]"
  xlims = range(dsamp[,xvar], na.rm=T) ; ylims = range(dsamp[,yvar], na.rm=T) 
  rect <- makerect("Still")
  myrectdf <- cbind(rbind(makerect("Still"),makerect("Cognac"),makerect("Sparkling")),
                      data.frame(Type=c("Still","Cognac","Sparkling")))
  names(myrectdf) <- c("X1","X2","X3","X4","Type")
  
  makeplot <- function()
  {
    gp <- ggplot(ds, aes_string(x=xvar, y=yvar, color="REG"))+
    annotate("rect",xmin = rect[1], xmax = rect[2], ymin = rect[3], ymax = rect[4], fill="lightgreen", alpha=0.3) +
    #geom_rect(data = myrectdf, mapping = aes_string(xmin = "X1", xmax = "X2", ymin = "X3", ymax = "X4"), color="black")+
    geom_point(data=dsamp, aes_string(x=xvar, y=yvar), color="grey", alpha=0.2)+
    geom_point(aes(shape=CN))+stat_ellipse()+facet_grid(scenario~Type)+
    xlab(xlab)+ylab(ylab)+xlim(xlims)+ylim(ylims)+ #scale_colour_brewer(palette="Paired")+
        ggtitle(mytitle)+theme(legend.position="bottom")+theme_bw()+
      scale_color_discrete(name="Region")+scale_shape_discrete(name="Country")
    return(gp)
  }
  
  gp <- makeplot()
  print(gp)
  
  mytitle <- "Cold damage related risks"
  xvar = "WCR_Avg" ; xlab = "Winter cold risk index [°C]"; 
  yvar = "SFR_Avg" ; ylab= "Spring frost risk index [°C]"
  xlims = range(dsamp[,xvar], na.rm=T) ; ylims = range(dsamp[,yvar], na.rm=T) 
  rect <- makerect()
  gp <- makeplot()
  print(gp)
  
  mytitle <- "Earliness and Water availibility"
  xvar = "HI_Avg" ; xlab = "Huglin Index [°C.days]"
  yvar = "DI_Avg" ; ylab= "Dryness Index [mm]"
  xlims = range(dsamp[,xvar], na.rm=T) ; ylims = range(dsamp[,yvar], na.rm=T) 
  rect <- makerect()
  gp <- makeplot()
  print(gp)
  
  mytitle <- "Earliness and Water availibility"
  xvar = "GST49_Avg" ; xlab = "Growing season temperature [°C]"
  yvar = "GSP49_Avg" ; ylab= "Growing season precipitation [mm]"
  xlims = range(dsamp[,xvar], na.rm=T) ; ylims = range(dsamp[,yvar], na.rm=T) 
  rect <- makerect()
  gp <- makeplot()
  print(gp)
  
  mytitle <-"Ripening period thermal conditions"
  xvar = "HSI_Avg" ; xlab =  "Heat Stress Index [°C]"
  yvar = "CI_Avg" ; ylab= "Cool night index [°C] "
  xlims = range(dsamp[,xvar], na.rm=T) ; ylims = range(dsamp[,yvar], na.rm=T) 
  rect <- makerect()
  gp <- makeplot()
  print(gp)
dev.off()


# 2 - Selected regions - Obs only - global criterions area ----
  # A subset of selected regions
ds <- d[d$CN_REG %in% myRegions$CN_REG,]
unique(ds$scenario)
ds <- ds[ds$scenario == "OBS",] # Only obs data
  # A 10% sample of the all wine producing regions
set.seed(55) ; ech <- sample(1:nrow(d),round(nrow(d)*0.1)) ; dsamp <- d[ech,]
dsamp <- dsamp[dsamp$scenario == "OBS",-grep("Type",names(dsamp))] # We don't want to discriminate the type of wine for the whole data set

pdf("Climagram_GeneralClimateFeatures_ObsOnly_RegionsSubSet.pdf", width = 8, height = 5)
  mytitle <- "Diseases risks and Fruit developement conditions"
  xvar = "DSP_Avg" ; xlab = "Disease period precipitation [mm]" ; 
  yvar = "FSP_Avg" ; ylab= "Fruit development period precipitation [mm]"
  xlims = range(dsamp[,xvar], na.rm=T) ; ylims = range(dsamp[,yvar], na.rm=T) 
  rect <- makerect()
  makeplot <- function()
  {
    gp <- ggplot(ds, aes_string(x=xvar, y=yvar, color="CN_REG"))+
      annotate("rect",xmin = rect[1], xmax = rect[2], ymin = rect[3], ymax = rect[4], fill="lightgreen", alpha=0.3) +
      geom_point(data=dsamp, aes_string(x=xvar, y=yvar), color="grey", alpha=0.2)+
      geom_point(aes(shape=Type),size=3, alpha=0.5)+stat_ellipse()+
      xlab(xlab)+ylab(ylab)+xlim(xlims)+ylim(ylims)+ 
      ggtitle(mytitle)+theme()+theme_bw()+
      #scale_colour_brewer(palette="Set1", name="Region", type = "colour")+
      #scale_colour_manual( name="Region", values=RColorBrewer::brewer.pal(n = 12, name = "Paired", alpha=0.5) )+
      scale_colour_discrete(name="Region" )+
      scale_shape_manual(name="Wine type", values=c(4,1,18))
    
    return(gp)
  }
  plot_allwinereg <- function()
  {
    gp <- ggplot(dsamp, aes_string(x=xvar, y=yvar))+
      geom_point(data=dsamp, aes_string(x=xvar, y=yvar), color="forestgreen", alpha=0.2)+xlab(xlab)+ylab(ylab)+xlim(xlims)+ylim(ylims)+ 
      ggtitle(mytitle)+theme(legend.position="none")+theme_bw()
    return(gp)
  } 
  gpall <- plot_allwinereg()
  print(gpall)

  gp <- makeplot()
  print(gp)
  
  mytitle <- "Cold damage related risks"
  xvar = "WCR_Avg" ; xlab = "Winter cold risk index [°C]"; 
  yvar = "SFR_Avg" ; ylab= "Spring frost risk index [°C]"
  xlims = range(dsamp[,xvar], na.rm=T) ; ylims = range(dsamp[,yvar], na.rm=T) 
  rect <- makerect()
  gpall <- plot_allwinereg()
  print(gpall)
  gp <- makeplot()
  print(gp)
  
  mytitle <- "Earliness and Water availibility"
  xvar = "HI_Avg" ; xlab = "Huglin Index [°C.days]"
  yvar = "DI_Avg" ; ylab= "Dryness Index [mm]"
  xlims = range(dsamp[,xvar], na.rm=T) ; ylims = range(dsamp[,yvar], na.rm=T) 
  rect <- makerect()
  gpall <- plot_allwinereg()
  print(gpall)
  gp <- makeplot()
  print(gp)
  
  mytitle <- "Earliness and Water availibility"
  xvar = "GST49_Avg" ; xlab = "Growing season temperature [°C]"
  yvar = "GSP49_Avg" ; ylab= "Growing season precipitation [mm]"
  xlims = range(dsamp[,xvar], na.rm=T) ; ylims = range(dsamp[,yvar], na.rm=T) 
  rect <- makerect()
  gpall <- plot_allwinereg()
  print(gpall)
  gp <- makeplot()
  print(gp)
  
  mytitle <-"Ripening period thermal conditions"
  xvar = "HSI_Avg" ; xlab =  "Heat Stress Index [°C]"
  yvar = "CI_Avg" ; ylab= "Cool night index [°C] "
  xlims = range(dsamp[,xvar], na.rm=T) ; ylims = range(dsamp[,yvar], na.rm=T) 
  rect <- makerect()
  gpall <- plot_allwinereg()
  print(gpall)
  gp <- makeplot()
  print(gp)
dev.off()



# 3 - Selected regions - Obs only - global criterions area OBS AND PLUS 2C ----
# A subset of selected regions
ds <- d[d$CN_REG %in% myRegions$CN_REG,]
unique(ds$scenario)

ds$fac_scenario <- factor(ds$scenario, labels=c("2001-2020", "+2°C (2042)"))
# A 10% sample of the all wine producing regions
set.seed(55) ; ech <- sample(1:nrow(d),round(nrow(d)*0.1)) ; dsamp <- d[ech,]
dsamp <- dsamp[dsamp$scenario == "OBS",-grep("Type",names(dsamp))] # We don't want to discriminate the type of wine for the whole data set

pdf("Climagram_GeneralClimateFeatures_ObsPlus2C_RegionsSubSet.pdf", width = 10, height = 5)
  mytitle <- "Diseases risks and Fruit developement conditions"
  xvar = "DSP_Avg" ; xlab = "Disease period precipitation [mm]" ; 
  yvar = "FSP_Avg" ; ylab= "Fruit development period precipitation [mm]"
  xlims = range(dsamp[,xvar], na.rm=T) ; ylims = range(dsamp[,yvar], na.rm=T) 
  rect <- makerect()
  makeplot <- function()
  {
    gp <- ggplot(ds, aes_string(x=xvar, y=yvar, color="CN_REG"))+
      annotate("rect",xmin = rect[1], xmax = rect[2], ymin = rect[3], ymax = rect[4], fill="lightgreen", alpha=0.3) +
      geom_point(data=dsamp, aes_string(x=xvar, y=yvar), color="grey", alpha=0.2)+
      geom_point(aes(shape=Type),size=3, alpha=0.5)+stat_ellipse()+
      xlab(xlab)+ylab(ylab)+xlim(xlims)+ylim(ylims)+ 
      ggtitle(mytitle)+theme()+theme_bw()+
      #scale_colour_brewer(palette="Set1", name="Region", type = "colour")+
      #scale_colour_manual( name="Region", values=RColorBrewer::brewer.pal(n = 12, name = "Paired", alpha=0.5) )+
      scale_colour_discrete(name="Region" )+
      scale_shape_manual(name="Wine type", values=c(4,1,18))+facet_wrap(fac_scenario~.,ncol=2)
    
    return(gp)
  }
  gp <- makeplot()
  print(gp)
  
  mytitle <- "Cold damage related risks"
  xvar = "WCR_Avg" ; xlab = "Winter cold risk index [°C]"; 
  yvar = "SFR_Avg" ; ylab= "Spring frost risk index [°C]"
  xlims = range(dsamp[,xvar], na.rm=T) ; ylims = range(dsamp[,yvar], na.rm=T) 
  rect <- makerect()
  gp <- makeplot()
  print(gp)
  
  mytitle <- "Earliness and Water availibility"
  xvar = "HI_Avg" ; xlab = "Huglin Index [°C.days]"
  yvar = "DI_Avg" ; ylab= "Dryness Index [mm]"
  xlims = range(dsamp[,xvar], na.rm=T) ; ylims = range(dsamp[,yvar], na.rm=T) 
  rect <- makerect()
  gp <- makeplot()
  print(gp)
  

  mytitle <- "Earliness and Water availibility"
  xvar = "GST49_Avg" ; xlab = "Growing season temperature [°C]"
  yvar = "GSP49_Avg" ; ylab= "Growing season precipitation [mm]"
  xlims = range(dsamp[,xvar], na.rm=T) ; ylims = range(dsamp[,yvar], na.rm=T) 
  rect <- makerect()
  gp <- makeplot()
  print(gp)
  
  mytitle <-"Ripening period thermal conditions"
  xvar = "HSI_Avg" ; xlab =  "Heat Stress Index [°C]"
  yvar = "CI_Avg" ; ylab= "Cool night index [°C] "
  xlims = range(dsamp[,xvar], na.rm=T) ; ylims = range(dsamp[,yvar], na.rm=T) 
  rect <- makerect()
  gp <- makeplot()
  print(gp)
  
  
dev.off()


# A subset of selected regions
ds <- d[d$CN_REG %in% myRegions$CN_REG,]
ds$fac_scenario <- factor(ds$scenario, labels=c("2001-2020", "+2°C (2042)"))
ds$Type <- factor(ds$Type, labels = c("Cognac","Champagne","Vins tranquilles"))
# A 100% sample of the all wine producing regions
set.seed(55) ; ech <- sample(1:nrow(d),round(nrow(d)*0.3)) ;
dsamp <- d[ech,]
dsamp <- d
dsamp <- dsamp[,-grep("Type",names(dsamp))] # We don't want to discriminate the type of wine for the whole data set

pdf("Climagram_IndHuglin-Secheresse.pdf", width = 8.5, height = 7)
  mytitle <- "2001-2020 : Précocité et disponibilité en eau "
  xvar = "HI_Avg" ; xlab = "Indice de Huglin [°C.jour]"
  yvar = "DI_Avg" ; ylab= "Indice de Sècheresse [mm]"
  xlims = range(dsamp[,xvar], na.rm=T) ; ylims = range(dsamp[,yvar], na.rm=T) 

  # Plot vierge
  ggplot(dsamp[dsamp$scenario=="OBS",], aes_string(x=xvar, y=yvar))+
    geom_hline(yintercept = di.lims, )+geom_vline(xintercept = hi.lims)+
    geom_point(color="forestgreen", alpha=0.2)+
    xlab(xlab)+ylab(ylab)+xlim(xlims)+ylim(ylims)+ 
    ggtitle(mytitle)+theme(legend.position = "none")+theme_bw()+
    annotate("text",label=di.names, x=rep(xlims[2]-100,length(di.names)), y=di.lab.y, angle=90, hjust=.5, fontface=2, size=4, vjust=0) +    
    annotate("text",label=hi.names, x=hi.lab.y, y=rep(ylims[1],length(hi.names)), angle=0, fontface=2, size=4, vjust=0)  
    
  ggplot(ds[ds$scenario=="OBS",], aes_string(x=xvar, y=yvar, color="CN_REG"))+
  geom_hline(yintercept = di.lims, )+geom_vline(xintercept = hi.lims)+
  annotate("rect",xmin = rect[1], xmax = rect[2], ymin = rect[3], ymax = rect[4], fill="lightgreen", alpha=0.3) +
  geom_point(data=dsamp[dsamp$scenario=="OBS",], aes_string(x=xvar, y=yvar), color="grey", alpha=0.2)+
  geom_point(aes(shape=Type),size=3, alpha=0.5)+stat_ellipse()+
  xlab(xlab)+ylab(ylab)+xlim(xlims)+ylim(ylims)+ 
  ggtitle(mytitle)+theme()+theme_bw()+
  #scale_colour_brewer(palette="Set1", name="Region", type = "colour")+
  #scale_colour_manual( name="Region", values=RColorBrewer::brewer.pal(n = 12, name = "Paired", alpha=0.5) )+
    annotate("text",label=di.names, x=rep(xlims[2]-100,length(di.names)), y=di.lab.y, angle=90, hjust=.5, fontface=2, size=4, vjust=0) +    
    annotate("text",label=hi.names, x=hi.lab.y, y=rep(ylims[1],length(hi.names)), angle=0, fontface=2, size=4, vjust=0)  +
    scale_colour_discrete(name="Region" )+
    scale_shape_manual(name="Type de vin", values=c(4,1,18))

  mytitle <- "+2°C : Précocité et disponibilité en eau "
  xlims = range(dsamp[,xvar], na.rm=T) ; ylims = range(dsamp[,yvar], na.rm=T) 
  ggplot(ds[ds$scenario=="plus2",], aes_string(x=xvar, y=yvar, color="CN_REG"))+
    geom_hline(yintercept = di.lims, )+geom_vline(xintercept = hi.lims)+
    annotate("rect",xmin = rect[1], xmax = rect[2], ymin = rect[3], ymax = rect[4], fill="lightgreen", alpha=0.3) +
    geom_point(data=dsamp[dsamp$scenario=="plus2",], aes_string(x=xvar, y=yvar), color="grey", alpha=0.2)+
    geom_point(aes(shape=Type),size=3, alpha=0.5)+stat_ellipse()+
    xlab(xlab)+ylab(ylab)+xlim(xlims)+ylim(ylims)+ 
    ggtitle(mytitle)+theme()+theme_bw()+
    #scale_colour_brewer(palette="Set1", name="Region", type = "colour")+
    #scale_colour_manual( name="Region", values=RColorBrewer::brewer.pal(n = 12, name = "Paired", alpha=0.5) )+
    annotate("text",label=di.names, x=rep(xlims[2]-100,length(di.names)), y=di.lab.y, angle=90, hjust=.5, fontface=2, size=4, vjust=0) +    
    annotate("text",label=hi.names, x=hi.lab.y, y=rep(ylims[1],length(hi.names)), angle=0, fontface=2, size=4, vjust=0)  +
    scale_colour_discrete(name="Region" )+
    scale_shape_manual(name="Wine type", values=c(4,1,18))
dev.off()
