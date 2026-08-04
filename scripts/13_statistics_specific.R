# 2023-12-20 B. Bois
# Statistics of potential areas for viticulture, specific (Cognac) criteria
# (uses the maps produced by script 11_filter_vineyard_maps_specific.R,
# written under the mydir subdirectory below -- must match script 11's own
# `mydir`. Shared reference tables (ColTable_PotArea.csv, ZoomTable.csv) and
# basemaps are assumed to stay at the top level, common to all criteria sets)
library(terra)
library(ggplot2)

source("R/config.R")
setwd(file.path(data_root, "2023_WWV_ANALYSIS"))
mydir <- "AAA_CognacCriterions"


####### LOAD FILES ###########
  # Potential areas
obs <- rast(file.path(mydir, "rasters_zonespotentielles/AAA_ZonesPotentielles_OBS2001-2020.tif"))
p2c <- rast(file.path(mydir, "rasters_zonespotentielles/AAA_ZonesPotentielles_plus2C.tif"))
p4c <- rast(file.path(mydir, "rasters_zonespotentielles/AAA_ZonesPotentielles_plus4C.tif"))
pot <- rast(list(obs=obs,plus2C=p2c,plus4C=p4c))
time(pot) <- NULL

  # Surface area of the potential zones
cellarea <- cellSize(obs,mask=T,unit="km")
plot(cellarea)
global(cellarea,"sum",na.rm=T) /10^6
potbin <- pot
pot <- pot*cellarea
plot(potbin)
plot(pot)

  # HYBRID RASTER COMBINING ALL OPTIONS
# Code: 1 = favourable climate, 0 = unfavourable climate
# THREE DIGITS: the first is for the +4C climate, the second for the +2C
# climate, and the third for the current (2001-2020) climate
rpot <- potbin$obs
rpot <- rpot+potbin$plus2C*10
rpot <- rpot+potbin$plus4C*100
coltable_pot <- read.csv("ColTable_PotArea.csv")
levels(rpot) <- coltable_pot[,1:2]
add.alpha <- function(col, alpha=1){
  apply(sapply(col, col2rgb)/255, 2, function(x)
    rgb(x[1], x[2], x[3], alpha=alpha)) }
coltable_pot$colhex <- apply(coltable_pot[,c("color","alpha")],1,function(x) add.alpha(x[1],x[2]))

# Countries of the world
wld <- vect(file.path(sig_root, "world/World_AdmBoundaries_Countries/world-administrative-boundaries.shp"))
wld$name <- sub("The former Yugoslav Republic of ", "", wld$name)
wld$name <- sub("U.K. of Great Britain and Northern Ireland","United Kingdom", wld$name)

  # Google Earth hybrid basemap
ge <- rast("GE_hybrid_0.1.tif")
crs(ge) <- crs(obs)
ge <- crop(ge, obs)
gena <- resample(ge, obs)
gena <- gena*(obs*0+1)

hrge <- rast("GE_hybrid_0.04.tif")
crs(hrge) <- crs(obs)
hrge <- crop(hrge, obs)
hrge <- resample(hrge, obs)
#hrge <- mask(hrge, obs)
plotRGB(hrge)

  # Lighten the raster (more brightness)
palege <- hrge+50
palege[palege > 255] <- 255
plotRGB(palege, ext=ext(c(-10,20,35,60)))
plotRGB(palege)

  # World DEM
dem <- rast(file.path(sig_root, "_MNT/SRTM_10minarc_WorldClim/alt.bil"))
dem <- resample(dem,obs, method="cubic")
plot(dem)

  # Areas currently planted with vineyards
vpts <-  vect(file.path(data_root, "VGDB_Points_LouisDelelee/pts_global_v4.0.shp"))
rvpts <- rasterize(x = vpts, y = obs)
rvpts <- cellSize(rvpts,mask=T,unit="km")
plot(rvpts)


####### COMPUTE AND PLOT POTENTIAL AREAS ###########
  # Total area of pixels where vineyards are identified -----
global(rvpts,  "sum", na.rm=T)*100/10^6 # 23 Mha... much larger than the actual vineyard area of 7.3 Mha (OIV 2022)
areavgbd <- zonal(rvpts,wld,fun="sum", na.rm=T)
areavgbd$CN <- wld$iso_3166_1_
areavgbd$Country <- wld$name
areavgbd$Continent <- wld$continent
areavgbd <- areavgbd[order(areavgbd$area, decreasing=T),]

  # VGDB (current vineyard) area under favourable climate ----
  # potential areas
vgdbpot <- rvpts*pot
vgdbarea <- global(rvpts,  "sum", na.rm=T)
vgdbarea  <- rbind(data.frame(area=vgdbarea), global(vgdbpot,  "sum", na.rm=T))
vgdbarea$prop <- vgdbarea$sum/vgdbarea$sum[1]
vgdbarea

  # Total by country and by type of potential -----

r2040 <- rpot == 110 | rpot == 111
r2000 <- rpot == 111
r20002040  <- rpot == 11 | rpot == 111
pot$all21st <- cellarea*r2000
pot$plus2plus4 <- cellarea*r2040
pot$obsplus2 <- cellarea*r20002040
time(pot) <- NULL
potbin <- pot > 0
plot(potbin)

  # Compute the area by country and by period
names(pot)
area <- zonal(pot,wld,fun="sum", na.rm=T)
names(area) <- paste0("areapot_",names(area))
area$CN <- wld$iso_3166_1_
area$Country <- wld$name
area$World_region <- wld$region
area$Continent <- wld$continent
wld$CountryArea_km2  <- expanse(wld, unit="km")
area$CountryArea_km2 <- wld$CountryArea_km2

head(area)
summary(area)
area <- area[!is.na(area$areapot_all21st),]

  # World total
globarea <- unlist(global(cellarea,  "sum", na.rm=T))
totarea <- apply(area[,grep("^area",names(area))],2,sum, na.rm=T)
totarea/10^6
globarea/10^6
totarea/globarea

  # Area sub-tables
wldreg_area <- aggregate(area[,grep("^area",names(area))],by=list(World_region=area$World_region),sum, na.rm=T)
wldreg_area <- wldreg_area[order(wldreg_area$areapot_all21st,decreasing=T),]
continent_area <- aggregate(area[,grep("^area",names(area))],by=list(Continent=area$Continent),sum, na.rm=T)
continent_area <- continent_area[order(continent_area$areapot_all21st,decreasing=T),]

  # Save tables
area <- area[order(area$areapot_all21st, decreasing=T),]
head(area)
write.csv(area,file.path(mydir, "out_tables/Superficies_ZonesPotentielles_km2_VignoblesActuelsInclus.csv"), row.names=T)

  

#### Graph of available areas by period
  clong <- tidyr::pivot_longer(continent_area, cols=grep("^area",names(continent_area)), names_to = "Period", values_to = "area_km2")
  clong$x1 <- rep(as.Date(c("2001-01-01","2033-01-01","2075-01-01", "2001-01-01","2033-01-01","2001-01-01")),length(unique(clong$Period)))
  clong$x2 <- rep(as.Date(c("2020-12-31","2052-12-31","2094-12-31", "2094-12-31","2094-12-31","2052-12-31")),length(unique(clong$Period)))
  clong <- clong[clong$area_km2>0,]
  head(clong)
  clong$PeriodFac <- factor(clong$Period, labels=c("All 21st Centuty", "Recent (2001-2020)", "Recent and +2°C",
                                             "+2°C (2040)", "+2°C and +4°C", "+4°C") )
pdf(file.path(mydir, "out_plots/Continent_areas_per_period.pdf"), height = 5, width=8)
  ggplot(clong)+geom_segment(aes(x=x1,xend=x2,y=PeriodFac, yend=PeriodFac,size=area_km2/1000,color=Continent))+facet_wrap(Continent~.)+scale_size_continuous(name="Area [1000 km²]",range=c(1,10), breaks = c(0,100,500,1000,2000,3000,3500))+xlab("Time")+ylab(NULL)
dev.off()

### Graph of areas by continent
names(area)
pdf(file.path(mydir, "out_plots/Bar_Area_Per_Continent.pdf"), width = 7, height = 5)
  area$myarea <- area$areapot_all21st
  ggplot(data=area[area$myarea > 0, ], aes(x=factor(Continent, levels=continent_area$Continent), y=myarea/1000,  fill=World_region))+geom_bar(stat="identity")+ylab("Area [1000 km²]")+xlab(NULL)+ggtitle("Area with climate features adapted for wine production\n throughout the 21st century")+ylim(c(0,3000))

  area$myarea <- area$areapot_plus2plus4
  ggplot(data=area[area$myarea > 0, ], aes(x=factor(Continent, levels=continent_area$Continent), y=myarea/1000,  fill=World_region))+geom_bar(stat="identity")+ylab("Area [1000 km²]")+xlab(NULL)+ggtitle("Area with climate features adapted for wine production\n from 2042 (+2C°)")+ylim(c(0,3000))
  
  area$myarea <- area$areapot_obsplus2
  ggplot(data=area[area$myarea > 0, ], aes(x=factor(Continent, levels=continent_area$Continent), y=myarea/1000,  fill=World_region))+geom_bar(stat="identity")+ylab("Area [1000 km²]")+xlab(NULL)+ggtitle("Area with climate features adapted for wine production\n until 2070's (not suitable for +4°C)")+ylim(c(0,3000))
  
  area$myarea <- area$areapot_plus4C
  ggplot(data=area[area$myarea > 0, ], aes(x=factor(Continent, levels=continent_area$Continent), y=myarea/1000,  fill=World_region))+geom_bar(stat="identity")+ylab("Area [1000 km²]")+xlab(NULL)+ggtitle("Area with climate features adapted for wine production\n from 2084 (+4°C)")+ylim(c(0,3000))
dev.off()


pdf(file.path(mydir, "out_plots/Bar_Area_Per_Country.pdf"), width = 10, height = 6)
  ggplot(data=area[area$areapot_all21st>0,], aes(x=reorder(Country,-areapot_all21st), y=areapot_all21st/1000))+  geom_bar(stat="identity", aes(fill=Continent))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10) )+ylab("Area [1000 km²]")+xlab(NULL)+ggtitle("Area with climate features adapted for wine production\n throughout the 21st century")
  ggplot(data=area[area$areapot_plus2plus4>0,], aes(x=reorder(Country,-areapot_plus2plus4), y=areapot_plus2plus4/1000))+  geom_bar(stat="identity", aes(fill=Continent))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10) )+ylab("Area [1000 km²]")+xlab(NULL)+ggtitle("Area with climate features adapted for wine production\n from 2042 (+2C°)")
  ggplot(data=area[area$areapot_obsplus2>0,], aes(x=reorder(Country,-areapot_obsplus2), y=areapot_obsplus2/1000))+  geom_bar(stat="identity", aes(fill=Continent))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10) )+ylab("Area [1000 km²]")+xlab(NULL)+ggtitle("Area with climate features adapted for wine production\n until 2070's (not suitable for +4°C)")
  ggplot(data=area[area$areapot_plus4C>0,], aes(x=reorder(Country,-areapot_plus4C), y=areapot_plus4C/1000))+  geom_bar(stat="identity", aes(fill=Continent))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10) )+ylab("Area [1000 km²]")+xlab(NULL)+ggtitle("Area with climate features adapted for wine production\n from 2084 (+4°C)")
dev.off()  



  # Remove areas currently planted with vineyards
narvpts <- is.na(rvpts)
#narvpts[!narvpts] <- NA
plot(narvpts, colNA="black")
pot2bin <- potbin * narvpts
pot2 <- pot * narvpts
plot(potbin)
  # Compute the area by country
area2 <- zonal(pot2,wld,fun="sum", na.rm=T)
area2$CN <- wld$iso_3166_1_
area2$Country <- wld$name
area2$World_region <- wld$region
area2$Continent <- wld$continent
area2$CountryArea_km2 <- wld$CountryArea_km2

  # World total
tot2area <- apply(area2[,c("obs","plus2C","plus4C")],2,sum, na.rm=T)
totarea/10^6 ; tot2area/10^6
totarea/10^6 - tot2area/10^6

totarea/globarea # Fraction of grapevine cultivable area other the world
tot2area/globarea # Fraction of grapevine cultivable area other the world excluding areas already cultivated with grapevines
(totarea - tot2area)/globarea




  # World total
globarea <- unlist(global(cellarea,"sum",na.rm=T))
totarea/10^6
globarea/10^6
sum(wld$CountryArea_km2)/10^6

totarea/globarea # Fraction of grapevine cultivable area other the world



##### "CURRENT VINEYARDS" ANALYSIS (VGDB POINTS L. DELELEE) ####
  # Maps of current vineyards meeting or not the selected climatic criteria
pdf(file.path(mydir, "outmaps/VignoblesActuels_HorsCriteres.pdf"), width = 7, height = 4)
  myfun <- function(){
    plot(vpts, cex=0.2,add=T, col="orangered")
    plot(vpotobs, cex=0.2,add=T)
  }
  vpotobs <- vgdbpot$lyr1 > 0
  vpotobs[!vpotobs] <- NA
  vpotobs <- vect(crds(vpotobs))
  plot(dem, main="2001-2020", col=terrain.colors(256, alpha=0.8), colNA="lightcyan2", legend=F, ylim=c(-60,85), fun=myfun)
  vpotobs <- vgdbpot$lyr2 > 0
  vpotobs[!vpotobs] <- NA
  vpotobs <- vect(crds(vpotobs))
  plot(dem, main="2042 (+2°C)", col=terrain.colors(256, alpha=0.8), colNA="lightcyan2", legend=F, ylim=c(-60,85), fun=myfun)
  vpotobs <- vgdbpot$lyr3 > 0
  vpotobs[!vpotobs] <- NA
  vpotobs <- vect(crds(vpotobs))
  plot(dem, main="2084 (+4°C)", col=terrain.colors(256, alpha=0.8), colNA="lightcyan2", legend=F, ylim=c(-60,85), fun=myfun)
dev.off()


##### POTENTIAL AREAS - "CURRENT" VINEYARDS INCLUDED ######
# Compute the areas
novpts <- is.na(rvpts)

  # World map of the potentially suitable wine-growing areas
pdf(file.path(mydir, "outmaps/CartesMondiale_ZonesPotentielles_VignoblesActuelsInclus.pdf"), width = 10,height = 6)
  plot(potbin$obs, legend=F, ylim=c(-60,80), colNA="lightskyblue", main="2001-2020", maxcell=40*10^6,
    fun=function(){ 
      plot(vpts, add=T, cex=0.05, col=rgb(0,0,0,0.1))
      lines(wld, col=rgb(0,0,0,0.3), lwd=0.1)
      })
  plot(potbin$plus2C, legend=F, ylim=c(-60,80), colNA="lightskyblue", main="Plus 2°C", maxcell=40*10^6,
       fun=function(){ 
         plot(vpts, add=T, cex=0.05, col=rgb(0,0,0,0.1))
         lines(wld, col=rgb(0,0,0,0.3), lwd=0.1)
       })
  plot(potbin$plus4C, legend=F, ylim=c(-60,80), colNA="lightskyblue", main="Plus 4°C", maxcell=40*10^6,
       fun=function(){ 
         plot(vpts, add=T, cex=0.05, col=rgb(0,0,0,0.1))
         lines(wld, col=rgb(0,0,0,0.3), lwd=0.1)
       })
dev.off()

  # Regional maps of the potential areas
pdf(file.path(mydir, "outmaps/Regional_ZonesPotentielles_VignoblesActuelsInclus.pdf"), width=12, height = 6)
zoomtab <- read.csv("ZoomTable.csv")
for(i in 1:nrow(zoomtab))
{
  myext <- ext(unlist(zoomtab[i,-1]))
  plot(potbin, legend=F, colNA="lightskyblue", ext=myext, nc=3,fun=function() lines(wld, col=rgb(0,0,0,0.3)))
}
dev.off()


# Vector of the areas not viable in the long run (for hatching)
r2080 <- potbin$plus4C == 0
r2080[!r2080] <- NA
plotRGB(gena)
plot(r2080, add=T)
v2080 <- as.polygons(r2080, aggregate=T)

# 
coltable_pot <- read.csv("ColTable_PotArea.csv")
add.alpha <- function(col, alpha=1){
  apply(sapply(col, col2rgb)/255, 2, function(x) 
    rgb(x[1], x[2], x[3], alpha=alpha)) }
coltable_pot$colhex <- apply(coltable_pot[,c("color2","alpha2")],1,function(x) add.alpha(x[1],x[2]))
rpot
levels(rpot) <- coltable_pot[,c("id","Suitability2")]
zoomtab <- read.csv("ZoomTable.csv")
myext <- ext(unlist(zoomtab[3,-1]))
myext
pdf(file.path(mydir, "outmaps/AAA_Global_and_Regional_ZonesPotentielles2040AllWines_VignoblesActuelsInclus.pdf"), width=12, height = 6)
  plot(rpot, legend=T, ylim=c(-60,80), colNA="lightskyblue", maxcell=40*10^6,
     fun=function(){ 
       plot(r2080, col=rgb(0.3,0,0,0.3), legend=F,add=T, maxcell=40*10^6)
       plot(vpts, add=T, cex=0.05, col=rgb(0,0,0,0.1))
       lines(wld, col=rgb(0,0,0,0.3), lwd=0.1)
     }, 
     mar=c(3.1, 3.1, 2.1, 10.1), 
     plg=list(title="Suitability"))

  zoomtab <- read.csv("ZoomTable.csv")
  for(i in 1:nrow(zoomtab))
  {  
    myext <- ext(unlist(zoomtab[i,-1]))
    plot(rpot, ext=myext, maxcell=40*10^6, colNA="lightskyblue",fun=function(){ 
      plot(r2080, col=rgb(0.3,0,0,0.3), legend=F, ext=myext,add=T, maxcell=40*10^6)
      plot(vpts, add=T, cex=0.05, col=rgb(0,0,0,0.1))
      lines(wld, col=rgb(0,0,0,0.3), lwd=0.1)
    }, mar=c(3.1, 3.1, 2.1, 10.1), 
    plg=list(title="Suitability"))
  }
dev.off()

r21st <-  rpot != 111 
r21st[!r21st] <- NA
rblanc <- hrge$GE_hybrid_0.04_1*0+1

pdf(file.path(mydir, "outmaps/AAA_Global_and_Regional_ZonesPotentielles2040AllWines_VignoblesActuelsInclus.pdf"), width=12, height = 6)
  plot(rpot, legend=T, ylim=c(-60,80), colNA="lightskyblue", maxcell=40*10^6,
       fun=function(){ 
         plot(r2080, col=rgb(0.3,0,0,0.3), legend=F,add=T, maxcell=40*10^6)
         plot(vpts, add=T, cex=0.05, col=rgb(0,0,0,0.1))
         lines(wld, col=rgb(0,0,0,0.3), lwd=0.1)
       }, 
       mar=c(3.1, 3.1, 2.1, 10.1), 
       plg=list(title="Suitability"))
  zoomtab <- read.csv("ZoomTable.csv")
  for(i in 1:nrow(zoomtab))
  {  
    myext <- ext(unlist(zoomtab[i,-1]))
    plot(rpot, ext=myext, maxcell=40*10^6, colNA="lightskyblue",fun=function(){ 
      plot(r2080, col=rgb(0.3,0,0,0.3), legend=F, ext=myext,add=T, maxcell=40*10^6)
      plot(vpts, add=T, cex=0.05, col=rgb(0,0,0,0.1))
      lines(wld, col=rgb(0,0,0,0.3), lwd=0.1)
    }, mar=c(3.1, 3.1, 2.1, 10.1), 
    plg=list(title="Suitability"))
  }
dev.off()
zoomtab <- read.csv("ZoomTable.csv")


plot(r21st)
v21st <- as.polygons(r21st, aggregate=T)


pdf(file.path(mydir, "outmaps/AAA_Global_and_Regional_ZonesPot-2001-2100-AllWines_VignoblesActuelsInclus.pdf"), width=12, height = 6)
  mymaxcell <- 40*10^6
  #mymaxcell <- 10^6
  #myext <- ext(unlist(zoomtab[3,-1]))
  myext <- ext(c(-180,180,-75,85))
  plotRGB(palege, colNA="lightskyblue", maxcell=mymaxcell, ext=myext)
  plot(r21st, add=T, col=rgb(1,1,1,0.7), maxcell=mymaxcell, ext=myext)
  plot(vpts, add=T, cex=0.05, col=rgb(0,0.8,0,0.2))
  sbar(d=5000,labels = c(0,2500,"5000 km"), halo = T,lonlat = T ,  ticks = T)
  for(i in 1:nrow(zoomtab))
  {  
    myext <- ext(unlist(zoomtab[i,-1]))
    plotRGB(hrge, ext=myext, colNA="lightskyblue", maxcell=mymaxcell)
    plot(r21st, add=T, col=rgb(1,1,1,0.7), maxcell=mymaxcell, ext=myext)
    plot(vpts, add=T, cex=0.05, col=rgb(0,0.8,0,0.2))
    sbar(d=1000,labels = c(0,500,"1000 km"), halo = T,lonlat = T ,  ticks = T)
  }
dev.off()


writeRaster(potbin$all21st,file.path(mydir, "AAA_ZonesPotentielles_2001-2100.tif"), overwrite=T)

#### INCLUDING CURRENTLY CULTIVATED AREAS IN THE ANALYSIS IS NOT ACTUALLY A
### PROBLEM, BUT TO ALSO HAVE FIGURES/MAPS FOR THE POTENTIAL AREAS THAT ARE
#### NOT ALREADY PLANTED WITH VINEYARDS, WE ALSO COMPUTE THOSE BELOW
##### POTENTIAL AREAS NOT ALREADY PLANTED WITH VINEYARDS ######

# World map of the potentially suitable wine-growing areas, excluding current vineyards
pdf(file.path(mydir, "outmaps/CartesMondiale_ZonesPotentielles_SansVignoblesActuels.pdf"), width = 10,height = 6)
plot(pot2bin$obs, legend=F, ylim=c(-60,80), colNA="lightskyblue", main="2001-2020", maxcell=40*10^6,
     fun=function(){ 
       plot(vpts, add=T, cex=0.05, col=rgb(0,0,0,0.1))
       lines(wld, col=rgb(0,0,0,0.3), lwd=0.1)
     })
plot(pot2bin$plus2C, legend=F, ylim=c(-60,80), colNA="lightskyblue", main="Plus 2°C", maxcell=40*10^6,
     fun=function(){ 
       plot(vpts, add=T, cex=0.05, col=rgb(0,0,0,0.1))
       lines(wld, col=rgb(0,0,0,0.3), lwd=0.1)
     })
plot(pot2bin$plus4C, legend=F, ylim=c(-60,80), colNA="lightskyblue", main="Plus 4°C", maxcell=40*10^6,
     fun=function(){ 
       plot(vpts, add=T, cex=0.05, col=rgb(0,0,0,0.1))
       lines(wld, col=rgb(0,0,0,0.3), lwd=0.1)
     })
dev.off()

# Regional maps of the potential areas
pdf(file.path(mydir, "outmaps/Regional_ZonesPotentielles_SansVignoblesActuelsInclus.pdf"), width=12, height = 6)
zoomtab <- read.csv("ZoomTable.csv")
for(i in 1:nrow(zoomtab))
{  
  myext <- ext(unlist(zoomtab[i,-1]))
  plot(pot2bin, legend=F, colNA="lightskyblue", ext=myext, nc=3,fun=function() lines(wld, col=rgb(0,0,0,0.3)))
}
dev.off()


