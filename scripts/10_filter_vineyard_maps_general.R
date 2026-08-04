
# 2023-11-10, B. Bois

# Mapping potentially suitable wine-growing areas in 2001-2020 and 2050 (+2C),
# based on general criteria defined with Kees in September 2023.

library(data.table)
library(terra)
library(ncdf4)
library(ggplot2)

source("R/config.R")
setwd(file.path(data_root, "2023_WWV_ANALYSIS"))

###### OBS one criterion map  #####
mync_dir_obs <- file.path(data_root, "AgroNetCDF")

#### Thresholds
th <- read.csv("Criterions_global.csv")
th <- th[!is.na(th$Rule),]
myrules <- paste(paste(th$Index, th$Statistic.over.time, sep="_"),th$Sign,th$Rule)

#### Identification of areas according to one criterion
i <- 1
myvar <- "HSI_Avg"
myrule  <- "r < 33"
r <- rast(paste0(mync_dir_obs, "/",myvar,"_OBS_2001-2020.nc"))
rrules <- eval(parse(text = myrule))  
rout <- r*rrules
rout[rout < 10] <- NA
plot(rout)
plot(rrules)
terra::plot(rrules, ext=ext(c(-20,40,30,70)))
wld <- vect(file.path(sig_root, "world/World_AdmBoundaries_Countries/world-administrative-boundaries.shp"))
vgdb <- vect(file.path(sig_root, "world/Vineyards_GeoDB_2015/v1.2.3/Vineyard_GeoDB_v1.2.3.shp"))
vgdb.pt <- vect(file.path(data_root, "VGDB_Points_LouisDelelee/pts_global_v4.0.shp"))
plot(rrules, ext=ext(c(-20,40,30,70)))
plot(wld, add=T)
#plot(vgdb, add=T, col=rgb(0,0,1,0.5))

# Countries of the world
wld <- vect(file.path(sig_root, "world/World_AdmBoundaries_Countries/world-administrative-boundaries.shp"))


#
#   # Google Earth hybrid basemap
# hrge <- rast("GE_hybrid_0.04.tif")
# crs(hrge) <- crs(rout)
# hrge <- crop(hrge, rout)
# hrge <- resample(hrge, rout)
# hrge <- mask(hrge, r)
#
#   # Lighten the raster (more brightness)
# palege <- hrge+50
# palege[palege > 255] <- 255
# plotRGB(palege, ext=ext(c(-10,20,35,60)))
#   # save the GE_hybrid basemap
# writeRaster(palege, "PaleGE_hybrid_TerraClimResoluion.tif", overwrite=T)
palege <- rast("PaleGE_hybrid_TerraClimResoluion.tif")

# Areas currently planted with vineyards
vpts <-  vect(file.path(data_root, "VGDB_Points_LouisDelelee/pts_global_v4.0.shp"))


###### OBS identification areas #####
mync_dir_obs <- file.path(data_root, "AgroNetCDF")

#### Thresholds
th <- read.csv("Criterions_global.csv")
th <- th[!is.na(th$Rule),]
paste(paste(th$Index, th$Statistic.over.time, sep="_"),th$Sign,th$Rule)

#### Identification of areas
i <- 1
myvar <- paste(th$Index[i], th$Statistic.over.time[i], sep="_")
rout <- rast(paste0(mync_dir_obs, "/",myvar,"_OBS_2001-2020.nc"))
rout <- rout*0+1
plot(rout)

pdf("outmaps/CarteCritereIndiv_OBS2001-2020.pdf", width = 8, height = 6)
  for(i in 1:nrow(th))
  {
    myvar <- paste(th$Index[i], th$Statistic.over.time[i], sep="_")
    r <- rast(paste0(mync_dir_obs, "/",myvar,"_OBS_2001-2020.nc"))
    myrule <-  paste("r",th$Sign,th$Rule)[i]
    r2 <- eval(parse(text = myrule))  
      mymaxcell <- 2*10^6
      myext <- ext(c(-180,180,-75,85))
      plotRGB(palege,  maxcell=mymaxcell, ext=myext, mar=c(3.1, 3.1, 2.1, 7.1),main=paste("Recent : 2001-2020\n\n",myvar,th$Sign[i],th$Rule[i]))
      plot(wld, add=T, lwd=0.1, col=rgb(0,0,0,0.5))
      plot(r2, add=T, alpha=0.5, maxcell=mymaxcell, ext=myext)
      plot(vpts, add=T, cex=0.05, col=rgb(0,0,0,0.2))
      sbar(d=5000,labels = c(0,2500,"5000 km"), halo = T,lonlat = T ,  ticks = T, cex=0.7)
    rout <- rout*r2
    r2[r2==0] <- NA
    routval <- r*r2
    writeRaster(routval, paste0("rasters_zonespotentielles/2001-2020/",myvar,"_ZonesPotentielles_OBS2001-2020.tif"), overwrite=T)
  }
dev.off()
plot(rout)
writeRaster(rout, "rasters_zonespotentielles/AAA_ZonesPotentielles_OBS2001-2020.tif", overwrite=T)


###### PLUS 2 C identification areas #####
mync_dir_obs <- file.path(data_root, "AgroNetCDF/TerraClimate_2c_climatologies20y")

#### Thresholds
th <- read.csv("Criterions_global.csv")
th <- th[!is.na(th$Rule),]
th$Statistic.over.time <- sub("Mode","Avg",th$Statistic.over.time)

#### Identification of areas
i <- 1
myvar <- paste(th$Index[i], th$Statistic.over.time[i], sep="_")
rout <- rast(paste0(mync_dir_obs, "/",myvar,"_plus2_2042.nc"))
rout <- rout*0+1

pdf("outmaps/CarteCritereIndiv_plus2C.pdf", width = 8, height = 6)
  for(i in 1:nrow(th))
  {
    myvar <- paste(th$Index[i], th$Statistic.over.time[i], sep="_")
    r <- rast(paste0(mync_dir_obs, "/",myvar,"_plus2_2042.nc"))
    myrule <-  paste("r",th$Sign,th$Rule)[i]
    r2 <- eval(parse(text = myrule))  
      mymaxcell <- 2*10^6
      mymaxcell <- 10^6
      myext <- ext(c(-180,180,-75,85))
      plotRGB(palege,  maxcell=mymaxcell, ext=myext, mar=c(3.1, 3.1, 2.1, 7.1),main=paste("+2°C\n\n",myvar,th$Sign[i],th$Rule[i]))
      plot(wld, add=T, lwd=0.1, col=rgb(0,0,0,0.5))
      plot(r2, add=T, alpha=0.5, maxcell=mymaxcell, ext=myext)
      plot(vpts, add=T, cex=0.05, col=rgb(0,0,0,0.2))
      sbar(d=5000,labels = c(0,2500,"5000 km"), halo = T,lonlat = T ,  ticks = T)
    rout <- rout*r2
    r2[r2==0] <- NA
    routval <- r*r2
    writeRaster(routval, paste0("rasters_zonespotentielles/plus2C/",myvar,"_ZonesPotentielles_plus2C.tif"), overwrite=T)
    #plot(routval):plot(wld,add=T)
  }
dev.off()
plot(rout)
writeRaster(rout, "rasters_zonespotentielles/AAA_ZonesPotentielles_plus2C.tif", overwrite=T)

###### PLUS 4 C identification areas #####
mync_dir_obs <- file.path(data_root, "AgroNetCDF/TerraClimate_4c_climatologies20y")

#### Thresholds
th <- read.csv("Criterions_global.csv")
th <- th[!is.na(th$Rule),]
th$Statistic.over.time <- sub("Mode","Avg",th$Statistic.over.time)

#### Identification of areas
i <- 1
myvar <- paste(th$Index[i], th$Statistic.over.time[i], sep="_")
rout <- rast(paste0(mync_dir_obs, "/",myvar,"_plus4_2084.nc"))
rout <- rout*0+1

pdf("outmaps/CarteCritereIndiv_plus4C.pdf")
  for(i in 1:nrow(th))
  {
    myvar <- paste(th$Index[i], th$Statistic.over.time[i], sep="_")
    r <- rast(paste0(mync_dir_obs, "/",myvar,"_plus4_2084.nc"))
    myrule <-  paste("r",th$Sign,th$Rule)[i]
    r2 <- eval(parse(text = myrule))  
      mymaxcell <- 2*10^6
      myext <- ext(c(-180,180,-75,85))
      plotRGB(palege, maxcell=mymaxcell, ext=myext, mar=c(3.1, 3.1, 2.1, 7.1),main=paste("+4°C\n\n",myvar,th$Sign[i],th$Rule[i]))
      plot(wld, add=T, lwd=0.1, col=rgb(0,0,0,0.5))
      plot(r2, add=T, alpha=0.5, maxcell=mymaxcell, ext=myext)
      plot(vpts, add=T, cex=0.05, col=rgb(0,0,0,0.2))
      sbar(d=5000,labels = c(0,2500,"5000 km"), halo = T,lonlat = T ,  ticks = T)
    rout <- rout*r2
    r2[r2==0] <- NA
    routval <- r*r2
    writeRaster(routval, paste0("rasters_zonespotentielles/plus4C/",myvar,"_ZonesPotentielles_plus4C.tif"), overwrite=T)
    #plot(routval):plot(wld,add=T)
  }
dev.off()
plot(rout)
writeRaster(rout, "rasters_zonespotentielles/AAA_ZonesPotentielles_plus4C.tif", overwrite=T)



##### CREATE BLUE WATER RESOURCE ESTIMATES
  # OBS
rQ <- rast(file.path(data_root, "AgroNetCDF/Q_Avg_OBS_2001-2020.nc"))
rWWB <- rast(file.path(data_root, "AgroNetCDF/WWB_Avg_OBS_2001-2020.nc"))
bwat <- rWWB-200+rQ
plot(bwat)
writeCDF(bwat, file.path(data_root, "AgroNetCDF/BlueWat_Avg_OBS_2001-2020.nc"), overwrite=T, varname="bwat", unit="mm", longname="Blue Water Estimate")
  # PLUS 2 C
rQ <- rast(file.path(data_root, "AgroNetCDF/TerraClimate_2c_climatologies20y/Q_Avg_plus2_2042.nc"))
rWWB <- rast(file.path(data_root, "AgroNetCDF/TerraClimate_2c_climatologies20y/WWB_Avg_plus2_2042.nc"))
bwat <- rWWB-200+rQ
plot(bwat)
writeCDF(bwat, file.path(data_root, "AgroNetCDF/TerraClimate_2c_climatologies20y/bwat_Avg_plus2_2042.nc"), overwrite=T, varname="bwat", unit="mm", longname="Blue Water Estimate")
  # PLUS 4 C
rQ <- rast(file.path(data_root, "AgroNetCDF/TerraClimate_4c_climatologies20y/Q_Avg_plus4_2084.nc"))
rWWB <- rast(file.path(data_root, "AgroNetCDF/TerraClimate_4c_climatologies20y/WWB_Avg_plus4_2084.nc"))
bwat <- rWWB-200+rQ
plot(bwat)
writeCDF(bwat, file.path(data_root, "AgroNetCDF/TerraClimate_4c_climatologies20y/bwat_Avg_plus4_2084.nc"), overwrite=T, varname="bwat", unit="mm", longname="Blue Water Estimate")
