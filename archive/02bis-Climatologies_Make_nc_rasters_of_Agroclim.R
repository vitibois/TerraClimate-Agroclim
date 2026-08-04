# 2023-11-10, B. Bois - extraction de tous les pixels de la vgdb
# par 100 lots 
library(data.table)
library(terra)
library(ncdf4)
source("E:/R_functions/koppen.r")
setwd("E:/Benjamin/_Recherches_Dijon/_AAA_TerraClimate")


# "Fake" Year for +2 and +4 °C. 
# In fact  the  years corresponding to +2 and +4°C from prin-industrial level 
# of the RCP8.5 ensemble temperature simualtion using 39 CMIP5 GCM 
# Wang, X., Jiang, D., & Lang, X. (2018). Climate Change of 4°C Global Warming above Pre-industrial Levels. Advances in Atmospheric Sciences, 35(7), 757‑770. https://doi.org/10.1007/s00376-018-7160-4
# Directories according to the type of data we want to select 
# mync_dir <- "TerraClimate_2c_climatologies20y"
# mync_prefix <-"TerraClimate2C_"
# scenario <- "plus2"
# 
# mync_dir <- "TerraClimate_4c_climatologies20y"
# mync_prefix <-"TerraClimate4C_"
# scenario <- "plus4"

mync_dir <- "TerraClimate_obs8110_climatologies20y"
mync_prefix <-"TerraClimate19812010_"
scenario <- "CLIMOBS19812010"

mync_dir <- "TerraClimate_obs0120_climatologies20y"
mync_prefix <-"TerraClimate20012020_"
scenario <- "CLIMOBS20012020"

mync_dir <- "TerraClimate_obs6190_climatologies30y"
mync_prefix <-"TerraClimate19611990_"
scenario <- "CLIMOBS19611990"


# Creation d'un raster
var <- "tmin"
rast_TC <- terra::rast(paste0("E:/Benjamin/_Recherches_Dijon/_AAA_TerraClimate/Download_TerraClimate/climatologies/",mync_dir,"/",mync_prefix,var,".nc"))
rNA <- rast_TC[[1]]
values(rNA) <- NA
rNA

myfun <- mean
funname <- "Avg"


if(scenario == "OBS") year2discard = 2000
mypath <- "E:/Benjamin/_Recherches_Dijon/_AAA_TerraClimate/AgroNetCDF/"
dir.create(paste0(mypath,mync_dir), showWarnings = F)

# List of tiles
tilesdf <- read.table(paste0("E:/Benjamin/_Recherches_Dijon/_AAA_TerraClimate/Download_TerraClimate/climatologies/",mync_dir,"/AAA_TilesLists_TerraClimate.txt"), header=T)

i_group <- tilesdf$ID[1]

  # A first loop to load all agrocliamtic indices for all pixels 
deb <- Sys.time()
#  i_group <- tilesdf$ID[1]
agrostat <- c()
for(i_group in tilesdf$ID)
{
  agro <- data.table::fread(paste0("CLIMATOLOGIES_TC_CSV_DATA/",mync_dir,"/",mync_prefix,"TC_csv_data_Tile_",i_group,".csv"))
  agrostat <- rbind(agrostat,agro)
  rm(agro);gc()
}
fin <- Sys.time()
fin-deb

##### 2check #####
period <- agrostat$year[1]

deb2 <- Sys.time()
  # A loop to create NetCDF
myv <- "KOPPEN"
myvs <- c('HI','CI','DI','KOPPEN','GST','GSP','GST49','GSP49','DSP','DST','RipeR','FSP','FST','DORT','SFR','WCR','HSI', 'WWB', 'Q')
for(myv in myvs)
{
  r <- rNA
  agrostat$KOPPEN
  r[agrostat$pixID] <- as.data.frame(agrostat)[,myv]
  #plot(r)
  # Create ncdf
    # path and file name, set dname
  ncpath <- paste0(mypath,mync_dir,"/")
  ncname <- paste0(myv,"_",funname,"_",scenario,"_",period)  
  ncfname <- paste(ncpath, ncname, ".nc", sep="")
  dname <- myv  
  longname <- paste0(myv," ",funname," - dataset: ",scenario," period ",period)  
    # Write NetCDF
  writeCDF(r, ncfname, dname, longname=longname, unit="", split=FALSE, overwrite=T)
  
  cat(paste0("\nFin de loop pour variable ",myv," - ",date(),"\n"))
}
fin2 <- Sys.time()
fin2-deb2

