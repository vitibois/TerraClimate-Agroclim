# 2023-11-10, B. Bois - extraction de tous les pixels de la vgdb
# par 100 lots 
library(data.table)
library(terra)
library(ncdf4)

setwd("E:/Benjamin/_Recherches_Dijon/_AAA_TerraClimate")

# List of tiles
tilesdf <- read.table("AAA_TilesLists_TerraClimate.txt", header=T)


# Creation d'un raster
MyVar <- "tmin"
year <- "2020"
rast_TC <- terra::rast(paste0("E:/Benjamin/_Recherches_Dijon/_AAA_TerraClimate/Download_TerraClimate/",MyVar,"_",year,".nc"))
rNA <- rast_TC[[1]]
values(rNA) <- NA
rNA


# The "mode" function --> for Koppen
myfun <-  function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
funname <- "Avg" # We use the Avg acronym for Koppen, but it's the Mode indeed that is used

myfun <- function(x) quantile(x, probs = 0.25)
funname <- "Q25"

myfun <- function(x) mean(x, na.rm=F)
funname <- "Avg"


scenario <- "OBS"
period <- "2001-2020"
if(scenario == "OBS") year2discard = 2000
mypath <- "E:/Benjamin/_Recherches_Dijon/_AAA_TerraClimate/AgroNetCDF/"

 myv <- "KOPPEN"
#myvs <- c('HI','CI','DI','GST','GSP','GST49','GSP49','DSP','DST','RipeR','FSP','FST','DORT','SFR','WCR','HSI')
myvs <- c('HI','CI','DI', 'Q','WWB','GST','GSP','GST49','GSP49','DSP','DST','RipeR','FSP','FST','DORT','SFR','WCR','HSI')
deb <- Sys.time()
  # Compter 5 to 8 min par indicateur, selon le calcul réalisé
for(myv in myvs)
{
#  i_group <- tilesdf$ID[1]
  agrostat <- c()
  for(i_group in tilesdf$ID)
  {
    df <- data.table::fread(paste0("./TC_CSV_DATA/TC_csv_data_TerraClimate_Obs_Tile_",i_group,".csv"))
    df <- df[year!=year2discard]
    agro <- df[,.(agrostat=myfun(get(myv))), by=pixID]
    agrostat <- rbind(agrostat,agro)
    rm(agro);gc()
  }
  r <- rNA
  r[agrostat$pixID] <- agrostat$agrostat

  # Create ncdf
    # path and file name, set dname
  ncpath <- mypath
  dir.create(ncpath, showWarnings = F)
  ncname <- paste0(myv,"_",funname,"_",scenario,"_",period)  
  ncfname <- paste(ncpath, ncname, ".nc", sep="")
  dname <- myv  
  longname <- paste0(myv," ",funname," - dataset: ",scenario," period ",period)  
    # Write NetCDF
  writeCDF(r, ncfname, dname, longname=longname, unit="", split=FALSE, overwrite=T)
  
  cat(paste0("\nFin de loop pour variable ",myv," - ",date(),"\n"))
}
fin <- Sys.time()
fin-deb

r <- rast(ncfname)
r
plot(r)
