# Liste de commande pour extraire les données des indices agroclimatiques 
# des pixels  TerraClimat (rasters "AgroNetCDF")aux points de la VGDB
# Crée par Benjamin Bois en 2023
# Mises à jour :
#   B. Bois, le 31 juillet 2025
#   - généralisation pour extraire directement les données des toutes les périodews
 
library(ncdf4)
library(MASS)
library(parallel)
library(rgdal)
library(maptools)
library(terra)

setwd("E:/Benjamin/_Recherches_Dijon/_AAA_TerraClimate/AgroNetCDF/")

# Winegrape Vineyard Geodatabase

vgdbpath <- "E:/SIG/world/Wine_Vineyards_GeoDB/v2/Points_Winegrape_Vineyards_Geodatabase_v2.1.gpkg"
vgdb.version <- sub("\\.gpkg","",sub("*.*Geodatabase_v","", vgdbpath))
mypts <- vect(vgdbpath)

mypts$CN_REG <- paste(mypts$CN, mypts$WINE_REGION, sep="_")

plot(mypts)
mypts
dim(geom(mypts))
x <- geom(mypts)[,"x"]
y <- geom(mypts)[,"y"]

####### Extraction Indices agroclimatiques actuels  #######
  # Liste des repertoires (donc des périodes)
mync_dirs <- sub("^\\./","",grep("TerraClimate_",list.dirs(),value=T))


# BOUCLE POUR TRAITER TOUS LES DOSSIERS N'AVOIR QU'UN GROS FICHIERS AVER TOUTES 
# LES DONNEES
mync_dir <- mync_dirs[1]

# Liste des variables / NetCDF à extraire
MyVars <-  c('CI_Avg','DI_Avg','DORT_Avg','DSP_Avg','DST_Avg','FSP_Avg','FST_Avg','GSP_Avg','GSP49_Avg','GST_Avg','GST49_Avg','HI_Avg','HSI_Avg','KOPPEN_Avg','Q_Avg','RipeR_Avg','SFR_Avg','WCR_Avg','WWB_Avg')


bigout <- c() # Un tableau pour compiler les données de tous les indices et tous les scenarii
for(mync_dir in mync_dirs)
{  
  myncs <- list.files(pattern=".nc$", path = mync_dir)
  mync1 <- grep(paste0("^",MyVars[1]), myncs, value=T)
  mysuff <- sub(MyVars[1],"",mync1)
  myncs <- paste0(MyVars, mysuff)
  myncinfo <- do.call(rbind,strsplit(myncs, "_"))
  myvar <- myncinfo[,1]
  mystat <- myncinfo[,2]
  myscen <- myncinfo[,3]
  
  i <- 1
  if(exists(x = "out")) rm(out) ; gc()
  for(i in 1:length(myncs))
  {
      mync <- rast(paste0(mync_dir,"/",myncs[i]))
      if(!exists(x = "out")){
        out <- terra::extract(mync,mypts, cells=T)
        out <- data.frame(pixID=out$cell, x=x,y=y, CN_REG=mypts$CN_REG, scenario=myscen[i])
      } 
      out <- cbind(out,terra::extract(mync,mypts)[,2])
      colnames(out)[ncol(out)] <- paste0(myvar[i],"_",mystat[i])
  }
  xy <- paste(out$x, out$y)
  
  Dupli <- duplicated(xy)
  sum(Dupli)
  out <- out[!Dupli,]
  
  bigout <- rbind(bigout, out)
}

DirOut <- paste0("../Extraction_TerraClimatPoints/VGDB_v",vgdb.version)
dir.create(DirOut)
NameOut <- paste0("AgroIndices_TerraClimate_VGBD_Pts_v",vgdb.version,"AllPeriods.csv")
data.table::fwrite(x = bigout, paste0(DirOut,"/",NameOut), row.names = F)  
