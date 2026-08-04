# Liste de commandes pour calculer des climatologies sur une période donnée à 
# partie de fichiers d'années individuelles
# Créé par Benjamin BOIS, 2025-05-14, à partir du code élaboré par Sebastien Nicolas (2023)
# MAJ : juil 2025, B. BOis

library(terra)
library(sf)
library(tidyverse)
library(tidyterra)

setwd("E:/Benjamin/_Recherches_Dijon/_AAA_TerraClimate/Download_TerraClimate")
yearmin <- 2001
yearmax <- 2020
avgyear <- round(mean(yearmin:yearmax))


# Répertoire de sortie : 
outdir <- paste0("TerraClimate_obs",substr(yearmin,3,4),substr(yearmax,3,4),"_climatologies20y")
ncpath <- paste0("climatologies/",outdir)
dir.create(ncpath)
prefix <- paste0("TerraClimate",yearmin,yearmax)

# Liste des fichiers
NClist <- list.files(pattern="nc$") ; NClist

  # Pas utile pour l'instant...mais si on veut gagner de l'espace on peut 
  # arrondir et écrire en Interger les pet et ppt
MyVarsDecPlaces <- c(tmin=2, tmax=2, vpd=2, pet=0, ppt=0)
MyVars <- names(MyVarsDecPlaces)

v <- MyVars[1]

for(v in MyVars[-1])
{
  MyFiles <- grep(paste0("_",v,"_"),NClist,value=T)
  MyR <- rast(MyFiles)
  crs(MyR) <- "epsg:4326"
  #plot(MyR,1)
  i <- 1
  MyRAvg <- c()
  for(i in 1:12)
  {
    MyIndex <- grep(paste0("_",i,"$"),names(MyR))
    MyRi <- subset(MyR,MyIndex)
    MyRAvg <- c(MyRAvg, mean(MyRi))
  }
  MyRAvg <- rast(MyRAvg)
  ncname <- paste0(prefix,"_",v)  
  ncfname <- paste0(ncpath,"/", ncname, ".nc")

    # Write NetCDF
  writeCDF(MyRAvg,ncfname, varname=v, overwrite=T)
}

