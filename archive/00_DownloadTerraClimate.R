# Liste de commandes pour telecharger les NCDF TerraClimate
# Créé par Benjamin BOIS, 2025-05-14, à partir du code élaboré par Sebastien Nicolas (2023)
# MAJ : juil 2025, B. Bois

setwd("E:/Benjamin/_Recherches_Dijon/_AAA_TerraClimate/Download_TerraClimate/")
# http://thredds.northwestknowledge.net:8080/thredds/catalog/TERRACLIMATE_ALL/data/TerraClimate_ws_2024.nc
# http://thredds.northwestknowledge.net:8080/thredds/catalog/TERRACLIMATE_ALL/data/catalog.html?dataset=TERRACLIMATE_ALL_SCAN/data/TerraClimate_pet_2015.nc

 ### -----  partie téléchargement ----
terraurl <- "https://climate.northwestknowledge.net/TERRACLIMATE-DATA/"
terraprefix <- "TerraClimate"

yearmin <-2001
yearmax <- 2020

TerraClimate_base_url <- function(year,var){
   paste0(terraurl, terraprefix,"_",var,"_",year,".nc")
 }

 # Variables à télécharger :
 MyVars <- c("tmin","tmax","vpd","pet","ppt")
 
i <- 2001
j <- 'tmin'
for(i in yearmin:yearmax){
 #for(j in c("tmin","tmax","ws","vpd","pet","srad","ppt","def","aet","PDSI"))
 deb <- Sys.time()

 for(j in MyVars)
      downloader::download(url = TerraClimate_base_url(year = i, var=j),
     destfile = paste0(terraprefix,"_",j,"_",i,".nc"),
     mode="wb")
 cat("Terminé pour année ",i," - ",date(),'\nTemps Ecoulé pour un loop: ',Sys.time()-deb,"\n")
}

