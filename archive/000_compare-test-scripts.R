setwd("E:/Benjamin/_Recherches_Dijon/_AAA_TerraClimate/2023_WWV_ANALYSIS/TC_CSV_DATA/AATEST/")
 do <- as.data.frame(data.table::fread("TC_csv_data_TerraClimate_Obs_Tile_58.csv"))
 dn <- as.data.frame(data.table::fread("ATC_csv_data_TerraClimate_Obs_Tile_58.csv"))

 do <- do[order(do$pixID, do$year),]
 dn <- dn[order(dn$pixID, dn$year),]
 ech <- sample(1:nrow(do), 2000)
 myvar <- "KOPPEN"
 plot(do[ech,myvar], dn[ech,myvar], col=rgb(0,0,0,0.2), pch=20);abline(a=0,b=1, col="red")

setwd("E:/Benjamin/_Recherches_Dijon/_AAA_TerraClimate/2023_WWV_ANALYSIS/CLIMATOLOGIES_TC_CSV_DATA/TerraClimate_2c_climatologies20y")
 do <- as.data.frame(data.table::fread("TerraClimate2C_ATC_csv_data_Tile_58.csv"))
 dn <- as.data.frame(data.table::fread("TerraClimate2C_TC_csv_data_Tile_58.csv"))
 names(do)
 do <- do[order(do$pixID, do$year),]
 dn <- dn[order(dn$pixID, dn$year),]
 ech <- sample(1:nrow(do), 2000)
 myvar <- "FSP"
 plot(do[ech,myvar], dn[ech,myvar], col=rgb(0,0,0,0.2), pch=20);abline(a=0,b=1, col="red")
 
 setwd("E:/Benjamin/_Recherches_Dijon/_AAA_TerraClimate/2023_WWV_ANALYSIS/TC_CSV_DATA/")
 d <- as.data.frame(data.table::fread("TC_csv_data_TerraClimate_Obs_Tile_1.csv"))
head(d) 
summary(d) 
  source("E:/R_functions/koppen.r") 
koppen(table.codes.only = T)
mylet <- koppen(table.codes.only = T)[d$KOPPEN,"letcode"]
sort(summary(factor(mylet)))
