# 2023-11-10, B. Bois - extraction de tous les pixels de TC er calcul d'indices
# par 100 lots 
setwd("E:/Benjamin/_Recherches_Dijon/_AAA_TerraClimate")

pacman::p_load(
  dplyr,
  tidyverse,
  ggplot2,
  sf,
  terra,
  tidyterra,
  geodata
)

# Sources
source("E:/R_functions/gccm.r")
source("E:/R_functions/koppen.r")




# Créations des lots
var <- "tmin"
year <- "2020"
rast_TC <- terra::rast(paste0("E:/Benjamin/_Recherches_Dijon/_AAA_TerraClimate/Download_TerraClimate/",var,"_",year,".nc"))
r <- rast_TC[[1]]
r.coords <- crds(r)
pixID <- terra::extract(r, r.coords, cells=T)[,1]
r.coords <- cbind(r.coords, pixID)
rm(pixID)
head(r.coords)
range(r.coords[,1])
range(r.coords[,2])
mybreaks <- c(seq.int(0,nrow(r.coords), by=round(nrow(r.coords)/99)),nrow(r.coords))
mygroups <- cut(1:nrow(r.coords),mybreaks, labels=F)
groups_ID <- unique(mygroups)
n_groups <- length(groups_ID)
write.table(data.frame(ID=groups_ID, Breaks=mybreaks[-1]), "AAA_TilesLists_TerraClimate.txt", row.names=F)
##### Extraction d'un groupe
i_group <- 1
i_group <- mygroups[which(r.coords[,"y"] < 0)[1]] # Test a group with both northern and southern hemisphere data


####### LOOOOOOOOOP #######
Deb <- Sys.time()
for(i_group in 1:n_groups)
{
  iDeb <- Sys.time()
  r.crds <- r.coords[mygroups==i_group,]
  mypixs <- r.crds[,'pixID']
  # var <- "tmin"
  # year <- 2020
  dat <- c()
    # Il faut une année de plus car pour l'hémisphère sur, on "perd" les 6 premiers mois
  for(year in 2000:2020)
  {
    rast_dat <- r.crds
    # extraction donnée .nc
    for(var in c("tmin","tmax","ppt","pet"))
    {
      rast_TC <- terra::rast(paste0("E:/Benjamin/_Recherches_Dijon/_AAA_TerraClimate/Download_TerraClimate/",var,"_",year,".nc"))
      rast_dat <- cbind(rast_dat, terra::extract(rast_TC,r.crds[,1:2]))
      ## Check the adress pixID works
      # r <- rast_TC$tmin_1
      # values(r) <- NA
      # dim(r)
      # r[r.crds[,"pixID"]] <- rast_dat[,"tmin_1"]
      # plot(r)
      # 
    }
    dat <- rbind(dat, cbind(year=year, rast_dat))
  }

    # Some data are missing for given pixels/years  :-( !!!!
    # We remove them from the data set
  MyColSums <- Rfast::colsums(as.matrix(dat))
  whichnas <- which(is.na(MyColSums))
  if(length(whichnas)>0)
  {
    rows2remove <- c()
    for(nas in whichnas) rows2remove <- c(rows2remove, which(is.na(dat[,nas])))
    dat <- dat[-rows2remove,]
  }  
  # Calcul des indices agroclimatiques
  TNInd <- grep("tmin", names(dat))
  TXInd <- grep("tmax", names(dat))
  TM <- (dat[,TNInd]+dat[,TXInd])/2
  names(TM) <- sub("tmin","tavg", names(TM))
  dat <- cbind(dat,TM)
  rm(TM) ; gc()

  isSH <- which(dat$y < 0)
  
  # GET Northern hemisphere data
  sdat <- dat[-isSH,] 
  if(length(isSH)==0) sdat <- dat
  TNInd <- grep("tmin", names(sdat))
  TXInd <- grep("tmax", names(sdat))
  RRInd <- grep("ppt", names(sdat))
  ETInd <- grep("pet", names(sdat))
  LATInd <- grep("^y$", names(sdat))
  TMInd <- grep("tavg", names(sdat))
  
  # Add southern hemisphere data, shifted by 6 months
  if(length(isSH) > 0)
  {
    SHsdat <- dat[isSH,] 
    # Shifting data  from July to Dec year N to N+1 so that calculation index 
    # uses data from the previous year in Sourthern hemisphere. 
    # i.e. Growing season temperature of year 2020 is the temperature average
    # from [Oct to Dec of year 2019 and Jan to March/April of year 2020]
    falseyears <- SHsdat$year-1
    falseyears[falseyears == min(falseyears)] <- min(SHsdat$year)
    yearpix <- paste(SHsdat$pixID,SHsdat$year)
    false_yearpix <- paste(SHsdat$pixID,falseyears)
    matchyearpix <- match(false_yearpix,yearpix)
    SHsdat[,c(TNInd[7:12], TXInd[7:12], RRInd[7:12], ETInd[7:12])] <- 
      SHsdat[matchyearpix,c(TNInd[7:12], TXInd[7:12], RRInd[7:12], ETInd[7:12])]
    
    # SHsdat <- SHsdat[order(SHsdat$pixID, SHsdat$year),]
    # SHsdat[,c(TNInd[7:12], TXInd[7:12], RRInd[7:12], ETInd[7:12])] <-
    #   apply(SHsdat[,c(TNInd[7:12], TXInd[7:12], RRInd[7:12], ETInd[7:12])], 2, function(x)
    #     unlist(tapply(x, INDEX = SHsdat$pixID, function(y) y <- y[c(1,1:(length(y)-1))])))

    # Changing indices 
    SHsdat[,TNInd] <- SHsdat[,TNInd[c(7:12,1:6)]]
    SHsdat[,TXInd] <- SHsdat[,TXInd[c(7:12,1:6)]]
    SHsdat[,RRInd] <- SHsdat[,RRInd[c(7:12,1:6)]]
    SHsdat[,ETInd] <- SHsdat[,ETInd[c(7:12,1:6)]]
    SHsdat[,TMInd] <- SHsdat[,TMInd[c(7:12,1:6)]]
    
    # Binding south and north hemisphere data
    sdat <- rbind(sdat, SHsdat)
#    rm(SHsdat) ; gc()
  }
  
  # Adding sept-dec precipitations from the previous year to calculate winter
  # water balance
  falseyears <- sdat$year-1
  range(falseyears)
  falseyears[falseyears == min(falseyears)] <- min(sdat$year)
  yearpix <- paste(sdat$pixID,sdat$year)
  false_yearpix <- paste(sdat$pixID,falseyears)
  matchyearpix <- match(false_yearpix,yearpix)
  sum(is.na(matchyearpix))
  WINTERHYDRO <- sdat[matchyearpix,c(RRInd[10:12], ETInd[10:12])]
  names(WINTERHYDRO) <- paste0("m1_",names(WINTERHYDRO))
  sdat <- cbind(sdat,WINTERHYDRO)
  M1RRInd <- grep("m1_ppt", colnames(sdat))
  M1ETInd <- grep("m1_pet", colnames(sdat))
  
  # As matrix ... for quicker calculation
  sdat <- as.matrix(sdat)
  
      # Each index of the Multicrieria classification system is calculated separately
  # trying to limitate as much as possible calculation time using 
  # matrix calculation
  
  # Matrix of number of days in the month
  ndaysmat <- t(t(sdat[,RRInd[4:9]])*0+1*c(30,31,30,31,31,30))
  if(length(isSH) > 0)
    ndaysmat[isSH,] <- t(t(ndaysmat[isSH,])*0+1*c(31,30,31,31,28,30))
  
  # HI 
  khi <- k.hi.hall(lat = sdat[,LATInd])
  HI <- sdat[,c(TMInd[4:9],TXInd[4:9])]-10
  HI[HI < 0] <- 0
  HI <- (HI[,1:6]+HI[,7:12])/2
  HI <- HI * ndaysmat
  HI <- rowSums(HI)*khi
  
  #CI
  CI <- sdat[,TNInd[9]]
  summary(CI)
  
  #DI
  # Soil evaporation coefficient
  jpm5 <- sdat[,RRInd[4:9]]/5
  jpm5test <- jpm5 > ndaysmat
  jpm5[jpm5test] <- ndaysmat[jpm5test]
  
  # Soil evaporation
  kdi <- c(0.1,0.3,0.5,0.5,0.5,0.5)
  es <- t(t(sdat[,ETInd[4:9]]/ndaysmat* jpm5) * (1-kdi))
  # grapevine transpiration
  tv <- t(t(sdat[,ETInd[4:9]]) * kdi)
  # Soil water balance
  rr <- sdat[,RRInd[4:9]]
  DI <- c()
  i <- 1
  for(i in 1:ncol(rr))
  {
    if(i==1)
    {
      DI <- 200 + rr[,i] - es[,i] - tv[,i]
      DI[DI > 200] <- 200
    } else {
      DI <- DI + rr[,i] - es[,i] - tv[,i] 
      DI[DI > 200] <- 200
    }
  }

  # Winter water balance (complementary to the dryness index)
  BHM0 <- DI
  BHM0[BHM0<0] <- 0
  rr <- sdat[,c(M1RRInd, RRInd[1:3])]
  et <- sdat[,c(M1ETInd,ETInd[1:3])]
  WWB <- c()
  for(i in 1:6)
  {
    if(i==1)
    {
      WWB <- DI + rr[,i] - et[,i] # water balance
      Q <- WWB-200 # Runoff (i.e. surplus)
      Q[Q < 0] <- 0  # Runoff (surplus) is not negative
      WWB[WWB > 200] <- 200 # Water balance can't be more than SWC
      WWB[WWB < 0] <- 0 # Water balance can't be negative
    } else {
      WWB <- WWB + rr[,i] - et[,i] # water balance
      Qm <- WWB-200 # Runoff (i.e. surplus)
      Qm[Qm < 0] <- 0  
      Q <- Q+Qm # Runoff (surplus) is not negative
      WWB[WWB > 200] <- 200 # Water balance can't be more than SWC 
      WWB[WWB < 0] <- 0 # Water balance can't be negative
    }
  }

  # HI DI and CI indices
  # GCCM <- t(apply(sdat[,c(TNInd[4:9], TXInd[4:9], RRInd[4:9], ETInd[4:9], LATInd, TMInd[4:9])],1, function(x){
  #   gccm(tn = x[1:6], tx = x[1:6+6], tm =  x[1:6+6*4+1],rr = x[1:6+6*2], 
  #        et0 = x[1:6+6*3],use.hall.khi.estimate = T, lat = x[6*4+1], val.only = T)    
  # }))
  
  # Test : 
  # ech <- sample(1:length(DI),1000)
  # summary(GCCM[,"DI"]) ; summary(DI)
  # summary(GCCM[,"HI"]) ; summary(HI)
  # summary(GCCM[,"CI"]) ; summary(CI)
  # plot(GCCM[ech,"DI"],DI[ech]);abline(a=0,b=1)
  # plot(GCCM[ech,"HI"],HI[ech]);abline(a=0,b=1)
  # plot(GCCM[ech,"CI"],CI[ech]);abline(a=0,b=1)
   
  # Koppen
  KOP <- koppen.fast(teta = sdat[,TMInd], rr = sdat[,RRInd])
  # Test : 
  # rr = sdat[1:20,RRInd]
  # teta = sdat[1:20,TMInd]
  # lats = abs(sdat[1:20,LATInd])
  # write.table(teta,"clipboard", sep="\t", row.names=F)
  
  
  # AVGST calculation
    # GST <- apply(sdat[,c(TMInd[4:10])],1, mean, na.rm=T)
    # GSR <- apply(sdat[,c(RRInd[4:10])],1, sum, na.rm=T)
  GST <- rowMeans(sdat[,c(TMInd[4:10])])
  GSP <- rowSums(sdat[,c(RRInd[4:10])])
  
  
  # AVGST49 calculation
    # GST49 <- apply(sdat[,c(TMInd[4:9])],1, mean, na.rm=T)
    # GSR49 <- apply(sdat[,c(RRInd[4:9])],1, sum, na.rm=T)
  GST49 <- rowMeans(sdat[,c(TMInd[4:9])])
  GSP49 <- rowSums(sdat[,c(RRInd[4:9])])
  
  #Dormant season temperature - Marco A.F. Conceicao tropical vineyard limit at 16°C
  DORT <- rowMeans(sdat[,c(TMInd[c(10:12,1:3)])])
  
  # Spring conditions (April-July ou SH : October-January)
  DST <- rowMeans(sdat[,c(TMInd[4:7])])
  DSP <- rowSums(sdat[,c(RRInd[4:7])])
  
  # Ripening conditions (September --> Botrytis risk / SH : March)
  RipeT <- sdat[,c(TMInd[9])]
  RipeR <- sdat[,c(RRInd[9])]
  
  # Fruit development conditions (July to September / SH: January to March)
  FST <- rowMeans(sdat[,c(TMInd[7:9])])
  FSP <- rowSums(sdat[,c(RRInd[7:9])])
  
  # Extreme indices
  SFR <- Rfast::rowMins(as.matrix(sdat[,c(TNInd[4:6])]), value=T)
  WCR <- Rfast::rowMins(as.matrix(sdat[,c(TNInd)]), value=T)
  HSI <- Rfast::rowMaxs(as.matrix(sdat[,c(TXInd)]), value=T)
  
  # 
  Agro <- cbind(sdat[,c("pixID","x","y","year")], HI, CI, DI, Q, WWB,  KOPPEN=KOP, GST, GSP, GST49, GSP49, DSP, DST, RipeR,  FSP, FST, DORT, SFR, WCR, HSI)
  #Agro <- cbind(sdat[,c("pixID","x","y","pixID","year")],GCCM,  KOPPEN=KOP, GST, GSP, GST49, GSP49, DSP, DST, RipeR, RipeT, FSP, FST, DORT,  SFR, WCR, HSI)
  names(Agro) <- sub("^CI.tmin..$","CI", names(Agro))
  
    # Writing data in csv file 
  data.table::fwrite(Agro, paste0("./TC_CSV_DATA/TC_csv_data_TerraClimate_Obs_Tile_",i_group,".csv"))
    
    # Fin de la boucle pour le tile i_group
  cat(paste0("\nFin de loop ",i_group," sur ",n_groups," - ",date(),"\n"))
  iFin <- Sys.time()
}
Fin <- Sys.time()
Fin-Deb
