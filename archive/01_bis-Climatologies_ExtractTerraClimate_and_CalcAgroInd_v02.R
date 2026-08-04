# 2023-11-10, B. Bois - extraction de tous les pixels de la vgdb
# par 100 lots et calcul d'indides agrocliamtiques
# POUR LES CLIMATOLOGIES (MOYENNES SUR 20 ou 30 ANS)

# Mises à jour
# BBois le 31 juillet 2027 : remise en forme du script
 
  setwd("E:/Benjamin/_Recherches_Dijon/_AAA_TerraClimate/CLIMATOLOGIES_TC_CSV_DATA/")
  
pacman::p_load(tidyverse,sf,terra,tidyterra,geodata)
  
  # Sources
  source("E:/R_functions/gccm.r")
  source("E:/R_functions/koppen.r")
  
  
  # "Fake" Year for +2 and +4 °C. 
  # In fact  the  years corresponding to +2 and +4°C from prin-industrial level 
  # of the RCP8.5 ensemble temperature simualtion using 39 CMIP5 GCM 
  # Wang, X., Jiang, D., & Lang, X. (2018). Climate Change of 4°C Global Warming above Pre-industrial Levels. Advances in Atmospheric Sciences, 35(7), 757‑770. https://doi.org/10.1007/s00376-018-7160-4
  # Directories according to the type of data we want to select 
  # mync_dir <- "TerraClimate_2c_climatologies20y"
  # mync_prefix <-"TerraClimate2C_"
  # year <- 2042

  # mync_dir <- "TerraClimate_4c_climatologies20y"
  # mync_prefix <-"TerraClimate4C_"
  # year <- 2084
  # 
  # mync_dir <- "TerraClimate_obs8110_climatologies20y"
  # mync_prefix <-"TerraClimate19812010_"
  # year <- 1995
  
  mync_dir <- "TerraClimate_obs0120_climatologies20y"
  mync_prefix <-"TerraClimate20012020_"
  year <- 2010
  
 # mync_dir <- "TerraClimate_obs6190_climatologies30y"
 # mync_prefix <-"TerraClimate19611990_"
 # year <- 1976
  
  

# Créations des lots de tableaux
var <- "tmin"
rast_TC <- terra::rast(paste0("E:/Benjamin/_Recherches_Dijon/_AAA_TerraClimate/Download_TerraClimate/climatologies/",mync_dir,"/",mync_prefix,var,".nc"))
r <- rast_TC[[1]]
r.coords <- crds(r, na.rm=T)
pixID <- terra::cellFromXY(r, r.coords)
r.coords <- cbind(r.coords, pixID)
head(r.coords)
range(r.coords[,1])
range(r.coords[,2])
mybreaks <- c(seq.int(0,nrow(r.coords), by=round(nrow(r.coords)/99)),nrow(r.coords))
mygroups <- cut(1:nrow(r.coords),mybreaks, labels=F)
groups_ID <- unique(mygroups)
n_groups <- length(groups_ID)
write.table(data.frame(ID=groups_ID, Breaks=mybreaks[-1]), paste0("E:/Benjamin/_Recherches_Dijon/_AAA_TerraClimate/Download_TerraClimate/climatologies/",mync_dir,"/AAA_TilesLists_TerraClimate.txt"), row.names=F)

##### Extraction d'un groupe
i_group <- 1
i_group <- mygroups[which(r.coords[,"y"] < 0)[1]] # Test a group with both northern and southern hemisphere data
  
Deb <- Sys.time()
for(i_group in 1:n_groups)
{
    r.crds <- r.coords[mygroups==i_group,]
    mypixs <- r.crds[,'pixID']
    # extraction donnée .nc
    var <- "tmin"
    rast_dat <- r.crds
    for(var in c("tmin","tmax","ppt","pet"))
    {
      rast_TC <-  terra::rast(paste0("E:/Benjamin/_Recherches_Dijon/_AAA_TerraClimate/Download_TerraClimate/climatologies/",mync_dir,"/",mync_prefix,var,".nc"))
      rast_dat <- cbind(rast_dat, terra::extract(rast_TC,r.crds[,1:2]))
      names(rast_dat) <- sub("Z1=","",names(rast_dat))
      ## Check the adress pixID works
      #  r <- rast_TC$`tmin_Z1=1`
      #  values(r) <- NA
      #  dim(r)
      #  r[r.crds[,"pixID"]] <- rast_dat[,"tmin_1"]
      #  plot(r)
      # 
    }
    dat <- cbind(year= year, rast_dat)
  
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
    
    #### NORTHERN HEMISPHERE ----
    sdat <- dat[-isSH,] 
    if(length(isSH)==0) sdat <- dat
    TNInd <- grep("tmin", names(sdat))
    TXInd <- grep("tmax", names(sdat))
    RRInd <- grep("ppt", names(sdat))
    ETInd <- grep("pet", names(sdat))
    LATInd <- grep("^y$", names(sdat))
    TMInd <- grep("tavg", names(sdat))
    
    if(length(isSH) > 0)
    {
      SHsdat <- dat[isSH,] 
      
      # Changing columns order (year start on July 1st in Southern Hemisphere)!
      SHsdat[,TNInd] <- SHsdat[,TNInd[c(7:12,1:6)]]
      SHsdat[,TXInd] <- SHsdat[,TXInd[c(7:12,1:6)]]
      SHsdat[,RRInd] <- SHsdat[,RRInd[c(7:12,1:6)]]
      SHsdat[,ETInd] <- SHsdat[,ETInd[c(7:12,1:6)]]
      SHsdat[,TMInd] <- SHsdat[,TMInd[c(7:12,1:6)]]
      
      # Binding south and north hemisphere data
      sdat <- rbind(sdat, SHsdat)
  #    rm(SHsdat) ; gc()
    }
    
    sdat <- as.matrix(sdat)
    
    # As Order of pixel might have been modified for subsets with both 
    # Northern hemisphere and Southern hemisphere data, we recalculate the location
    # of Southern Hemisphere lines
    isSH <- which(sdat[,"y"] < 0)
  
    # Each index of the Multicrieria classification system is calculated separately
    # trying to limitate as much as possible calculation time using 
    # matrix calculation
    
    # Matrix of number of days in the month 
    ndaysmat <- t(t(sdat[,RRInd])*0+1*c(31,28,31,30,31,30,31,31,30,31,30,31))
    if(length(isSH) > 0)
      ndaysmat[isSH,] <- ndaysmat[isSH,c(7:12,1:6)]
    
    # WI
    WI <- sdat[,c(TMInd[4:10])]-10 
    WI[WI < 0] <- 0
    WI <- WI * ndaysmat[,4:10]
    WI <- rowSums(WI)
    

    # GLADSTONES BEDD - PLACEHOLDER OFR GLADSTONES BEDD
    # BEDD <- ...
  
    # HI 
    khi <- k.hi.hall(lat = sdat[,LATInd])
    HI <- sdat[,c(TMInd[4:9],TXInd[4:9])]-10
    HI[HI < 0] <- 0
    HI <- (HI[,1:6]+HI[,7:12])/2
    HI <- HI * ndaysmat[,4:9]
    HI <- rowSums(HI)*khi
    
    #CI
    CI <- sdat[,TNInd[9]]
    summary(CI)
    
    #DI
    # Soil evaporation coefficient
    jpm5 <- sdat[,RRInd[4:9]]/5   # NOTE THAT WE DO NOT ROUND THE DAYS WITH SOIL EVAPORATION HERE
    jpm5test <- jpm5 > ndaysmat[,4:9] # Identify monthes for which jpm5 are higher than ndays in the month
    jpm5[jpm5test] <- ndaysmat[,4:9][jpm5test] # jpm5 > ndays are set to ndays (high cutoff)
    
    # Soil evaporation
    kdi <- c(0.1,0.3,0.5,0.5,0.5,0.5)
    es <- t(t(sdat[,ETInd[4:9]]/ndaysmat[,4:9]* jpm5) * (1-kdi))
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
    rr <- sdat[, RRInd[c(10:12,1:3)]]
    et <- sdat[,ETInd[c(10:12,1:3)]]
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

    # Koppen
    KOP <- koppen.fast(teta = sdat[,TMInd], rr = sdat[,RRInd])
    
    # Mean Annual Temperature
    ANT <- rowMeans(sdat[,c(TMInd)])
    ANP <- rowSums(sdat[,c(RRInd)])
    
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
    RipeP <- sdat[,c(RRInd[9])]
    
    # Fruit development conditions (July to September / SH: January to March)
    FST <- rowMeans(sdat[,c(TMInd[7:9])])
    FSP <- rowSums(sdat[,c(RRInd[7:9])])
    
    # Extreme indices
    SFR <- Rfast::rowMins(as.matrix(sdat[,c(TNInd[4:6])]), value=T)
    WCR <- Rfast::rowMins(as.matrix(sdat[,c(TNInd)]), value=T)
    HSI <- Rfast::rowMaxs(as.matrix(sdat[,c(TXInd)]), value=T)
    
    # Table of agroclimatic indices
    Agro <- cbind(sdat[,c("pixID","x","y","year")], KOPPEN=KOP, 
                                                    HI, CI, 
                                                    DI, Q, WWB,  
                                                    WI, 
                                                    ANT, ANP,
                                                    GST, GSP, GST49, GSP49, 
                                                    DST, DSP, RipeT, RipeP,
                                                    FST, FSP, DORT,
                                                    SFR, WCR, HSI)
    #Agro <- cbind(sdat[,c("pixID","x","y","pixID","year")],GCCM,  KOPPEN=KOP, GST, GSP, GST49, GSP49, DSP, DST, RipeR, RipeT, FSP, FST, DORT,  SFR, WCR, HSI)
    names(Agro) <- sub("^CI.tmin..$","CI", names(Agro))
    
      # Writing data in csv file 
    dir.create(mync_dir, showWarnings = F)
    data.table::fwrite(Agro, paste0(mync_dir,"/",mync_prefix,"TC_csv_data_Tile_",i_group,".csv"))
      
      # Fin de la boucle pour le tile i_group
    cat(paste0("\nFin de loop ",i_group," sur ",n_groups," - ",date(),"\n"))
}
Fin <- Sys.time()
Fin-Deb
