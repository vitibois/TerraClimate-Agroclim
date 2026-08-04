# Descriptive statistics and climagram plots of the agroclimatic indices at
# vineyard points (VGDB), comparing the observed (2001-2020) and +2C (2042)
# scenarios, for a curated subset of wine regions.
#
# Reads the combined per-scenario point extraction produced by script 07
# (scripts/07_extract_vgdb_point_data.R) and the following reference tables,
# expected in the analysis working directory below (not produced by this
# pipeline, maintained separately by the author):
#   - ReferencesRegions_v02_Mendoza.csv: which CN_REG regions to include,
#     and their wine Type/Style.
#   - Criterions.csv: climatic thresholds ("rules") delimiting the region
#     highlighted on each plot.
#
# Originally written by Benjamin Bois.

library(terra)
library(ggplot2)
library(data.table)

source("R/config.R")

setwd(file.path(data_root, "2023_WWV_ANALYSIS"))

## ---- read data sets ----------------------------------------------------------

vgdb_csv <- file.path(data_root, "Extraction_TerraClimatPoints", paste0("VGDB_v", vgdb_version),
  paste0("AgroIndices_TerraClimate_VGDB_Pts_v", vgdb_version, "AllPeriods.csv"))
vgdb_data <- as.data.frame(fread(vgdb_csv))

dobs <- cbind(period = "2001-2020", vgdb_data[vgdb_data$scenario == "CLIMOBS20012020", ])
dfut <- cbind(period = "2042", vgdb_data[vgdb_data$scenario == "plus2", ])
# Merge obs and future data sets
same_names <- names(dobs)[names(dobs) %in% names(dfut)]
d <- rbind(dobs[, same_names], dfut[, same_names])

length(unique(dobs$CN_REG))

# Change region names
# The regional AVA "Napa Valley" is identified as "AVA Napa Valley"
d$CN_REG <- sub("US_Napa Valley", "US_AVA Napa Valley", d$CN_REG)
# Other "finer" AVAs within Napa are renamed as "Napa Valley"
napa <- c("US_Calistoga", "US_St Helena", "US_Rutherford", "US_Oakville",
          "US_Yountville", "US_Oak Knool District")
d$CN_REG[d$CN_REG %in% napa] <- "US_Napa Valley"
d$CN_REG[d$CN_REG == "NZ_Wairau Valley"] <- "NZ_Marlborough"

# Selected regions
myRegions <- read.csv("ReferencesRegions_v02_Mendoza.csv")
myRegions <- myRegions[myRegions$Include == TRUE, ]
mymatch <- match(d$CN_REG, myRegions$CN_REG)
d$Type <- myRegions$Type[mymatch]
d$Style <- myRegions$Style[mymatch]
d$Type[is.na(d$Type)] <- "Undertermined"
d$Style[is.na(d$Style)] <- "Undertermined"
summary(factor(d$Type)); summary(factor(d$Style))

myRegions$CN <- unlist(lapply(strsplit(myRegions$CN_REG, split = "_"), function(x) x[1]))
myRegions$REG <- unlist(lapply(strsplit(myRegions$CN_REG, split = "_"), function(x) x[2]))
myRegions$CN <- factor(myRegions$CN, levels = unique(myRegions$CN))
myRegions$REG <- factor(myRegions$REG, levels = unique(myRegions$REG))

# Identify Region and Country from CN_REG ("<country>_<region>")
d$CN <- substr(d$CN_REG, 1, 2)
d$REG <- substr(d$CN_REG, 4, nchar(d$CN_REG))

# We remove Koppen A (Tropical) and E (Polar) climates
d <- d[d$KOPPEN_Avg %in% 4:27, ]

# Thresholds selected
th <- read.csv("Criterions.csv")
th <- th[!is.na(th$Rule), ]
th$var <- paste0(th$Index, "_", th$Statistic.over.time)
myrules <- paste(paste(th$Index, th$Statistic.over.time, sep = "_"), th$Sign, th$Rule)

# Rainfall classes matching the DI classes
di.classes <- c(-1000, -100, 50, 150, 200)
rr49.classes <- c(0, 200, 350, 475)
di.lims <- c(-100, 50, 150)
di.lab.y <- c(-150, -25, 100, 200)
hi.lims <- c(1500, 1800, 2100, 2400, 3000, 4000)
hi.lab.y <- c(1200, 1650, 1950, 2250, 2700, 3400)
# hi.names <- c("Very\ncool","Cool","Temp.","Temp.\nwarm", "Warm", "Hot")
# di.names <- c("Very dry","Mod. dry","Sub-humid","Humid")
# myxlab <- "Huglin index [°C.days]"
# myylab <- "Dryness index [mm]"

hi.names <- c("Très\nfrais", "Frais", "Temp.", "Temp.\nchaud", "Chaud", "Très chaud")
di.names <- c("Très Sec", "Sec", "Sub-humide", "Humide")
myxlab <- "Indice de Huglin [°C.jours]"
myylab <- "Indice de sècheresse [mm]"
myzlab <- "Indice de stress thermique [°C]"

## ---- plots ---------------------------------------------------------------------

# Adds the rectangle mapping the criterion space retained for a given wine Type
makerect <- function(Type = "Still") {
  rxmax <- th[th$Type == Type & th$var == xvar & th$Sign == "<", "Rule"]
  if (length(rxmax) == 0) rxmax <- xlims[2]
  rxmin <- th[th$Type == Type & th$var == xvar & th$Sign == ">", "Rule"]
  if (length(rxmin) == 0) rxmin <- xlims[1]
  rymax <- th[th$Type == Type & th$var == yvar & th$Sign == "<", "Rule"]
  if (length(rymax) == 0) rymax <- ylims[2]
  rymin <- th[th$Type == Type & th$var == yvar & th$Sign == ">", "Rule"]
  if (length(rymin) == 0) rymin <- ylims[1]
  c(rxmin, rxmax, rymin, rymax)
}

# 1 - Selected regions - Obs and +2C facets + global criterion area ----
ds <- d[d$CN_REG %in% myRegions$CN_REG, ]
ds$CN <- factor(as.character(ds$CN), levels = myRegions$CN)
ds$REG <- factor(as.character(ds$REG), levels = myRegions$REG)
# A 10% sample of all wine-producing regions
set.seed(55); ech <- sample(1:nrow(d), round(nrow(d) * 0.1)); dsamp <- d[ech, ]
dsamp <- dsamp[, -grep("Type", names(dsamp))] # don't discriminate wine type for the whole data set

mywidth <- 6
pdf("./out_plots/Mendoza_Climagram_GeneralClimateFeatures_Obs_and_Plus2C_RegionsSubSet.pdf",
    width = mywidth, height = 6)
  mytitle <- "Diseases risks and Fruit developement conditions"
  xvar <- "DSP_Avg"; xlab <- "Disease period precipitation [mm]"
  yvar <- "FSP_Avg"; ylab <- "Fruit development period precipitation [mm]"
  xlims <- range(dsamp[, xvar], na.rm = TRUE); ylims <- range(dsamp[, yvar], na.rm = TRUE)
  rect <- makerect("Still")
  myrectdf <- cbind(rbind(makerect("Still"), makerect("Cognac"), makerect("Sparkling")),
                     data.frame(Type = c("Still", "Cognac", "Sparkling")))
  names(myrectdf) <- c("X1", "X2", "X3", "X4", "Type")

  makeplot <- function() {
    ggplot(ds, aes_string(x = xvar, y = yvar, color = "REG")) +
      annotate("rect", xmin = rect[1], xmax = rect[2], ymin = rect[3], ymax = rect[4],
               fill = "lightgreen", alpha = 0.3) +
      geom_point(data = dsamp, aes_string(x = xvar, y = yvar), color = "grey", alpha = 0.2) +
      geom_point(aes(shape = CN)) + stat_ellipse() + facet_grid(scenario ~ Type) +
      xlab(xlab) + ylab(ylab) + xlim(xlims) + ylim(ylims) +
      ggtitle(mytitle) + theme(legend.position = "bottom") + theme_bw() +
      scale_color_discrete(name = "Region") + scale_shape_discrete(name = "Country")
  }
  print(makeplot())

  mytitle <- "Cold damage related risks"
  xvar <- "WCR_Avg"; xlab <- "Winter cold risk index [°C]"
  yvar <- "SFR_Avg"; ylab <- "Spring frost risk index [°C]"
  xlims <- range(dsamp[, xvar], na.rm = TRUE); ylims <- range(dsamp[, yvar], na.rm = TRUE)
  rect <- makerect()
  print(makeplot())

  mytitle <- "Earliness and Water availability"
  xvar <- "HI_Avg"; xlab <- "Huglin Index [°C.days]"
  yvar <- "DI_Avg"; ylab <- "Dryness Index [mm]"
  xlims <- range(dsamp[, xvar], na.rm = TRUE); ylims <- range(dsamp[, yvar], na.rm = TRUE)
  rect <- makerect()
  print(makeplot())

  mytitle <- "Earliness and Water availability"
  xvar <- "GST49_Avg"; xlab <- "Growing season temperature [°C]"
  yvar <- "GSP49_Avg"; ylab <- "Growing season precipitation [mm]"
  xlims <- range(dsamp[, xvar], na.rm = TRUE); ylims <- range(dsamp[, yvar], na.rm = TRUE)
  rect <- makerect()
  print(makeplot())

  mytitle <- "Ripening period thermal conditions"
  xvar <- "HSI_Avg"; xlab <- "Heat Stress Index [°C]"
  yvar <- "CI_Avg"; ylab <- "Cool night index [°C]"
  xlims <- range(dsamp[, xvar], na.rm = TRUE); ylims <- range(dsamp[, yvar], na.rm = TRUE)
  rect <- makerect()
  print(makeplot())
dev.off()

# 2 - Selected regions - Obs only - global criterion area ----
ds <- d[d$CN_REG %in% myRegions$CN_REG, ]
ds <- ds[ds$scenario == "OBS", ]
set.seed(55); ech <- sample(1:nrow(d), round(nrow(d) * 0.1)); dsamp <- d[ech, ]
dsamp <- dsamp[dsamp$scenario == "OBS", -grep("Type", names(dsamp))]

pdf("Climagram_GeneralClimateFeatures_ObsOnly_RegionsSubSet.pdf", width = 8, height = 5)
  mytitle <- "Diseases risks and Fruit developement conditions"
  xvar <- "DSP_Avg"; xlab <- "Disease period precipitation [mm]"
  yvar <- "FSP_Avg"; ylab <- "Fruit development period precipitation [mm]"
  xlims <- range(dsamp[, xvar], na.rm = TRUE); ylims <- range(dsamp[, yvar], na.rm = TRUE)
  rect <- makerect()

  makeplot <- function() {
    ggplot(ds, aes_string(x = xvar, y = yvar, color = "CN_REG")) +
      annotate("rect", xmin = rect[1], xmax = rect[2], ymin = rect[3], ymax = rect[4],
               fill = "lightgreen", alpha = 0.3) +
      geom_point(data = dsamp, aes_string(x = xvar, y = yvar), color = "grey", alpha = 0.2) +
      geom_point(aes(shape = Type), size = 3, alpha = 0.5) + stat_ellipse() +
      xlab(xlab) + ylab(ylab) + xlim(xlims) + ylim(ylims) +
      ggtitle(mytitle) + theme() + theme_bw() +
      scale_colour_discrete(name = "Region") +
      scale_shape_manual(name = "Wine type", values = c(4, 1, 18))
  }
  plot_allwinereg <- function() {
    ggplot(dsamp, aes_string(x = xvar, y = yvar)) +
      geom_point(data = dsamp, aes_string(x = xvar, y = yvar), color = "forestgreen", alpha = 0.2) +
      xlab(xlab) + ylab(ylab) + xlim(xlims) + ylim(ylims) +
      ggtitle(mytitle) + theme(legend.position = "none") + theme_bw()
  }
  print(plot_allwinereg())
  print(makeplot())

  mytitle <- "Cold damage related risks"
  xvar <- "WCR_Avg"; xlab <- "Winter cold risk index [°C]"
  yvar <- "SFR_Avg"; ylab <- "Spring frost risk index [°C]"
  xlims <- range(dsamp[, xvar], na.rm = TRUE); ylims <- range(dsamp[, yvar], na.rm = TRUE)
  rect <- makerect()
  print(plot_allwinereg())
  print(makeplot())

  mytitle <- "Earliness and Water availability"
  xvar <- "HI_Avg"; xlab <- "Huglin Index [°C.days]"
  yvar <- "DI_Avg"; ylab <- "Dryness Index [mm]"
  xlims <- range(dsamp[, xvar], na.rm = TRUE); ylims <- range(dsamp[, yvar], na.rm = TRUE)
  rect <- makerect()
  print(plot_allwinereg())
  print(makeplot())

  mytitle <- "Earliness and Water availability"
  xvar <- "GST49_Avg"; xlab <- "Growing season temperature [°C]"
  yvar <- "GSP49_Avg"; ylab <- "Growing season precipitation [mm]"
  xlims <- range(dsamp[, xvar], na.rm = TRUE); ylims <- range(dsamp[, yvar], na.rm = TRUE)
  rect <- makerect()
  print(plot_allwinereg())
  print(makeplot())

  mytitle <- "Ripening period thermal conditions"
  xvar <- "HSI_Avg"; xlab <- "Heat Stress Index [°C]"
  yvar <- "CI_Avg"; ylab <- "Cool night index [°C]"
  xlims <- range(dsamp[, xvar], na.rm = TRUE); ylims <- range(dsamp[, yvar], na.rm = TRUE)
  rect <- makerect()
  print(plot_allwinereg())
  print(makeplot())
dev.off()

# 3 - Selected regions - Obs and +2C - global criterion area ----
ds <- d[d$CN_REG %in% myRegions$CN_REG, ]
ds$fac_scenario <- factor(ds$scenario, labels = c("2001-2020", "+2°C (2042)"))
set.seed(55); ech <- sample(1:nrow(d), round(nrow(d) * 0.1)); dsamp <- d[ech, ]
dsamp <- dsamp[dsamp$scenario == "OBS", -grep("Type", names(dsamp))]

pdf("Climagram_GeneralClimateFeatures_ObsPlus2C_RegionsSubSet.pdf", width = 10, height = 5)
  mytitle <- "Diseases risks and Fruit developement conditions"
  xvar <- "DSP_Avg"; xlab <- "Disease period precipitation [mm]"
  yvar <- "FSP_Avg"; ylab <- "Fruit development period precipitation [mm]"
  xlims <- range(dsamp[, xvar], na.rm = TRUE); ylims <- range(dsamp[, yvar], na.rm = TRUE)
  rect <- makerect()

  makeplot <- function() {
    ggplot(ds, aes_string(x = xvar, y = yvar, color = "CN_REG")) +
      annotate("rect", xmin = rect[1], xmax = rect[2], ymin = rect[3], ymax = rect[4],
               fill = "lightgreen", alpha = 0.3) +
      geom_point(data = dsamp, aes_string(x = xvar, y = yvar), color = "grey", alpha = 0.2) +
      geom_point(aes(shape = Type), size = 3, alpha = 0.5) + stat_ellipse() +
      xlab(xlab) + ylab(ylab) + xlim(xlims) + ylim(ylims) +
      ggtitle(mytitle) + theme() + theme_bw() +
      scale_colour_discrete(name = "Region") +
      scale_shape_manual(name = "Wine type", values = c(4, 1, 18)) +
      facet_wrap(fac_scenario ~ ., ncol = 2)
  }
  print(makeplot())

  mytitle <- "Cold damage related risks"
  xvar <- "WCR_Avg"; xlab <- "Winter cold risk index [°C]"
  yvar <- "SFR_Avg"; ylab <- "Spring frost risk index [°C]"
  xlims <- range(dsamp[, xvar], na.rm = TRUE); ylims <- range(dsamp[, yvar], na.rm = TRUE)
  rect <- makerect()
  print(makeplot())

  mytitle <- "Earliness and Water availability"
  xvar <- "HI_Avg"; xlab <- "Huglin Index [°C.days]"
  yvar <- "DI_Avg"; ylab <- "Dryness Index [mm]"
  xlims <- range(dsamp[, xvar], na.rm = TRUE); ylims <- range(dsamp[, yvar], na.rm = TRUE)
  rect <- makerect()
  print(makeplot())

  mytitle <- "Earliness and Water availability"
  xvar <- "GST49_Avg"; xlab <- "Growing season temperature [°C]"
  yvar <- "GSP49_Avg"; ylab <- "Growing season precipitation [mm]"
  xlims <- range(dsamp[, xvar], na.rm = TRUE); ylims <- range(dsamp[, yvar], na.rm = TRUE)
  rect <- makerect()
  print(makeplot())

  mytitle <- "Ripening period thermal conditions"
  xvar <- "HSI_Avg"; xlab <- "Heat Stress Index [°C]"
  yvar <- "CI_Avg"; ylab <- "Cool night index [°C]"
  xlims <- range(dsamp[, xvar], na.rm = TRUE); ylims <- range(dsamp[, yvar], na.rm = TRUE)
  rect <- makerect()
  print(makeplot())
dev.off()

# 4 - Huglin x Dryness climagram, Obs and +2C ----
ds <- d[d$CN_REG %in% myRegions$CN_REG, ]
ds$fac_scenario <- factor(ds$scenario, labels = c("2001-2020", "+2°C (2042)"))
ds$Type <- factor(ds$Type, labels = c("Cognac", "Champagne", "Vins tranquilles"))
dsamp <- d
dsamp <- dsamp[, -grep("Type", names(dsamp))]

pdf("Climagram_IndHuglin-Secheresse.pdf", width = 8.5, height = 7)
  mytitle <- "2001-2020 : Précocité et disponibilité en eau "
  xvar <- "HI_Avg"; xlab <- "Indice de Huglin [°C.jour]"
  yvar <- "DI_Avg"; ylab <- "Indice de Sècheresse [mm]"
  xlims <- range(dsamp[, xvar], na.rm = TRUE); ylims <- range(dsamp[, yvar], na.rm = TRUE)

  print(ggplot(dsamp[dsamp$scenario == "OBS", ], aes_string(x = xvar, y = yvar)) +
    geom_hline(yintercept = di.lims) + geom_vline(xintercept = hi.lims) +
    geom_point(color = "forestgreen", alpha = 0.2) +
    xlab(xlab) + ylab(ylab) + xlim(xlims) + ylim(ylims) +
    ggtitle(mytitle) + theme(legend.position = "none") + theme_bw() +
    annotate("text", label = di.names, x = rep(xlims[2] - 100, length(di.names)),
             y = di.lab.y, angle = 90, hjust = .5, fontface = 2, size = 4, vjust = 0) +
    annotate("text", label = hi.names, x = hi.lab.y, y = rep(ylims[1], length(hi.names)),
             angle = 0, fontface = 2, size = 4, vjust = 0))

  rect <- makerect()
  print(ggplot(ds[ds$scenario == "OBS", ], aes_string(x = xvar, y = yvar, color = "CN_REG")) +
    geom_hline(yintercept = di.lims) + geom_vline(xintercept = hi.lims) +
    annotate("rect", xmin = rect[1], xmax = rect[2], ymin = rect[3], ymax = rect[4],
             fill = "lightgreen", alpha = 0.3) +
    geom_point(data = dsamp[dsamp$scenario == "OBS", ], aes_string(x = xvar, y = yvar),
               color = "grey", alpha = 0.2) +
    geom_point(aes(shape = Type), size = 3, alpha = 0.5) + stat_ellipse() +
    xlab(xlab) + ylab(ylab) + xlim(xlims) + ylim(ylims) +
    ggtitle(mytitle) + theme() + theme_bw() +
    annotate("text", label = di.names, x = rep(xlims[2] - 100, length(di.names)),
             y = di.lab.y, angle = 90, hjust = .5, fontface = 2, size = 4, vjust = 0) +
    annotate("text", label = hi.names, x = hi.lab.y, y = rep(ylims[1], length(hi.names)),
             angle = 0, fontface = 2, size = 4, vjust = 0) +
    scale_colour_discrete(name = "Region") +
    scale_shape_manual(name = "Type de vin", values = c(4, 1, 18)))

  mytitle <- "+2°C : Précocité et disponibilité en eau "
  xlims <- range(dsamp[, xvar], na.rm = TRUE); ylims <- range(dsamp[, yvar], na.rm = TRUE)
  print(ggplot(ds[ds$scenario == "plus2", ], aes_string(x = xvar, y = yvar, color = "CN_REG")) +
    geom_hline(yintercept = di.lims) + geom_vline(xintercept = hi.lims) +
    annotate("rect", xmin = rect[1], xmax = rect[2], ymin = rect[3], ymax = rect[4],
             fill = "lightgreen", alpha = 0.3) +
    geom_point(data = dsamp[dsamp$scenario == "plus2", ], aes_string(x = xvar, y = yvar),
               color = "grey", alpha = 0.2) +
    geom_point(aes(shape = Type), size = 3, alpha = 0.5) + stat_ellipse() +
    xlab(xlab) + ylab(ylab) + xlim(xlims) + ylim(ylims) +
    ggtitle(mytitle) + theme() + theme_bw() +
    annotate("text", label = di.names, x = rep(xlims[2] - 100, length(di.names)),
             y = di.lab.y, angle = 90, hjust = .5, fontface = 2, size = 4, vjust = 0) +
    annotate("text", label = hi.names, x = hi.lab.y, y = rep(ylims[1], length(hi.names)),
             angle = 0, fontface = 2, size = 4, vjust = 0) +
    scale_colour_discrete(name = "Region") +
    scale_shape_manual(name = "Wine type", values = c(4, 1, 18)))
dev.off()
