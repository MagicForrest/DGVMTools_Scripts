#!/usr/bin/R

### TODO
# 1. Loop across biomes


library(DGVMTools)
library(raster)
library(Cairo)
source("~/Projects/DGVMTools/Additional/plotUtils.v1.0.R")

t1 <- Sys.time()


##########################################################################################################
################ HERE SUPPLY VARIOUS RUN SETTINGS 

# Analysis label and plot directory
analysis.label <- "r8498"
plot.dir <- "/home/forrest/Projects/SPITFIRE/Results/FireMIP2019/DevPlots/r8498"
if(!file.exists(plot.dir)){dir.create(plot.dir)}

# save run-specific plots to run directory
savePlotsToRunDir <- FALSE

# Overlay
map.overlay <- "world"

# location of statistic metrics/summaries on plots (if selected)
metrics.xy= c(x=-140, y=-25)
summaries.xy= c(x=-130, y=-40)

# File reading options 
read.full <- FALSE
write <- TRUE
verbose <- TRUE

# resolution
resolution <- "HD"

# which plots to make
doIndividualPlots <- TRUE # plot individual run benckmarks
doMultiPlots <- TRUE # make multipanel plots with a panel for each run (groups of runs to plot together should be defined below)
doDiffPlots <- FALSE

# DGVMDDirectory
DGVMData.dir <- "/home/forrest/DGVMData/"

##########################################################################################################
################ HERE CHOOSE REGIONS TO SUM OVER

regions <- list(
  Global = list(id = "Global", name = "Global", extent = extent(-180,180,-90,90)),
  NorthernHemisphere = list(id = "NorthernHemisphere", name = "Northern Hemisphere", extent = extent(-180,180,0,90)),
  SouthernHemisphere = list(id = "SouthernHemisphere", name = "Southern Hemisphere", extent = extent(-180,180,-90,0)),
  Africa = list( id = "Africa", name = "Africa", extent =  extent(-20, 55, -30, 36)),
  Europe = list( id = "Europe", name = "Europe", extent =  extent(-30, 40, 36, 70)),
  Asia = list( id = "Asia", name = "Asia", extent =  extent(40, 180, -10, 80)),
  NorthAmerica = list( id = "NorthAmerica", name = "North America", extent =  extent(-170, -70, 25, 75)),
  SouthAmerica = list( id = "SouthAmerica", name = "South America", extent = extent(-180, -50, -60, 25)),
  Australia = list( id = "Australia", name = "Australia", extent = extent(110, 160, -45 ,10)),
  Mediterranean = list( id = "Med", name = "Mediterranean", extent = extent(10, 40, 28 ,48)),
  CentralAsia = list( id = "CentralAsia", name = "Central Asia", extent = extent(25, 140, 40, 55)),
  SouthEastAsia = list( id = "SouthEastAsia", name = "South East Asia", extent = extent(90, 140, 10, 40)),
  CentralNorthAmerica = list( id = "CentralNorthAmerica", name = "Central North America", extent = extent(-110, -85, 30, 50)),
  Boreal = list( id = "Boreal", name = "Boreal", extent = extent(-180, 180, 60, 90)),
  NHAfrica = list( id = "NHAfrica", name = "Northern Hemisphere Africa", extent = extent(-20, 50, 0, 25)),
  SHAfrica = list( id = "SHAfrica", name = "Southern Hemisphere Africa", extent = extent(5, 50, -30, 0))
)

##########################################################################################################
################ HERE CHOOSE WHICH BENCHMARKS TO DO

# biomes to plot
biomes <- list(dataset = "HandPBiomes", classification = Smith2014BiomeScheme)

# Define a list of benchmarks
benchmark.instruction.list <- list(
  
  GFED4_Annual = list("dataset" = "GFED4",
                      "name" = "GFED4",
                      "correct" = FALSE,
                      "variable" = "burntfraction_std",
                      "subannual.resolution" = "Year",
                      "subannual.aggregate.method" = "sum",
                      "spatial.aggregate.method" = "w.sum",
                      "first.year" = 1996,
                      "last.year" = 2013,
                      "layer.name" = "Total",
                      "unit" = "%",
                      "cuts" = c(0,0.002,0.005,0.01,0.02,0.05,0.10,0.2,0.50,1),
                      #"show.summary" = NULL,
                      show.metrics = c("NME", "NME_2", "r2", "m", "c"))
  
  # "Simard2011" = list("dataset" = "Simard2011",
  #                     "name" = "Simard et al. 2011",
  #                     "correct" = FALSE,
  #                     "variable" = "canopyheight_std",
  #                     "layer.name" = "CanopyHeight",
  #                     "unit" = "m"),
  # 
  # "Avitabile+Thurner_with_correction" = list("dataset" = "AvitabileThurner",
  #                                            "name" = "Avitabile + Thurner",
  #                                            "variable" = "vegC_std",
  #                                            "correct" = TRUE,
  #                                            "layer.name" = "Tree",
  #                                            "unit" =  bquote("kgC" ~ m^{"-2"})),
  # 
  # "Avitabile+Thurner" = list("dataset" = "AvitabileThurner",
  #                            "name" = "Avitabile + Thurner",
  #                            "variable" = "vegC_std",
  #                            "correct" = FALSE,
  #                            "layer.name" = "Tree",
  #                            "unit" =  bquote("kgC" ~ m^{"-2"})),
  # 
  # "MODISTreecover" = list("dataset" = "MOD44B.006",
  #                         "name" = "MODIS MOD44B",
  #                         "variable" = "vegcover_std",
  #                         "correct" = FALSE,
  #                         "layer.name" = "Tree",
  #                         "unit" =  "%"),
  # 
  # "MODISTreecover_with_correction" = list("dataset" = "MOD44B.006",
  #                                         "name" = "MODIS MOD44B",
  #                                         "variable" = "vegcover_std",
  #                                         "correct" = TRUE,
  #                                         "layer.name" = "Tree",
  #                                         "unit" =  "%"),
  # 
  # "Beer2010" = list("dataset" = "Beer2010",
  #                   "name" = "Beer et al. 2010 GPP",
  #                   "variable" = "aGPP_std",
  #                   "correct" = FALSE,
  #                   "layer.name" = "Total",
  #                   "unit" = bquote("kgC" ~ m^{"-2"} ~ y^{"-1"})) 
  
)






###################################################################
################ HERE DEFINE THE RUNS

# r8498

PNV_SPITFIRE <- defineSource(id = "PNV_SPITFIRE",
                             name = "SPITFIRE (PNV)",
                             dir = "/home/forrest/GuessRuns/FireMIP/Dev/r8498/PNV_SPITFIRE",
                             format = GUESS, 
                             forcing.data = "CRUJRA")

# duration sensitivty tests
PNV_SPITFIRE_8hr <- defineSource(id = "PNV_SPITFIRE_8hr",
                                 name = "SPITFIRE (PNV, 8hr)",
                                 dir = "/home/forrest/GuessRuns/FireMIP/Dev/r8498/PNV_SPITFIRE_8hr",
                                 format = GUESS, 
                                 forcing.data = "CRUJRA")

PNV_SPITFIRE_12hr <- defineSource(id = "PNV_SPITFIRE_12hr",
                                  name = "SPITFIRE (PNV, 12hr)",
                                  dir = "/home/forrest/GuessRuns/FireMIP/Dev/r8498/PNV_SPITFIRE_12hr",
                                  format = GUESS, 
                                  forcing.data = "CRUJRA")

PNV_SPITFIRE_daylength <- defineSource(id = "PNV_SPITFIRE_daylength",
                                       name = "SPITFIRE (PNV, daylength)",
                                       dir = "/home/forrest/GuessRuns/FireMIP/Dev/r8498/PNV_SPITFIRE_daylength",
                                       format = GUESS, 
                                       forcing.data = "CRUJRA")

# fuel moisture sensitivity tests
PNV_SPITFIRE_VPD <- defineSource(id = "PNV_SPITFIRE_VPD",
                                 name = "SPITFIRE (PNV, VPD)",
                                 dir = "/home/forrest/GuessRuns/FireMIP/Dev/r8498/PNV_SPITFIRE_VPD",
                                 format = GUESS, 
                                 forcing.data = "CRUJRA")

PNV_SPITFIRE_NoSoilMoist <- defineSource(id = "PNV_SPITFIRE_NoSoilMoist",
                                         name = "SPITFIRE (PNV, NoSoilMoist)",
                                         dir = "/home/forrest/GuessRuns/FireMIP/Dev/r8498/PNV_SPITFIRE_NoSoilMoist",
                                         format = GUESS, 
                                         forcing.data = "CRUJRA")

PNV_SPITFIRE_VPD_NoSoilMoist <- defineSource(id = "PNV_SPITFIRE_VPD_NoSoilMoist",
                                             name = "SPITFIRE (PNV, VPD, NoSoilMoist)",
                                             dir = "/home/forrest/GuessRuns/FireMIP/Dev/r8498/PNV_SPITFIRE_VPD_NoSoilMoist",
                                             format = GUESS, 
                                             forcing.data = "CRUJRA")

# ignitions tests
PNV_SPITFIRE_Lightning <- defineSource(id = "PNV_SPITFIRE_Lightning",
                                       name = "SPITFIRE (PNV, Lightning only)",
                                       dir = "/home/forrest/GuessRuns/FireMIP/Dev/r8498/PNV_SPITFIRE_Lightning",
                                       format = GUESS, 
                                       forcing.data = "CRUJRA")

PNV_SPITFIRE_Human <- defineSource(id = "PNV_SPITFIRE_Human",
                                   name = "SPITFIRE (PNV, Human only)",
                                   dir = "/home/forrest/GuessRuns/FireMIP/Dev/r8498/PNV_SPITFIRE_Human",
                                   format = GUESS, 
                                   forcing.data = "CRUJRA")
# misc
PNV_SPITFIRE_1hrSigma <- defineSource(id = "PNV_SPITFIRE_1hrSigma",
                                      name = "SPITFIRE (PNV, 1hr sigma)",
                                      dir = "/home/forrest/GuessRuns/FireMIP/Dev/r8498/PNV_SPITFIRE_1hrSigma",
                                      format = GUESS, 
                                      forcing.data = "CRUJRA")

PNV_SPITFIRE_SapSizeDecrease <- defineSource(id = "PNV_SPITFIRE_SapSizeDecrease",
                                             name = "SPITFIRE (PNV, sap size decreased)",
                                             dir = "/home/forrest/GuessRuns/FireMIP/Dev/r8498/PNV_SPITFIRE_SapSizeDecrease/",
                                             format = GUESS, 
                                             forcing.data = "CRUJRA")

PNV_SPITFIRE_OriginalFBD <- defineSource(id = "PNV_SPITFIRE_OriginalFBD",
                                         name = "SPITFIRE (PNV, Original FBD)",
                                         dir = "/home/forrest/GuessRuns/FireMIP/Dev/r8498/PNV_SPITFIRE_OriginalFBD/",
                                         format = GUESS, 
                                         forcing.data = "CRUJRA")

# wind limit (r8511)
PNV_SPITFIRE_NoWindLimit <- defineSource(id = "PNV_SPITFIRE_NoWindLimit",
                                         name = "SPITFIRE (PNV, No Wind Limit)",
                                         dir = "/home/forrest/GuessRuns/FireMIP/Dev/r8511/PNV_SPITFIRE_NoWindLimit/",
                                         format = GUESS, 
                                         forcing.data = "CRUJRA")

PNV_SPITFIRE_LasslopWindLimit <- defineSource(id = "PNV_SPITFIRE_LasslopWindLimit",
                                              name = "SPITFIRE (PNV, Lasslop Wind Limit)",
                                              dir = "/home/forrest/GuessRuns/FireMIP/Dev/r8511/PNV_SPITFIRE_LasslopWindLimit/",
                                              format = GUESS, 
                                              forcing.data = "CRUJRA")

#########################################################################################
################ HERE SELECT WHICH RUNS SHOULD BE PLOTTED INDIVIDUALLY

runs <- list(
  
  PNV_SPITFIRE,
  PNV_SPITFIRE_8hr,
  PNV_SPITFIRE_12hr,
  PNV_SPITFIRE_daylength,
  PNV_SPITFIRE_VPD,
  PNV_SPITFIRE_NoSoilMoist,
  PNV_SPITFIRE_VPD_NoSoilMoist,
  PNV_SPITFIRE_Human,
  PNV_SPITFIRE_Lightning,
  PNV_SPITFIRE_1hrSigma,
  PNV_SPITFIRE_SapSizeDecrease,
  PNV_SPITFIRE_OriginalFBD,
  PNV_SPITFIRE_NoWindLimit,
  PNV_SPITFIRE_LasslopWindLimit
)

#########################################################################################
################ HERE SELECT WHICH RUNS SHOULD BE PLOTTED TOGETHER AS A GROUP
################ (NOTE THE RUNS MUST BE INCLUDED IN THE 'runs' LIST ABOVE)

comparison.groups <- list(
  
  list(runs = list(PNV_SPITFIRE,
                   PNV_SPITFIRE_8hr,
                   PNV_SPITFIRE_12hr,
                   PNV_SPITFIRE_daylength),
       name = "Fire Duration",
       id = "FireDuration"
  ),
  list(runs = list(PNV_SPITFIRE,
                   PNV_SPITFIRE_VPD,
                   PNV_SPITFIRE_NoSoilMoist,
                   PNV_SPITFIRE_VPD_NoSoilMoist),
       name = "Fuel Moisture",
       id = "FuelMoisture"
  ),
  list(runs = list(PNV_SPITFIRE,
                   PNV_SPITFIRE_Human,
                   PNV_SPITFIRE_Lightning),
       name = "Ignitions",
       id = "Ignitions"
  ), 
  list(runs = list(PNV_SPITFIRE,
                   PNV_SPITFIRE_1hrSigma,
                   PNV_SPITFIRE_SapSizeDecrease,
                   PNV_SPITFIRE_OriginalFBD),
       name = "Misc",
       id = "Misc"
  ),
  list(runs = list(PNV_SPITFIRE,
                   PNV_SPITFIRE_NoWindLimit,
                   PNV_SPITFIRE_LasslopWindLimit),
       name = "Wind Limits",
       id = "WindLimit"
  )
)



##########################################################################################################
################ NORMALLY NOTHING TO CHANGE AFTER HERE!

# actually not currently used
BAFracToArea <- function(x, area.unit = "ha"){
  
  x <- addArea(x, unit = area.unit, verbose = FALSE)
  layers.to.calculate <- layers(x)
  for(this.layer in layers.to.calculate) x <- layerOp(x, "*", c(this.layer, "Area"), this.layer)
  x <- layerOp(x, NULL, "Area")
  return(x)
  
}

# for each benchmark
for(this.benchmark in benchmark.instruction.list) {
  
  # First read the benchmarking data
  dataset <- defineSource(id = paste(this.benchmark$dataset),
                          name = this.benchmark$name, 
                          dir = file.path(DGVMData.dir, this.benchmark$dataset, resolution), 
                          format = DGVMData)
  
  
  for(this.region in regions) {
    
    this.Data.Field <- getField(source = dataset, 
                                var = this.benchmark$variable,
                                subannual.resolution = this.benchmark$subannual.resolution,
                                subannual.aggregate.method = this.benchmark$subannual.aggregate.method,
                                first.year = this.benchmark$first.year,
                                last.year = this.benchmark$last.year,
                                spatial.extent = this.region$extent,
                                spatial.extent.id = this.region$id,
                                spatial.aggregate.method = this.benchmark$spatial.aggregate.method,
                                write = write,
                                read.full = read.full,
                                verbose = verbose)
    
    this.Data.Field@data <- setnames(this.Data.Field@data, layers(this.Data.Field), this.benchmark$variable)
    
    print(paste0("### Read benchmarking data: ", this.benchmark$name))
    
    # show.summary
    units <- character(0)
    if(length(this.benchmark$show.summary) > 0 ) {
      
      # calculate sums
      Data.sum <- aggregateSpatial(this.Data.Field, method = this.benchmark$show.summary)
      
      # sort out units
      if(this.benchmark$variable == "burntfraction_std") {
        Data.sum <- Data.sum@data / (10^4 * 10^6)
        units <- "Mha"
      }
      
      if(this.benchmark$show.summary == "w.sum") summary.string = "Sum"
      
      # make the labels data.frame
      sum.labels.vec <- c(paste(summary.string, "=", round(Data.sum,3) , units))
      field.names.vec <- c(this.Data.Field@source@name)
      
    }
    
    # data frame for stats
    metrics.df <- data.frame(stringsAsFactors = FALSE, row.names = NULL)
    
    
    #################################################################################
    ########  BENCHMARK EACH RUN
    
    # list of Comparison objects for plotting together later
    all.Comparisons <- list()
    all.Fields <- list()
    all.Fields[[this.benchmark$dataset]] <- this.Data.Field
    
    for(run in runs){
      
      if(!savePlotsToRunDir) {
        run.plot.dir <- file.path(plot.dir, "Runs", run@id)
        if(!file.exists(run.plot.dir)){dir.create(run.plot.dir, recursive = TRUE)}
      }
      else {
        run.plot.dir <- run@dir
      }
      
      
      # read the model data
      quant <- lookupQuantity(this.benchmark$variable)
      this.Model.Field <- getField(source = run, 
                                   var = this.benchmark$variable,
                                   subannual.resolution = this.benchmark$subannual.resolution,
                                   subannual.aggregate.method = this.benchmark$subannual.aggregate.method,
                                   first.year = this.benchmark$first.year,
                                   last.year = this.benchmark$last.year,
                                   spatial.extent = this.region$extent,
                                   spatial.extent.id = this.region$id,
                                   spatial.aggregate.method = this.benchmark$spatial.aggregate.method,
                                   write = write,
                                   read.full = read.full,
                                   verbose = verbose)
      
      print(paste0("* Read data: ", run@name, " vs. ", this.benchmark$name))
      
      # do the comparison
      this.comparison <- compareLayers(field1 = this.Model.Field, 
                                       field2 = this.Data.Field,
                                       layers1 = this.benchmark$variable,
                                       layers2 = this.benchmark$variable,
                                       keepall1 = FALSE,
                                       keepall2 = FALSE,
                                       tolerance = 0.1,#tolerance,
                                       verbose = FALSE,
                                       show.stats = FALSE,
                                       override.quantity = TRUE)
      print(paste0("** Done comparison: ", run@name, " vs. ", this.benchmark$name))
      
      
      # add stats to data.frame
      stats.row <- append(this.Model.Field@source@name, this.comparison@stats)
      names(stats.row) <- append("Source", names(stats.row)[2:length(names(stats.row))])
      metrics.df <- rbind(metrics.df, stats.row, stringsAsFactors=FALSE)
      
      # show.summary
      if(length(this.benchmark$show.summary)>0) {
        
        # calculate sums
        Model.sum <- aggregateSpatial(this.Model.Field, method = this.benchmark$show.summary)
        
        # sort out units
        if(this.benchmark$variable == "burntfraction_std") {    Model.sum <- Model.sum@data / (10^4 * 10^6)}
        
        # make the labels data.frame
        sum.labels.vec <- append(sum.labels.vec, paste(summary.string, "=", round(Model.sum,3) , units))
        field.names.vec <- append(field.names.vec, this.Model.Field@source@name)
        
        # for this plot
        labels.wanted <- which(field.names.vec %in% c(this.Model.Field@source@name, this.Data.Field@source@name))
        sum.labels.df <-  data.frame(label = sum.labels.vec[labels.wanted], 
                                     Facet = field.names.vec[labels.wanted], 
                                     x = summaries.xy["x"], y = summaries.xy["y"])
      }
      
      
      if(doIndividualPlots) {
        
        if(doDiffPlots) {
          # plot the difference
          this.diff.plot <- plotTemporal(list(this.Data.Field, this.Model.Field),
                                         text.multiplier = 2)
          
          # add metric stats to plot
          if(length(this.benchmark$show.metrics) > 0) {
            # select the stats we want
            selected.metrics.df <- metrics.df[which(metrics.df$Source %in% c(this.Data.Field@source@name, this.Model.Field@source@name)), append("Source", this.benchmark$show.metrics)]
            # build a data.frame of the required stats (starting from and empty data.frame)
            metrics.label.df <- data.frame(stringsAsFactors = FALSE)
            for(row.i in 1:nrow(selected.metrics.df))  metrics.label.df <- rbind(metrics.label.df,
                                                                                 list("label" = paste(names(selected.metrics.df[row.i,this.benchmark$show.metrics]), round(selected.metrics.df[row.i,this.benchmark$show.metrics],3), sep = "=", collapse = "\n"),
                                                                                      Facet = selected.metrics.df$Source[row.i],
                                                                                      x = stats.xy["x"],
                                                                                      y = stats.xy["y"]),
                                                                                 stringsAsFactors = FALSE)
            # add the stats to the plot
            this.diff.plot  <- this.diff.plot  + geom_text(data = metrics.label.df,  mapping = aes(x = x, y = y, label = label))
          }
          
          magicPlot(this.diff.plot,
                    filename = file.path(run.plot.dir , paste(this.benchmark$dataset, "TS", "Diff", this.region$id, analysis.label, paste0(this.benchmark$first.year, "-", this.benchmark$last.year),  sep = ".")),
                    height = 900,
                    width = 1800)
        }
        
        # plot the values
        this.values.plot <- plotTemporal(list(this.Data.Field, this.Model.Field), 
                                         text.multiplier = 2,
                                         nrow = 2,
                                         col.by = "Source")
        
        if(length(this.benchmark$show.summary) > 0) {   
          this.values.plot  <- this.values.plot  + geom_text(data = sum.labels.df,  mapping = aes(x = x, y = y, label = label))
        }
         
        magicPlot(this.values.plot, 
                  filename = file.path(run.plot.dir , paste(this.benchmark$dataset, "TS", "2-up", this.region$id, analysis.label, paste0(this.benchmark$first.year, "-", this.benchmark$last.year), sep = ".")),
                  height = 900,
                  width = 1800)
        
        print(paste0("*** Done plots: ", run@name, " vs. ", this.benchmark$name))
        
      }
      
      # save for group plots
      all.Comparisons[[run@id]] <- this.comparison
      all.Fields[[run@id]] <- this.Model.Field
      
    } # for each run
    
    
    #################################################################################
    ########  IF REQUESTED MAKE MULTIPANEL PLOTS - one panel per run
    
    if(doMultiPlots) {
      
      ### do each comparison group in turn
      for(group in comparison.groups) {
        
        # make a directory for the comparison
        local.group.dir <- file.path(plot.dir, "Groups", group$id)
        dir.create(local.group.dir, showWarnings = FALSE, recursive = TRUE)
        
        # extract the runs we want to plot together for this group
        local.Comparisons <- list()
        local.Fields <- list()
        local.Fields[[this.benchmark$dataset]] <- all.Fields[[this.benchmark$dataset]]
        local.names <- list()
        local.names[[this.benchmark$dataset]] <- all.Fields[[this.benchmark$dataset]]@source@name
        for(run in group$runs){
          local.Comparisons[[run@id]] <- all.Comparisons[[run@id]]
          local.Fields[[run@id]] <- all.Fields[[run@id]]
          local.names[[run@id]] <- local.Comparisons[[run@id]]@source1@name
          
        }
        
        
        # plot the difference
        if(doDiffPlots) {
          # calculate optimal number of columns, based on having one more row than columns
          opt.ncols  <- 1
          while(opt.ncols * (opt.ncols + 1) < length(local.Comparisons)) { opt.ncols <- opt.ncols +1 }
          
          # make and save
          this.diff.plot <- plotSpatialComparison(local.Comparisons, 
                                                  map.overlay = "world",
                                                  text.multiplier = 2,
                                                  ncol = opt.ncols)
          
          # add stats to plot
          if(length(this.benchmark$show.metrics) > 0) {
            # select the stats we want
            selected.metrics.df <- metrics.df[which(metrics.df$Source %in% local.names), append("Source", this.benchmark$show.metrics)]
            # build a data.frame of the required stats
            metrics.label.df <- data.frame(stringsAsFactors = FALSE)
            for(row.i in 1:nrow(selected.metrics.df))  metrics.label.df <- rbind(metrics.label.df,
                                                                                 list("label" = paste(names(selected.metrics.df[row.i,this.benchmark$show.metrics]), round(selected.metrics.df[row.i,this.benchmark$show.metrics],3), sep = "=", collapse = "\n"), 
                                                                                      Facet = paste(selected.metrics.df$Source[row.i], "-", this.Data.Field@source@name, sep = " "), 
                                                                                      x = stats.xy["x"], 
                                                                                      y = stats.xy["y"]),
                                                                                 stringsAsFactors = FALSE)
            # add the stats to the plot
            this.diff.plot  <- this.diff.plot  + geom_text(data = metrics.label.df,  mapping = aes(x = x, y = y, label = label))
          }
          
          magicPlot(this.diff.plot, 
                    filename = file.path(local.group.dir, paste(this.benchmark$dataset, "TS", "Diff", this.region$id, analysis.label, paste0(this.benchmark$first.year, "-", this.benchmark$last.year),  sep = ".")),
                    height = 900,
                    width = 1800)
        }
        
        # plot the values
        
        # calculate optimal number of columns, based on having one more row than columns
        opt.ncols  <- 1
        while(opt.ncols * (opt.ncols + 1) < length(local.Comparisons)) {
          opt.ncols <- opt.ncols +1
        }
        
        this.values.plot <- plotTemporal(local.Fields, 
                                         text.multiplier = 2,
                                         ncol = 1,
                                         col.by = "Source")
        
        if(length(this.benchmark$show.summary) > 0) {
          # extract labe
          labels.wanted <- which(field.names.vec %in% append(local.names, this.Data.Field@source@name))
          sum.labels.df <-  data.frame(label = sum.labels.vec[labels.wanted], 
                                       Facet = field.names.vec[labels.wanted], 
                                       x = summaries.xy["x"], y = summaries.xy["y"])
          this.values.plot  <- this.values.plot  + geom_text(data = sum.labels.df,  mapping = aes(x = x, y = y, label = label))
        }
        
        magicPlot(this.values.plot, 
                  height = 900,
                  width = 1800,
                  filename = file.path(local.group.dir, paste(this.benchmark$dataset, "TS", "2-up", this.region$id, analysis.label, paste0(this.benchmark$first.year, "-", this.benchmark$last.year), sep = ".")))
        
        # also write (all) the metrics to a text file in the group dir
        write.table(x = metrics.df, file = file.path(local.group.dir, paste("ComparisonMetrics", this.benchmark$dataset, "txt", sep = ".")), append = FALSE)
        
        print(paste0("### Done group: ", group$name, " vs. ", this.benchmark$name))
        
      } # for each comparison group
      
    } # if do multiplots
    
  } # for each region 
  
} # for each benchmark variable 
