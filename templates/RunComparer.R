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

# Years
first.year = 1976
last.year = 2000

# Re-read 
reread <- FALSE


##########################################################################################################
################ HERE SELECT WHICH TYPES OF PLOTS TO MAKE

# group comparisons and difference plots
doMultiPlots <- TRUE # make multipanel plots for the groups defined above (with a panel for each run)
doDifferencePlots <- TRUE # do difference plots for run pairs (pairs defined above)

# only for plots with multiple layers
doLayerPlots <- TRUE # plot individual layers (usuually) PFTs
doAggregates <- TRUE # plot some aggregates 

# only for monthly variables
doMonthly <- FALSE 
doSeasonal <- TRUE 

# if biome classification - maybe move this to another script
doBiomes <- TRUE 



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
                                             name = "SPITFIRE (PNV, spa size decreased)",
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
################ HERE SELECT WHICH RUNS SHOULD BE PLOTTED (INDIVIDUALLY)

runs <- list(
  
  PNV_SPITFIRE,
  # PNV_SPITFIRE_8hr,
  # PNV_SPITFIRE_12hr,
  # PNV_SPITFIRE_daylength,
  # PNV_SPITFIRE_VPD,
  # PNV_SPITFIRE_NoSoilMoist,
  # PNV_SPITFIRE_VPD_NoSoilMoist,
  # PNV_SPITFIRE_Human,
  # PNV_SPITFIRE_Lightning,
  # PNV_SPITFIRE_1hrSigma,
  # PNV_SPITFIRE_SapSizeDecrease,
  # PNV_SPITFIRE_OriginalFBD,
  PNV_SPITFIRE_NoWindLimit,
  PNV_SPITFIRE_LasslopWindLimit
  
)

#########################################################################################
################ HERE SELECT WHICH RUNS SHOULD BE PLOTTED TOGETHER AS A GROUP

comparison.groups <- list(
  
  # list(runs = list(PNV_SPITFIRE,
  #                  PNV_SPITFIRE_8hr,
  #                  PNV_SPITFIRE_12hr,
  #                  PNV_SPITFIRE_daylength),
  #      name = "Fire Duration",
  #      id = "FireDuration"
  # ),
  # list(runs = list(PNV_SPITFIRE,
  #                  PNV_SPITFIRE_VPD,
  #                  PNV_SPITFIRE_NoSoilMoist,
  #                  PNV_SPITFIRE_VPD_NoSoilMoist),
  #      name = "Fuel Moisture",
  #      id = "FuelMoisture"
  # ),
  # list(runs = list(PNV_SPITFIRE,
  #                  PNV_SPITFIRE_Human,
  #                  PNV_SPITFIRE_Lightning),
  #      name = "Ignitions",
  #      id = "Ignitions"
  # ),
  # list(runs = list(PNV_SPITFIRE,
  #                  PNV_SPITFIRE_1hrSigma,
  #                  PNV_SPITFIRE_SapSizeDecrease,
  #                  PNV_SPITFIRE_OriginalFBD),
  #      name = "Misc",
  #      id = "Misc"
  # ),
  list(runs = list(PNV_SPITFIRE,
                   PNV_SPITFIRE_NoWindLimit,
                   PNV_SPITFIRE_LasslopWindLimit),
       name = "Wind Limits",
       id = "WindLimit"
  )
  
)


##########################################################################################################
################ HERE SELECT PAIRS OF RUNS WHICH SHOULD BE COMPARED DIRECTLY AGAINST EACH OTHER  

difference.pairs <- list(
  # list("base" =  PNV_SPITFIRE, "new" = PNV_SPITFIRE_8hr),
  # list("base" =  PNV_SPITFIRE, "new" = PNV_SPITFIRE_12hr),
  # list("base" =  PNV_SPITFIRE, "new" = PNV_SPITFIRE_daylength),
  # list("base" =  PNV_SPITFIRE, "new" = PNV_SPITFIRE_VPD),
  # list("base" =  PNV_SPITFIRE, "new" = PNV_SPITFIRE_NoSoilMoist),
  # list("base" =  PNV_SPITFIRE, "new" = PNV_SPITFIRE_VPD_NoSoilMoist),
  # list("base" =  PNV_SPITFIRE, "new" = PNV_SPITFIRE_Lightning),
  # list("base" =  PNV_SPITFIRE, "new" = PNV_SPITFIRE_Human),
  # list("base" =  PNV_SPITFIRE, "new" = PNV_SPITFIRE_1hrSigma),
  # list("base" =  PNV_SPITFIRE, "new" = PNV_SPITFIRE_SapSizeDecrease),
  # list("base" =  PNV_SPITFIRE, "new" = PNV_SPITFIRE_OriginalFBD),
  list("base" =  PNV_SPITFIRE, "new" = PNV_SPITFIRE_NoWindLimit),
  list("base" =  PNV_SPITFIRE, "new" = PNV_SPITFIRE_LasslopWindLimit)
)



##########################################################################################################
################ HERE CHOOSE WHICH VARIABLES AND PFT AGGREGATES TO PLOT

# plot some variables for each run
vars.to.plot <- list(
  
  #lai = list(var = "lai", cuts = seq(0, 10, 0.5), PFT = TRUE)#,
  #fpc = list(var = "fpc", cuts = seq(0, 1.3, 0.1), PFT = TRUE),
  mfirefrac = list(var = "mfirefrac", cuts = c(0,0.002,0.005,0.01,0.02,0.05,0.10,0.2,0.50,1), PFT = FALSE)
  #tot_runof = list(var = "tot_runoff", cuts = seq(0,1000,50), PFT = FALSE)
  
)


# PFT aggregates to plot - these are each a list of arguments to layerOp (which are used by a do.call)
PFT.aggs <- list(
  
  # Maximums
  MaxPFT = list(operator = "max.layer", layers = c(".PFT"), new.layer = "MaxPFT"),
  MaxTree = list(operator = "max.layer", layers = c(".Tree"), new.layer = "MaxTree"),
  MaxWoody = list(operator = "max.layer", layers = c(".Tree", ".Shrub"), new.layer = "MaxWoody"),
  
  # Totals
  Tree = list(operator = "+", layers = c(".Tree"), new.layer = "Tree"),
  Grass = list(operator = "+", layers = c(".Grass"), new.layer = "Grass"),
  Woody = list(operator = "+", layers = c(".Tree", ".Shrub"), new.layer = "Woody")
  
  # Could also add ".Boreal" or ".Broadleaved" (for example)
  # ...
  
)

##########################################################################################################
################ HERE CHOOSE WHICH BIOMES TO PLOT

# biomes to plot
biomes <- list(dataset = "HandPBiomes", classification = Smith2014BiomeScheme)

# resolution
resolution <- "HD"


##########################################################################################################
################ NORMALLY NOTHING TO CHANGE AFTER HERE!




#biome.scheme
if(doBiomes) {
  
  # get the biome dataset
  biome.scheme <- biomes$classification
  
  biome.src <- defineSource(id = paste(biomes$dataset, resolution, sep = "."),
                            name = "PNV Vegetation Types", 
                            dir = file.path("/home/forrest/DGVMData/HandP_PNV/HD"),
                            format = DGVMData) 
  
  biome.data.field <- getField(source = biome.src, 
                               var = biomes$classification,
                               verbose = FALSE)
  
}


for(var.details in vars.to.plot) {
  
  
  var.str <- var.details$var
  var.cuts <- var.details$cuts
  var.isPFT <- var.details$PFT
  print(var.str)
  
  field.list <- list()
  field.seasonal.list <- list()
  field.annual.list <- list()
  
  #################################################################################
  ########  FOR EACH RUN MAKE SOME STANDARD PLOTS OF THIS VARIABLE
  
  for(run in runs){
    
    if(!savePlotsToRunDir) {
      run.plot.dir <- file.path(plot.dir, "Runs", run@id)
      if(!file.exists(run.plot.dir)){dir.create(run.plot.dir, recursive = TRUE)}
    }
    else {
      run.plot.dir <- run@dir
    }
    
    
    print(run.plot.dir)
    print(run@id)
    
    quant <- lookupQuantity(var.str, run@format)
    
    this.model.field <- getField(source = run, 
                                 var = quant,
                                 first.year = first.year,
                                 last.year = last.year,
                                 year.aggregate.method = "mean",
                                 read.full = reread,
                                 write = TRUE,
                                 verbose = TRUE)
    
    # determine if it is monthly
    isMonthly <- "Month" %in% getDimInfo(this.model.field)
    
    # calculate optimal number of columns, based on having one more row than columns
    opt.ncols  <- 1
    if(isMonthly) npanels <- length(layers(this.model.field)) * getDimInfo(this.model.field, "size")$Month
    else npanels <- length(layers(this.model.field)) 
    
    while(opt.ncols * (opt.ncols + 1) < npanels) {
      opt.ncols <- opt.ncols +1
    }
    
    ### always make a summary plot of all layers 
    summary.plot <- plotSpatial(this.model.field,
                                cuts = var.cuts,
                                map.overlay = map.overlay,
                                text.multiplier = 2, 
                                ncol = opt.ncols)
    magicPlot(summary.plot, 
              filename = file.path(run.plot.dir, paste(var.str, "Summary", paste0(first.year, "-", last.year), analysis.label, sep = ".")))
    
    ### If variable is do sums monthly, also 
    if(isMonthly) {
      
      # sum to annual and plot that
      this.model.field.annual <- aggregateSubannual(this.model.field, "sum")
      annual.plot <- plotSpatial(this.model.field.annual,
                                 cuts = var.cuts,
                                 map.overlay = map.overlay,
                                 text.multiplier = 2)
      magicPlot(annual.plot, 
                filename = file.path(run.plot.dir, paste(var.str, "AnnualSum", paste0(first.year, "-", last.year), analysis.label, sep = ".")))
      
      if(doMonthly) {
        
        # plot months individually
        for(month in all.months) {
          month.plot <- plotSpatial(this.model.field, 
                                    month = month@index,
                                    cuts = var.cuts,
                                    map.overlay = map.overlay,
                                    text.multiplier = 2)
          magicPlot(month.plot, 
                    filename = file.path(run.plot.dir, paste(var.str, month@id, paste0(first.year, "-", last.year), analysis.label, sep = ".")))
        }
        
      }
      
      if(doSeasonal) {
        
        # calculate seasonal averages
        this.model.field.seasonal <- aggregateSubannual(this.model.field, "sum", target = "Season")
        
        # plot all seasons together
        season.plot <- plotSpatial(this.model.field.seasonal,
                                   cuts = var.cuts,
                                   map.overlay = map.overlay,
                                   text.multiplier = 2)
        magicPlot(season.plot, 
                  filename = file.path(run.plot.dir, paste(var.str, "Seasons", "Sum", paste0(first.year, "-", last.year), analysis.label, sep = ".")))
        
        # plot seasons individually
        for(season in all.seasons) {
          season.plot <- plotSpatial(this.model.field.seasonal, 
                                     seasons = season@id,
                                     cuts = var.cuts,
                                     map.overlay = map.overlay,
                                     text.multiplier = 2)
          magicPlot(season.plot, 
                    filename = file.path(run.plot.dir, paste(var.str, season@id, "Sum", paste0(first.year, "-", last.year), analysis.label, sep = ".")))
        }
        
      }
      
    }
    
    
    ### individual plots
    if(doLayerPlots) {
      
      for(this.layer in names(this.model.field)) {
        individual.plot <- plotSpatial(this.model.field,
                                       layers = this.layer,
                                       map.overlay = map.overlay,
                                       text.multiplier = 2,
                                       cuts = var.cuts)
        magicPlot(individual.plot, 
                  filename = file.path(run.plot.dir, paste(var.str, this.layer, paste0(first.year, "-", last.year), analysis.label, sep = ".")))
      }
      message("*** Done initial layer plots")  
    }
    
    ### aggregate plots  
    if(doAggregates && var.isPFT) {
      
      for(this.aggregate in PFT.aggs){
        
        
        this.arguments <- list(x = this.model.field)
        this.model.field <- do.call(layerOp, append(this.arguments, this.aggregate))
        
        # don't do cuts is max.laeyers
        local.cuts <- var.cuts
        if(this.aggregate$operator == "max.layer") local.cuts <- waiver() 
        
        if(this.aggregate$new.layer %in% names(this.model.field)) {
          aggregate.plot <- plotSpatial(this.model.field,
                                        layers = this.aggregate$new.layer,
                                        map.overlay = map.overlay, 
                                        text.multiplier = 2,
                                        cuts = local.cuts)
          magicPlot(aggregate.plot, 
                    filename = file.path(run.plot.dir, paste(var.str, this.aggregate$new.layer, paste0(first.year, "-", last.year), analysis.label, sep = ".")))
        }
      }
      message("*** Done aggregate layer plots")  
      
    } # for each aggregate
    
    # save for plotting all together later and comparison plots
    field.list[[run@id]] <- this.model.field
    if(isMonthly && doSeasonal) field.seasonal.list[[run@id]] <- this.model.field.seasonal
    if(isMonthly) field.annual.list[[run@id]] <- this.model.field.annual
    
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
      local.field.list <- list()
      local.field.annual.list <- list()
      local.field.seasonal.list <- list()
      for(run in group$runs){
        local.field.list[[run@id]] <- field.list[[run@id]]
        local.field.annual.list[[run@id]] <- field.annual.list[[run@id]]
        local.field.seasonal.list[[run@id]] <- field.seasonal.list[[run@id]]
      }
      
      
      ### individual layers plots
      if(doLayerPlots) {
        
        all.layers <- c()
        for(run in group$runs){
          all.layers <- append(all.layers, names(local.field.list[[run@id]]))
        }
        all.layers <- unique(all.layers)
        
        for(this.layer in all.layers) {
          group.layer.plot <- plotSpatial(local.field.list,
                                          layers = this.layer,
                                          map.overlay = map.overlay,
                                          text.multiplier = 2,
                                          cuts = var.cuts)
          if(length(local.field.list) > 1) group.layer.plot <- group.layer.plot + facet_wrap(~Facet, nrow = 2)
          magicPlot(group.layer.plot, 
                    filename = file.path(local.group.dir, paste(var.str, group$id, this.layer, paste0(first.year, "-", last.year), analysis.label, sep = ".")))
        } # for each layer
        
      } # if doIndividual plots
      
      ### If variable is do sums monthly, also 
      if(isMonthly) {
        
        # plot the annual sums
        annual.plot <- plotSpatial(local.field.annual.list,
                                   cuts = var.cuts,
                                   map.overlay = map.overlay,
                                   text.multiplier = 2)
        magicPlot(annual.plot, 
                  filename = file.path(local.group.dir, paste(var.str, "AnnualSum", paste0(first.year, "-", last.year), analysis.label, sep = ".")))
        
        if(doMonthly) {
          
          # plot months individually
          for(month in all.months) {
            month.plot <- plotSpatial(local.field.list, 
                                      month = month@index,
                                      cuts = var.cuts,
                                      map.overlay = map.overlay,
                                      text.multiplier = 2)
            magicPlot(month.plot, 
                      filename = file.path(local.group.dir, paste(var.str, month@id, paste0(first.year, "-", last.year), analysis.label, sep = ".")))
          }
          
        }
        
        if(doSeasonal) {
          
          # plot all seasons together
          season.plot <- plotSpatial(local.field.seasonal.list, 
                                     cuts = var.cuts,
                                     map.overlay = map.overlay,
                                     text.multiplier = 2)
          magicPlot(season.plot, 
                    filename = file.path(local.group.dir, paste(var.str, "Seasons", "Sum", paste0(first.year, "-", last.year), analysis.label, sep = ".")))
          
          # plot seasons individually
          for(season in all.seasons) {
            season.plot <- plotSpatial(local.field.seasonal.list, 
                                       seasons = season@id,
                                       cuts = var.cuts,
                                       map.overlay = map.overlay,
                                       text.multiplier = 2)
            print(season.plot)
            stop()
            magicPlot(season.plot, 
                      filename = file.path(local.group.dir, paste(var.str, season@id, "Sum", paste0(first.year, "-", last.year), analysis.label, sep = ".")))
          
          
            
          }
          
        }
        
      }
      
      
      ### layer aggregates
      if(doAggregates && var.isPFT) {
        
        for(this.aggregate in PFT.aggs){
          
          # don't cut for maxes 
          local.cuts <- var.cuts
          if(this.aggregate$operator == "max.layer") local.cuts <- waiver() 
          
          group.aggregate.plot <- plotSpatial(local.field.list,
                                              layers = this.aggregate$new.layer,
                                              map.overlay = map.overlay, 
                                              text.multiplier = 2,
                                              cuts = local.cuts)
          if(length(local.field.list) > 1) group.aggregate.plot <- group.aggregate.plot + facet_wrap(~Facet, nrow = 2)
          magicPlot(group.aggregate.plot, 
                    filename = file.path(local.group.dir, paste(var.str, group$id, this.aggregate$new.layer, paste0(first.year, "-", last.year), analysis.label, sep = ".")))
        } # for each aggregate
        
      } # if doAggregates
      
      print(paste0("*** Done group: ", group$name))
      
    } # for each comparison group 
    
  } # if do multiplots
  
  
  ######################################################
  #### NOW DIFFERENCE PLOTS
  
  if(doDifferencePlots) {
    
    for(this.pair in difference.pairs) {
      
      base <- this.pair$base
      new <- this.pair$new
      base.field <- field.list[[base@id]]
      new.field <- field.list[[new@id]]
      
      # make a directory for the comparison
      local.comparison.dir <- file.path(plot.dir, "Comparisons", paste(new@id, base@id, sep ="-"))
      dir.create(local.comparison.dir, showWarnings = FALSE, recursive = TRUE)
      
      # make zero-layers for missing layers
      for(this.layer in names(base.field)) {
        if(!this.layer %in% names(new.field)) {
          layerOp(new.field, 0, layers = this.layer)
        }
      }
      
      for(this.layer in names(new.field)) {
        if(!this.layer %in% names(base.field)) {
          layerOp(base.field, 0, layers = this.layer)
        }
      }
      
      for(layer in names(new.field)){
        
        # make the comparison layer
        this.comparison <- compareLayers(new.field, base.field, layers1 = layer, show.stats = FALSE)
        
        
        # plot the difference
        diff.layer.plot<- plotSpatialComparison(this.comparison, 
                                                map.overlay = "world",
                                                text.multiplier = 2)
        magicPlot(diff.layer.plot, 
                  filename = file.path(local.comparison.dir, paste(layer, quant@id, "Diff", analysis.label, paste0(first.year, "-", last.year), sep = ".")))
        
        # plot the values
        abs.layer.plot <- plotSpatialComparison(this.comparison, 
                                                type = "values",
                                                map.overlay = "world",
                                                text.multiplier = 2)
        if(length(local.field.list) > 1) abs.layer.plot <- abs.layer.plot + facet_wrap(~Facet, nrow = 2)
        magicPlot(abs.layer.plot, 
                  filename = file.path(local.comparison.dir, paste(layer, quant@id, "2-up", analysis.label, paste0(first.year, "-", last.year), sep = ".")))
        
      } # for each layer 
      
      print(paste0("*** Done difference: ", paste(new@id, base@id, sep ="-")))
      
    } # end for each difference pair
    
  } # end if doDifference plots
  
} # for each variable 





#######################################################
#### NOW BIOMES

if(doBiomes) {
  
  biome.objects.list <- list()
  
  
  # calculate biomes for each run and save the plots
  for(run in runs){
    
    if(!savePlotsToRunDir) {
      run.plot.dir <- file.path(plot.dir, "Runs", run@id)
      if(!file.exists(run.plot.dir)){dir.create(run.plot.dir, recursive = TRUE)}
    }
    else {
      run.plot.dir <- run@dir
    }
    
    this.model.field <- getField(source = run, 
                                 var = "lai",
                                 first.year = first.year,
                                 last.year = last.year,
                                 year.aggregate.method = "mean",
                                 read.full = reread,
                                 write = TRUE,
                                 verbose = TRUE)
    
    # combine shade intolerant PFT layers
    layerOp(this.model.field, "+", c("BNE", "BINE"), "BNE")
    layerOp(this.model.field, 0, "BINE")
    
    if("TeIBS" %in% names(this.model.field)) {
      layerOp(this.model.field, "+", c("TeBS", "TeIBS"), "TeBS")
      layerOp(this.model.field, 0, "TeIBS")
    }
    
    if("IBS" %in% names(this.model.field)) {
      layerOp(this.model.field, "+", c("TeBS", "IBS"), "TeBS")
      layerOp(this.model.field, 0, "IBS")
    }
    
    if("BIBS" %in% names(this.model.field)) {
      layerOp(this.model.field, "+", c("BIBS", "BNE"), "BNE")
      layerOp(this.model.field, 0, "BIBS")
    }
    
    # calculate the biomes from the model output
    model.biomes <- getScheme(source = run, scheme = biome.scheme, first.year = 1961, last.year = 1990, year.aggregate.method = "mean")
    biome.objects.list[[run@id]] <- model.biomes
    
    
    # comapare biomes
    comparison.layer <- compareLayers(model.biomes, biome.data.field, layers1 = biome.scheme@id, tolerance = 0.1, show.stats = FALSE)
    
    # plot without data
    biome.plot <- plotSpatial(model.biomes,
                              text.multiplier = 2,
                              map.overlay = map.overlay)
    magicPlot(biome.plot, 
              filename = file.path(run.plot.dir, paste(biome.scheme@id, "Biomes",paste0(first.year, "-", last.year), analysis.label, sep = ".")))
    
    # plot with data
    biome.plot.with.data <- plotSpatialComparison(comparison.layer, 
                                                  type = "values",
                                                  map.overlay = map.overlay,
                                                  text.multiplier = 2,
                                                  facet.order = c(run@name, biome.data.field@source@name),
                                                  nrow = 2)
    magicPlot(biome.plot.with.data, 
              filename = file.path(run.plot.dir, paste(biome.scheme@id, "BiomesWithData", paste0(first.year, "-", last.year), analysis.label, sep = ".")))
    
  } # end for each run
  
  
  if(doMultiPlots){
    
    ### do each comparison group in turn
    for(group in comparison.groups) {
      
      # make a directory for the comparison
      local.group.dir <- file.path(plot.dir, "Groups", group$id)
      dir.create(local.group.dir, showWarnings = FALSE, recursive = TRUE)
      
      # extract the runs we want to plot together
      local.biomes.list <- list()
      for(run in group$runs){
        local.biomes.list[[run@id]] <- biome.objects.list[[run@id]]
      }
      
      # add the biome data
      local.biomes.list <- append(local.biomes.list, biome.data.field)
      
      # plot the biomes
      biome.group.plot <- plotSpatial(local.biomes.list,
                                      layers = biome.scheme@id,
                                      map.overlay = map.overlay,
                                      text.multiplier = 2,
                                      title = NULL)
      #p <- p + facet_wrap(~Facet, nrow = 4)
      biome.group.plot <- biome.group.plot + theme(legend.position='bottom')
      biome.group.plot <- biome.group.plot + theme(legend.title=element_blank())
      biome.group.plot <- biome.group.plot + guides(fill=guide_legend(ncol=2))
      
      magicPlot(biome.group.plot, 
                filename = file.path(local.group.dir, paste(biome.scheme@id, group$id, analysis.label, paste0(first.year, "-", last.year), sep = ".")))
      
      print(paste0("*** Done biome group: ", group$id))
      
    } # end for each group
    
  } # end if do multiplots 
  
  
  
  ######################################################
  #### NOW DIFFERENCE PLOTS
  
  if(doDifferencePlots) {
    
    for(this.pair in difference.pairs) {
      
      base <- this.pair$base
      new <- this.pair$new
      base.field <- biome.objects.list[[base@id]]
      new.field <- biome.objects.list[[new@id]]
      
      # make a directory for the comparison
      local.comparison.dir <- file.path(plot.dir, "Comparisons", paste(new@id, base@id, sep ="-"))
      dir.create(local.comparison.dir, showWarnings = FALSE, recursive = TRUE)
      
      # make the comparison layer
      this.comparison <- compareLayers(new.field, base.field, layers1 = biome.scheme@id, show.stats = FALSE)
      
      # plot the difference
      biome.diff.plot <- plotSpatialComparison(this.comparison, 
                                               map.overlay = "world",
                                               text.multiplier = 2)
      magicPlot(biome.diff.plot, 
                filename = file.path(local.comparison.dir, paste(biome.scheme@id, "Diff", analysis.label, paste0(first.year, "-", last.year),  sep = ".")))
      
      
      # plot the values
      biome.values.plot <- plotSpatialComparison(this.comparison, 
                                                 type = "values",
                                                 map.overlay = "world",
                                                 text.multiplier = 2)
      if(length(local.field.list) > 1) biome.values.plot <- biome.values.plot + facet_wrap(~Facet, nrow = 2)
      magicPlot(biome.values.plot, 
                filename = file.path(local.comparison.dir, paste(biome.scheme@id, "2-up", analysis.label, paste0(first.year, "-", last.year), sep = ".")))
      
      print(paste0("*** Done biome difference: ", paste(new@id, base@id, sep ="-")))
      
    } # end for each difference.pair
    
  } # end if doDifference plots
  
  print("*** Done biomes")
  
} # end if do biomes



t2 <- Sys.time()
print(t2-t1)

