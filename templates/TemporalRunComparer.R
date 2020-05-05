#!/usr/bin/R

### TODO
# 1. Loop across biomes

##### LIBRARIES ETC #####

library(DGVMTools)
library(raster)
library(Cairo)
source("~/Tools/DGVMTools_Scripts/utils/PlotUtils.R")

t1 <- Sys.time()


##### VARIOUS RUN SETTINGS #####

# Analysis label and plot directory
analysis.label <- "r8630"
plot.dir <- file.path("/home/mforrest/Projects/FireMIP/plots/global", analysis.label)
if(!file.exists(plot.dir)){dir.create(plot.dir)}

# save run-specific plots to run directory
savePlotsToRunDir <- FALSE


# Years
first.year = 1998
last.year = 2015

# Re-read 
reread <- FALSE

# preferred difference in between number of rows and columns (can be adjusted depending on
# the dimensions (Lon/Lat) of the plot extent), but is somewhat arivially asthetic
preferred.rows.more.than.cols <- 1

# plot type string
plot.type.string = "TS"

##### SELECT PLOT TYPES #####

# group comparisons and difference plots
doRunPlots <- TRUE # make plots for each run in it's own directory
doMultiPlots <- TRUE # make multipanel plots for the groups defined above (with a panel for each run)
doDifferencePlots <- TRUE # do difference plots for run pairs (pairs defined above)

# only for plots with multiple layers
doLayerPlots <- FALSE # plot individual layers (usually) PFTs
doAggregates <- TRUE # plot some aggregates 

# only for monthly variables
doMonthly <- FALSE
doSeasonal <- TRUE 


#### DEFINE THE RUNS #####

# Reference runs
Old_PNV <- defineSource(id = "Old_PNV",
                        name = "Old PNV",
                        dir = "/home/mforrest/Projects/FireMIP/runs/global/r8617/Base",
                        format = GUESS,
                        forcing.data = "CRUJRA")

# r8630
PNV <- defineSource(id = "PNV",
                    name = "PNV",
                    dir = "/home/mforrest/Projects/FireMIP/runs/global/r8630/PNV",
                    format = GUESS,
                    forcing.data = "CRUJRA")

# r8630
PNV_openmpi <- defineSource(id = "PNV_openmpi",
                            name = "PNV OpenMPI",
                            dir = "/home/mforrest/Projects/FireMIP/runs/global/r8630_openmpi/PNV",
                            format = GUESS,
                            forcing.data = "CRUJRA")

LUH2 <- defineSource(id = "LUH2",
                     name = "LUH2 (gross)",
                     dir = "/home/mforrest/Projects/FireMIP/runs/global/r8630_openmpi/LUH2",
                     format = GUESS,
                     forcing.data = "CRUJRA")

LUH2_net <- defineSource(id = "LUH2_net",
                         name = "LUH2 (net)",
                         dir = "/home/mforrest/Projects/FireMIP/runs/global/r8630_openmpi/LUH2_net",
                         format = GUESS,
                         forcing.data = "CRUJRA")


##### MAKE RUN LIST #####
runs <- list(
  
  Old_PNV,
  PNV,
  PNV_openmpi,
  # SF0,
  # SF3,
  LUH2,
  LUH2_net
  
)

##### DEFINE PLOT GROUPS #####
## Note that the runs must be included in th "runs" list above
run.groups <- list(
  
  list(runs = list(Old_PNV,
                   PNV,
                   PNV_openmpi),
       name = "PNV",
       id = "PNV"
  ),
  
  list(runs = list(PNV,
                   LUH2,
                   LUH2_net),
       name = "LandUse",
       id = "LandUse"
  )
  
)

##### DEFINE RUN PAIRS #####
## These runs will be compared directly against each otherS

run.pairs <- list(
  
  list("base" = PNV, "new" = LUH2),
  list("base" = LUH2, "new" = LUH2_net)
)

##### DEFINE AND SELECT REGIONS #####

regions <- list(
  # Global = list(id = "Global", name = "Global", extent = extent(-180,180,-90,90)),
  # NorthernHemisphere = list(id = "NorthernHemisphere", name = "Northern Hemisphere", extent = extent(-180,180,0,90)),
  # SouthernHemisphere = list(id = "SouthernHemisphere", name = "Southern Hemisphere", extent = extent(-180,180,-90,0)),
  # Africa = list( id = "Africa", name = "Africa", extent =  extent(-20, 55, -30, 36)),
  # Europe = list( id = "Europe", name = "Europe", extent =  extent(-30, 40, 36, 70)),
  # Asia = list( id = "Asia", name = "Asia", extent =  extent(40, 180, -10, 80)),
  # NorthAmerica = list( id = "NorthAmerica", name = "North America", extent =  extent(-170, -70, 25, 75)),
  # SouthAmerica = list( id = "SouthAmerica", name = "South America", extent = extent(-180, -50, -60, 25)),
  # Australia = list( id = "Australia", name = "Australia", extent = extent(110, 160, -45 ,10)),
  # Mediterranean = list( id = "Med", name = "Mediterranean", extent = extent(10, 40, 28 ,48)),
  # CentralAsia = list( id = "CentralAsia", name = "Central Asia", extent = extent(25, 140, 40, 55)),
  # SouthEastAsia = list( id = "SouthEastAsia", name = "South East Asia", extent = extent(90, 140, 10, 40)),
  # CentralNorthAmerica = list( id = "CentralNorthAmerica", name = "Central North America", extent = extent(-110, -85, 30, 50)),
  # Boreal = list( id = "Boreal", name = "Boreal", extent = extent(-180, 180, 60, 90)),
  # NHAfrica = list( id = "NHAfrica", name = "Northern Hemisphere Africa", extent = extent(-20, 50, 0, 25)),
  SHAfrica = list( id = "SHAfrica", name = "Southern Hemisphere Africa", extent = extent(5, 50, -30, 0))
)


#### VARIABLES AND PFT AGGREGATE SELECTION ####

# plot some variables for each run
vars.to.plot <- list(
  
  #list(var = "lai", cuts = seq(0, 10, 0.5), PFT = TRUE)
  #flist(var = "fpc", cuts = seq(0, 1.3, 0.1), PFT = TRUE),
  #list(var = "mfirefrac", cuts = c(0,0.002,0.005,0.01,0.02,0.05,0.10,0.2,0.50,1), PFT = FALSE, agg.method = "sum")#,
  list(var = "monthly_burned_area", 
       cuts = c(0,0.002,0.005,0.01,0.02,0.05,0.10,0.2,0.50,1), 
       PFT = FALSE, 
       spatial.agg.method = "w.sum",
       subannual.agg.method = "sum")
  #,
  #list(var = "monthly_burned_area", cuts = c(0,0.002,0.005,0.01,0.02,0.05,0.10,0.2,0.50,2), PFT = FALSE, agg.method = "sum")#,
  #list(var = "mSAV", PFT = FALSE, agg.method = "mean"),
  #list(var = "mFBD", PFT = FALSE, agg.method = "mean")#,
  #list(var = "mRoS", PFT = FALSE, agg.method = "mean"),
  #list(var = "meff_wind", PFT = FALSE, agg.method = "mean"),
  #list(var = "mfire_size", PFT = FALSE, agg.method = "mean"),
  #list(var = "mfireintens", PFT = FALSE, agg.method = "mean")#,
  #list(var = "mMoE", PFT = FALSE, agg.method = "mean")
  
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


##### MAIN LOOP #####
## Normlly nothing to chaneg after here!


for(var.details in vars.to.plot) {
  
  ##### READ RUNS AND DO RUN PLOTS #####
  
  var.str <- var.details$var
  var.isPFT <- var.details$PFT
  var.spatial.agg.method <- var.details$spatial.agg.method
  print(var.str)
  
  field.list <- list()
  field.seasonal.list <- list()
  field.annual.list <- list()
  
  for(this.region in regions) {
    
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
                                   spatial.extent = this.region$extent,
                                   spatial.extent.id = this.region$id,
                                   spatial.aggregate.method = var.spatial.agg.method,
                                   read.full = reread,
                                   write = TRUE,
                                   verbose = TRUE)
      
      # determine if it is monthly
      isMonthly <- "Month" %in% getDimInfo(this.model.field)
      
      if(isMonthly) {
        # sum to annual and plot that
        this.model.field.annual <- aggregateSubannual(this.model.field, var.details$subannual.agg.method)
        # calculate seasonal aggregates if required
        if(doSeasonal) this.model.field.seasonal <- aggregateSubannual(this.model.field,  method = var.details$subannual.agg.method, target = "Season")
      }
      
      ##### RUN PLOTS #####
      if(doRunPlots) {
        
        summary.plot <- plotTemporal(this.model.field,
                                     text.multiplier = 2, 
                                     ncol = 1#nOptCols(this.model.field, preferred.rows.more.than.cols, "Temporal")    # calculate optimal number of columns, based on having one more row than columns
        )
        
        magicPlot(summary.plot, 
                  height = 900,
                  width = 1800,
                  filename = file.path(run.plot.dir, paste(var.str, plot.type.string, this.region$id, "Summary", paste0(first.year, "-", last.year), analysis.label, sep = ".")))
        
        ### If variable is do sums monthly, also 
        if(isMonthly) {
          
          annual.plot <- plotTemporal(this.model.field.annual,
                                      text.multiplier = 2,
                                      ncol = 1#nOptCols(this.model.field.annual, preferred.rows.more.than.cols, "Temporal")    # calculate optimal number of columns, based on having one more row than columns
          )
          
          magicPlot(annual.plot, 
                    height = 900,
                    width = 1800,
                    filename = file.path(run.plot.dir, paste(var.str, plot.type.string, this.region$id, "AnnualSum", paste0(first.year, "-", last.year), analysis.label, sep = ".")))
          
          if(doMonthly) {
            
            # plot months individually
            for(month in all.months) {
              month.plot <- plotTemporal(selectMonths(this.model.field, month@index), 
                                         text.multiplier = 2,
                                         ncol = 1#nOptCols(this.model.field, preferred.rows.more.than.cols, "Temporal")    # calculate optimal number of columns, based on having one more row than columns
              )
              magicPlot(month.plot, 
                        height = 900,
                        width = 1800,
                        filename = file.path(run.plot.dir, paste(var.str, plot.type.string, this.region$id, month@id, paste0(first.year, "-", last.year), analysis.label, sep = ".")))
            } # for each month
            
          } # if doMonthly
          
          if(doSeasonal) {
            
            # plot all seasons together
            season.plot <- plotTemporal(this.model.field.seasonal,
                                        text.multiplier = 2,
                                        ncol = 1#nOptCols(this.model.field.seasonal, preferred.rows.more.than.cols, "Temporal")    # calculate optimal number of columns, based on having one more row than columns
            )
            magicPlot(season.plot, 
                      height = 900,
                      width = 1800,
                      filename = file.path(run.plot.dir, paste(var.str, plot.type.string, this.region$id, "Seasonal", var.details$subannual.agg.method, paste0(first.year, "-", last.year), analysis.label, sep = ".")))
            
            # plot seasons individually
            for(season in all.seasons) {
              season.plot <- plotTemporal(selectSeasons(this.model.field.seasonal, season@id), 
                                          text.multiplier = 2)
              magicPlot(season.plot, 
                        height = 900,
                        width = 1800,
                        filename = file.path(run.plot.dir, paste(var.str, plot.type.string, this.region$id, season@id, var.details$subannual.agg.method, paste0(first.year, "-", last.year), analysis.label, sep = ".")))
            }
            
          }
          
        }
        
        
        ### individual plots
        if(doLayerPlots) {
          
          for(this.layer in names(this.model.field)) {
            individual.plot <- plotTemporal(this.model.field,
                                            layers = this.layer,
                                            text.multiplier = 2,
                                            ncol = 1#nOptCols(this.model.field, preferred.rows.more.than.cols, "Temporal", nlayers = 1)
            )
            
            magicPlot(individual.plot, 
                      height = 900,
                      width = 1800,
                      filename = file.path(run.plot.dir, paste(var.str, plot.type.string, this.region$id, this.layer, paste0(first.year, "-", last.year), analysis.label, sep = ".")))
          }
          message("*** Done initial layer plots")  
        }
        
        ### aggregate plots  
        if(doAggregates && var.isPFT) {
          
          for(this.aggregate in PFT.aggs){
            
            
            this.arguments <- list(x = this.model.field)
            this.model.field <- do.call(layerOp, append(this.arguments, this.aggregate))
            
            
            if(this.aggregate$new.layer %in% names(this.model.field)) {
              aggregate.plot <- plotTemporal(this.model.field,
                                             layers = this.aggregate$new.layer,
                                             text.multiplier = 2)
              magicPlot(aggregate.plot, 
                        height = 900,
                        width = 1800,
                        filename = file.path(run.plot.dir, paste(var.str, plot.type.string, this.region$id, this.aggregate$new.layer, paste0(first.year, "-", last.year), analysis.label, sep = ".")))
            }
          }
          message("*** Done aggregate layer plots")  
          
        } # for each aggregate
        
      } # end if doRunPlots
      
      # save for plotting all together later and comparison plots
      field.list[[run@id]] <- this.model.field
      if(isMonthly && doSeasonal) field.seasonal.list[[run@id]] <- this.model.field.seasonal
      if(isMonthly) field.annual.list[[run@id]] <- this.model.field.annual
      
      rm(this.model.field, this.model.field.seasonal, this.model.field.annual)
      gc()
      
      
    } # for each run
    
    
    #####  GROUP (MULTIPANEL) PLOTS  ######
    ##  IF REQUESTED MAKE MULTIPANEL PLOTS - one panel per run
    
    if(doMultiPlots) {
      
      ### do each comparison group in turn
      for(group in run.groups) {
        
        # make a directory for the comparison groups
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
            group.layer.plot <- plotTemporal(local.field.list,
                                             layers = this.layer,
                                             col.by = "Source",
                                             text.multiplier = 2)
            magicPlot(group.layer.plot,  
                      height = 900,
                      width = 1800,
                      filename = file.path(local.group.dir, paste(var.str, plot.type.string, this.region$id, group$id, this.layer, paste0(first.year, "-", last.year), analysis.label, sep = ".")))
          } # for each layer
          
        } # if doIndividual plots
        
        ### If variable is do sums monthly, also 
        if(isMonthly) {
          
          # plot the annual sums
          annual.plot <- plotTemporal(local.field.annual.list,
                                      text.multiplier = 2,
                                      col.by = "Source",
                                      ncol = 1#nOptCols(length(local.field.annual.list), preferred.rows.more.than.cols)
          )
          magicPlot(annual.plot, 
                    height = 900,
                    width = 1800,
                    filename = file.path(local.group.dir, paste(var.str, plot.type.string, this.region$id, "AnnualSum", paste0(first.year, "-", last.year), analysis.label, sep = ".")))
          
          if(doMonthly) {
            
            # plot months individually
            for(month in all.months) {
              
              temp.field.list <- list()
              for(temp.Field in local.field.list) {
                temp.field.list[[temp.Field@id]] <-selectMonths(temp.Field, month@index)
              }
              
              month.plot <- plotTemporal(temp.field.list,
                                         text.multiplier = 2,
                                         col.by = "Source",
                                         ncol = 1#nOptCols(length(local.field.list), preferred.rows.more.than.cols)
              )
              magicPlot(month.plot, 
                        height = 900,
                        width = 1800,
                        filename = file.path(local.group.dir, paste(var.str, plot.type.string, this.region$id, month@id, paste0(first.year, "-", last.year), analysis.label, sep = ".")))
              rm(temp.field.list)
              
            }
            
          }
          
          if(doSeasonal) {
            
            # plot all seasons together
            season.plot <- plotTemporal(local.field.seasonal.list, 
                                        col.by = "Source",
                                        text.multiplier = 2)
            
            magicPlot(season.plot, 
                      height = 900,
                      width = 1800,
                      filename = file.path(local.group.dir, paste(var.str, plot.type.string, this.region$id, "Seasons", var.details$subannual.agg.method, paste0(first.year, "-", last.year), analysis.label, sep = ".")))
            
            # plot seasons individually
            for(season in all.seasons) {
              
              temp.field.list <- list()
              for(temp.Field in local.field.seasonal.list) {
                temp.field.list[[temp.Field@id]] <-selectSeasons(temp.Field, season@id)
              }
              
              season.plot <- plotTemporal(temp.field.list, 
                                          text.multiplier = 2,
                                          col.by = "Source",
                                          ncol = 1#nOptCols(length(local.field.seasonal.list), preferred.rows.more.than.cols, ntimes = 1)
              )
              magicPlot(season.plot, 
                        height = 900,
                        width = 1800,
                        filename = file.path(local.group.dir, paste(var.str, plot.type.string, this.region$id, season@id, var.details$subannual.agg.method, paste0(first.year, "-", last.year), analysis.label, sep = ".")))
              
              rm(temp.field.list)
              
            } # for each season
            
          } # if do seasonal plots 
          
        } # if variable is monthly
        
        
        ### layer aggregates
        if(doAggregates && var.isPFT) {
          
          for(this.aggregate in PFT.aggs){
            
            
            group.aggregate.plot <- plotTemporal(local.field.list,
                                                 layers = this.aggregate$new.layer,
                                                 text.multiplier = 2)
            if(length(local.field.list) > 1) group.aggregate.plot <- group.aggregate.plot + facet_wrap(~Facet, nrow = 2)
            magicPlot(group.aggregate.plot, 
                      height = 900,
                      width = 1800,
                      filename = file.path(local.group.dir, paste(var.str, plot.type.string, this.region$id, group$id, this.aggregate$new.layer, paste0(first.year, "-", last.year), analysis.label, sep = ".")))
          } # for each aggregate
          
        } # if doAggregates
        
        print(paste0("*** Done group: ", group$name))
        
      } # for each comparison group 
      
    } # if do multiplots
    
    
    #####  DIFFERENCE PLOTS  #####
    
    if(doDifferencePlots) {
      
      for(this.pair in run.pairs) {
        
        base <- this.pair$base
        new <- this.pair$new
        base.field <- field.list[[base@id]]
        new.field <- field.list[[new@id]]
        
        # make a directory for the comparison
        local.comparison.dir <- file.path(plot.dir, "Pairs", paste(new@id, base@id, sep ="-"))
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
          nMonths <- 1
          if("Month" %in% getDimInfo(this.comparison)) nMonths <- getDimInfo(this.comparison, "size")$Month
          diff.layer.plot<- plotTemporalComparison(this.comparison, 
                                                   text.multiplier = 2,
                                                   col.by = "Source",
                                                   ncol = 1 #nOptCols(nMonths, rows.more.than.cols = preferred.rows.more.than.cols)
          )
          magicPlot(diff.layer.plot, 
                    height = 900,
                    width = 1800,
                    filename = file.path(local.comparison.dir, paste(layer, quant@id, plot.type.string, this.region$id, "Diff", analysis.label, paste0(first.year, "-", last.year), sep = ".")))
          
          # plot the values
          abs.layer.plot <- plotTemporalComparison(this.comparison, 
                                                   type = "values",
                                                   text.multiplier = 2,
                                                   col.by = "Source",
                                                   ncol = 1# nOptCols(nMonths * 2, rows.more.than.cols = preferred.rows.more.than.cols)
          )
          magicPlot(abs.layer.plot, 
                    height = 900,
                    width = 1800,
                    filename = file.path(local.comparison.dir, paste(layer, quant@id, plot.type.string, this.region$id, "2-up", analysis.label, paste0(first.year, "-", last.year), sep = ".")))
          
          # do Difference plots for Monthly variables
          ### If variable is do sums monthly, also 
          if(isMonthly) {
            
            # make the annual comparison layer
            this.comparison <- compareLayers(field.annual.list[[new@id]], field.annual.list[[base@id]], layers1 = layer, show.stats = FALSE)
            
            # plot difference
            diff.layer.plot <- plotTemporalComparison(this.comparison, 
                                                      text.multiplier = 2,
                                                      col.by = "Source",
                                                      ncol = 1)
            magicPlot(diff.layer.plot, 
                      height = 900,
                      width = 1800,
                      filename = file.path(local.comparison.dir, paste(layer, quant@id, plot.type.string, this.region$id, "AnnualSum", "Diff", analysis.label, paste0(first.year, "-", last.year), sep = ".")))
            
            # plot the values
            abs.layer.plot <- plotTemporalComparison(this.comparison, 
                                                     type = "values",
                                                     text.multiplier = 2,
                                                     col.by = "Source",
                                                     ncol = 1)
            magicPlot(abs.layer.plot, 
                      height = 900,
                      width = 1800,
                      filename = file.path(local.comparison.dir, paste(layer, quant@id, plot.type.string, this.region$id, "AnnualSum", "2-up", analysis.label, paste0(first.year, "-", last.year), sep = ".")))
            

            if(doMonthly) {
              
              # plot months individually
              for(month in all.months) {
                
                # make the annual comparison layer
                this.comparison <- compareLayers(selectMonths(new.field, month@index), selectMonths(base.field, month@index), layers1 = layer, show.stats = FALSE)
                
                diff.layer.plot <- plotTemporalComparison(this.comparison, 
                                                          text.multiplier = 2,
                                                          col.by = "Source",
                                                          ncol = 1)
                magicPlot(diff.layer.plot, 
                          height = 900,
                          width = 1800,
                          filename = file.path(local.comparison.dir, paste(layer, quant@id, plot.type.string, this.region$id, month@abbreviation, "Diff", analysis.label, paste0(first.year, "-", last.year), sep = ".")))
                
                # plot the values
                abs.layer.plot <- plotTemporalComparison(this.comparison, 
                                                         type = "values",
                                                         text.multiplier = 2,
                                                         col.by = "Source",
                                                         ncol = 1)
                magicPlot(abs.layer.plot, 
                          height = 900,
                          width = 1800,
                          filename = file.path(local.comparison.dir, paste(layer, quant@id, plot.type.string ,this.region$id, month@abbreviation, "2-up", analysis.label, paste0(first.year, "-", last.year), sep = ".")))
              } # for each month
              
            } # end if doMonthly
            
            if(doSeasonal) {
              
              base.field <- field.seasonal.list[[base@id]]
              new.field <- field.seasonal.list[[new@id]]
              
              seasonal.comparison <- compareLayers(new.field, base.field, layers1 = layer, show.stats = FALSE)
              
              # plot the difference
              nSeasons <- 1
              if("Season" %in% getDimInfo(seasonal.comparison)) nSeasons <- getDimInfo(seasonal.comparison, "size")$Season
              
              
              # plot difference
              diff.layer.plot <- plotTemporalComparison(seasonal.comparison, 
                                                        text.multiplier = 2,
                                                        col.by = "Source",
                                                        ncol = 1#nOptCols(nSeasons, rows.more.than.cols = preferred.rows.more.than.cols)
              )
              
              magicPlot(diff.layer.plot, 
                        height = 900,
                        width = 1800,
                        filename = file.path(local.comparison.dir, paste(layer, quant@id, plot.type.string, this.region$id, "SeasonalSum", "Diff", analysis.label, paste0(first.year, "-", last.year), sep = ".")))
              
              # plot the values
              abs.layer.plot <- plotTemporalComparison(seasonal.comparison, 
                                                       type = "values",
                                                       col.by = "Source",
                                                       text.multiplier = 2)
              
              magicPlot(abs.layer.plot, 
                        height = 900,
                        width = 1800,
                        filename = file.path(local.comparison.dir, paste(layer, quant@id, plot.type.string, this.region$id, "SeasonalSum", "2-up", analysis.label, paste0(first.year, "-", last.year), sep = ".")))
              
              # plot seasons individually
              for(season in all.seasons) {
                
                # make the annual comparison layer
                diff.layer.plot <- plotTemporalComparison(seasonal.comparison, 
                                                          seasons = season@abbreviation,
                                                          text.multiplier = 2,
                                                          col.by = "Source",
                                                          ncol = 1)
                magicPlot(diff.layer.plot, 
                          height = 900,
                          width = 1800,
                          filename = file.path(local.comparison.dir, paste(layer, quant@id, plot.type.string, this.region$id, season@abbreviation, "Diff", analysis.label, paste0(first.year, "-", last.year), sep = ".")))
                
                # plot the values
                abs.layer.plot <- plotTemporalComparison(this.comparison, 
                                                         type = "values",
                                                         text.multiplier = 2,
                                                         col.by = "Source",
                                                         ncol = 1)
                magicPlot(abs.layer.plot, 
                          height = 900,
                          width = 1800,
                          filename = file.path(local.comparison.dir, paste(layer, quant@id, plot.type.string, this.region$id, season@abbreviation, "2-up", analysis.label, paste0(first.year, "-", last.year), sep = ".")))
                
              } # for each season
              
            } # if do seasonal plots 
            
          } # if variable is monthly
          
        } # for each layer 
        
        print(paste0("*** Done difference: ", paste(new@id, base@id, sep ="-")))
        
      } # end for each difference pair
      
    } # end if doDifference plots
    
  } # for each region
  
  
} # for each variable 


t2 <- Sys.time()
print(t2-t1)

