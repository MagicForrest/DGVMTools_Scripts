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
analysis.label <- "Example"
plot.dir <- "/home/forrest/TestPlots"
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

# preferred difference in between number of rows and columns (can be adjusted depending on
# the dimensions (Lon/Lat) of the plot extent), but is somewhat arivially asthetic
preferred.rows.more.than.cols <- 1

##### SELECT PLOT TYPES #####

# group comparisons and difference plots
doRunPlots <- FALSE # make plots for each run in it's own directory
doMultiPlots <- TRUE # make multipanel plots for the groups defined above (with a panel for each run)
doDifferencePlots <- TRUE # do difference plots for run pairs (pairs defined above)

# only for plots with multiple layers
doLayerPlots <- TRUE # plot individual layers (usually) PFTs
doAggregates <- TRUE # plot some aggregates 

# only for monthly variables
doMonthly <- TRUE 
doSeasonal <- TRUE 


##### DEFINE THE RUNS #####
## But not that these runs won't be read or plotted if the are not includes in "runs",
## "run.groups" or "run.pairs"

# wind limit (r8511)
PNV_SPITFIRE_NoWindLimit <- defineSource(id = "PNV_SPITFIRE_NoWindLimit",
                                         name = "SPITFIRE (PNV, No Wind Limit)",
                                         dir = "/home/forrest/GuessRuns/FireMIP/PNV_SPITFIRE_NoWindLimit/",
                                         format = GUESS, 
                                         forcing.data = "CRUJRA")


PNV_SPITFIRE_LasslopWindLimit <- defineSource(id = "PNV_SPITFIRE_LasslopWindLimit",
                                              name = "SPITFIRE (PNV, Lasslop Wind Limit)",
                                              dir = "/home/forrest/GuessRuns/FireMIP/PNV_SPITFIRE_LasslopWindLimit/",
                                              format = GUESS, 
                                              forcing.data = "CRUJRA")





##### RUNS FOR PLOTTING #####

runs <- list(
  
  PNV_SPITFIRE_NoWindLimit,
  PNV_SPITFIRE_LasslopWindLimit
  
)

###### DEFINE GROUPS #####

run.groups <- list(
  
  list(runs = list(PNV_SPITFIRE_NoWindLimit,
                   PNV_SPITFIRE_LasslopWindLimit),
       name = "Wind Limits",
       id = "WindLimit"
  )
  
)


##### DEFINE RUN PAIRS #####
## These runs will be compared directly against each otherS

run.pairs <- list(
  list("base" =  PNV_SPITFIRE_NoWindLimit, "new" = PNV_SPITFIRE_LasslopWindLimit)
)



#### VARIABLES AND PFT AGGREGATE SELECTION ####

# plot some variables for each run
vars.to.plot <- list(
  
  #lai = list(var = "lai", cuts = seq(0, 10, 0.5), PFT = TRUE)#,
  #fpc = list(var = "fpc", cuts = seq(0, 1.3, 0.1), PFT = TRUE),
  mfirefrac = list(var = "mfirefrac", cuts = c(0,0.002,0.005,0.01,0.02,0.05,0.10,0.2,0.50,1), PFT = FALSE, agg.method = "sum")
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
  var.cuts <- var.details$cuts
  var.isPFT <- var.details$PFT
  print(var.str)
  
  field.list <- list()
  field.seasonal.list <- list()
  field.annual.list <- list()
  
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
    
    if(isMonthly) {
      # sum to annual and plot that
      this.model.field.annual <- aggregateSubannual(this.model.field, var.details$agg.method)
      # calculate seasonal aggregates if required
      if(doSeasonal) this.model.field.seasonal <- aggregateSubannual(this.model.field, var.details$agg.method, target = "Season")
    }
    
    ##### RUN PLOTS #####
    if(doRunPlots) {
      
      summary.plot <- plotSpatial(this.model.field,
                                  cuts = var.cuts,
                                  map.overlay = map.overlay,
                                  text.multiplier = 2, 
                                  ncol = nOptCols(this.model.field, preferred.rows.more.than.cols, "Spatial")    # calculate optimal number of columns, based on having one more row than columns
      )
      
      magicPlot(summary.plot, 
                filename = file.path(run.plot.dir, paste(var.str, "Summary", paste0(first.year, "-", last.year), analysis.label, sep = ".")))
      
      ### If variable is do sums monthly, also 
      if(isMonthly) {
        
          annual.plot <- plotSpatial(this.model.field.annual,
                                   cuts = var.cuts,
                                   map.overlay = map.overlay,
                                   text.multiplier = 2,
                                   ncol = nOptCols(this.model.field.annual, preferred.rows.more.than.cols, "Spatial")    # calculate optimal number of columns, based on having one more row than columns
        )
        
        magicPlot(annual.plot, 
                  filename = file.path(run.plot.dir, paste(var.str, "AnnualSum", paste0(first.year, "-", last.year), analysis.label, sep = ".")))
        
        if(doMonthly) {
          
          # plot months individually
          for(month in all.months) {
            month.plot <- plotSpatial(this.model.field, 
                                      month = month@index,
                                      cuts = var.cuts,
                                      map.overlay = map.overlay,
                                      text.multiplier = 2,
                                      ncol = nOptCols(this.model.field, preferred.rows.more.than.cols, "Spatial")    # calculate optimal number of columns, based on having one more row than columns
            )
            magicPlot(month.plot, 
                      filename = file.path(run.plot.dir, paste(var.str, month@id, paste0(first.year, "-", last.year), analysis.label, sep = ".")))
          } # for each month
          
        } # if doMonthly
        
        if(doSeasonal) {
          
          # plot all seasons together
          season.plot <- plotSpatial(this.model.field.seasonal,
                                     cuts = var.cuts,
                                     map.overlay = map.overlay,
                                     text.multiplier = 2,
                                     ncol = nOptCols(this.model.field.seasonal, preferred.rows.more.than.cols, "Spatial")    # calculate optimal number of columns, based on having one more row than columns
          )
          magicPlot(season.plot, 
                    filename = file.path(run.plot.dir, paste(var.str, "Seasons", var.details$agg.method, paste0(first.year, "-", last.year), analysis.label, sep = ".")))
          
          # plot seasons individually
          for(season in all.seasons) {
            season.plot <- plotSpatial(this.model.field.seasonal, 
                                       seasons = season@id,
                                       cuts = var.cuts,
                                       map.overlay = map.overlay,
                                       text.multiplier = 2)
            magicPlot(season.plot, 
                      filename = file.path(run.plot.dir, paste(var.str, season@id, var.details$agg.method, paste0(first.year, "-", last.year), analysis.label, sep = ".")))
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
                                         cuts = var.cuts, 
                                         ncol = nOptCols(this.model.field, preferred.rows.more.than.cols, "Spatial", nlayers = 1))
          
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
      
    } # end if doRunPlots
    
    # save for plotting all together later and comparison plots
    field.list[[run@id]] <- this.model.field
    if(isMonthly && doSeasonal) field.seasonal.list[[run@id]] <- this.model.field.seasonal
    if(isMonthly) field.annual.list[[run@id]] <- this.model.field.annual
    
    
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
          group.layer.plot <- plotSpatial(local.field.list,
                                          layers = this.layer,
                                          map.overlay = map.overlay,
                                          text.multiplier = 2,
                                          cuts = var.cuts)
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
                                   text.multiplier = 2,
                                   ncol = nOptCols(length(local.field.annual.list), preferred.rows.more.than.cols))
        magicPlot(annual.plot, 
                  filename = file.path(local.group.dir, paste(var.str, "AnnualSum", paste0(first.year, "-", last.year), analysis.label, sep = ".")))
        
        if(doMonthly) {
          
          # plot months individually
          for(month in all.months) {
            month.plot <- plotSpatial(local.field.list, 
                                      month = month@index,
                                      cuts = var.cuts,
                                      map.overlay = map.overlay,
                                      text.multiplier = 2,
                                      ncol = nOptCols(length(local.field.list), preferred.rows.more.than.cols))
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
          season.plot <- season.plot + facet_grid(cols = vars(Season), rows = vars(Source), switch = "y")
          magicPlot(season.plot, 
                    filename = file.path(local.group.dir, paste(var.str, "Seasons", var.details$agg.method, paste0(first.year, "-", last.year), analysis.label, sep = ".")))
          
          # plot seasons individually
          for(season in all.seasons) {
            
            season.plot <- plotSpatial(local.field.seasonal.list, 
                                       seasons = season@id,
                                       cuts = var.cuts,
                                       map.overlay = map.overlay,
                                       text.multiplier = 2                                       ,
                                       ncol = nOptCols(length(local.field.seasonal.list), preferred.rows.more.than.cols, ntimes = 1))
            magicPlot(season.plot, 
                      filename = file.path(local.group.dir, paste(var.str, season@id, var.details$agg.method, paste0(first.year, "-", last.year), analysis.label, sep = ".")))
            
          } # for each season
          
        } # if do seasonal plots 
        
      } # if variable is monthly
      
      
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
        diff.layer.plot<- plotSpatialComparison(this.comparison, 
                                                map.overlay = "world",
                                                text.multiplier = 2,
                                                ncol = nOptCols(nMonths, rows.more.than.cols = preferred.rows.more.than.cols))
        magicPlot(diff.layer.plot, 
                  filename = file.path(local.comparison.dir, paste(layer, quant@id, "Diff", analysis.label, paste0(first.year, "-", last.year), sep = ".")))
        
        # plot the values
        abs.layer.plot <- plotSpatialComparison(this.comparison, 
                                                type = "values",
                                                map.overlay = "world",
                                                text.multiplier = 2,
                                                
                                                ncol = nOptCols(nMonths * 2, rows.more.than.cols = preferred.rows.more.than.cols))
        magicPlot(abs.layer.plot, 
                  filename = file.path(local.comparison.dir, paste(layer, quant@id, "2-up", analysis.label, paste0(first.year, "-", last.year), sep = ".")))
        
        # do Difference plots for Monthly variables
        ### If variable is do sums monthly, also 
        if(isMonthly) {
          
          # make the annual comparison layer
          this.comparison <- compareLayers(field.annual.list[[new@id]], field.annual.list[[base@id]], layers1 = layer, show.stats = FALSE)
          
          # plot difference
          diff.layer.plot <- plotSpatialComparison(this.comparison, 
                                                   map.overlay = "world",
                                                   text.multiplier = 2,
                                                   ncol = 1)
          magicPlot(diff.layer.plot, 
                    filename = file.path(local.comparison.dir, paste(layer, quant@id, "AnnualSum", "Diff", analysis.label, paste0(first.year, "-", last.year), sep = ".")))
          
          # plot the values
          abs.layer.plot <- plotSpatialComparison(this.comparison, 
                                                  type = "values",
                                                  map.overlay = "world",
                                                  text.multiplier = 2,
                                                  ncol = 1)
          magicPlot(abs.layer.plot, 
                    filename = file.path(local.comparison.dir, paste(layer, quant@id, "AnnualSum", "2-up", analysis.label, paste0(first.year, "-", last.year), sep = ".")))
          
          if(doMonthly) {
            
            # plot months individually
            for(month in all.months) {
              
              # make the annual comparison layer
              this.comparison <- compareLayers(selectMonths(new.field, month@index), selectMonths(base.field, month@index), layers1 = layer, show.stats = FALSE)
              
              diff.layer.plot <- plotSpatialComparison(this.comparison, 
                                                       map.overlay = "world",
                                                       text.multiplier = 2,
                                                       ncol = 1)
              magicPlot(diff.layer.plot, 
                        filename = file.path(local.comparison.dir, paste(layer, quant@id, month@abbreviation, "Diff", analysis.label, paste0(first.year, "-", last.year), sep = ".")))
              
              # plot the values
              abs.layer.plot <- plotSpatialComparison(this.comparison, 
                                                      type = "values",
                                                      map.overlay = "world",
                                                      text.multiplier = 2,
                                                      cuts = var.cuts,
                                                      ncol = 1)
              magicPlot(abs.layer.plot, 
                        filename = file.path(local.comparison.dir, paste(layer, quant@id, month@abbreviation, "2-up", analysis.label, paste0(first.year, "-", last.year), sep = ".")))
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
            diff.layer.plot <- plotSpatialComparison(seasonal.comparison, 
                                                     map.overlay = "world",
                                                     text.multiplier = 2,
                                                     ncol = nOptCols(nSeasons, rows.more.than.cols = preferred.rows.more.than.cols))
            
            magicPlot(diff.layer.plot, 
                      filename = file.path(local.comparison.dir, paste(layer, quant@id, "SeasonalSum", "Diff", analysis.label, paste0(first.year, "-", last.year), sep = ".")))
            
            # plot the values
            abs.layer.plot <- plotSpatialComparison(seasonal.comparison, 
                                                    type = "values",
                                                    map.overlay = "world",
                                                    cuts = var.cuts,
                                                    text.multiplier = 2)
            abs.layer.plot <- abs.layer.plot + facet_grid(cols = vars(Season), rows = vars(Source), switch = "y")
            
            magicPlot(abs.layer.plot, 
                      filename = file.path(local.comparison.dir, paste(layer, quant@id, "SeasonalSum", "2-up", analysis.label, paste0(first.year, "-", last.year), sep = ".")))
            
            # plot seasons individually
            for(season in all.seasons) {
              
              # make the annual comparison layer
              diff.layer.plot <- plotSpatialComparison(seasonal.comparison, 
                                                       map.overlay = "world",
                                                       seasons = season@abbreviation,
                                                       text.multiplier = 2,
                                                       ncol = 1)
              magicPlot(diff.layer.plot, 
                        filename = file.path(local.comparison.dir, paste(layer, quant@id, season@abbreviation, "Diff", analysis.label, paste0(first.year, "-", last.year), sep = ".")))
              
              # plot the values
              abs.layer.plot <- plotSpatialComparison(this.comparison, 
                                                      type = "values",
                                                      map.overlay = "world",
                                                      text.multiplier = 2,
                                                      cuts = var.cuts,
                                                      ncol = 1)
              magicPlot(abs.layer.plot, 
                        filename = file.path(local.comparison.dir, paste(layer, quant@id, season@abbreviation, "2-up", analysis.label, paste0(first.year, "-", last.year), sep = ".")))
              
            } # for each season
            
          } # if do seasonal plots 
          
        } # if variable is monthly
        
      } # for each layer 
      
      print(paste0("*** Done difference: ", paste(new@id, base@id, sep ="-")))
      
    } # end for each difference pair
    
  } # end if doDifference plots
  
  
} # for each variable 


t2 <- Sys.time()
print(t2-t1)

