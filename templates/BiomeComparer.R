#!/usr/bin/R

### TODO
# 1. Loop across biomes


library(DGVMTools)
library(raster)
library(Cairo)
source("~/Tools/DGVMTools_Scripts/utils/PlotUtils.R")

t1 <- Sys.time()


##########################################################################################################
################ HERE SUPPLY VARIOUS RUN SETTINGS 

# Analysis label and plot directory
analysis.label <- "r8572"
plot.dir <- "/home/matthew/Projects/FireMIP/plots/spatial/r8572"
if(!file.exists(plot.dir)){dir.create(plot.dir)}

DGVMData.dir <- "/home/matthew/DGVMData/"


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

##########################################################################################################
################ HERE SELECT WHICH TYPES OF PLOTS TO MAKE

# group comparisons and difference plots
doMultiPlots <- TRUE # make multipanel plots for the groups defined above (with a panel for each run)
doDifferencePlots <- TRUE # do difference plots for run pairs (pairs defined above)

# only for plots with multiple layers
doLayerPlots <- FALSE # plot individual layers (usually) PFTs
doAggregates <- FALSE # plot some aggregates 

# only for monthly variables
doMonthly <- FALSE 
doSeasonal <- FALSE 

# if biome classification - maybe move this to another script
doBiomes <- TRUE 



##### DEFINE THE RUNS #####

# r8572
Base <- defineSource(id = "Base",
                     name = "Base",
                     dir = "/spare/GuessRuns/r8572/Base",
                     format = GUESS, 
                     forcing.data = "CRUJRA")

Daylength<- defineSource(id = "Daylength",
                         name = "Daylength",
                         dir = "/spare/GuessRuns/r8572/Daylength",
                         format = GUESS, 
                         forcing.data = "CRUJRA")

NoWindLimit<- defineSource(id = "NoWindLimit",
                           name = "NoWindLimit",
                           dir = "/spare/GuessRuns/r8572/NoWindLimit",
                           format = GUESS, 
                           forcing.data = "CRUJRA")

HumIgn1.0 <- defineSource(id = "HumIgn1.0",
                          name = "HumIgn1.0",
                          dir = "/spare/GuessRuns/r8572/HumIgn1.0",
                          format = GUESS, 
                          forcing.data = "CRUJRA")

# r5874
Base_r8574 <- defineSource(id = "Base_r8574",
                           name = "Base (r8574)",
                           dir = "/spare/GuessRuns/r8574/Base",
                           format = GUESS, 
                           forcing.data = "CRUJRA")


# r5875
HoffmanFBD <- defineSource(id = "HoffmanFBD",
                           name = "Hoffman FBD",
                           dir = "/spare/GuessRuns/r8575/HoffmanFBD",
                           format = GUESS, 
                           forcing.data = "CRUJRA")

AllHoffmanFBD <- defineSource(id = "AllHoffmanFBD",
                              name = "AllHoffman FBD",
                              dir = "/spare/GuessRuns/r8575/AllHoffmanFBD",
                              format = GUESS, 
                              forcing.data = "CRUJRA")

HoffmanFBD_Daylength <- defineSource(id = "HoffmanFBD_Daylength",
                                     name = "Hoffman FBD, Daylength",
                                     dir = "/spare/GuessRuns/r8575/HoffmanFBD_Daylength",
                                     format = GUESS, 
                                     forcing.data = "CRUJRA")

AllHoffmanFBD_Daylength <- defineSource(id = "AllHoffmanFBD_Daylength",
                                        name = "AllHoffman FBD, Daylength",
                                        dir = "/spare/GuessRuns/r8575/AllHoffmanFBD_Daylength",
                                        format = GUESS, 
                                        forcing.data = "CRUJRA")


# r8578
HoffmanFBD_Daylength_Sigma <- defineSource(id = "HoffmanFBD_Daylength_Sigma",
                                           name = "Hoffman FBD, Daylength, Sigma",
                                           dir = "/spare/GuessRuns/r8578/HoffmanFBD_Daylength_Sigma",
                                           format = GUESS, 
                                           forcing.data = "CRUJRA")

AllHoffmanFBD_Daylength_Sigma <- defineSource(id = "AllHoffmanFBD_Daylength_Sigma",
                                              name = "AllHoffman FBD, Daylength, Sigma",
                                              dir = "/spare/GuessRuns/r8578/AllHoffmanFBD_Daylength_Sigma",
                                              format = GUESS, 
                                              forcing.data = "CRUJRA")

HoffmanFBD_Sigma <- defineSource(id = "HoffmanFBD_Sigma",
                                 name = "Hoffman FBD, Sigma",
                                 dir = "/spare/GuessRuns/r8578/HoffmanFBD_Sigma",
                                 format = GUESS, 
                                 forcing.data = "CRUJRA")

AllHoffmanFBD_Sigma <- defineSource(id = "AllHoffmanFBD_Sigma",
                                    name = "AllHoffman FBD, Sigma",
                                    dir = "/spare/GuessRuns/r8578/AllHoffmanFBD_Sigma",
                                    format = GUESS, 
                                    forcing.data = "CRUJRA")

AllHoffmanFBD_Sigma_NoSoilMoist <- defineSource(id = "AllHoffmanFBD_Sigma_NoSoilMoist",
                                                name = "AllHoffman FBD, Sigma, NoSoilMoist",
                                                dir = "/spare/GuessRuns/r8578/AllHoffmanFBD_Sigma_NoSoilMoist",
                                                format = GUESS, 
                                                forcing.data = "CRUJRA")

AllHoffmanFBD_Sigma_VPD <- defineSource(id = "AllHoffmanFBD_Sigma_VPD",
                                                name = "AllHoffman FBD, Sigma, VPD",
                                                dir = "/spare/GuessRuns/r8578/AllHoffmanFBD_Sigma_VPD",
                                                format = GUESS, 
                                                forcing.data = "CRUJRA")

##### MAKE RUN LIST #####
runs <- list(
  
  # Base,
  # Daylength,
  # NoWindLimit,
  # HumIgn1.0,
  #Base_r8574,
  #HoffmanFBD,
  #AllHoffmanFBD,
  #HoffmanFBD_Sigma,
  AllHoffmanFBD_Sigma,
  # HoffmanFBD_Daylength,
  # AllHoffmanFBD_Daylength,
  # HoffmanFBD_Daylength_Sigma,
  #AllHoffmanFBD_Daylength_Sigma
  AllHoffmanFBD_Sigma_NoSoilMoist,
  AllHoffmanFBD_Sigma_VPD
  
)

##### DEFINE PLOT GROUPS #####
## Note that the runs must be included in th "runs" list above

plot.groups <- list(
  
  # list(runs = list(Base,
  #                  Daylength,
  #                  NoWindLimit,
  #                  HumIgn1.0,
  
)

##### DEFINE PLOT GROUPS #####
## Note that the runs must be included in th "runs" list above

comparison.groups <- list(
  
  # list(runs = list(Base,
  #                  Daylength,
  #                  NoWindLimit,
  #                  HumIgn1.0,
  #                  Base_r8574),
  #      name = "r8572",
  #      id = "r8572"
  # ),
  # list(runs = list(Base_r8574,
  #                  HoffmanFBD,
  #                  AllHoffmanFBD),
  #      name = "FBD",
  #      id = "FBD"
  # ),
  #list(runs = list(Base_r8574,
  #                 Daylength,
  #                 HoffmanFBD,
  #                 HoffmanFBD_Daylength,
  #                 AllHoffmanFBD,
  #                 AllHoffmanFBD_Daylength),
  #     name = "FBD_Daylength",
  #     id = "FBD_Daylength"
  #)
  # list(runs = list(#Base_r8574,
  #   # Daylength,
  #   HoffmanFBD,
  #   HoffmanFBD_Sigma,
  #   # HoffmanFBD_Daylength,
  #   # HoffmanFBD_Daylength_Sigma,
  #   AllHoffmanFBD,
  #   AllHoffmanFBD_Sigma#,
  #   # AllHoffmanFBD_Daylength,
  #   # AllHoffmanFBD_Daylength_Sigma
  # ),
  # name = "Sigma2",
  # id = "Sigma2"
  # )
  # list(runs = list(
  #   AllHoffmanFBD_Sigma,
  #   AllHoffmanFBD_Sigma_NoSoilMoist
  # ),
  # name = "No Soil Moisure",
  # id = "NoSoilMoist",
  list(runs = list(
    AllHoffmanFBD_Sigma,
    AllHoffmanFBD_Sigma_NoSoilMoist,
    AllHoffmanFBD_Sigma_VPD
  ),
  name = "Alt. Fuel Moist.",
  id = "VPD"
  )
)



##### DEFINE RUN PAIRS #####
## These runs will be compared directly against each otherS

difference.pairs <- list(
  #list("base" = Base, "new" = Daylength),
  # list("base" = Base, "new" = NoWindLimit),
  #list("base" = Base, "new" = HumIgn1.0),
  #list("base" = Base_r8574, "new" = HoffmanFBD),
  #list("base" = Base_r8574, "new" = AllHoffmanFBD)
)


##########################################################################################################
################ HERE CHOOSE WHICH VARIABLES AND PFT AGGREGATES TO PLOT

# plot some variables for each run
vars.to.plot <- list(
  
  #lai = list(var = "lai", cuts = seq(0, 10, 0.5), PFT = TRUE)#,
  #fpc = list(var = "fpc", cuts = seq(0, 1.3, 0.1), PFT = TRUE),
  #mfirefrac = list(var = "mfirefrac", cuts = c(0,0.002,0.005,0.01,0.02,0.05,0.10,0.2,0.50,1), PFT = FALSE)
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
                            dir = file.path(DGVMData.dir, "HandP_PNV/HD"),
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
    
    
    ### always make a summary plot of all layers 
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
      
      # sum to annual and plot that
      this.model.field.annual <- aggregateSubannual(this.model.field, "sum")
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
        }
        
      }
      
      if(doSeasonal) {
        
        # calculate seasonal averages
        this.model.field.seasonal <- aggregateSubannual(this.model.field, "sum", target = "Season")
        
        # plot all seasons together
        season.plot <- plotSpatial(this.model.field.seasonal,
                                   cuts = var.cuts,
                                   map.overlay = map.overlay,
                                   text.multiplier = 2,
                                   ncol = nOptCols(this.model.field.seasonal, preferred.rows.more.than.cols, "Spatial")    # calculate optimal number of columns, based on having one more row than columns
        )
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
          season.plot <- season.plot + facet_grid(cols = vars(Source), rows = vars(Season, switch = "y"))
          magicPlot(season.plot, 
                    filename = file.path(local.group.dir, paste(var.str, "Seasons", "Sum", paste0(first.year, "-", last.year), analysis.label, sep = ".")))
          
          # plot seasons individually
          for(season in all.seasons) {
            
            
            
            season.plot <- plotSpatial(local.field.seasonal.list, 
                                       seasons = season@id,
                                       cuts = var.cuts,
                                       map.overlay = map.overlay,
                                       text.multiplier = 2                                       ,
                                       ncol = nOptCols(length(local.field.seasonal.list), preferred.rows.more.than.cols, ntimes = 1))
            magicPlot(season.plot, 
                      filename = file.path(local.group.dir, paste(var.str, season@id, "Sum", paste0(first.year, "-", last.year), analysis.label, sep = ".")))
            
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
  
  
  ######################################################
  #### NOW DIFFERENCE PLOTS
  
  if(doDifferencePlots) {
    
    for(this.pair in difference.pairs) {
      
      base <- this.pair$base
      new <- this.pair$new
      base.field <- field.list[[base@id]]
      new.field <- field.list[[new@id]]
      
      # make a directory for the comparison
      local.comparison.dir <- file.path(plot.dir, "Differences", paste(new@id, base@id, sep ="-"))
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
        this.comparison <- compareLayers(new.field, base.field, layers1 = layer, show.stats = FALSE, override.quantity = TRUE)
        
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
          this.comparison <- compareLayers(field.annual.list[[new@id]], field.annual.list[[base@id]], layers1 = layer, show.stats = FALSE, override.quantity = TRUE)
          
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
                                                    text.multiplier = 2)
            abs.layer.plot <- abs.layer.plot + facet_grid(cols = vars(Season), rows = vars(Source), switch = y)
            
            magicPlot(abs.layer.plot, 
                      filename = file.path(local.comparison.dir, paste(layer, quant@id, "SeasonalSum", "2-up", analysis.label, paste0(first.year, "-", last.year), sep = ".")))
            
            # plot seasons individually
            for(season in all.seasons) {
              
              # make the annual comparison layer
              this.comparison <- compareLayers(selectSeasons(new.field, season@abbreviation), selectSeasons(base.field, season@abbreviation), layers1 = layer, show.stats = FALSE)
              
              diff.layer.plot <- plotSpatialComparison(this.comparison, 
                                                       map.overlay = "world",
                                                       text.multiplier = 2,
                                                       ncol = 1)
              magicPlot(diff.layer.plot, 
                        filename = file.path(local.comparison.dir, paste(layer, quant@id, season@abbreviation, "Diff", analysis.label, paste0(first.year, "-", last.year), sep = ".")))
              
              # plot the values
              abs.layer.plot <- plotSpatialComparison(this.comparison, 
                                                      type = "values",
                                                      map.overlay = "world",
                                                      text.multiplier = 2,
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
    
    
    # calculate the biomes from the model output
    model.biomes <- getScheme(source = run, scheme = biome.scheme, first.year = 1961, last.year = 1990, year.aggregate.method = "mean", read.full = reread)
    biome.objects.list[[run@id]] <- model.biomes
    
    
    # comapare biomes
    comparison.layer <- compareLayers(model.biomes, biome.data.field, layers1 = biome.scheme@id, tolerance = 0.1, show.stats = FALSE, override.quantity = TRUE)
    
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
                                      title = NULL,
                                      ncol = nOptCols(x = length(local.biomes.list), rows.more.than.cols = 2))

      biome.group.plot <- biome.group.plot + theme(legend.position='bottom')
      biome.group.plot <- biome.group.plot + theme(legend.title=element_blank())
      biome.group.plot <- biome.group.plot + guides(fill=guide_legend(ncol=2))
      
      magicPlot(biome.group.plot, 
                height = 1100,
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
      magicPlot(biome.values.plot, 
                filename = file.path(local.comparison.dir, paste(biome.scheme@id, "2-up", analysis.label, paste0(first.year, "-", last.year), sep = ".")))
      
      print(paste0("*** Done biome difference: ", paste(new@id, base@id, sep ="-")))
      
    } # end for each difference.pair
    
  } # end if doDifference plots
  
  print("*** Done biomes")
  
} # end if do biomes



t2 <- Sys.time()
print(t2-t1)

