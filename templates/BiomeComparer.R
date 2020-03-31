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

# biomes to plot
biomes <- list(dataset = "HandPBiomes", classification = Smith2014BiomeScheme)


##########################################################################################################
################ HERE CHOOSE WHICH BIOMES TO PLOT

biomes <- list(dataset = "HandPBiomes", classification = Smith2014BiomeScheme)

# resolution
resolution <- "HD"

# Analysis label and plot directory
analysis.label <- "Test"
plot.dir <- "/home/matthew/Projects/FireMIP/plots/spatial/Test"
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


##### MAKE RUN LIST #####
runs <- list(
  
  Base,
  Daylength,
 
  
)

##### DEFINE PLOT GROUPS #####
## Note that the runs must be included in th "runs" list above

run.groups <- list(
  
  list(runs = list(Base,
                    Daylength),
       name = "Example Plot Group",
       id = "Example")
  
)


##### DEFINE RUN PAIRS #####
## These runs will be compared directly against each otherS

run.pairs <- list(
  list("base" = Base, "new" = Daylength),
)






##########################################################################################################
################ NORMALLY NOTHING TO CHANGE AFTER HERE!




# get the biome dataset
biome.scheme <- biomes$classification

biome.src <- defineSource(id = paste(biomes$dataset, resolution, sep = "."),
                          name = "PNV Vegetation Types", 
                          dir = file.path(DGVMData.dir, "HandP_PNV/HD"),
                          format = DGVMData) 

biome.data.field <- getField(source = biome.src, 
                             var = biomes$classification,
                             verbose = FALSE)


#######################################################
#### NOW BIOMES


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
            height = 1200,
            filename = file.path(run.plot.dir, paste(biome.scheme@id, "BiomesWithData", paste0(first.year, "-", last.year), analysis.label, sep = ".")))
  
} # end for each run


if(doMultiPlots){
  
  ### do each comparison group in turn
  for(group in run.groups) {
    
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
              height = 1200,
              filename = file.path(local.group.dir, paste(biome.scheme@id, group$id, analysis.label, paste0(first.year, "-", last.year), sep = ".")))
    
    print(paste0("*** Done biome group: ", group$id))
    
  } # end for each group
  
} # end if do multiplots 



######################################################
#### NOW DIFFERENCE PLOTS

if(doDifferencePlots) {
  
  for(this.pair in run.pairs) {
    
    base <- this.pair$base
    new <- this.pair$new
    base.field <- biome.objects.list[[base@id]]
    new.field <- biome.objects.list[[new@id]]
    
    # make a directory for the comparison
    local.comparison.dir <- file.path(plot.dir, "Pairs", paste(new@id, base@id, sep ="-"))
    dir.create(local.comparison.dir, showWarnings = FALSE, recursive = TRUE)
    
    # make the comparison layer
    this.comparison <- compareLayers(new.field, base.field, layers1 = biome.scheme@id, show.stats = FALSE)
    
    # plot the difference
    biome.diff.plot <- plotSpatialComparison(this.comparison, 
                                             map.overlay = "world",
                                             text.multiplier = 2)
    magicPlot(biome.diff.plot, 
              height = 1200,
              filename = file.path(local.comparison.dir, paste(biome.scheme@id, "Diff", analysis.label, paste0(first.year, "-", last.year),  sep = ".")))
    
    
    # plot the values
    biome.values.plot <- plotSpatialComparison(this.comparison, 
                                               type = "values",
                                               map.overlay = "world",
                                               text.multiplier = 2)
    magicPlot(biome.values.plot, 
              height = 1200,
              filename = file.path(local.comparison.dir, paste(biome.scheme@id, "2-up", analysis.label, paste0(first.year, "-", last.year), sep = ".")))
    
    print(paste0("*** Done biome difference: ", paste(new@id, base@id, sep ="-")))
    
  } # end for each difference.pair
  
} # end if doDifference plots

print("*** Done biomes")




t2 <- Sys.time()
print(t2-t1)

