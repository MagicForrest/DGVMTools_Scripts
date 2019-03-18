library(DGVMTools)

### DEFINE A SOURCE

GUESS.run <- defineSource(id = "LPJ-GUESS_Example",
                          dir = "/home/forrest/Workshops/2019-03-13_Lund/LPJ-GUESS",
                          format = GUESS,
                          name = "LPJ-GUESS Example Run")

###  GET AND PLOT SPATIAL MEAN (AS A TIME SERIES)

LAI.spatial.mean <- getField(source = GUESS.run, 
                             var = "lai",
                             spatial.aggregate.method = "mean")

print(plotTemporal(LAI.spatial.mean))



###  AGGREGATION - LAYERS, SPACE AND YEARS

# read the full dataset
LAI.full <- getField(source = GUESS.run, 
                     var = "lai")

# make tree and grass layers
LAI.full <- layerOp(LAI.full, "+", ".Tree", "Tree")
LAI.full <- layerOp(LAI.full, "+", ".Grass", "Grass")

# aggregate years and plot Tree and Grass
LAI.year.mean <- aggregateYears(LAI.full, method = "mean")
print(plotSpatial(LAI.year.mean, layers = c("Tree", "Grass")))

# aggregate over space and plot Tree and Grass
LAI.spatial.mean <- aggregateSpatial(LAI.full, method = "mean")
print(plotTemporal(LAI.spatial.mean, layers = c("Tree", "Grass")))

# look at the variability through years
LAI.year.var <- aggregateYears(LAI.full, method = "var")
print(plotSpatial(LAI.year.var, layers = c("Tree", "Grass")))


# free up memory
rm(LAI.full)


### BIOME CLASSIFICATION

Global.biomes <- getBiomes(GUESS.run, 
                           Smith2014BiomeScheme, 
                           year.aggregate.method = "mean",
                           first.year = 1961,
                           last.year = 1990)

print(plotSpatial(Global.biomes))





### COMPARISON TO DATA - SAATCHI BIOMASS


# get standard variable vegC_std
GUESS.vegC <- getField(GUESS.run,
                       "vegC_std",
                       first.year = 2000,
                       last.year = 2010,
                       year.aggregate.method = "mean")

# calculate Tree total and have a look
GUESS.vegC <- layerOp(GUESS.vegC, "+", ".Tree", "Tree")
print(plotSpatial(GUESS.vegC, "Tree"))

# define the data Source (Saatchi data @ HD)
Saatchi.dataset <- defineSource(id = "Saatch2011",
                                name = "Saatchi et al. 2011 tropical biomass",
                                format = DGVMData,
                                dir = system.file("extdata", "DGVMData", "Saatchi2011", "HD", package = "DGVMTools"))

# get the data field
Saatchi.vegC <- getField(source = Saatchi.dataset,
                         var = "vegC_std")

# have a look and plot
print(plotSpatial(Saatchi.vegC))

# compare layers to produce a Comparison object
vegC.comparison <- compareLayers(field1 = GUESS.vegC, field2 = Saatchi.vegC, layers1 = "Tree", layers2 = "Tree")

# plot difference map
print(plotSpatialComparison(vegC.comparison))

# plot side-by-side
print(plotSpatialComparison(vegC.comparison, type = "values"))

# make scatter plot
print(plotScatterComparison(vegC.comparison))



### ONE MORE PLOT TYPE: SEASONAL CYCLE 
# London gridcell


# define a gridecll and the get a Field of the monthly LAI for the gridcell
gridcell <- data.frame(Lon = c(0.25), Lat = c(51.25))
London.mlai <- getField(source = GUESS.run, 
                        var = "mlai",
                        spatial.extent = gridcell,
                        spatial.extent.id = "London baby")


# plot subannual cycle of LAI
print(plotSubannual(field = London.mlai,
                    year.col.gradient = TRUE,
                    alpha = 0.5))




### EXTRA
biomes.source <- defineSource(id = "DataBiomes",
                              dir = system.file("extdata", "DGVMData", "HandP_PNV", "HD", package = "DGVMTools"), # this would normally just be a character string containing a path
                              format = DGVMData,
                              name = "Haxeltine and Prentice Biomes")
biomes.data <- getField(source = biomes.source, var = "Smith2014")
print(plotSpatial(biomes.data))

biome.comparison <- compareLayers(field1 = Global.biomes, field2 = biomes.data, layers1 = "Smith2014", layers2 = "Smith2014")

# have a look at this comparison object
print(biome.comparison )

# plot difference map
print(plotSpatialComparison(biome.comparison, type = "difference"), map.overlay = "world")



# plot side-by-side
print(plotSpatialComparison(biome.comparison, type = "values", map.overlay = "world"))

# and make nicer with a bit of ggplot2 action ()
biome.plot <- plotSpatialComparison(biome.comparison, type = "values", map.overlay = "world")

# change panel layout
biome.plot <- biome.plot + facet_wrap(~Facet, nrow = 2)
print(biome.plot)

# move legend and make two columns
biome.plot <- biome.plot + theme(legend.position = "bottom")
biome.plot <- biome.plot + guides(fill = guide_legend(ncol = 3))
print(biome.plot)


rm(LAI.year.var, LAI.spatial.mean, LAI.year.mean)


