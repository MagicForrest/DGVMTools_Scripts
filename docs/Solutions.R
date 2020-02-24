### Load DGVMTools and define the LPJ-GUESS Run that we will use for all tasks
library(DGVMTools) 
LPJ.source <- defineSource(  name = "LPJ-GUESS Example Run",
                      dir = "/home/forrest/LPJ-GUESS_Global_Run",
                      format = "LPJ-GUESS")

### Task 1 - plot global evergreen NPP for 2005 to 2015

LPJ.anpp <- getField(source = LPJ.source, 
                     var = "anpp",
                     first.year = 2005, 
                     last.year = 2015, 
                     year.aggregate.method = "mean")

# lets have a look
print(plotSpatial(LPJ.anpp))

# calculate evergreen total 
LPJ.anpp <- layerOp(x = LPJ.anpp, operator = "+", layers = ".Evergreen", new.layer = "Evergreen")

# plot the evergreen layer
print(plotSpatial(LPJ.anpp, layers = "Evergreen"))

# free up memory
rm(LPJ.anpp)
gc()


### Task 2 - plot the subannual (monthly) cycle of LAI in Frankfurt

LPJ.mlai <- getField(source = LPJ.source, 
                     var = "mlai")

LPJ.mlai.Frankfurt <- selectGridcells(x = LPJ.mlai, 
                                      gridcells = c(8.25, 50.75),
                                      spatial.extent.id = "Frankfurt")

print(plotSubannual(LPJ.mlai.Frankfurt))

print(plotSubannual(LPJ.mlai.Frankfurt, year.col.gradient = TRUE))

# clean up
rm(LPJ.mlai)
gc()


### Task 3 - plot of LPJ-GUESS tree biomass (mean 1990-2010) vs Avitabile biomass

# first get the LPJ-GUESS data (vegC_std)
LPJ.vegC <- getField(source = LPJ.source, 
                     var = "vegC_std", 
                     first.year = 1990, 
                     last.year = 2010, 
                     year.aggregate.method = "mean")

# calculate the Tree sum
LPJ.vegC <- layerOp(x = LPJ.vegC,
                    operator = "+", 
                    layers = ".Tree", 
                    new.layer = "Tree")

# second get the Avitabile biomass
Avitabile.src <- defineSource(dir = "/media/forrest/Data/DGVMData/Avitabile2016/HD",
                              name = "Avitabile Tropical Biomass",
                              format = DGVMData)

Avitabile.vegC <- getField(source = Avitabile.src, 
                     var = "vegC_std")

# lets check it
print(Avitabile.vegC)
print(plotSpatial(Avitabile.vegC))

# now compare 
vegC.comparison <- compareLayers(field1 = LPJ.vegC,
                                 field2 = Avitabile.vegC,
                                 layers1 = "Tree",
                                 layers2 = "Tree")

# also plot difference                      
print(plotSpatialComparison(vegC.comparison))

# and values 
print(plotSpatialComparison(vegC.comparison, type = "values"))

# use ggplot2 to lay it out better
comparison.plot <- plotSpatialComparison(vegC.comparison, type = "values")
comparison.plot <- comparison.plot + facet_wrap(~Facet, ncol = 1)
print(comparison.plot)
