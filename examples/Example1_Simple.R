#### STEP 0 - Load package
library(DGVMTools)
library(raster)


##### STEP 1 - Define the ModelRun
run <- defineSource(name = "Example Run",
                    dir = "/home/forrest/LPJ-GUESS_Global_Run",
                    format = GUESS)


##### STEP 2 - Get the model object - note yearly averaging is also being done

model.field <- getField(source = run,
                               var = "lai",
                               year.aggregate.method = "mean")



##### STEP 3 - Plot the Field
print(plotSpatial(model.field))

##### TADAA.







############# ALSO... #############

# take a look at the Field by printing it
print(model.field)

# convert Field to Raster* object
model.raster <- as.Raster(model.field)

# check it out
print(model.raster)
plot(model.raster)

# also covert to a data.frame
model.df <- as.data.frame(model.field)
print(head(model.df))

rm(run, model.field, model.raster, model.df)
