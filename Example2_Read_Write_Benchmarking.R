#### Load package
library(DGVMTools)


##### Define the ModelRun
run <- defineSource(  name = "Example Run",
                      dir = "/home/forrest/Workshops/2019-03-13_Lund/LPJ-GUESS",
                      format = "LPJ-GUESS")


##### Benchmark first reading of the full file
##### File is gzipped, 119 Mb gzipped, 836 Mb uncompressed

t1 <- Sys.time()
field1 <- getField(source = run,
                        var = "lai",
                        year.aggregate.method = "mean",
                        read.full = TRUE,
                        write = TRUE)
t2 <- Sys.time()
print(t2-t1)

print(plotSpatial(field1))



##### Benchmark re-read

t1 <- Sys.time()
field2 <- getField(source = run,
                        var = "lai",
                        year.aggregate.method = "mean",
                        read.full = FALSE)
t2 <- Sys.time()
print(t2-t1)

print(plotSpatial(field2))



##### Save as a netCDF 

field.test <- getField(source = run,
                   var = "lai",
                   read.full = TRUE,
                   first.year = 1950,
                   last.year = 1970,
                   write = TRUE)

writeNetCDF(x = field.test, filename = "~/example.lai.nc")

rm(field1, field2, field.test)
