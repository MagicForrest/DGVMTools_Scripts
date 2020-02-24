library(DGVMTools)
library(Cairo)
source("~/Projects/DGVMTools/Additional/plotUtils.v1.0.R")


# All model runs
model.runs <- list (
  LasslopWind = defineSource(id = "LasslopWind",
                             name = "Lasslop wind limit",
                             dir = "/home/forrest/LocalRuns/FireMIP/PNV_SPITFIRE_WindLasslop/",
                             format = GUESS)
  ,
  NoWind = defineSource(id = "NoWind",
                        name = "No wind limit",
                        dir = "/home/forrest/LocalRuns/FireMIP/PNV_SPITFIRE_WindNoLimit/",
                        format = GUESS)
  ,
  AndrewsWind = defineSource(id = "AndrewsWind",
                             name = "Andrews wind limit",
                             dir = "/home/forrest/LocalRuns/FireMIP/PNV_SPITFIRE_WindAndrews/",
                             format = GUESS)
  ,
  RothermelWind = defineSource(id = "RothermelWind",
                               name = "Rothermel wind limit",
                               dir = "/home/forrest/LocalRuns/FireMIP/PNV_SPITFIRE_WindRothermel/",
                               format = GUESS)
  # ,
  # Dev = defineSource(id = "Dev",
  #                           name = "Dev (Lasslop wind limit, original FBD )",
  #                           dir = "/home/forrest/LocalRuns/FireMIP/PNV_SPITFIRE_WindLasslop_OriginalFBD//",
  #                           format = GUESS)
)

# data sets
GFED4.Source <- defineSource(id = "GFED4",
                             name = "GFED4",
                             dir  = "/data/forrest/DGVMData/GFED4/HD/",
                             DGVMData)


# analysis name 
analysis.name <- "WindLimits"

# plot base directory
plot.base.dir <- "/home/forrest/Projects/SPITFIRE/Results/FireMIP2019/LocalRuns"
if(!file.exists(plot.base.dir)){dir.create(plot.base.dir)}

# make a sub directory for the plots and also a pdf
plot.dir <- file.path(plot.base.dir, analysis.name)
if(!file.exists(plot.dir)){dir.create(plot.dir)}
pdf.booklet <- file.path(plot.base.dir, paste(analysis.name, "pdf", sep = "."))
pdf(file = pdf.booklet , onefile = TRUE, width = 10, height = 5)


# read the gridcells list
gridcells <- read.table("/home/forrest/ModelSetups/FireMIP/gridlist_SPITFIRE_devsites.txt", stringsAsFactors = FALSE)
names(gridcells) <- c("Lon", "Lat", "Site")

# set the colours
preferred.plot.cols <- c("red3", "deepskyblue", "green2", "gold", "magenta")
source.cols <- c()

seasonal.vars <-c("mfirefrac") 

for(this.var in seasonal.vars) {
  
  
  
  # read the model runs
  all.Fields <- list()
  nFields <- 0
  for(this.run in model.runs) {
    all.Fields[[length(all.Fields)+1]] <- getField(source = this.run,
                                                   var = this.var, 
                                                   first.year = 1995,
                                                   last.year = 2013, 
                                                   write = TRUE, 
                                                   read.full = TRUE)
    nFields <- nFields +1
    source.cols[[this.run@name]] <- preferred.plot.cols[[nFields]]
  }
  
  
  # also include GFED4 if appropriate
  if(this.var == "mfirefrac") {
    all.Fields[[length(all.Fields)+1]] <- getField(source = GFED4.Source,
                                                   var = "burntfraction_std",
                                                   spatial.extent = gridcells[, c(1,2)], 
                                                   spatial.extent.id = "SPITFIRE-Dev_Gridcells", 
                                                   write = TRUE, 
                                                   read.full = TRUE)
    
    all.Fields[[length(all.Fields)]]@quant <- lookupQuantity("mfirefrac", context = GUESS@quantities)
    all.Fields[[length(all.Fields)]]@data <- setnames(all.Fields[[length(all.Fields)]]@data, "Total", "mfirefrac")
    nFields <- nFields +1
    source.cols[["GFED4"]] <- "grey25" 
  }
  
  this.var.Fields <- list()
  for(this.gridcell.index in 1:NROW(gridcells)) {
    
    this.gridcell.Fields <- list()
    for(this.Field in all.Fields) {
      
      this.gridcell.Fields[[length(this.gridcell.Fields)+1]] <- selectGridcells(this.Field,
                                                                                c(gridcells$Lon[this.gridcell.index], gridcells$Lat[this.gridcell.index]),
                                                                                spatial.extent.id = gridcells$Site[this.gridcell.index])
    }
    
    # plot the last three
    single.gridcell.plot <- plotSubannual(this.gridcell.Fields,
                        col.by = "Source",
                        cols = source.cols,
                        scales = "free")
    print(single.gridcell.plot)
    single.gridcell.plot <-  single.gridcell.plot + theme(text = element_text(size = theme_get()$text$size * 4))
    magicPlot(single.gridcell.plot, filename = file.path(plot.dir, paste(this.var, gsub(" ", "_", gridcells$Site[this.gridcell.index]), analysis.name, sep = ".")), height = 900, width = 1800)
   
    this.var.Fields <- append(this.var.Fields, this.gridcell.Fields)
  }
  
  summary.plot <- plotSubannual(this.var.Fields, 
                      col.by = "Source", 
                      cols = source.cols, 
                      scales = "free")
  print(summary.plot)
  magicPlot(single.gridcell.plot, filename = file.path(plot.dir, paste(this.var, "Summary", analysis.name, sep = ".")), height = 900, width = 1800)
  
}

dev.off(); dev.off()