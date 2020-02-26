library(Cairo)
library(DGVMTools)

magicPlot <- function(p, filename, type = "png", height = NULL, width = NULL, dpi = "auto", bg = "transparent"){
  
  if(!is.null(p))  {ar <- getGGAR(p)}
  else {ar <- 1}
  
  # set a basic height for scaling if nothing provided
  if(missing(height) && missing(width)){
    if(type == "png") {height <- 700 }
    else if(type == "pdf") {height <- 7}
  }
  
  # scale to height
  if(is.null(width)) width <- height * ar
  
  # or scale to width
  if(is.null(height)) height <- width * ar
  
  # make and write the plot
  Cairo(file = paste(filename, type, sep ="."), type = type, width=width, height = height, dpi = dpi, bg = bg )
  print(p)
  dev.off()
  
}

getGGAR <- function(p){
  
  
  # Get the axes on the panels in terms of plot units
  y.range <- layer_scales(p)$y$range$range
  x.range <- layer_scales(p)$x$range$range
  y.size <- y.range[2] - y.range[1]
  x.size <- x.range[2] - x.range[1]
  
  #print(y.size)
  #print(x.size)
  
  # Get the numer of facet
  p.built <- ggplot_build(p) # do this once to for efficiency
  all.facets <- c()
  for(table.counter in 1:length(p.built$data)){
    all.facets <- append(all.facets, unique(p.built$data[[table.counter]]$PANEL) )
  }
  n.facets <- length(unique(all.facets))
  #print(n.facets)
  
  # Given the number of facets and the nrow and ncol parameters, determine the numbers of rows and columns
  par <- p.built$layout$facet$params
  rows_cols <- wrap_dims(n.facets, par$nrow, par$ncol)
  #print(rows_cols)
  
  # And combine the rowsxcols with the panel sizes to make an approximate aspect ratio
  aspect.ratio <- (x.size * rows_cols[2]) / (y.size * rows_cols[1])
  
  return(aspect.ratio)
  
}

#' @param x Either the number of panels or something from which the number of panels can be derived
#' currently possible options are: a simple numeric (givien the number of panels) or a Field (from 
#' from which the number of panels are derived from the layers and additional possible facets)
#' @param rows.more.than.cols Preffered number of rows more than columns as an integer (can be negative).  For
#' example a value of "1" will optimise for more row than columns, "2" for two more rows, "0" for the same 
#' number of rows as columns, "-1" for one less row than columns.  Note, the function doesn't guarantee 
#' that it will work out like this, just that it will try to respect this preference if there is a choice.
#' @param type The type of plot as a character string (used for determing panels from dimensions).  Currently only
#' implemented for "Spatial".
#' @param nlayers Numeric for number of layers (in case only subset a subset of layers in x are plotted)
#' @param ntimes Numeric for number of time periods (Years and/or Days/Months/Seasons) (in case only subset a subset of
#' time periods in x are plotted)
#' 

nOptCols <- function(x, rows.more.than.cols = 0, type = "Spatial", nlayers, ntimes) {
  
  # simple numeric
  if(is.numeric(x)) {
    npanels <- x
  }
  #  a single field (assume plotting all layers and Days/Months/Seasons/Years)
  else if(is.Field(x)){
    
    if(tolower(type) == "spatial") {
      if(missing(nlayers)) nlayers <- length(layers(x))
      if(missing(ntimes)) {
        ntimes <- 1
        if("Day" %in% getDimInfo(x)) ntimes <- ntimes * getDimInfo(x, "size")$Day
        else if("Month" %in% getDimInfo(x)) ntimes <- ntimes * getDimInfo(x, "size")$Month
        else if("Season" %in% getDimInfo(x)) ntimes <- ntimes * getDimInfo(x, "size")$Season
        if("Year" %in% getDimInfo(x)) ntimes <- ntimes * getDimInfo(x, "size")$Year
      }
      npanels <- ntimes * nlayers
      
    }
    
    else {
      stop(paste0("nOptCols not implemented for type = ", type))
    }
    
  }
  else {
    stop(paste0("nOptCols not implemented for class = ", class(x)[1]))
  }
  
  opt.ncols <- 1
  while(opt.ncols * (opt.ncols + rows.more.than.cols) < npanels) {
    opt.ncols <- opt.ncols + 1
  }
  if(rows.more.than.cols < 0) opt.ncols <- opt.ncols + rows.more.than.cols
  
  if(opt.ncols >= 10) warning("More than 10 columns for optimal facetting, that is a lot")
  return(opt.ncols)
  
} 



