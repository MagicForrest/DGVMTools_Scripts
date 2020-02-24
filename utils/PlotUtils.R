library(Cairo)

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
