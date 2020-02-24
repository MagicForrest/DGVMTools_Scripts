# DGVMTools_Scripts

This repostory comprises some additional scripts and documentation for learning and using DGVMTools, as well as some convnience utilities for plotting using Cairo graphics (requires external code (easy to install, https://www.cairographics.org/).

Also check out vignette included in the package, this is the standard R way to offer detailed documentation for packegs. Also look at the function-by-function examples available in R with `example("<function_name>")`.  At time of writing the vignettes are fairly complete, but some functions are missing examples.

## Contents

* **examples** contains the three standard example scripts that I use for workshops, starting from the simple open-read-plot example to more complicated thing
* **templates** contains scripts that can be the basis of fairly complicated analysis.  Although complex, they are highly configurable with all the options to change at the top.  These are the scripts I commonly use during model development/analysis.
* **utils** currently just contains a single function "magicPlot" which is a convenience wrapper to the Cairo graphics backend.  The idea is that you throw the function a plot object (as is returned by ggplot2 for examples) and a file name, and the fucntios writes the plot to the file.  It attempts to guess an appropriate aspect ratio for the plot, which works well for spatial plots but *not at all* for temporal plots.  For temporal plots use the height and width arguments.
* **docs** A few miscellaneous resources for learning DGVMTools.
