# BiomasaFP
Functions to analyse data from www.ForestPlots.net.

**BiomasaFP** can:

* Apply inbuilt and user generated allometric equations to estimate aboveground biomass of individual trees

* Calculate plot-level aboveground biomass, woody productivity and stem dynamics

* Fit height-diameter models

Please see https://martin-sullivan.shinyapps.io/BiomasaFP/ to explore how to use the package.

Visit www.ForestPlots.net for information about the ForestPlots database and how to get involved.

## Getting started
BiomasaFP is designed to be used with outputs from the ForestPlots.net dataset. The functions assume that your data has certain columns as found in the output of ForestPlots.net. This means you can do analyses with relatively little user input if your data has been output from ForestPlots.net. If your data are not in ForestPlots (why???!), you might find the BIOMASS R package (https://cran.r-project.org/web/packages/BIOMASS/index.html) easier to use, but note that BIOMASS does not give you woody productivity and stem dynamics (see https://rpubs.com/martinsulli/BiomassAndBiomasa for using the packages together).

The BiomasaFP R package needs three inputs: (1) tree by tree data from the advanced search, (2) metadata from the query library (Basic plot informaton>Plot Information for R Package V1.1) and (3) wood density of individual trees from the query library (Wood Density>Wood Density of Individual trees).

Once you have read these in, the *mergefp* function combines them into a common format that is used by the other functions in the package.

*SummaryAGWP* is arguably the main function in the package - it gives you AGB, AGWP and stem dynamics.

So while there is more you can do with the package (see https://martin-sullivan.shinyapps.io/BiomasaFP/), the basic route from data input to analysis is just a few lines of code.


`trees<-read.csv("TreeData.csv")`

`md <-read.csv("MetaData.csv")`

`wd <- read.csv("WoodDensity.csv")`

`# Merge data`

`dat<-mergefp(trees,md,wd)`

`# Calculate AGB, AGWP etc`

`SummaryAGWP(dat)`



