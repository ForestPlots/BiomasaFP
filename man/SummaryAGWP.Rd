% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SummaryAGWP.r
\name{SummaryAGWP}
\alias{SummaryAGWP}
\title{SummaryAGWP}
\usage{
SummaryAGWP(
  xdataset,
  AGBEquation,
  dbh = "D4",
  rec.meth = 0,
  height.data = NULL,
  palm.eq = TRUE,
  treefern = TRUE
)
}
\arguments{
\item{xdataset}{Object returned by \code{\link{mergefp}}.}

\item{AGBEquation}{Allometric equation function to use when estimating AGB, for example \code{\link{AGBChv14}}.}

\item{dbh}{Name of column containing diameter data. Default is "D4".}

\item{rec.meth}{Method used to estimate AGWP of recruits. If 0 (default), estimates growth from starting diameter of 0mm. If another value is provided, then growth is estimated from a starting diameter of 100mm.}

\item{height.data}{Dataframe containing model type and parameters for height diameter model. These are matched on to the xdataset object at tree level. If NULL (default), then regional height-diameter equations (from Felpaush et al. 2012) are used.}

\item{palm.eq}{Logical. If TRUE the family level diameter based equation from Goodman et al 2013 is used to estimate AGB of monocots. Otherwise AGB of monocots is estimated using the allometric equation supplied to AGBEquation.}

\item{treefern}{Logical. If TRUE the height based equation from Tiepolo et al. 2002 is used to estimate the biomass of treeferns. If FALSE biomass is esimated using the allometric equation supplied to AGBEquation.}
}
\value{
A data frame with PlotViewID, CensusNo, and observed and unobserved elements of AGWP, stem dynamics and AGB mortality.
}
\description{
Function to estimate AGWP by plotview and census interval, including estimating unobserved recruitment and growth of trees that died between census periods using the Talbot et al. (2014) census interval correction.
}
\references{
Talbot et al. 2014. Methods to estimate aboveground wood productivity from long-term forest inventory plots. Forest Ecology and Management 320: 30-38.
}
\author{
Martin Sullivan, Gabriela Lopez-Gonzalez
}
