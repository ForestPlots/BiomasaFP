% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Mortality_rate_functions.r
\name{life.and.death}
\alias{life.and.death}
\alias{unobs.recs}
\title{Calculation of instantanious mortality and recruitment rates - used in SummaryAGWP2}
\usage{
life.and.death(time, mort.pc, rec.area)

unobs.recs(dat)
}
\arguments{
\item{time}{Census interval}

\item{mort.pc}{Per capita mortaility rate}

\item{rec.area}{Area-based recruitment rate}

\item{dat}{Dataframe returned by \code{CalcAGB}.}
}
\description{
Accessory functions used by SummaryAGWP2.
\code{life.and.death} is called within \code{unobs.recs}.
}
